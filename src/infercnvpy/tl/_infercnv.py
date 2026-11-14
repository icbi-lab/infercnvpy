import re
import sys
import warnings
from collections import deque
from collections.abc import Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from typing import Any, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from numpy.lib.stride_tricks import sliding_window_view

from infercnvpy._logging import SKIP_OWN_FRAMES, log, log_result
from infercnvpy._numba import n_kernel_chunks, numba_threads
from infercnvpy._util import _as_matrix, _is_backed, _is_dask, _params_for_uns
from infercnvpy.tl._infercnv_kernels import (
    _cnv_csr,
    _cnv_dense,
    _count_kept,
    _fill_kept,
    _gene_values,
    _group_sums_csr,
    _group_sums_dense,
    _n_unsorted_rows,
)


def infercnv(
    adata: AnnData,
    *,
    reference_key: str | None = None,
    reference_cat: None | str | Sequence[str] = None,
    reference: np.ndarray | None = None,
    lfc_clip: float = 3,
    window_size: int = 100,
    step: int = 10,
    dynamic_threshold: float | None = 1.5,
    noise_reference: Literal["reference", "all"] = "reference",
    exclude_chromosomes: Sequence[str] | None = ("chrX", "chrY"),
    chunksize: int = 5000,
    n_jobs: int | None = None,
    lazy: bool | None = None,
    inplace: bool = True,
    layer: str | None = None,
    key_added: str = "cnv",
    calculate_gene_values: bool = False,
) -> None | tuple[dict, sp.csr_matrix, np.ndarray | None]:
    """Infer Copy Number Variation (CNV) by averaging gene expression over genomic regions.

    This method is heavily inspired by `infercnv <https://github.com/broadinstitute/inferCNV/>`_
    but more computationally efficient. The method is described in more detail
    in on the :ref:`infercnv-method` page.

    There, you can also find instructions on how to :ref:`prepare input data <input-data>`.

    The expression matrix can be in memory (numpy array or scipy sparse matrix), a
    :class:`dask.array.Array` (also with sparse blocks), or the matrix of a backed AnnData
    (`anndata.read_h5ad(..., backed="r")`). Lazy matrices are read block by block and never
    loaded entirely, see :ref:`large-data`. By default, the blocks are read ahead in background
    threads and each is processed with `n_jobs` threads; the result is in memory. With `lazy`, the
    result of a dask array stays a dask array instead.

    .. versionchanged:: 0.7.0
        By default, the noise threshold is computed from the reference cells (`noise_reference`),
        not separately for each chunk of `chunksize` cells, so the result no longer depends on
        `chunksize`. The CNV matrix and the gene values are float32.

    Parameters
    ----------
    adata
        annotated data matrix
    reference_key
        Column name in adata.obs that contains tumor/normal annotations.
        If this is set to None, the average of all cells is used as reference.
    reference_cat
        One or multiple values in `adata.obs[reference_key]` that annotate
        normal cells.
    reference
        Directly supply an array of average normal gene expression. Overrides
        `reference_key` and `reference_cat`.
    lfc_clip
        Clip log fold changes at this value
    window_size
        size of the running window (number of genes in to include in the window)
    step
        only compute every nth running window where n = `step`. Set to 1 to compute
        all windows.
    dynamic_threshold
        Remove noise: values within `dynamic_threshold` standard deviations of the
        center of the smoothed, cell-centered gene expression are set to 0, see
        `noise_reference`. Set to `None` to disable this step.
    noise_reference
        Which cells define the noise band of `dynamic_threshold`.

        `"reference"`
            The band follows R inferCNV. It is `mean ± dynamic_threshold * sd`, where
            `mean` is the mean over the reference cells (`reference_cat`) and `sd` is
            the mean of the standard deviations of the individual reference cells.
            Without reference cells (`reference_key` not set, or an explicit
            `reference`), all cells are used.
        `"all"`
            The band is `0 ± dynamic_threshold * sd`, where `sd` is the standard
            deviation over all cells. It grows with the fraction of tumor cells.
    exclude_chromosomes
        List of chromosomes to exclude. The default is to exclude genosomes.
    chunksize
        Number of cells processed at a time for in-memory and backed data. It bounds
        the memory of the intermediate dense matrix (`chunksize` x windows) and does not
        change the result. For a dask array, its blocks are processed instead.
    n_jobs
        Number of threads. Default: all cores (numba's setting, see
        :func:`numba.set_num_threads`). Ignored with `lazy`, where the dask scheduler
        decides and each block runs in a single thread.
    lazy
        Only for a dask array: return the CNV matrix (and the gene values) as lazy dask
        arrays, computed by :meth:`dask.array.Array.map_blocks` with the current dask
        scheduler. The noise band is computed up front with one dask computation.
        Default: `True` if a :class:`distributed.Client` is active, otherwise `False`:
        the blocks are read ahead in background threads with the synchronous scheduler and
        the result is in memory.
    inplace
        If True, save the results in adata.obsm, otherwise return the CNV matrix.
    layer
        Layer from adata to use. If `None`, use `X`.
    key_added
        Key under which the cnv matrix will be stored in adata if `inplace=True`.
        Will store the matrix in `adata.obsm["X_{key_added}"] and additional information
        in `adata.uns[key_added]`.
    calculate_gene_values
        If True per gene CNVs will be calculated and stored in `adata.layers["gene_values_{key_added}"]`.
        The value of a gene is the mean of all windows that contain it, so many genes share
        similar values. Not all genes are in a window, because the window size and step
        size are not always a multiple of the number of genes; these and the dropped genes
        are NaN. This needs a dense float32 cells x genes matrix in memory.

    Returns
    -------
    If `inplace=True`, returns `None` and sets the fields below. Otherwise
    returns a tuple `(chr_pos, cnv_matrix, gene_values)` with the same values;
    `gene_values` is `None` unless `calculate_gene_values=True`.
    Sets the following fields:

    `adata.obsm["X_{key_added}"]` : :class:`scipy.sparse.csr_matrix`
        Smoothed and denoised CNV values (float32), cells x genomic windows, ordered by
        chromosome and genomic position. With `lazy`, a :class:`dask.array.Array` of
        CSR blocks.
    `adata.uns[key_added]["chr_pos"]` : :class:`dict`
        Maps each chromosome to the column of `X_{key_added}` where it starts.
    `adata.uns[key_added]["window_start"]`, `adata.uns[key_added]["window_end"]` : :class:`numpy.ndarray`
        Genomic span of each window (int64, base pairs): the start of its first gene and the
        largest end of its genes. Used to select genomic regions, e.g. by :func:`infercnvpy.pl.spatial_cnv`.
    `adata.uns[key_added]["params"]` : :class:`dict`
        Parameters of this call, stored as given. `reference_key` and `reference_cat`
        are `None` when an explicit `reference` is passed.
    `adata.layers["gene_values_{key_added}"]` : :class:`numpy.ndarray`
        Only if `calculate_gene_values=True`: CNV values per gene (float32), `NaN` for
        genes that are in no window. With `lazy`, a :class:`dask.array.Array`.
    """
    if not adata.var_names.is_unique:
        raise ValueError("Ensure your var_names are unique!")
    if {"chromosome", "start", "end"} - set(adata.var.columns) != set():
        raise ValueError(
            "Genomic positions not found. There need to be `chromosome`, `start`, and `end` columns in `adata.var`. "
        )
    if noise_reference not in ("reference", "all"):
        raise ValueError(f'`noise_reference` must be "reference" or "all", got {noise_reference!r}.')
    if window_size < 1 or step < 1 or chunksize < 1:
        raise ValueError("`window_size`, `step` and `chunksize` must be at least 1.")
    if adata.n_obs == 0:
        raise ValueError("infercnv needs at least one cell.")

    params = _params_for_uns(
        # an explicit `reference` overrides `reference_key` and `reference_cat`
        reference_key=reference_key if reference is None else None,
        reference_cat=reference_cat if reference is None else None,
        lfc_clip=lfc_clip,
        window_size=window_size,
        step=step,
        dynamic_threshold=dynamic_threshold,
        noise_reference=noise_reference,
        exclude_chromosomes=exclude_chromosomes,
        chunksize=chunksize,
        layer=layer,
        calculate_gene_values=calculate_gene_values,
    )

    positions: Any = adata.var[["chromosome", "start", "end"]]
    # the lazy `var` of `anndata.experimental.read_lazy` has to be loaded
    var: pd.DataFrame = positions.to_memory() if hasattr(positions, "to_memory") else positions
    var_mask = var["chromosome"].isnull()
    n_unplaced = int(var_mask.sum())
    if exclude_chromosomes is not None:
        var_mask = var_mask | var["chromosome"].isin(exclude_chromosomes)
    if not var.loc[~var_mask, "chromosome"].astype(str).str.startswith("chr").any():
        raise ValueError(
            'None of the chromosome names in `adata.var["chromosome"]` start with "chr". '
            'infercnvpy only uses chromosomes named like "chr1". Add the prefix, e.g. with '
            '`adata.var["chromosome"] = "chr" + adata.var["chromosome"].astype(str)`.'
        )
    expr = _as_matrix(adata.X if layer is None else adata.layers[layer], lazy=True)
    if lazy is None:
        lazy = _is_dask(expr) and _distributed_client_active()
    elif lazy and not _is_dask(expr):
        raise ValueError("`lazy=True` needs the expression matrix as a dask array.")

    with numba_threads(n_jobs):
        source = _Rows(expr, chunksize)
        start = log.debug(
            "infercnv: %d cells x %d genes in %d blocks (%s%s)",
            *(adata.n_obs, adata.n_vars, source.n_blocks, source.kind, ", lazy" if lazy else ""),
        )
        layout = _Layout(var, var_mask.to_numpy(), window_size, step)
        ref = _get_reference(adata, reference_key, reference_cat, reference, layer, source=source, lazy=lazy)
        plan = _Plan.build(layout, ref, lfc_clip, dynamic_threshold, noise_reference)
        band_rows = None
        if (
            noise_reference == "reference"
            and reference is None
            and reference_key is not None
            and reference_cat is not None
        ):
            band_rows = np.flatnonzero(_obs_column(adata, reference_key).isin(_categories(reference_cat)))
        res, per_gene_mtx, band = (_run_lazy if lazy else _run)(plan, source, band_rows, calculate_gene_values)

    if dynamic_threshold is not None:
        log.debug(
            "infercnv: removed values within %.4g ± %.4g (from %s)",
            band[1],
            band[0],
            "all cells" if band_rows is None else f"{len(band_rows)} reference cells",
        )
    msg, args = "%d cells x %d windows on %d chromosomes", [res.shape[0], res.shape[1], len(layout.chr_pos)]
    if n_unplaced:
        msg += ", dropped %d of %d genes without a genomic position"
        args += [n_unplaced, adata.n_vars]
    log_result(
        "infercnv",
        msg,
        *args,
        dropped=n_unplaced,
        total=adata.n_vars,
        remedy="add positions with io.genomic_position_from_gtf or io.genomic_position_from_biomart",
        time=start,
    )

    if inplace:
        adata.obsm[f"X_{key_added}"] = res
        adata.uns[key_added] = {
            "chr_pos": layout.chr_pos,
            "window_start": layout.window_start,
            "window_end": layout.window_end,
            "params": params,
        }
        if calculate_gene_values:
            adata.layers[f"gene_values_{key_added}"] = per_gene_mtx
        return None
    return layout.chr_pos, res, per_gene_mtx


def _natural_sort(l: Sequence):
    """Natural sort without third party libraries.

    Adapted from: https://stackoverflow.com/a/4836734/2340703
    """

    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


def _obs_column(adata: AnnData, key: str) -> pd.Series:
    """`adata.obs[key]` as a pandas Series (the `obs` of `anndata.experimental.read_lazy` holds xarray arrays)."""
    column = adata.obs[key]
    return column if isinstance(column, pd.Series) else pd.Series(np.asarray(column))


def _categories(reference_cat: str | Sequence[str]) -> np.ndarray:
    return np.array([reference_cat] if isinstance(reference_cat, str) else reference_cat)


class _Layout:
    """Genes in genomic order and the geometry of the running windows, built once per call.

    Genes are ordered by chromosome (natural sort; only names starting with "chr", without "chrM") and start
    position. Per chromosome `c`, with `g` genes and window length `n = min(window_size, g)`, the pyramid weights are
    the convolution of boxes of `a` and `b` genes, and `K` windows starting every `step` genes are kept. A chromosome
    with at most `window_size` genes is a single unweighted window (`a = g`, `b = 1`).
    """

    def __init__(self, var: pd.DataFrame, var_mask: np.ndarray, window_size: int, step: int):
        used = var.loc[~var_mask, ["chromosome", "start", "end"]]
        chromosomes = _natural_sort([x for x in used["chromosome"].unique() if x.startswith("chr") and x != "chrM"])
        order, geometry, scatter, starts, ends = [], [], [], [], []
        self.chr_pos: dict[str, int] = {}
        n_genes = n_windows = self.n_covered = 0
        for chromosome in chromosomes:
            # `sort_values` is not stable: genes with the same start keep the order infercnvpy has always used
            genes = used.loc[used["chromosome"] == chromosome].sort_values("start").index.to_numpy()
            g = len(genes)
            if window_size < g:
                n, a = window_size, (window_size + 1) // 2
            else:  # a single unweighted window: box(g) * box(1)
                n = a = g
            b = n + 1 - a
            K = len(range(0, g - n + 1, step))
            pos = used.loc[genes]
            starts.append(pos["start"].to_numpy(dtype=np.int64)[: (K - 1) * step + 1 : step])
            ends.append(sliding_window_view(pos["end"].to_numpy(dtype=np.int64), n)[::step][:K].max(axis=1))
            order.append(var.index.get_indexer(genes))
            geometry.append((n_genes, g, a, b, K, n_windows, n))
            # per gene: the windows `first .. last` that contain it, and its pyramid weight in each
            i = np.arange(g)
            first = np.maximum(0, (i - n + step) // step)
            last = np.minimum(K - 1, i // step)
            scatter.append(
                (n_windows + first, last - first + 1, i - first * step, np.full(g, n), np.full(g, min(a, b)))
            )
            self.chr_pos[chromosome] = n_windows
            self.n_covered += min(g, (K - 1) * step + n)
            n_genes += g
            n_windows += K
        #: column in `var` of each gene in genomic order
        self.order = np.concatenate(order)
        #: genomic start of the first gene and largest end of any gene in each window
        self.window_start = np.concatenate(starts)
        self.window_end = np.concatenate(ends)
        #: position in genomic order of each column of `var`, -1 for genes that are not used
        self.col_map = np.full(len(var), -1, dtype=np.int64)
        self.col_map[self.order] = np.arange(n_genes)
        self.n_vars = len(var)
        self.n_windows = n_windows
        g0s, gs, As, Bs, Ks, w0s, ns = (np.array(x, dtype=np.int64) for x in zip(*geometry, strict=True))
        self.geometry = (g0s, gs, As, Bs, Ks, w0s, step)
        self.step = step
        first, n_in, offset, length, ab = (np.concatenate(x).astype(np.int64) for x in zip(*scatter, strict=True))
        #: largest number of windows that contain one gene
        self.span = int(n_in.max())
        m = np.arange(self.span)
        k = offset[:, np.newaxis] - step * m  # position of the gene in its windows `first`, `first + 1`, ...
        weights = np.minimum(np.minimum(k + 1, length[:, np.newaxis] - k), ab[:, np.newaxis]).astype(np.float64)
        weights[m >= n_in[:, np.newaxis]] = 0
        #: per gene: first window that contains it, and its pyramid weights from there on (0 past its last window)
        self.scatter = (first, weights)
        #: per window: the sum of its pyramid weights
        self.window_weight = np.repeat((As * Bs).astype(np.float64), Ks)
        self.gene_geometry = (gs, ns, Ks, w0s, g0s, step, self.order, self.n_covered)


@dataclass(frozen=True)
class _Plan:
    """Everything the kernels need besides the cells: layout, reference bounds and noise settings."""

    layout: _Layout
    lo: np.ndarray
    hi: np.ndarray
    base: np.ndarray  # clipped log fold change of an expression of 0
    cap: np.floating
    base_sums: np.ndarray  # weighted window sums of `base`
    max_nnz: int  # cells with at most this many values are smoothed as a difference to `base`
    dynamic_threshold: float | None
    noise_reference: str

    @classmethod
    def build(cls, layout: _Layout, reference: np.ndarray, lfc_clip: float, dynamic_threshold, noise_reference):
        # With one reference `lo == hi`: the bounded difference is the plain difference.
        ref = np.asarray(reference)[:, layout.order].astype(np.float32)
        lo, hi = ref.min(axis=0), ref.max(axis=0)
        cap = np.float32(lfc_clip)
        zero = np.float32(0)
        base = np.asarray(np.clip(zero - np.clip(zero, lo, hi), -cap, cap), dtype=np.float32)
        first, weights = layout.scatter
        base_sums = np.zeros(layout.n_windows + weights.shape[1])
        np.add.at(base_sums, first[:, np.newaxis] + np.arange(weights.shape[1]), base[:, np.newaxis] * weights)
        # a value costs `span` additions in the sparse way, the prefix sums about 2 per gene
        n_genes = layout.order.size
        max_nnz = min(int(_SPARSE_COST * n_genes / layout.span), n_genes)
        return cls(
            layout, lo, hi, base, cap, base_sums[: layout.n_windows], max_nnz, dynamic_threshold, noise_reference
        )


def _windows(
    plan: _Plan, X, rows: np.ndarray, band=None, *, out: bool = True, serial: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median-centred windows (float32) of the rows `rows` of `X`, their (sum, sum of squares), and their counts.

    The counts are the numbers of values outside the noise `band`, all 0 without a band. With `out=False`, the windows
    are not kept (the first array has no rows).
    """
    lay = plan.layout
    n_rows = len(rows)
    windows = np.empty((n_rows if out else 0, lay.n_windows), dtype=np.float32)
    rowstats = np.empty((n_rows, 2), dtype=np.float64)
    counts = np.empty(n_rows, dtype=np.int64)
    thr, centre = (np.inf, 0.0) if band is None else band
    n_chunks = n_kernel_chunks(n_rows, serial=serial)
    shared = (
        rows, lay.col_map, plan.lo, plan.hi, plan.base, plan.cap, *lay.geometry, *lay.scatter, lay.window_weight,
        plan.base_sums, plan.max_nnz, thr, centre, windows, rowstats, counts, n_chunks,
    )  # fmt: skip
    if sp.issparse(X):
        _cnv_csr(X.data, X.indices, X.indptr, *shared, serial=serial)
    else:
        _cnv_dense(X, *shared, serial=serial)
    return windows, rowstats, counts


#: cells with fewer values than `_SPARSE_COST * genes / span` take the sparse way (see `_infercnv_kernels`)
_SPARSE_COST = 0.9

#: number of blocks read ahead, each in its own thread (measured on lazily read h5ad files: reads overlap)
_READ_AHEAD = {"dask": 4, "backed": 2}

#: `(threshold, centre)` that removes nothing
_NO_BAND = (-1.0, 0.0)


def _noise_band(plan: _Plan, rowstats: np.ndarray) -> tuple[float, float]:
    """`(threshold, centre)` of the noise band from the per-cell (sum, sum of squares) of the windows."""
    if plan.dynamic_threshold is None:
        return _NO_BAND
    n_windows = plan.layout.n_windows
    sums, squares = rowstats[:, 0], rowstats[:, 1]
    mean = sums.sum() / (sums.size * n_windows)
    if plan.noise_reference == "all":  # standard deviation over all values, around 0
        sd = np.sqrt(max(squares.sum() / (sums.size * n_windows) - mean**2, 0.0))
        return plan.dynamic_threshold * float(sd), 0.0
    # R inferCNV (`clear_noise_via_ref_mean_sd`): mean of the per-cell standard deviations, around the mean
    cell_var = np.maximum(squares - sums**2 / n_windows, 0.0) / max(n_windows - 1, 1)
    return plan.dynamic_threshold * float(np.sqrt(cell_var).mean()), float(mean)


class _CsrBuilder:
    """A CSR matrix filled block by block of rows with the windows outside the noise band.

    The data and index arrays are allocated for a dense result: the operating system only maps the pages that are
    written. `finish` shrinks them in place (`ndarray.resize`, a `realloc`). If that allocation fails, they grow in
    place as needed. The blocks are never stacked.
    """

    def __init__(self, n_rows: int, n_cols: int):
        self.shape = (n_rows, n_cols)
        self.indptr = np.zeros(n_rows + 1, dtype=np.int64)
        idx_dtype = np.int32 if n_cols < np.iinfo(np.int32).max else np.int64
        try:
            self.data = np.empty(n_rows * n_cols, dtype=np.float32)
            self.indices = np.empty(n_rows * n_cols, dtype=idx_dtype)
        except MemoryError:
            self.data = np.empty(0, dtype=np.float32)
            self.indices = np.empty(0, dtype=idx_dtype)
        self.n = 0

    def add(self, windows: np.ndarray, counts: np.ndarray, band: tuple[float, float], serial: bool = False):
        """Append the rows `windows`, with `counts` values outside the noise `band` per row."""
        r0, r1 = self.n, self.n + len(windows)
        indptr = self.indptr[r0 : r1 + 1]
        np.cumsum(counts, out=indptr[1:])
        indptr[1:] += indptr[0]
        nnz = int(indptr[-1])
        if nnz > self.data.size:  # room for the remaining rows at the density so far, and at least 25 % more
            size = max(nnz, int(nnz / r1 * self.shape[0] * 1.05), self.data.size * 5 // 4)
            self.data.resize(size, refcheck=False)
            self.indices.resize(size, refcheck=False)
        _fill_kept(windows, *band, indptr, self.data, self.indices, serial=serial)
        self.n = r1

    def finish(self) -> sp.csr_matrix:
        nnz = int(self.indptr[-1])
        self.data.resize(nnz, refcheck=False)
        self.indices.resize(nnz, refcheck=False)
        idx_dtype = np.int32 if max(nnz, self.shape[1]) < np.iinfo(np.int32).max else np.int64
        indices = self.indices.astype(idx_dtype, copy=False)
        return sp.csr_matrix((self.data, indices, self.indptr.astype(idx_dtype)), shape=self.shape, copy=False)


def _to_csr(windows: np.ndarray, counts: np.ndarray, band: tuple[float, float], serial: bool = False) -> sp.csr_matrix:
    """The entries of `windows` outside the noise band (`counts` per row) as a CSR matrix."""
    builder = _CsrBuilder(*windows.shape)
    builder.add(windows, counts, band, serial=serial)
    return builder.finish()


def _gene_values_of(plan: _Plan, out: np.ndarray, band: tuple[float, float], serial: bool = False) -> np.ndarray:
    gv = np.full((out.shape[0], plan.layout.n_vars), np.nan, dtype=np.float32)
    n_chunks = n_kernel_chunks(out.shape[0], serial=serial)
    _gene_values(out, *plan.layout.gene_geometry, *band, gv, n_chunks, serial=serial)
    return gv


def _run(plan: _Plan, source: "_Rows", band_rows: np.ndarray | None, gene_values: bool):
    """Noise band (a pass over `band_rows`, or all cells, unless there is a single block), then the CNV matrix.

    Blocks are read ahead in threads and each is processed with all numba threads. In memory, the band comes from a
    single kernel call over its rows.
    """
    band: tuple[float, float] | None
    if plan.dynamic_threshold is not None and source.n_blocks > 1:
        if source.kind == "in memory":
            rows = np.arange(source.n_obs) if band_rows is None else band_rows
            rowstats = _windows(plan, source.X, rows, out=False)[1]
        else:
            rowstats = np.concatenate([_windows(plan, X, rows, out=False)[1] for X, rows in source.rows_of(band_rows)])
        band = _noise_band(plan, rowstats)
    else:
        band = None if plan.dynamic_threshold is not None else _NO_BAND
    builder = _CsrBuilder(source.n_obs, plan.layout.n_windows)
    gv = np.full((source.n_obs, plan.layout.n_vars), np.nan, dtype=np.float32) if gene_values else None
    for X, rows in source.rows_of(None):
        windows, rowstats, counts = _windows(plan, X, rows, band)
        if band is None:  # a single block defines the band
            band = _noise_band(plan, rowstats if band_rows is None else rowstats[band_rows])
            _count_kept(windows, *band, counts)
        if gv is not None:
            gv[builder.n : builder.n + len(rows)] = _gene_values_of(plan, windows, band)
        builder.add(windows, counts, band)
    assert band is not None
    return builder.finish(), gv, band


def _run_lazy(plan: _Plan, source: "_Rows", band_rows: np.ndarray | None, gene_values: bool):
    """The CNV matrix (and gene values) as lazy dask arrays; the serial kernels run inside dask tasks.

    The noise band is computed first, by one dask computation over the blocks that hold `band_rows` (or all blocks).
    """
    import dask

    X = source.X
    band = _NO_BAND
    if plan.dynamic_threshold is not None:
        blocks = X.to_delayed().ravel()
        tasks = [dask.delayed(_dask_stats_task)(blocks[i], plan, local) for i, local in source.blocks_of(band_rows)]
        band = _noise_band(plan, np.concatenate(dask.compute(*tasks)))
    lay = plan.layout
    res = X.map_blocks(
        partial(_lazy_cnv_block, plan, band),
        dtype=np.float32,
        chunks=(X.chunks[0], (lay.n_windows,)),
        meta=sp.csr_matrix((0, lay.n_windows), dtype=np.float32),
    )
    gv = None
    if gene_values:
        gv = X.map_blocks(
            partial(_lazy_gene_values_block, plan, band),
            dtype=np.float32,
            chunks=(X.chunks[0], (lay.n_vars,)),
            meta=np.empty((0, lay.n_vars), dtype=np.float32),
        )
    return res, gv, band


def _dask_stats_task(block, plan: _Plan, local: np.ndarray | None) -> np.ndarray:
    X = _kernel_ready(block)
    return _windows(plan, X, np.arange(X.shape[0]) if local is None else local, out=False, serial=True)[1]


def _lazy_cnv_block(plan: _Plan, band, block) -> sp.csr_matrix:
    X = _kernel_ready(block)
    windows, _, counts = _windows(plan, X, np.arange(X.shape[0]), band, serial=True)
    return _to_csr(windows, counts, band, serial=True)


def _lazy_gene_values_block(plan: _Plan, band, block) -> np.ndarray:
    X = _kernel_ready(block)
    return _gene_values_of(plan, _windows(plan, X, np.arange(X.shape[0]), serial=True)[0], band, serial=True)


class _Rows:
    """Row blocks of the expression matrix: in memory, backed (h5ad), or a dask array.

    In memory and backed, blocks have `chunksize` rows; a dask array keeps its own blocks. `rows_of` reads the blocks
    ahead in threads (dask: with the synchronous scheduler) and yields them on this thread, so the kernels can use all
    numba threads.
    """

    def __init__(self, X, chunksize: int):
        self.n_obs = X.shape[0]
        if _is_dask(X):
            self.kind = "dask"
            if len(X.chunks[1]) > 1:
                X = X.rechunk({1: -1})
            sizes = list(X.chunks[0])
        else:
            self.kind = "backed" if _is_backed(X) else "in memory"
            if self.kind == "in memory":
                X = _kernel_ready(X)
            sizes = [chunksize] * (self.n_obs // chunksize) + [self.n_obs % chunksize] * bool(self.n_obs % chunksize)
        self.X = X
        self.bounds = np.concatenate([[0], np.cumsum(sizes, dtype=np.int64)])

    @property
    def n_blocks(self) -> int:
        return len(self.bounds) - 1

    def blocks_of(self, rows: np.ndarray | None) -> Iterator[tuple[int, np.ndarray | None]]:
        """`(block, local row indices)` of the blocks that hold any of `rows` (sorted), or `(block, None)` for all."""
        for i in range(self.n_blocks):
            if rows is None:
                yield i, None
                continue
            s, e = self.bounds[i], self.bounds[i + 1]
            local = rows[np.searchsorted(rows, s) : np.searchsorted(rows, e)] - s
            if len(local):
                yield i, local

    def rows_of(self, rows: np.ndarray | None) -> Iterator[tuple[Any, np.ndarray]]:
        """`(X, row indices into X)` of every block that holds any of `rows` (sorted; `None`: all), in order.

        `X` is an in-memory, kernel-ready matrix that holds the block.
        """
        blocks = list(self.blocks_of(rows))
        for (i, local), (X, row0) in zip(blocks, self.read_ahead([i for i, _ in blocks]), strict=True):
            s, e = int(self.bounds[i]), int(self.bounds[i + 1])
            yield X, row0 + (np.arange(e - s) if local is None else local)

    def read_ahead(self, blocks: list[int]) -> Iterator[tuple[Any, int]]:
        """`read` of each of `blocks`. Backed and dask blocks are read ahead in threads while earlier ones are used."""
        if self.kind == "in memory" or len(blocks) < 2:
            yield from map(self.read, blocks)
            return
        n_ahead = _READ_AHEAD[self.kind]
        with ThreadPoolExecutor(n_ahead) as pool:
            ahead: deque = deque()
            for i in blocks:
                ahead.append(pool.submit(self.read, i))
                if len(ahead) > n_ahead:
                    yield ahead.popleft().result()
            while ahead:
                yield ahead.popleft().result()

    def read(self, i: int) -> tuple[Any, int]:
        """Block `i` as an in-memory, kernel-ready matrix, and the row of that matrix where the block starts."""
        s, e = int(self.bounds[i]), int(self.bounds[i + 1])
        if self.kind == "in memory":
            return self.X, s
        block = self.X.blocks[i].compute(scheduler="synchronous") if self.kind == "dask" else self.X[s:e]
        return _kernel_ready(block), 0


def _kernel_ready(X):
    """`X` as a CSR matrix or C-contiguous array of float32 or float64, without copying when it already is one."""
    if hasattr(X, "to_scipy_sparse"):  # pydata/sparse blocks
        X = X.to_scipy_sparse()
    if sp.issparse(X):
        if X.format != "csr":
            X = X.tocsr()
        if X.dtype not in (np.float32, np.float64):
            X = X.astype(np.float32)
        if _canonical(X):
            X.has_canonical_format = True
        else:  # duplicate entries are summed, as by `toarray()`
            X = X.copy()
            X.sum_duplicates()
        return X
    X = np.asarray(X)
    return np.ascontiguousarray(X if X.dtype in (np.float32, np.float64) else X.astype(np.float32))


def _canonical(X) -> bool:
    """Whether the CSR matrix `X` has sorted indices without duplicates (scipy's own check is serial)."""
    return _n_unsorted_rows(X.indices, X.indptr) == 0


def _distributed_client_active() -> bool:
    """Whether a `dask.distributed` client is the default scheduler (without importing distributed)."""
    distributed = sys.modules.get("distributed")
    if distributed is None:
        return False
    try:
        distributed.default_client()
    except ValueError:
        return False
    return True


def _get_reference(
    adata: AnnData,
    reference_key: str | None,
    reference_cat: None | str | Sequence[str],
    reference: np.ndarray | None,
    layer: str | None,
    *,
    source: "_Rows | None" = None,
    lazy: bool = False,
) -> np.ndarray:
    """Parameter validation extraction of reference gene expression.

    If multiple reference categories are given, compute the mean per
    category. The means are float64 sums of the rows divided by their number, computed in
    parallel (`source` is the expression matrix, if it is already wrapped). For a dask array or a
    backed matrix, all means are computed in a single pass over the blocks that hold reference
    cells (with `lazy`, as one dask computation).

    Returns a 2D array with reference categories in rows, genes in columns.
    If there's just one category, it's still a 2D array.
    """
    X: Any = adata.X if layer is None else adata.layers[layer]

    if reference is None:
        cats, codes = None, np.zeros(adata.n_obs, dtype=np.int64)
        if reference_key is None or reference_cat is None:
            warnings.warn(
                "Using the mean of all cells as reference. For better results, "
                "provide either `reference`, or both `reference_key` and `reference_cat`.",
                UserWarning,
                skip_file_prefixes=SKIP_OWN_FRAMES,
            )
        else:
            labels = _obs_column(adata, reference_key)
            cats = _categories(reference_cat)
            cats_in_obs = pd.Index(cats).isin(labels.unique())
            if not np.all(cats_in_obs):
                raise ValueError(
                    "The following reference categories were not found in "
                    f"adata.obs[reference_key]: {cats[~cats_in_obs]}"
                )
            unique = pd.Index(pd.unique(cats))
            codes = unique.get_indexer(labels)
        if source is None:
            source = _Rows(_as_matrix(X, lazy=True), chunksize=5000)
        reference = _means(source, codes, 1 if cats is None else len(unique), lazy=lazy)
        if cats is not None:
            reference = reference[unique.get_indexer(cats)]  # a repeated category repeats its mean

    reference = np.asarray(reference)
    if reference.ndim == 1:
        reference = reference[np.newaxis, :]

    if reference.shape[1] != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference


def _means(source: _Rows, codes: np.ndarray, n_groups: int, *, lazy: bool) -> np.ndarray:
    """Per group, the column means of the rows with `codes == group` (-1: no group), in one pass over the data.

    Blocks without any of these rows are not read.
    """
    rows = None if (codes >= 0).all() else np.flatnonzero(codes >= 0)
    blocks = list(source.blocks_of(rows))
    if source.kind == "in memory":
        rows = np.arange(len(codes)) if rows is None else rows
        sums = _group_sums(source.X, rows, codes[rows], n_groups)
    elif lazy:  # one dask computation, in parallel with the current scheduler
        import dask

        delayed = source.X.to_delayed().ravel()
        tasks = [
            dask.delayed(_dask_block_sums)(delayed[i], *_block_groups(source, i, local, codes), n_groups)
            for i, local in blocks
        ]
        sums = np.sum(dask.compute(*tasks), axis=0)
    else:
        sums = np.zeros((n_groups, source.X.shape[1]))
        for (i, local), (X, row0) in zip(blocks, source.read_ahead([i for i, _ in blocks]), strict=True):
            local, groups = _block_groups(source, i, local, codes)
            sums += _group_sums(X, row0 + local, groups, n_groups)
    counts = np.bincount(codes[codes >= 0], minlength=n_groups)
    return sums / counts[:, np.newaxis]


def _block_groups(source: _Rows, i: int, local: np.ndarray | None, codes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rows of block `i` (local indices) that belong to a group, and their groups."""
    s, e = int(source.bounds[i]), int(source.bounds[i + 1])
    if local is None:
        local = np.arange(e - s)
    return local, codes[s + local]


def _dask_block_sums(block, local: np.ndarray, groups: np.ndarray, n_groups: int) -> np.ndarray:
    return _group_sums(_kernel_ready(block), local, groups, n_groups)


def _group_sums(X, rows: np.ndarray, groups: np.ndarray, n_groups: int) -> np.ndarray:
    """Per group (float64): the sum of the rows `rows` of the kernel-ready matrix `X` whose group is in `groups`."""
    sums = np.zeros((n_kernel_chunks(len(rows), per_thread=1), n_groups, X.shape[1]))
    if sp.issparse(X):
        _group_sums_csr(X.data, X.indices, X.indptr, rows, groups, sums)
    else:
        _group_sums_dense(X, rows, groups, sums)
    return sums.sum(axis=0)
