import os
import re
import threading
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse
from anndata import AnnData
from scanpy import logging
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import thread_map

from infercnvpy._util import _fetch_rows, _is_lazy, _read_units


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
    noise_sd: Literal["global", "chunk"] = "global",
    noise_sample_size: int = 50_000,
    exclude_chromosomes: Sequence[str] | None = ("chrX", "chrY"),
    chunksize: int = 5000,
    n_jobs: int | None = None,
    inplace: bool = True,
    layer: str | None = None,
    key_added: str = "cnv",
    calculate_gene_values: bool = False,
) -> None | tuple[dict, scipy.sparse.csr_matrix, np.ndarray | None]:
    """Infer Copy Number Variation (CNV) by averaging gene expression over genomic regions.

    This method is heavily inspired by `infercnv <https://github.com/broadinstitute/inferCNV/>`_
    but more computationally efficient. The method is described in more detail
    in on the :ref:`infercnv-method` page.

    There, you can also find instructions on how to :ref:`prepare input data <input-data>`.

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
        Values `< dynamic threshold * STDDEV` will be set to 0, where STDDEV is
        the stadard deviation of the smoothed gene expression. Set to `None` to disable
        this step.
    noise_sd
        How the standard deviation used for `dynamic_threshold` is determined.
        `"global"` (default) estimates it once from cells spread evenly over the dataset (see `noise_sample_size`),
        so the result does not depend on `chunksize` and is exact for datasets with at most `noise_sample_size`
        cells. `"chunk"` computes it separately for every chunk of cells (the behaviour of
        infercnvpy <= 0.6, where the result depends on `chunksize`).
    noise_sample_size
        Number of cells used to estimate the noise standard deviation if `noise_sd="global"`. Only whole read units
        (the blocks of a dask array, or `chunksize` cells otherwise) are used, because a block is read as a whole
        anyway, and at least 5 of them: for a dask array with 50,000-cell blocks that is 250,000 cells. Their
        results are kept in memory (up to 2 GB), so they are not computed twice.
    exclude_chromosomes
        List of chromosomes to exclude. The default is to exclude genosomes.
    chunksize
        Process dataset in chunks of cells. This allows to run infercnv on
        datasets with many cells, where the dense matrix would not fit into memory.
        `adata.X` (or the layer) may be a dask array or a backed matrix: it is then read chunk by chunk and
        never loaded entirely. Peak memory is about `n_jobs * chunksize * n_genes` values for dense data,
        and much less for sparse data.
    n_jobs
        Number of jobs for parallel processing. Default: use all cores.
        Data will be submitted to workers in chunks, see `chunksize`.
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
        As many genes will be included in each segment the resultant per gene value will be an average of the genes included in the segment.
        Additionally not all genes will be included in the per gene CNV, due to the window size and step size not always being a multiple of
        the number of genes. Any genes not included in the per gene CNV will be filled with NaN.
        Note this will significantly increase the memory and computation time, it is recommended to decrease the chunksize to ~100 if this is set to True.


    Returns
    -------
    Depending on inplace, either return the smoothed and denoised gene expression
    matrix sorted by genomic position, or add it to adata.
    """
    if not adata.var_names.is_unique:
        raise ValueError("Ensure your var_names are unique!")
    if {"chromosome", "start", "end"} - set(adata.var.columns) != set():
        raise ValueError(
            "Genomic positions not found. There need to be `chromosome`, `start`, and `end` columns in `adata.var`. "
        )

    # `var` can be lazy (`anndata.experimental.read_lazy`): load just the three columns that are needed
    pos = adata.var[["chromosome", "start", "end"]]
    pos = pos.to_memory() if hasattr(pos, "to_memory") else pos
    pos.index = adata.var_names

    var_mask = pos["chromosome"].isnull()
    if np.sum(var_mask):
        logging.warning(f"Skipped {np.sum(var_mask)} genes because they don't have a genomic position annotated. ")  # type: ignore
    if exclude_chromosomes is not None:
        var_mask = var_mask | pos["chromosome"].isin(exclude_chromosomes)

    keep = np.flatnonzero(~np.asarray(var_mask))
    expr = adata.X if layer is None else adata.layers[layer]
    if scipy.sparse.issparse(expr):
        expr = expr.tocsr()
    # Lazy matrices are read in large units (whole dask blocks) of which several are in memory at the same time:
    # use few workers by default.
    max_workers = (min(_default_n_jobs(), 4) if _is_lazy(expr) else _default_n_jobs()) if n_jobs is None else n_jobs

    reference = _get_reference(
        adata, reference_key, reference_cat, reference, layer, chunksize=chunksize, max_workers=max_workers
    )[:, keep]

    var = pos.iloc[keep]
    layout = _build_layout(var, window_size, step, with_gene_map=calculate_gene_values)
    # Genes of `adata` in the order the kernel needs them: by chromosome, then genomic position.
    # Genes that are not used are dropped here, so the full matrix is never copied or densified.
    cols = keep[layout.order]
    reference = np.asarray(reference)[:, layout.order]
    n_obs = adata.shape[0]
    units = _read_units(expr, chunksize)

    def _load(start: int, stop: int):
        # the only read of the rows of this unit (a whole dask block, if lazy)
        return _fetch_rows(expr, start, min(stop, n_obs))

    def _pieces(raw):
        """The rows of a loaded unit, restricted to the used genes, as pieces of at most `chunksize` rows."""
        for a in range(0, raw.shape[0], chunksize):
            yield raw[a : a + chunksize][:, cols]

    def _centered(unit):
        return [
            _chunk_centered(piece, reference, layout, lfc_clip, calculate_gene_values=calculate_gene_values)
            for piece in _pieces(_load(*unit))
        ]

    cache: dict = {}
    if dynamic_threshold is not None and noise_sd == "global":
        noise_thres = dynamic_threshold * _estimate_noise_sd(
            _load,
            _pieces,
            units,
            noise_sample_size,
            reference,
            layout,
            lfc_clip,
            max_workers,
            cache,
            use_cache=not calculate_gene_values,
        )

        def _process(unit):
            parts = cache.pop(unit[0], None) or _centered(unit)
            return [_apply_threshold(x_res, gene_res, noise_thres) for x_res, gene_res in parts]

    else:

        def _process(unit):
            return [
                _apply_threshold(x_res, gene_res, dynamic_threshold * np.std(x_res) if dynamic_threshold else None)
                for x_res, gene_res in _centered(unit)
            ]

    unit_results = thread_map(_process, units, tqdm_class=tqdm, max_workers=max_workers)
    chunks = [c for parts in unit_results for c, _ in parts]
    gene_chunks = [g for parts in unit_results for _, g in parts]

    res = scipy.sparse.vstack(chunks)
    chr_pos = layout.chr_pos

    if calculate_gene_values:
        # Genes that are not part of any window (or were excluded) are filled with NaN.
        per_gene_mtx = np.full((adata.shape[0], adata.shape[1]), np.nan, dtype=gene_chunks[0].dtype)
        per_gene_mtx[:, cols[layout.covered_genes]] = np.vstack(gene_chunks)
    else:
        per_gene_mtx = None

    if inplace:
        adata.obsm[f"X_{key_added}"] = res
        adata.uns[key_added] = {"chr_pos": chr_pos}

        if calculate_gene_values:
            adata.layers[f"gene_values_{key_added}"] = per_gene_mtx

    else:
        return chr_pos, res, per_gene_mtx


def _default_n_jobs() -> int:
    """Number of CPUs this process may use (respects affinity masks, e.g. `taskset`/cgroups)."""
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)


def _natural_sort(l: Sequence):
    """Natural sort without third party libraries.

    Adapted from: https://stackoverflow.com/a/4836734/2340703
    """

    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


@dataclass
class _ChromBlock:
    """Smoothing of one chromosome, expressed on the genes sorted by genomic position.

    `W` (`n_genes x n_windows`) holds the pyramid weights, so `x[:, g0:g1] @ W` are the (strided) running means.
    `G` (`n_windows x n_covered`) averages, for each covered gene, all windows it is part of.
    """

    g0: int
    g1: int
    w0: int
    W: np.ndarray
    cov: np.ndarray  # positions (within the chromosome) of genes in >=1 window
    G: np.ndarray | None


@dataclass
class _Layout:
    """Precomputed, chunk-independent description of the genomic smoothing.

    Attributes
    ----------
    chr_pos
        Maps each chromosome to the column of the smoothed matrix where it starts.
    order
        Positions of the genes (in the `var` the layout was built from) sorted by chromosome and genomic start.
        Chunks have to be restricted to these columns, in this order, before they are passed to the kernel.
    blocks
        One `_ChromBlock` per chromosome.
    n_windows
        Total number of columns of the smoothed matrix.
    W_sparse
        All chromosome blocks as one sparse block-diagonal `n_genes x n_windows` matrix. Each gene is only part of
        about `window_size / step` windows, so this is much cheaper than the dense blocks for sparse input.
    covered_genes
        Positions in `order` of the genes that are part of at least one window (the columns of the per-gene values).
    """

    chr_pos: dict
    order: np.ndarray
    blocks: list
    n_windows: int
    W_sparse: scipy.sparse.csr_matrix
    covered_genes: np.ndarray


def _chromosome_weights(m: int, window_size: int, step: int) -> np.ndarray:
    """Dense `m x n_windows` weight matrix of the running means over `m` position-sorted genes.

    Windows use a pyramidal weighting and only every `step`-th window is kept. If the chromosome has no more genes
    than `window_size`, a single window averaging all genes is used.
    """
    if window_size < m:
        r = np.arange(1, window_size + 1)
        pyramid = np.minimum(r, r[::-1]).astype(np.float64)
        starts = np.arange(0, m - window_size + 1, step)
        W = np.zeros((m, len(starts)))
        idx = starts[None, :] + np.arange(window_size)[:, None]
        W[idx, np.arange(len(starts))[None, :]] = (pyramid / pyramid.sum())[:, None]
    else:
        W = np.full((m, 1), 1.0 / m)
    return W


def _build_layout(var: pd.DataFrame, window_size: int, step: int, with_gene_map: bool = False) -> _Layout:
    """Build the smoothing layout for the genes in `var` (columns `chromosome` and `start`)."""
    chrom = np.asarray(var["chromosome"])
    start = np.asarray(var["start"])
    chromosomes = _natural_sort([x for x in pd.unique(chrom) if x.startswith("chr") and x != "chrM"])

    order, blocks, chr_pos, covered = [], [], {}, []
    g0 = w0 = 0
    for chr_ in chromosomes:
        pos = np.flatnonzero(chrom == chr_)
        order.append(pos[np.argsort(start[pos], kind="stable")])
        W = _chromosome_weights(len(pos), window_size, step)
        counts = (W > 0).sum(axis=1)
        cov = np.flatnonzero(counts > 0)
        G = ((W[cov] > 0) / counts[cov, None]).T if with_gene_map else None
        blocks.append(_ChromBlock(g0=g0, g1=g0 + len(pos), w0=w0, W=W, cov=cov, G=G))
        covered.append(g0 + cov)
        chr_pos[chr_] = w0
        g0 += len(pos)
        w0 += W.shape[1]

    empty = np.zeros(0, dtype=int)
    return _Layout(
        chr_pos=chr_pos,
        order=np.concatenate(order) if order else empty,
        blocks=blocks,
        n_windows=w0,
        W_sparse=scipy.sparse.block_diag([scipy.sparse.csr_matrix(b.W) for b in blocks], format="csr")
        if blocks
        else scipy.sparse.csr_matrix((0, 0)),
        covered_genes=np.concatenate(covered) if covered else empty,
    )


def _dense_to_csr(a: np.ndarray) -> scipy.sparse.csr_matrix:
    """Convert a dense (row-major) array to CSR; much faster than `csr_matrix(a)` for large arrays."""
    mask = a != 0
    indptr = np.zeros(a.shape[0] + 1, dtype=np.int64)
    np.cumsum(mask.sum(axis=1), out=indptr[1:])
    indices = np.nonzero(mask)[1].astype(np.int32 if a.shape[1] < 2**31 else np.int64)
    return scipy.sparse.csr_matrix((a[mask], indices, indptr), shape=a.shape)


def _running_mean(
    x: np.ndarray | scipy.sparse.spmatrix,
    n: int = 50,
    step: int = 10,
    gene_list: np.ndarray | list | None = None,
    calculate_gene_values: bool = False,
) -> tuple[np.ndarray, pd.DataFrame | None]:
    """Pyramid-weighted running mean over the columns of `x` (one chromosome, genes in positional order).

    Kept with its original signature for backwards compatibility and for tests of the smoothing itself;
    :func:`infercnv` uses the same weight matrices (:func:`_build_layout`) directly.

    Returns the smoothed matrix (every `step`-th window) and, if `calculate_gene_values`, the per-gene values
    (mean of all windows that contain a gene) as a data frame.
    """
    x = x.toarray() if scipy.sparse.issparse(x) else np.asarray(x)
    genes = np.asarray(gene_list) if gene_list is not None else np.arange(x.shape[1]).astype(str)
    var = pd.DataFrame({"chromosome": "chr1", "start": np.arange(len(genes))}, index=genes)
    layout = _build_layout(var, n, step, with_gene_map=calculate_gene_values)
    (block,) = layout.blocks  # a single chromosome
    smoothed = x @ block.W
    if not calculate_gene_values:
        return smoothed, None
    return smoothed, pd.DataFrame(smoothed @ block.G, columns=var.index[layout.order[layout.covered_genes]])


def _calculate_gene_averages(convolved_gene_names: np.ndarray, smoothed_x: np.ndarray) -> pd.DataFrame:
    """Per-gene value: the mean of all windows that contain the gene.

    `convolved_gene_names` has one row per window (the names of the genes in it), `smoothed_x` one column per window.
    Kept with its original signature for backwards compatibility; :func:`infercnv` uses a precomputed matrix instead.
    """
    names = np.asarray(convolved_gene_names)
    window_values = np.repeat(np.asarray(smoothed_x, dtype=float), names.shape[1], axis=1)  # one value per gene slot
    averaged = pd.DataFrame(window_values.T).groupby(pd.Index(names.ravel()), sort=False).mean().T
    return averaged.reset_index(drop=True)


# Sparse input with fewer non-zeros than this fraction is smoothed without densifying (see `_smooth_sparse`).
_SPARSE_PATH_MAX_DENSITY = 0.1


def _log_fold_change(x, reference, lfc_cap):
    """Log fold change vs. the reference, clipped at `lfc_cap`. Elementwise, so also valid on sparse `.data`.

    `reference` is `n_ref x n_genes` and has to be broadcastable against `x`. With multiple references the
    "bounded" difference is used: entries between the lowest and highest reference are considered having a
    logFC of 0. Note that `x - clip(x, lo, hi)` equals `max(x - hi, 0) + min(x - lo, 0)`.
    """
    if reference.shape[0] == 1:
        lfc = x - reference[0]
    else:
        lfc = x - np.clip(x, reference.min(axis=0), reference.max(axis=0))
    return np.clip(lfc, -lfc_cap, lfc_cap, out=lfc)


def _smooth_dense(tmp_x, reference, layout: _Layout, lfc_cap, dtype) -> np.ndarray:
    """Densify the chunk, compute clipped logFC and smooth with one BLAS matmul per chromosome."""
    x = tmp_x.toarray() if scipy.sparse.issparse(tmp_x) else np.array(tmp_x)  # `x` is a private copy either way
    x = x.astype(dtype, copy=False)
    if reference.shape[0] == 1:
        x -= reference[0]
    else:
        x -= np.clip(x, reference.min(axis=0), reference.max(axis=0))
    np.clip(x, -lfc_cap, lfc_cap, out=x)
    x_smoothed = np.empty((x.shape[0], layout.n_windows), dtype=dtype)
    for b in layout.blocks:
        np.matmul(x[:, b.g0 : b.g1], b.W.astype(dtype, copy=False), out=x_smoothed[:, b.w0 : b.w0 + b.W.shape[1]])
    return x_smoothed


def _smooth_sparse(tmp_x: scipy.sparse.csr_matrix, reference, layout: _Layout, lfc_cap, dtype) -> np.ndarray:
    """Same as `_smooth_dense`, but without ever densifying the chunk.

    Every zero in the chunk takes the same clipped logFC `c[gene]` (what the reference makes of an expression of 0).
    As smoothing is linear, `smooth(logFC) = c @ W + smooth(delta)`, where `delta = logFC - c` is non-zero only
    where the chunk has non-zeros. So the baseline is smoothed once and only the stored entries are touched.
    """
    base = _log_fold_change(np.zeros(1, dtype=dtype), reference, lfc_cap)  # `c`, shape (n_genes,)
    idx = tmp_x.indices
    if reference.shape[0] == 1:
        lfc = tmp_x.data.astype(dtype, copy=False) - reference[0][idx]
    else:
        lfc = tmp_x.data.astype(dtype, copy=False)
        lfc = lfc - np.clip(lfc, reference.min(axis=0)[idx], reference.max(axis=0)[idx])
    np.clip(lfc, -lfc_cap, lfc_cap, out=lfc)
    lfc -= base[idx]
    delta = scipy.sparse.csr_matrix((lfc, idx, tmp_x.indptr), shape=tmp_x.shape)
    W = layout.W_sparse.astype(dtype, copy=False)
    x_smoothed = (delta @ W).toarray()
    x_smoothed += W.T @ base
    return x_smoothed


def _chunk_centered(
    tmp_x, reference, layout: _Layout, lfc_cap, calculate_gene_values=False
) -> tuple[np.ndarray, np.ndarray | None]:
    """Steps 1-4 of infercnv for one chunk of cells: logFC, clipping, smoothing and per-cell centering.

    `tmp_x` and `reference` must be restricted to the genes in `layout.order`, in that order.
    Returns the dense `n_cells x n_windows` matrix and, if requested, the per-gene values.
    """
    dtype = tmp_x.dtype if tmp_x.dtype.kind == "f" else np.dtype(np.float64)  # keep precision; promote integers
    reference = np.asarray(reference).astype(dtype, copy=False)  # sparse means come back as np.matrix
    if scipy.sparse.issparse(tmp_x) and tmp_x.nnz < _SPARSE_PATH_MAX_DENSITY * tmp_x.shape[0] * tmp_x.shape[1]:
        x_smoothed = _smooth_sparse(tmp_x.tocsr(), reference, layout, lfc_cap, dtype)
    else:
        x_smoothed = _smooth_dense(tmp_x, reference, layout, lfc_cap, dtype)
    # Step 4 - center by cell
    x_res = x_smoothed - np.median(x_smoothed, axis=1)[:, np.newaxis]
    if calculate_gene_values:
        gene_res = np.hstack(
            [x_smoothed[:, b.w0 : b.w0 + b.W.shape[1]] @ b.G.astype(dtype, copy=False) for b in layout.blocks]
        )
        gene_res -= np.median(gene_res, axis=1)[:, np.newaxis]
    else:
        gene_res = None
    return x_res, gene_res


def _apply_threshold(x_res, gene_res, noise_thres) -> tuple[scipy.sparse.csr_matrix, np.ndarray | None]:
    """Step 5 - set values with an absolute value below `noise_thres` to 0 (in place) and make `x_res` sparse."""
    if noise_thres is not None:
        x_res[np.abs(x_res) < noise_thres] = 0
        if gene_res is not None:
            gene_res[np.abs(gene_res) < noise_thres] = 0
    return _dense_to_csr(x_res), gene_res


def _infercnv_chunk_kernel(
    tmp_x, reference, layout: _Layout, lfc_cap, dynamic_threshold, calculate_gene_values=False
) -> tuple[scipy.sparse.csr_matrix, np.ndarray | None]:
    """The complete infercnv work for one chunk, with the noise threshold determined from this chunk only."""
    x_res, gene_res = _chunk_centered(tmp_x, reference, layout, lfc_cap, calculate_gene_values)
    noise_thres = None if dynamic_threshold is None else dynamic_threshold * np.std(x_res)
    return _apply_threshold(x_res, gene_res, noise_thres)


# Keep the centered chunks computed while estimating the noise sd (so their rows are not read and computed a second
# time in the main pass), up to this many bytes.
_NOISE_CACHE_MAX_BYTES = 2e9


def _estimate_noise_sd(
    load, pieces, units, sample_size, reference, layout: _Layout, lfc_cap, max_workers, cache: dict, use_cache=True
) -> float:
    """Standard deviation of the centered, smoothed CNV values over (a sample of) all cells.

    `load(start, stop)` reads a unit of rows once and `pieces(raw)` splits it into chunks restricted to the used
    genes (see `infercnv`).

    With at most `sample_size` cells this is exactly the standard deviation over all cells. Otherwise, evenly spaced
    units (at least 5, if there are that many, and at least `sample_size` cells) are read and used as a whole: reading
    a unit loads all of its rows anyway, so using only some of them would waste the read.

    Unless `use_cache` is False (per-gene values are needed), the centered chunks of the units are stored in `cache`,
    keyed by the start of the unit, as long as they fit in `_NOISE_CACHE_MAX_BYTES`; the main pass then skips these
    units entirely.
    """
    n_obs = units[-1][1]
    if n_obs <= sample_size:
        selected = list(units)
    else:
        avg_rows = n_obs / len(units)
        k = min(len(units), max(5, int(np.ceil(sample_size / avg_rows))))
        selected = [units[i] for i in np.unique(np.linspace(0, len(units) - 1, k).round().astype(int))]

    lock = threading.Lock()
    cached_bytes = [0]

    def _stats(unit):
        start, stop = unit
        raw = load(start, stop)
        n = s1 = s2 = 0.0
        parts = [] if use_cache else None
        for piece in pieces(raw):
            x_res, _ = _chunk_centered(piece, reference, layout, lfc_cap)
            x64 = x_res.astype(np.float64, copy=False)
            n, s1, s2 = n + x64.size, s1 + x64.sum(), s2 + np.square(x64).sum()
            if parts is not None:
                parts.append((x_res, None))
        if parts is not None:
            nbytes = sum(x.nbytes for x, _ in parts)
            with lock:
                if cached_bytes[0] + nbytes <= _NOISE_CACHE_MAX_BYTES:
                    cached_bytes[0] += nbytes
                    cache[start] = parts
        return n, s1, s2

    n = s1 = s2 = 0.0
    for n_, a, b in thread_map(_stats, selected, tqdm_class=tqdm, max_workers=max_workers, desc="noise sd"):
        n, s1, s2 = n + n_, s1 + a, s2 + b
    return float(np.sqrt(max(s2 / n - (s1 / n) ** 2, 0.0)))


def _infercnv_chunk(tmp_x, var, reference, lfc_cap, window_size, step, dynamic_threshold, calculate_gene_values=False):
    """Run infercnv on a single chunk of cells. Convenience wrapper around `_infercnv_chunk_kernel`.

    Returns the chromosome start positions, the (sparse) CNV matrix and, if requested, a data frame with
    the per-gene values (columns are the genes that are part of at least one window).
    """
    layout = _build_layout(var, window_size, step, with_gene_map=calculate_gene_values)
    x_res, gene_res = _infercnv_chunk_kernel(
        tmp_x[:, layout.order],
        np.asarray(reference)[:, layout.order],
        layout,
        lfc_cap,
        dynamic_threshold,
        calculate_gene_values=calculate_gene_values,
    )
    if calculate_gene_values:
        gene_res = pd.DataFrame(gene_res, columns=var.index[layout.order[layout.covered_genes]].to_numpy())
    return layout.chr_pos, x_res, gene_res


def _get_reference(
    adata: AnnData,
    reference_key: str | None,
    reference_cat: None | str | Sequence[str],
    reference: np.ndarray | None,
    layer: str | None,
    chunksize: int = 5000,
    max_workers: int = 1,
) -> np.ndarray:
    """Parameter validation extraction of reference gene expression.

    If multiple reference categories are given, compute the mean per
    category. All means are computed in a single streaming pass over `X` (in chunks of `chunksize` cells), so
    this works for dask and backed matrices without loading them.

    Returns a 2D array with reference categories in rows, genes in columns.
    If there's just one category, it's still a 2D array.
    """
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if reference is None:
        if reference_key is None or reference_cat is None:
            logging.warning(
                "Using mean of all cells as reference. For better results, "
                "provide either `reference`, or both `reference_key` and `reference_cat`. "
            )  # type: ignore
            categories, obs_values = None, None
        else:
            obs_col = adata.obs[reference_key]
            obs_col = obs_col.to_pandas() if hasattr(obs_col, "to_pandas") else obs_col  # lazy obs
            if isinstance(reference_cat, str):
                reference_cat = [reference_cat]
            reference_cat = np.array(reference_cat)
            reference_cat_in_obs = np.isin(reference_cat, obs_col)
            if not np.all(reference_cat_in_obs):
                raise ValueError(
                    "The following reference categories were not found in "
                    "adata.obs[reference_key]: "
                    f"{reference_cat[~reference_cat_in_obs]}"
                )
            categories, obs_values = reference_cat, np.asarray(obs_col)

        reference = _streaming_means(X, categories, obs_values, chunksize, max_workers)

    if reference.ndim == 1:
        reference = reference[np.newaxis, :]

    if reference.shape[1] != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference


def _streaming_means(X, categories, obs_values, chunksize: int, max_workers: int = 1) -> np.ndarray:
    """Column means of `X` over all rows (`categories=None`), or per category of `obs_values`, in one pass.

    Per unit of rows (see `_read_units`), the sums of all categories are a single (sparse) indicator-matrix product,
    so no row subsetting of `X` is needed. Every row is read once.
    """
    n_groups = 1 if categories is None else len(categories)

    def _partial(unit):
        start, stop = unit
        chunk = _fetch_rows(X, start, stop)
        n = chunk.shape[0]
        # the indicator has the dtype of the data: scipy would otherwise upcast (copy) the whole chunk to float64
        dtype = chunk.dtype if chunk.dtype.kind == "f" else np.float64
        if categories is None:
            indicator = scipy.sparse.csr_matrix(np.ones((1, n), dtype=dtype))
        else:
            labels = obs_values[start : start + n]
            group, row = np.nonzero(categories[:, None] == labels[None, :])
            indicator = scipy.sparse.csr_matrix((np.ones(len(row), dtype=dtype), (group, row)), shape=(n_groups, n))
        part = indicator @ chunk
        part = part.toarray() if scipy.sparse.issparse(part) else np.asarray(part)
        return part.astype(np.float64, copy=False), np.asarray(indicator.sum(axis=1)).ravel().astype(np.int64)

    sums = np.zeros((n_groups, X.shape[1]), dtype=np.float64)
    counts = np.zeros(n_groups, dtype=np.int64)
    units = _read_units(X, chunksize)
    for part, count in thread_map(_partial, units, max_workers=max_workers, disable=len(units) < 3, desc="reference"):
        sums += part
        counts += count
    return sums / counts[:, None]
