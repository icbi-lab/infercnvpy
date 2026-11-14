"""The numba kernels of `tl.infercnv` against a plain numpy implementation of the algorithm."""

import threading

import numba
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import infercnvpy as cnv
from infercnvpy.tl._infercnv import _natural_sort
from infercnvpy.tl._infercnv_kernels import _median


def _naive_infercnv(adata, reference, *, lfc_clip, window_size, step, dynamic_threshold, noise_rows, mode):
    """infercnv with np.convolve, as infercnvpy <= 0.6 computed it (but with one noise band for all cells)."""
    var = adata.var
    mask = (var["chromosome"].isnull() | var["chromosome"].isin(["chrX", "chrY"])).to_numpy()
    used = var.loc[~mask]
    x = np.asarray(adata.X.toarray() if sp.issparse(adata.X) else adata.X, dtype=np.float64)
    if reference.shape[0] == 1:
        lfc = x - reference[0]
    else:
        lo, hi = reference.min(axis=0), reference.max(axis=0)
        lfc = np.where(x > hi, x - hi, np.where(x < lo, x - lo, 0))
    lfc = np.clip(lfc, -lfc_clip, lfc_clip)

    windows, gene_parts, gene_cols = [], [], []
    for chromosome in _natural_sort([c for c in used["chromosome"].unique() if c.startswith("chr") and c != "chrM"]):
        genes = used.loc[used["chromosome"] == chromosome].sort_values("start").index
        xs = lfc[:, var.index.get_indexer(genes)]
        g = len(genes)
        n = window_size if window_size < g else g
        r = np.arange(1, n + 1)
        pyramid = np.minimum(r, r[::-1]) if window_size < g else np.ones(n)
        smoothed = np.array([np.convolve(row, pyramid, mode="valid") for row in xs])[:, ::step] / pyramid.sum()
        member = np.zeros((smoothed.shape[1], g))
        for k in range(smoothed.shape[1]):
            member[k, k * step : k * step + n] = 1
        covered = member.sum(axis=0) > 0
        windows.append(smoothed)
        gene_parts.append((smoothed @ member[:, covered]) / member.sum(axis=0)[covered])
        gene_cols.append(var.index.get_indexer(genes[covered]))
    res = np.hstack(windows)
    res -= np.median(res, axis=1, keepdims=True)
    gene = np.hstack(gene_parts)
    gene -= np.median(gene, axis=1, keepdims=True)

    centre, thr = 0.0, -1.0
    if dynamic_threshold is not None:
        cells = res if noise_rows is None else res[noise_rows]
        if mode == "all":
            thr = dynamic_threshold * res.std()
        else:
            centre, thr = cells.mean(), dynamic_threshold * cells.std(axis=1, ddof=1).mean()
    res[np.abs(res - centre) < thr] = 0
    gene[np.abs(gene - centre) < thr] = 0
    gene_values = np.full(x.shape, np.nan)
    gene_values[:, np.concatenate(gene_cols)] = gene
    return res, gene_values


@pytest.fixture(scope="module")
def random_adata():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 40, 400
    chromosome = rng.choice(
        ["chr1", "chr2", "chr10", "chr3", "chrM", "chrX", "1", None],
        n_genes,
        p=[0.35, 0.25, 0.2, 0.03, 0.02, 0.05, 0.02, 0.08],
    )
    var = pd.DataFrame(
        {"chromosome": chromosome, "start": rng.integers(0, 400, n_genes)},  # many tied starts
        index=[f"g{i}" for i in range(n_genes)],
    )
    var["end"] = var["start"] + 10
    X = (rng.poisson(0.5, (n_cells, n_genes)) * rng.random((n_cells, n_genes))).astype(np.float32)
    cell_type = rng.choice(["normal", "immune", "tumor"], n_cells)
    obs = pd.DataFrame({"cell_type": cell_type}, index=[f"c{i}" for i in range(n_cells)])
    return AnnData(X=X, obs=obs, var=var)


FORMATS = {
    "dense": lambda X: X,
    "dense64": lambda X: X.astype(np.float64),
    "dense_int": lambda X: np.round(X * 3).astype(np.int64),
    "fortran": np.asfortranarray,
    "csr": sp.csr_matrix,
    "csc": sp.csc_matrix,
    "csr_array": sp.csr_array,
    "csr64": lambda X: sp.csr_matrix(X.astype(np.float64)),
}


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize(("window_size", "step"), [(30, 5), (7, 1), (1000, 10)])
@pytest.mark.parametrize("n_ref", [1, 2])
@pytest.mark.parametrize("noise_reference", ["reference", "all"])
@pytest.mark.parametrize("sparse_cost", [0.0, 1e9], ids=["prefix_sums", "sparse"])
def test_matches_naive(random_adata, fmt, window_size, step, n_ref, noise_reference, sparse_cost, monkeypatch):
    """Both ways of smoothing a cell (see `_infercnv_kernels`) match the naive implementation."""
    monkeypatch.setattr(cnv.tl._infercnv, "_SPARSE_COST", sparse_cost)
    adata = random_adata.copy()
    adata.X = FORMATS[fmt](adata.X)
    cats = ["normal", "immune"][:n_ref]
    reference = np.vstack([np.asarray(adata.X[(adata.obs["cell_type"] == c).to_numpy()].mean(axis=0)) for c in cats])
    noise_rows = np.flatnonzero(adata.obs["cell_type"].isin(cats))
    kw = {"lfc_clip": 1.0, "window_size": window_size, "step": step, "dynamic_threshold": 1.0}
    expected, expected_gv = _naive_infercnv(
        adata, reference.astype(np.float32), noise_rows=noise_rows, mode=noise_reference, **kw
    )
    _, res, gv = cnv.tl.infercnv(
        adata,
        reference_key="cell_type",
        reference_cat=cats,
        noise_reference=noise_reference,
        calculate_gene_values=True,
        inplace=False,
        **kw,
    )
    assert isinstance(res, sp.csr_matrix)
    assert res.dtype == gv.dtype == np.float32
    np.testing.assert_array_equal(res.toarray() != 0, expected != 0)
    np.testing.assert_allclose(res.toarray(), expected, atol=2e-7)
    np.testing.assert_array_equal(np.isnan(gv), np.isnan(expected_gv))
    np.testing.assert_allclose(gv, expected_gv, atol=5e-7)


def test_running_mean_values():
    """Smoothing of one chromosome: a pyramid of 3 genes; a chromosome of fewer genes than the window is its mean."""
    X = np.array([[1, 2, 3, 4, 5, 1, 2, 3], [6, 7, 8, 9, 10, 3, 5, 7]], dtype=np.float32)
    var = pd.DataFrame(
        {"chromosome": ["chr1"] * 5 + ["chr2"] * 3, "start": [5, 1, 2, 3, 4, 0, 1, 2], "end": 0},
        index=[f"g{i}" for i in range(8)],
    )
    adata = AnnData(X=X, var=var)
    kw = {"reference": np.zeros(8), "lfc_clip": 100, "step": 1, "dynamic_threshold": None, "inplace": False}
    chr_pos, res, gv = cnv.tl.infercnv(adata, window_size=3, calculate_gene_values=True, **kw)
    assert chr_pos == {"chr1": 0, "chr2": 3}
    # chr1 in genomic order is 2, 3, 4, 5, 1: windows 3, 4, 15/4 (pyramid 1, 2, 1); chr2 is a single window
    windows = np.array([[3, 4, 15 / 4, 2], [8, 9, 35 / 4, 5]])
    np.testing.assert_allclose(res.toarray(), windows - np.median(windows, axis=1, keepdims=True), rtol=1e-6)
    # a gene is the mean of the windows that contain it
    genes = np.array([[15 / 4, 3, 7 / 2, 43 / 12, 31 / 8, 2, 2, 2], [35 / 4, 8, 17 / 2, 103 / 12, 71 / 8, 5, 5, 5]])
    np.testing.assert_allclose(gv, genes - np.median(genes, axis=1, keepdims=True), rtol=1e-6)


@pytest.mark.parametrize("n", [1, 2, 3, 10, 11, 100, 101])
def test_median(n):
    rng = np.random.default_rng(n)
    for values in (rng.normal(size=n), rng.integers(0, 3, n).astype(np.float64), np.zeros(n)):
        assert _median(values, np.empty(n)) == np.median(values)


@pytest.fixture(scope="module")
def oligo():
    return cnv.datasets.oligodendroglioma()


REF = {"reference_key": "cell_type", "reference_cat": ["Microglia/Macrophage", "Oligodendrocytes (non-malignant)"]}


def _run(adata, **kw):
    return cnv.tl.infercnv(adata, **REF, inplace=False, **kw)


def _assert_identical(a, b):
    assert a[0] == b[0]
    assert (a[1] != b[1]).nnz == 0
    if a[2] is not None:
        np.testing.assert_array_equal(a[2], b[2])


@pytest.mark.parametrize("noise_reference", ["reference", "all"])
def test_chunksize_does_not_change_the_result(oligo, noise_reference):
    results = [
        _run(oligo, chunksize=chunksize, noise_reference=noise_reference, calculate_gene_values=True)
        for chunksize in (7, 50, 10_000)
    ]
    _assert_identical(results[0], results[2])
    _assert_identical(results[1], results[2])


def test_noise_reference(oligo):
    """The band comes from the reference cells; without reference cells from all cells; or from the std of all."""
    ref = _run(oligo)[1]
    all_cells = _run(oligo, noise_reference="all")[1]
    assert ref.nnz > all_cells.nnz  # the reference cells are less noisy than the tumor cells
    reference = cnv.tl._infercnv._get_reference(oligo, REF["reference_key"], REF["reference_cat"], None, None)
    # without reference cells, both the band of the reference cells and the global std use all cells
    explicit = cnv.tl.infercnv(oligo, reference=reference, inplace=False)[1]
    unthresholded = _run(oligo, dynamic_threshold=None)[1].toarray().astype(np.float64)
    centre, thr = unthresholded.mean(), 1.5 * unthresholded.std(axis=1, ddof=1).mean()
    expected = np.where(np.abs(unthresholded - centre) < thr, 0, unthresholded)
    np.testing.assert_array_equal(explicit.toarray(), expected)
    assert (explicit != ref).nnz > 0


def test_n_jobs(oligo):
    before = numba.get_num_threads()
    results = [_run(oligo, n_jobs=n_jobs, chunksize=50) for n_jobs in (1, 2, -1, 10_000, None)]
    assert numba.get_num_threads() == before
    for res in results[1:]:
        _assert_identical(res, results[0])
    for n_jobs in (0, -2):
        with pytest.raises(ValueError, match="n_jobs"):
            _run(oligo, n_jobs=n_jobs)


def test_runs_serially_outside_the_main_thread(oligo):
    """In another thread, the serial kernels run (a parallel kernel would abort with numba's workqueue layer)."""
    results = []
    thread = threading.Thread(target=lambda: results.append(_run(oligo, chunksize=50, calculate_gene_values=True)))
    thread.start()
    thread.join()
    _assert_identical(results[0], _run(oligo, chunksize=50, calculate_gene_values=True))


def test_records_noise_reference(oligo):
    adata = oligo.copy()
    cnv.tl.infercnv(adata, **REF, noise_reference="all")
    assert adata.uns["cnv"]["params"]["noise_reference"] == "all"
    with pytest.raises(ValueError, match="noise_reference"):
        cnv.tl.infercnv(adata, **REF, noise_reference="chunk")
    with pytest.raises(ValueError, match="dask array"):
        cnv.tl.infercnv(adata, **REF, lazy=True)


@pytest.mark.parametrize("kw", [{"window_size": 0}, {"step": 0}, {"chunksize": 0}])
def test_invalid_sizes(oligo, kw):
    with pytest.raises(ValueError, match="at least 1"):
        _run(oligo, **kw)


def test_no_cells(oligo):
    with pytest.raises(ValueError, match="at least one cell"):
        _run(oligo[:0])


def test_duplicate_sparse_entries(random_adata):
    """Duplicate entries of a CSR matrix are summed, as by `toarray()`; the input is left as it is."""
    X = sp.csr_matrix(random_adata.X)
    data, indices, indptr = [], [], [0]
    for r in range(X.shape[0]):  # every entry twice, with half the value, in reverse order first
        row = slice(X.indptr[r], X.indptr[r + 1])
        data += [X.data[row][::-1] / 2, X.data[row] / 2]
        indices += [X.indices[row][::-1], X.indices[row]]
        indptr.append(indptr[-1] + 2 * (row.stop - row.start))
    dup = sp.csr_matrix((np.concatenate(data), np.concatenate(indices), indptr), shape=X.shape)
    assert not dup.has_canonical_format
    before = dup.data.copy()
    adata = random_adata.copy()
    kw = {"reference": np.zeros(adata.n_vars), "window_size": 7, "inplace": False}
    want = cnv.tl.infercnv(adata, **kw)[1]
    adata.X = dup
    np.testing.assert_allclose(cnv.tl.infercnv(adata, **kw)[1].toarray(), want.toarray(), atol=1e-6)
    np.testing.assert_array_equal(dup.data, before)


@pytest.mark.parametrize("preallocated", [True, False])
def test_csr_builder(preallocated):
    """Blocks of rows go into one CSR matrix, also when its arrays have to grow (the dense allocation failed)."""
    from infercnvpy.tl._infercnv import _CsrBuilder

    rng = np.random.default_rng(0)
    windows = (rng.normal(size=(50, 30)) * (rng.random((50, 30)) < 0.3)).astype(np.float32)
    band = (0.5, 0.1)
    keep = (windows != 0) & (np.abs(windows.astype(np.float64) - band[1]) >= band[0])
    builder = _CsrBuilder(*windows.shape)
    if not preallocated:
        builder.data, builder.indices = builder.data[:0].copy(), builder.indices[:0].copy()
    for rows in (slice(0, 3), slice(3, 40), slice(40, 50)):
        builder.add(windows[rows], keep[rows].sum(axis=1), band)
    res = builder.finish()
    np.testing.assert_array_equal(res.toarray(), np.where(keep, windows, 0))
    assert res.indptr.dtype == res.indices.dtype == np.int32


@pytest.fixture(scope="module")
def oligo_log(oligo):
    """The oligodendroglioma counts, log-normalised as infercnv expects (and as `min_mean` assumes)."""
    import scanpy as sc

    adata = oligo.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def _gene_means(adata):
    """Mean of expm1 of the expression over all cells, per gene (the scale of `min_mean`)."""
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    return np.expm1(X.astype(np.float64)).mean(axis=0)


def test_min_mean_default_keeps_all_genes(oligo_log):
    """`min_mean=None` (default) changes nothing, and a threshold that removes no gene gives the same result."""
    default = _run(oligo_log, calculate_gene_values=True)
    _assert_identical(_run(oligo_log, min_mean=None, calculate_gene_values=True), default)
    _assert_identical(_run(oligo_log, min_mean=0, calculate_gene_values=True), default)


@pytest.mark.parametrize("fmt", ["csr", "dense", "dask", "backed"])
def test_min_mean_removes_lowly_expressed_genes(oligo_log, fmt, tmp_path, caplog):
    """The genes below `min_mean` are left out before windowing: the same as running on the other genes only."""
    import logging

    oligo = oligo_log
    means = np.unique(_gene_means(oligo))
    i = len(means) * 3 // 10
    threshold = float((means[i] + means[i + 1]) / 2)  # halfway between two genes, so rounding cannot move a gene
    keep = _gene_means(oligo) >= threshold
    expected = _run(oligo[:, keep].copy(), calculate_gene_values=True)

    adata = oligo.copy()
    if fmt == "dense":
        adata.X = adata.X.toarray()
    elif fmt == "dask":
        import dask.array as da

        adata.X = da.from_array(adata.X.toarray(), chunks=(50, -1))
    elif fmt == "backed":
        import anndata as ad

        adata.write_h5ad(tmp_path / "a.h5ad")
        adata = ad.read_h5ad(tmp_path / "a.h5ad", backed="r")
    with caplog.at_level(logging.INFO, logger="scverse.infercnvpy"):
        chr_pos, res, gv = cnv.tl.infercnv(
            adata, **REF, inplace=False, calculate_gene_values=True, min_mean=threshold, chunksize=60
        )
    assert chr_pos == expected[0]
    if fmt in ("csr", "dense"):
        assert (res != expected[1]).nnz == 0
    else:  # the reference means of lazy data are sums over blocks
        np.testing.assert_allclose(res.toarray(), expected[1].toarray(), rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(gv[:, keep], expected[2], rtol=1e-5, atol=1e-6)
    assert np.isnan(gv[:, ~keep]).all()
    excluded = oligo.var["chromosome"].isna() | oligo.var["chromosome"].isin(["chrX", "chrY"])
    n_low = int((~keep & ~excluded.to_numpy()).sum())
    [msg] = [r.getMessage() for r in caplog.records if r.name == "scverse.infercnvpy"]
    assert f"left out {n_low} genes with a mean expression below min_mean={threshold:g}" in msg


def test_min_mean_recorded_and_all_genes_removed(oligo_log):
    adata = oligo_log.copy()
    cnv.tl.infercnv(adata, **REF, min_mean=0.5)
    assert adata.uns["cnv"]["params"]["min_mean"] == 0.5
    with pytest.raises(ValueError, match="No genes left"):
        _run(oligo_log, min_mean=1e30)
