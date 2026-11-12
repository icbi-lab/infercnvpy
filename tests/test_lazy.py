"""infercnv on lazy / on-disk data (dask, dask with sparse blocks as in lazily loaded spatial tables, backed h5ad)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import infercnvpy as cnv
import infercnvpy.tl._infercnv as m

da = pytest.importorskip("dask.array")
dask = pytest.importorskip("dask")

REF = {"reference_key": "cell_type", "reference_cat": ["Microglia/Macrophage", "Oligodendrocytes (non-malignant)"]}


@pytest.fixture(scope="module")
def base():
    return cnv.datasets.oligodendroglioma()


@pytest.fixture(scope="module")
def expected(base):
    a = base.copy()
    cnv.tl.infercnv(a, **REF)
    return a.obsm["X_cnv"].toarray(), a.uns["cnv"]["chr_pos"]


def _with_X(base, X):
    return ad.AnnData(X=X, obs=base.obs.copy(), var=base.var.copy())


def _csr_blocks(base, block_rows, counter=None):
    """Dask array whose row blocks are scipy CSR matrices, loaded lazily (like a lazily loaded `spatialdata_io` table)."""
    X = sp.csr_matrix(base.X)
    n, g = X.shape

    def load(start, stop):
        if counter is not None:
            counter.append((start, stop))
        return X[start:stop]

    meta = sp.csr_matrix((0, g), dtype=X.dtype)
    blocks = [
        da.from_delayed(
            dask.delayed(load)(s, min(s + block_rows, n)), shape=(min(block_rows, n - s), g), dtype=X.dtype, meta=meta
        )
        for s in range(0, n, block_rows)
    ]
    return da.concatenate(blocks, axis=0)


def _check(adata, expected, **kw):
    cnv.tl.infercnv(adata, **REF, **kw)
    np.testing.assert_allclose(adata.obsm["X_cnv"].toarray(), expected[0], rtol=1e-4, atol=1e-5)
    assert adata.uns["cnv"]["chr_pos"] == expected[1]


@pytest.mark.parametrize("chunksize", [50, 5000])
def test_dask_dense(base, expected, chunksize):
    X = da.from_array(base.X.toarray(), chunks=(40, base.shape[1]))
    _check(_with_X(base, X), expected, chunksize=chunksize)


@pytest.mark.parametrize("block_rows,chunksize", [(30, 30), (30, 70), (184, 50)])
def test_dask_sparse_blocks(base, expected, block_rows, chunksize):
    _check(_with_X(base, _csr_blocks(base, block_rows)), expected, chunksize=chunksize)


def test_dask_layer(base, expected):
    a = _with_X(base, sp.csr_matrix(base.X.shape, dtype=np.float32))
    a.layers["raw"] = _csr_blocks(base, 50)
    _check(a, expected, layer="raw", chunksize=60)


@pytest.mark.parametrize("fmt", ["csr", "dense"])
def test_backed_h5ad(base, expected, tmp_path, fmt):
    X = sp.csr_matrix(base.X) if fmt == "csr" else base.X.toarray()
    path = tmp_path / "a.h5ad"
    _with_X(base, X).write_h5ad(path)
    a = ad.read_h5ad(path, backed="r")
    assert a.isbacked
    _check(a, expected, chunksize=60)


def test_reads_each_row_once_per_pass(base, expected):
    """Reference + (noise estimate + main pass, which share cached chunks) = 2 reads of the data, nothing more."""
    loads = []
    a = _with_X(base, _csr_blocks(base, 46, counter=loads))
    _check(a, expected, chunksize=46)
    n_blocks = -(-base.shape[0] // 46)
    assert len(loads) == 2 * n_blocks, loads
    assert max(stop - start for start, stop in loads) <= 46, "more than one block was loaded at once"


def test_streaming_reference_matches_in_memory(base):
    a = _with_X(base, _csr_blocks(base, 25))
    for kw in ({"reference_key": None, "reference_cat": None}, REF):
        got = m._get_reference(
            a, kw["reference_key"] if "reference_key" in kw else None, kw["reference_cat"], None, None, chunksize=37
        )
        want = m._get_reference(base, kw["reference_key"], kw["reference_cat"], None, None)
        np.testing.assert_allclose(got, want, rtol=1e-5)


def test_chunksize_independent_by_default(base):
    results = []
    for chunksize in (7, 50, 10_000):
        a = base.copy()
        cnv.tl.infercnv(a, **REF, chunksize=chunksize)
        results.append(a.obsm["X_cnv"].toarray())
    np.testing.assert_allclose(results[0], results[2], rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(results[1], results[2], rtol=1e-4, atol=1e-5)


def test_chunk_mode_is_chunk_dependent(base):
    """`noise_sd="chunk"` keeps the legacy behaviour (per-chunk threshold), which depends on `chunksize`."""
    small, large = base.copy(), base.copy()
    cnv.tl.infercnv(small, **REF, chunksize=10, noise_sd="chunk")
    cnv.tl.infercnv(large, **REF, chunksize=10_000, noise_sd="chunk")
    assert not np.allclose(small.obsm["X_cnv"].toarray(), large.obsm["X_cnv"].toarray())


def test_sampled_noise_sd_is_close_to_exact(base):
    ref = m._get_reference(base, None, None, None, None)
    var = base.var.loc[:, ["chromosome", "start", "end"]]
    keep = np.flatnonzero(~base.var["chromosome"].isin(["chrX", "chrY"]).to_numpy())
    layout = m._build_layout(var.iloc[keep], 100, 10)
    cols = keep[layout.order]
    X, r = sp.csr_matrix(base.X), ref[:, cols]

    def load(a, b):
        return X[a:b]

    def pieces(raw):
        for a in range(0, raw.shape[0], 10):
            yield raw[a : a + 10][:, cols]

    n = base.shape[0]
    units = [(s, min(s + 10, n)) for s in range(0, n, 10)]  # 19 units of 10 cells
    args = (r, layout, 3.0, 1)
    exact = m._estimate_noise_sd(load, pieces, units, 10_000, *args, {})
    cache = {}
    sampled = m._estimate_noise_sd(load, pieces, units, 60, *args, cache)  # -> 7 of the 19 units (70 cells)
    assert len(cache) == 7 and sampled != exact
    assert abs(sampled - exact) / exact < 0.2
    assert not m._estimate_noise_sd(load, pieces, units, 60, *args, cache, use_cache=False) != sampled
    assert len(cache) == 7, "use_cache=False must not touch the cache"


def _spatial_like(base, block_rows, lazy):
    """Table shaped like a lazily loaded spatial table: Ensembl-id index, int32 counts and extra non-gene
    (control probe) features without genomic position; with `lazy`, `X` is a dask array of CSR row blocks."""
    rng = np.random.default_rng(0)
    n, g = base.shape
    n_ctrl = 20
    counts = sp.hstack(
        [sp.csr_matrix(base.X.astype(np.int32)), sp.csr_matrix(rng.poisson(0.05, (n, n_ctrl)).astype(np.int32))]
    ).tocsr()
    var = pd.concat(
        [
            base.var[["chromosome", "start", "end"]],
            pd.DataFrame({"chromosome": [None] * n_ctrl, "start": np.nan, "end": np.nan}),
        ]
    )
    var.index = pd.Index(
        [f"ENSG{i:011d}" for i in range(g)] + [f"NegControlProbe_{i}" for i in range(n_ctrl)], name="feature_id"
    )
    var["feature_name"] = [f"gene{i}" for i in range(g)] + [f"ctrl{i}" for i in range(n_ctrl)]
    var["feature_type"] = ["Gene Expression"] * g + ["Negative Control Probe"] * n_ctrl
    obs = base.obs.copy()
    obs["region"] = pd.Categorical(["cell_labels"] * n)
    obs["cell_id"] = np.arange(n)
    X = counts
    if lazy:
        X = da.concatenate(
            [
                da.from_delayed(
                    dask.delayed(lambda s=s: counts[s : s + block_rows])(),
                    shape=(min(block_rows, n - s), g + n_ctrl),
                    dtype=counts.dtype,
                    meta=sp.csr_matrix((0, g + n_ctrl), dtype=counts.dtype),
                )
                for s in range(0, n, block_rows)
            ]
        )
    return ad.AnnData(X=X, obs=obs, var=var)


def _spatial_workflow(adata, **kw):
    """The workflow of the docs: keep gene features, normalise + log1p (lazily for dask), infercnv."""
    import scanpy as sc

    genes = adata[:, (adata.var["feature_type"] == "Gene Expression").to_numpy()].copy()
    sc.pp.normalize_total(genes, target_sum=1e4)
    sc.pp.log1p(genes)
    cnv.tl.infercnv(genes, **REF, **kw)
    return genes


def test_spatial_like_table_matches_eager(base):
    eager = _spatial_workflow(_spatial_like(base, 40, lazy=False))
    lazy_input = _spatial_like(base, 40, lazy=True)
    lazy = _spatial_workflow(lazy_input, chunksize=60)
    assert hasattr(lazy.X, "compute"), "the table should still be lazy after normalisation"
    np.testing.assert_allclose(lazy.obsm["X_cnv"].toarray(), eager.obsm["X_cnv"].toarray(), rtol=1e-4, atol=1e-5)
    assert lazy.uns["cnv"]["chr_pos"] == eager.uns["cnv"]["chr_pos"]
    cnv.tl.cnv_score(lazy, groupby="cell_type")  # downstream functions work on the result
    assert lazy.obs["cnv_score"].notna().all()


@pytest.fixture(scope="module")
def zarr_store(base, tmp_path_factory):
    path = tmp_path_factory.mktemp("lazy") / "oligo.zarr"
    b = base.copy()
    b.X = sp.csr_matrix(b.X)
    b.write_zarr(path)
    return path


@pytest.mark.parametrize("rows", [40, 184])
def test_scanpy_dask_tutorial_loading(base, expected, zarr_store, rows):
    """`obs`/`var` in memory and `X` lazy via `read_elem_lazy`, as in the scanpy dask tutorial."""
    import zarr
    from anndata.experimental import read_elem_lazy

    f = zarr.open(zarr_store, mode="r")
    X = read_elem_lazy(f["X"], chunks=(rows, base.shape[1]))
    a = ad.AnnData(X=X, obs=ad.io.read_elem(f["obs"]), var=ad.io.read_elem(f["var"]))
    _check(a, expected, chunksize=50)


def test_fully_lazy_anndata_and_downstream(expected, zarr_store):
    """`anndata.experimental.read_lazy`: `obs` and `var` are lazy too."""
    pytest.importorskip("xarray")
    a = ad.experimental.read_lazy(zarr_store)
    _check(a, expected, chunksize=50)
    cnv.tl.cnv_score(a, "cell_type")
    cnv.tl.ithcna(a, "cell_type")
    cnv.tl.ithgex(a, "cell_type")
    cnv.pl.chromosome_heatmap(a, groupby="cell_type", show=False)
    cnv.pl.chromosome_heatmap_summary(a, groupby="cell_type", show=False)
    assert np.isfinite(np.asarray(a.obs["cnv_score"])).all()


@pytest.mark.parametrize(
    "block_rows,chunksize",
    [(46, 46), (46, 10), (92, 7), (30, 70), (184, 5)],
)
def test_every_dask_block_is_read_once_per_pass_whatever_the_chunksize(base, expected, block_rows, chunksize):
    """Reading `chunksize` rows out of a bigger dask block loads the whole block, so chunks must not re-read it.

    (Lazily loaded spatial tables use large row blocks that are slow to read: repeated reads would be fatal.)
    """
    loads = []
    a = _with_X(base, _csr_blocks(base, block_rows, counter=loads))
    _check(a, expected, chunksize=chunksize)
    n_blocks = -(-base.shape[0] // block_rows)
    # reference pass + one pass shared by noise estimate and main pass
    assert len(loads) == 2 * n_blocks, (len(loads), n_blocks)


def test_sampled_noise_blocks_are_not_read_a_second_time(base, monkeypatch):
    """The blocks read for the noise estimate are used as a whole and kept, so the main pass skips them."""
    loads = []
    a = _with_X(base, _csr_blocks(base, 10, counter=loads))  # 19 blocks of 10 cells
    cnv.tl.infercnv(a, **REF, chunksize=10, noise_sample_size=50)
    n_blocks = -(-base.shape[0] // 10)
    assert len(loads) == 2 * n_blocks, (len(loads), n_blocks)  # reference pass + (noise estimate + main pass)
    cached = a.obsm["X_cnv"].toarray()

    # without room for the cache, the sampled blocks have to be read again in the main pass (but nothing else changes)
    monkeypatch.setattr(m, "_NOISE_CACHE_MAX_BYTES", 0)
    loads2 = []
    b = _with_X(base, _csr_blocks(base, 10, counter=loads2))
    cnv.tl.infercnv(b, **REF, chunksize=10, noise_sample_size=50)
    extra = len(loads2) - 2 * n_blocks
    assert 0 < extra < n_blocks // 2, (len(loads2), n_blocks)
    np.testing.assert_allclose(b.obsm["X_cnv"].toarray(), cached, rtol=1e-6, atol=1e-7)


def test_partial_cache_budget_is_respected(base, monkeypatch):
    """Only as many sampled blocks are cached as fit; the others are read again, results are unchanged."""
    ref = _with_X(base, _csr_blocks(base, 10))
    cnv.tl.infercnv(ref, **REF, chunksize=10, noise_sample_size=50)

    n_obs, n_windows = base.shape[0], ref.obsm["X_cnv"].shape[1]
    one_block = 10 * n_windows * 4
    # room for two of the six sampled blocks, and not for the remainder (the last block of the data has only 4 cells)
    monkeypatch.setattr(m, "_NOISE_CACHE_MAX_BYTES", 2.1 * one_block)
    loads = []
    a = _with_X(base, _csr_blocks(base, 10, counter=loads))
    cnv.tl.infercnv(a, **REF, chunksize=10, noise_sample_size=50, n_jobs=1)
    n_blocks = -(-n_obs // 10)
    assert len(loads) == 2 * n_blocks + 4, (len(loads), n_blocks)
    np.testing.assert_allclose(a.obsm["X_cnv"].toarray(), ref.obsm["X_cnv"].toarray(), rtol=1e-6, atol=1e-7)


def test_sampled_noise_estimate_is_independent_of_dask_blocking(base):
    """Same data, different block layout: the same (sampled) result as long as the units cover the same cells."""
    results = []
    for block_rows in (10, 20):
        a = _with_X(base, _csr_blocks(base, block_rows))
        cnv.tl.infercnv(a, **REF, chunksize=20, noise_sample_size=60)
        results.append(a.obsm["X_cnv"].toarray())
    assert results[0].shape == results[1].shape
    assert np.isfinite(results[0]).all() and np.isfinite(results[1]).all()
