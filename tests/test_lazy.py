"""infercnv on lazy and on-disk data: dask arrays (dense, or sparse blocks as in lazily loaded spatial tables) and
backed h5ad. Adapted from the tests of scverse/infercnvpy#192 by Stephen Williams."""

from contextlib import nullcontext

import anndata as ad
import dask
import dask.array as da
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

import infercnvpy as cnv
from infercnvpy.tl._infercnv import _get_reference

REF = {"reference_key": "cell_type", "reference_cat": ["Microglia/Macrophage", "Oligodendrocytes (non-malignant)"]}


@pytest.fixture(scope="module")
def base():
    return cnv.datasets.oligodendroglioma()


@pytest.fixture(scope="module")
def expected(base):
    """In-memory result, with gene values."""
    return cnv.tl.infercnv(base, **REF, calculate_gene_values=True, inplace=False)


def _with_X(base, X):
    return ad.AnnData(X=X, obs=base.obs.copy(), var=base.var.copy())


def _csr_blocks(base, block_rows, counter=None):
    """Dask array whose row blocks are CSR matrices loaded lazily, like a lazily loaded `spatialdata_io` table."""
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


def _check(adata, expected, gene_values=True, **kw):
    """The result equals the in-memory one up to the precision of the reference (float64 sums in another order)."""
    chr_pos, res, gv = cnv.tl.infercnv(adata, **REF, calculate_gene_values=gene_values, inplace=False, **kw)
    if hasattr(res, "compute"):
        res, gv = dask.compute(res, gv)
    assert isinstance(res, sp.csr_matrix)
    assert res.dtype == np.float32
    assert chr_pos == expected[0]
    np.testing.assert_allclose(res.toarray(), expected[1].toarray(), rtol=1e-5, atol=1e-6)
    if gene_values:
        np.testing.assert_allclose(gv, expected[2], rtol=1e-5, atol=1e-6)
    return res


@pytest.mark.parametrize("block_rows", [40, 500])
def test_dask_dense(base, expected, block_rows):
    X = da.from_array(base.X.toarray(), chunks=(block_rows, 1000))  # also split along genes
    _check(_with_X(base, X), expected)


@pytest.mark.parametrize("block_rows", [30, 184])
def test_dask_sparse_blocks(base, expected, block_rows):
    _check(_with_X(base, _csr_blocks(base, block_rows)), expected)


def test_dask_layer(base, expected):
    a = _with_X(base, sp.csr_matrix(base.X.shape, dtype=np.float32))
    a.layers["raw"] = _csr_blocks(base, 50)
    _check(a, expected, layer="raw")


@pytest.mark.parametrize("fmt", ["csr", "dense"])
def test_backed_h5ad(base, expected, tmp_path, fmt):
    X = sp.csr_matrix(base.X) if fmt == "csr" else base.X.toarray()
    path = tmp_path / "a.h5ad"
    _with_X(base, X).write_h5ad(path)
    a = ad.read_h5ad(path, backed="r")
    assert a.isbacked
    _check(a, expected, chunksize=60)
    cnv.tl.infercnv(a, **REF, chunksize=60)
    assert isinstance(a.obsm["X_cnv"], sp.csr_matrix)


@pytest.mark.parametrize("dynamic_threshold", [None, 1.5])
def test_reads_each_block_once_per_pass(base, expected, dynamic_threshold):
    """Reference means and the noise band (blocks with reference cells), and the result: one read of a block each."""
    loads = []
    a = _with_X(base, _csr_blocks(base, 46, counter=loads))
    cnv.tl.infercnv(a, **REF, dynamic_threshold=dynamic_threshold)
    n_blocks = -(-base.n_obs // 46)
    is_ref = base.obs["cell_type"].isin(REF["reference_cat"]).to_numpy()
    with_ref = sum(is_ref[s : s + 46].any() for s in range(0, base.n_obs, 46))
    assert 0 < with_ref < n_blocks
    assert len(loads) == n_blocks + with_ref * (1 if dynamic_threshold is None else 2)
    assert max(stop - start for start, stop in loads) <= 46, "more than one block was loaded at once"


def test_streaming_reference_matches_in_memory(base):
    lazy = _with_X(base, _csr_blocks(base, 25))
    for key, cats in ((None, None), (REF["reference_key"], REF["reference_cat"])):
        with pytest.warns(UserWarning) if key is None else nullcontext():
            got = _get_reference(lazy, key, cats, None, None)
        want = np.asarray(_get_reference(base, key, cats, None, None))
        np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize("scheduler", ["threads", "synchronous"])
def test_lazy_output(base, expected, scheduler):
    """With `lazy=True`, the result stays a dask array of CSR blocks and equals the eager result exactly."""
    X = _csr_blocks(base, 50)
    eager = cnv.tl.infercnv(_with_X(base, X), **REF, calculate_gene_values=True, inplace=False)
    a = _with_X(base, X)
    with dask.config.set(scheduler=scheduler):
        cnv.tl.infercnv(a, **REF, lazy=True, calculate_gene_values=True)
        assert isinstance(a.obsm["X_cnv"], da.Array)
        assert isinstance(a.layers["gene_values_cnv"], da.Array)
        res, gv = dask.compute(a.obsm["X_cnv"], a.layers["gene_values_cnv"])
    assert (sp.csr_matrix(res) != eager[1]).nnz == 0
    np.testing.assert_array_equal(gv, eager[2])
    _check(_with_X(base, X), expected, lazy=True)


def test_lazy_is_opt_in(base):
    """Without `lazy=True`, the result of a dask array is in memory, whatever the active dask scheduler."""
    a = _with_X(base, _csr_blocks(base, 100))
    with dask.config.set(scheduler="threads"):
        cnv.tl.infercnv(a, **REF)
    assert isinstance(a.obsm["X_cnv"], sp.csr_matrix)


def test_tools_that_need_cnv_in_memory(base):
    """Tools without block-wise reading tell the user to load a lazy CNV matrix; cnv_score reads it block by block."""
    a = _with_X(base, _csr_blocks(base, 100))
    cnv.tl.infercnv(a, **REF, lazy=True)
    eager = a.copy()
    eager.obsm["X_cnv"] = eager.obsm["X_cnv"].compute()
    a.obsm["spatial"] = np.random.default_rng(0).random((a.n_obs, 2))
    a.obsp["spatial_connectivities"] = sp.eye(a.n_obs, format="csr")
    calls = [
        lambda: cnv.tl.ithcna(a, "cell_type", use_rep="X_cnv"),
        lambda: cnv.tl.smooth_cnv(a),
        lambda: cnv.pl.spatial_cnv(a, region="chr1", show=False),
    ]
    for call in calls:
        with pytest.raises(TypeError, match=r"needs `obsm\['X_cnv'\]` in memory, got a dask array.*\.compute\(\)"):
            call()
    assert cnv.tl.cnv_score(a, "cell_type", inplace=False) == pytest.approx(
        cnv.tl.cnv_score(eager, "cell_type", inplace=False)
    )


def test_lazy_without_threshold(base, expected):
    res = cnv.tl.infercnv(_with_X(base, _csr_blocks(base, 50)), **REF, dynamic_threshold=None, lazy=True, inplace=False)
    want = cnv.tl.infercnv(base, **REF, dynamic_threshold=None, inplace=False)[1]
    np.testing.assert_allclose(res[1].compute().toarray(), want.toarray(), rtol=1e-5, atol=1e-5)


def test_block_size_does_not_change_the_result(base):
    results = [
        cnv.tl.infercnv(_with_X(base, _csr_blocks(base, rows)), **REF, inplace=False)[1].toarray() for rows in (10, 97)
    ]
    np.testing.assert_allclose(results[0], results[1], rtol=1e-5, atol=1e-6)


def test_lazily_normalised_table(base):
    """A dask array of CSR blocks, normalised lazily by scanpy, gives the eager result."""
    counts = sp.csr_matrix(base.X)

    def normalised(X):
        a = _with_X(base, X)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        return a

    lazy_input = _with_X(base, counts)
    lazy_input.X = da.from_array(counts, chunks=(40, -1), asarray=False)
    lazy = normalised(lazy_input.X)
    assert isinstance(lazy.X, da.Array)
    eager = normalised(counts)
    want = cnv.tl.infercnv(eager, **REF, inplace=False)
    _check(lazy, want, gene_values=False)


@pytest.fixture(scope="module")
def zarr_store(base, tmp_path_factory):
    path = tmp_path_factory.mktemp("lazy") / "oligo.zarr"
    b = base.copy()
    b.X = sp.csr_matrix(b.X)
    b.write_zarr(path)
    return path


def test_read_elem_lazy(base, expected, zarr_store):
    """`obs`/`var` in memory and `X` lazy via `read_elem_lazy`, as in the scanpy dask tutorial."""
    import zarr
    from anndata.experimental import read_elem_lazy

    f = zarr.open(zarr_store, mode="r")
    X = read_elem_lazy(f["X"], chunks=(50, base.n_vars))
    _check(ad.AnnData(X=X, obs=ad.io.read_elem(f["obs"]), var=ad.io.read_elem(f["var"])), expected)


def test_read_lazy(expected, zarr_store):
    """`anndata.experimental.read_lazy`: `obs` and `var` are lazy too."""
    a = ad.experimental.read_lazy(zarr_store)
    _check(a, expected)
    cnv.tl.infercnv(a, **REF)
    assert isinstance(a.obsm["X_cnv"], sp.csr_matrix)


def test_unplaced_and_extra_features(base, expected):
    """Features without a position (e.g. control probes of a spatial panel) are skipped in lazy input too."""
    n_ctrl = 20
    ctrl = sp.csr_matrix(np.random.default_rng(0).poisson(0.5, (base.n_obs, n_ctrl)).astype(np.float32))
    X = sp.hstack([sp.csr_matrix(base.X), ctrl]).tocsr()
    var = pd.concat(
        [base.var[["chromosome", "start", "end"]], pd.DataFrame({"chromosome": [None] * n_ctrl, "start": np.nan})]
    )
    var.index = list(base.var_names) + [f"NegControlProbe_{i}" for i in range(n_ctrl)]
    a = ad.AnnData(X=da.from_array(X, chunks=(60, -1), asarray=False), obs=base.obs.copy(), var=var)
    chr_pos, res, gv = cnv.tl.infercnv(a, **REF, calculate_gene_values=True, inplace=False)
    np.testing.assert_allclose(res.toarray(), expected[1].toarray(), rtol=1e-5, atol=1e-6)
    assert np.isnan(gv[:, -n_ctrl:]).all()
