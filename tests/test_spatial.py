import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from scipy.spatial import cKDTree

import infercnvpy as cnv
from infercnvpy.tl._spatial import _hop_csr, _hop_dense


def _knn_graph(coords: np.ndarray, k: int = 6) -> sp.csr_matrix:
    """Binary kNN graph on the coordinates, symmetrised (like squidpy's `spatial_connectivities`)."""
    _, idx = cKDTree(coords).query(coords, k + 1)
    n = len(coords)
    G = sp.csr_matrix((np.ones(n * k), (np.repeat(np.arange(n), k), idx[:, 1:].ravel())), shape=(n, n))
    return ((G + G.T) > 0).astype(np.float64).tocsr()


def _reference(G: sp.csr_matrix, X, alpha, n_hops):
    """The smoothing operator built with scipy: diag(a) + diag((1 - a) / d) @ A, applied n_hops times."""
    A = G.tolil()
    A.setdiag(0)
    A = A.tocsr()
    d = np.asarray(A.sum(axis=1)).ravel()
    deg = np.diff(A.indptr)
    a = np.full(len(d), alpha) if alpha is not None else 1 / (deg + 1)
    a[d == 0] = 1
    S = sp.diags(a) + sp.diags(np.divide(1 - a, d, out=np.zeros_like(d), where=d > 0)) @ A
    Y = X.toarray() if sp.issparse(X) else X
    for _ in range(n_hops):
        Y = S @ Y
    return Y


@pytest.fixture
def graph_and_cnv():
    rng = np.random.default_rng(0)
    n, m = 300, 40
    G = sp.random(n, n, density=0.03, random_state=1, format="csr") + sp.eye(n)  # weighted, with a diagonal
    G = sp.lil_matrix(G)
    G[5, :] = 0  # a cell without neighbours
    X = rng.normal(size=(n, m)) * (rng.random((n, m)) < 0.3)
    return sp.csr_matrix(G), X.astype(np.float32)


@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix, sp.csc_matrix])
@pytest.mark.parametrize("alpha", [None, 0.0, 0.3, 1.0])
@pytest.mark.parametrize("n_hops", [1, 3])
def test_smooth_cnv_matches_scipy(graph_and_cnv, fmt, alpha, n_hops):
    G, X = graph_and_cnv
    adata = AnnData(obs=pd.DataFrame(index=[str(i) for i in range(X.shape[0])]))
    adata.obsm["X_cnv"] = fmt(X)
    adata.obsp["spatial_connectivities"] = G
    res = cnv.tl.smooth_cnv(adata, alpha=alpha, n_hops=n_hops, inplace=False)
    assert res.dtype == np.float32
    np.testing.assert_allclose(res, _reference(G, X, alpha, n_hops), rtol=1e-5, atol=1e-6)
    np.testing.assert_array_equal(res[5], X[5])  # no neighbours: unchanged


def test_serial_twin_matches(graph_and_cnv):
    G, X = graph_and_cnv
    args = (G.indptr, G.indices, G.data, -1.0)
    Xc = sp.csr_matrix(X)
    for fn, x in ((_hop_csr, (Xc.indptr, Xc.indices, Xc.data)), (_hop_dense, (X,))):
        par, ser = np.empty_like(X), np.empty_like(X)
        fn(*args, *x, par)
        fn(*args, *x, ser, serial=True)
        np.testing.assert_allclose(par, ser, rtol=1e-6)


def test_smooth_cnv_heatmap_round_trip(adata_infercnv):
    rng = np.random.default_rng(0)
    adata_infercnv.obsm["spatial"] = rng.random((adata_infercnv.n_obs, 2))
    adata_infercnv.obsp["spatial_connectivities"] = _knn_graph(adata_infercnv.obsm["spatial"])
    cnv.tl.smooth_cnv(adata_infercnv, n_hops=2)
    smoothed = adata_infercnv.obsm["X_cnv_smoothed"]
    assert smoothed.shape == adata_infercnv.obsm["X_cnv"].shape
    assert adata_infercnv.uns["cnv_smoothed"]["chr_pos"] == adata_infercnv.uns["cnv"]["chr_pos"]
    assert adata_infercnv.uns["cnv_smoothed"]["params"]["n_hops"] == 2
    cnv.pl.chromosome_heatmap(adata_infercnv, use_rep="cnv_smoothed", show=False)
    cnv.pl.spatial_cnv(
        adata_infercnv,
        region=["chr1", "chr7:0-60_000_000"],
        color="cnv_leiden",
        use_rep="cnv_smoothed",
        show=False,
    )


def test_window_positions(_adata_infercnv):
    uns = _adata_infercnv.uns["cnv"]
    start, end = uns["window_start"], uns["window_end"]
    assert len(start) == len(end) == _adata_infercnv.obsm["X_cnv"].shape[1]
    assert (end > start).all()
    bounds = [*sorted(uns["chr_pos"].values()), len(start)]
    for a, b in zip(bounds, bounds[1:], strict=False):
        assert (np.diff(start[a:b]) >= 0).all()


def test_smooth_cnv_errors(graph_and_cnv):
    G, X = graph_and_cnv
    adata = AnnData(obs=pd.DataFrame(index=[str(i) for i in range(X.shape[0])]), obsm={"X_cnv": X})
    with pytest.raises(KeyError, match="spatial_neighbors"):
        cnv.tl.smooth_cnv(adata)
    adata.obsp["spatial_connectivities"] = G
    with pytest.raises(ValueError, match="alpha"):
        cnv.tl.smooth_cnv(adata, alpha=2)
    with pytest.raises(ValueError, match="n_hops"):
        cnv.tl.smooth_cnv(adata, n_hops=0)


def _grid(side: int = 40, flip: float = 0.15, seed: int = 0):
    """Two halves of a grid (tumour on the left), with a fraction `flip` of the labels flipped."""
    rng = np.random.default_rng(seed)
    xy = np.stack(np.meshgrid(np.arange(side), np.arange(side)), -1).reshape(-1, 2).astype(float)
    truth = np.where(xy[:, 0] < side / 2, "tumour", "normal")
    noisy = truth.copy()
    flipped = rng.random(len(truth)) < flip
    noisy[flipped] = np.where(truth[flipped] == "tumour", "normal", "tumour")
    adata = AnnData(obs=pd.DataFrame({"call": pd.Categorical(noisy)}, index=[str(i) for i in range(len(xy))]))
    adata.obsm["spatial"] = xy
    adata.obsp["spatial_connectivities"] = _knn_graph(xy, k=8)
    return adata, truth


def test_spatial_regions_reduces_label_noise():
    adata, truth = _grid()
    before = np.mean(adata.obs["call"].to_numpy() == truth)
    cnv.tl.spatial_regions(adata, "call", beta=4)
    after = np.mean(adata.obs["call_spatial"].to_numpy() == truth)
    assert before < 0.9
    assert after > 0.98
    assert list(adata.obs["call_spatial"].cat.categories) == list(adata.obs["call"].cat.categories)
    assert adata.uns["call_spatial"]["params"]["n_sweeps"] < 20


def test_spatial_regions_beta_zero_and_determinism():
    adata, _ = _grid()
    res0 = cnv.tl.spatial_regions(adata, "call", beta=0, inplace=False)
    pd.testing.assert_series_equal(res0, adata.obs["call"], check_names=False)
    res1 = cnv.tl.spatial_regions(adata, "call", beta=4, inplace=False)
    res2 = cnv.tl.spatial_regions(adata, "call", beta=4, inplace=False)
    pd.testing.assert_series_equal(res1, res2)


def test_spatial_regions_probabilities_and_missing():
    adata, truth = _grid(flip=0.0)
    rng = np.random.default_rng(1)
    p_tumour = np.clip((truth == "tumour") * 0.4 + 0.3 + rng.normal(scale=0.2, size=len(truth)), 0.01, 0.99)
    adata.obsm["p"] = pd.DataFrame({"tumour": p_tumour, "normal": 1 - p_tumour}, index=adata.obs_names)
    calls = np.where(p_tumour > 0.5, "tumour", "normal").astype(object)
    calls[::7] = np.nan  # missing calls are filled from the probabilities and the neighbours
    adata.obs["call"] = pd.Categorical(calls)
    res = cnv.tl.spatial_regions(adata, "call", probabilities="p", beta=4, inplace=False)
    assert res.notna().all()
    assert np.mean(res.to_numpy() == truth) > np.mean(calls == truth) + 0.1
    assert np.mean(res.to_numpy() == truth) > 0.97


def test_spatial_regions_fills_missing_labels():
    adata, truth = _grid(flip=0.0)
    calls = truth.astype(object)
    calls[::5] = np.nan
    adata.obs["call"] = pd.Categorical(calls)
    res = cnv.tl.spatial_regions(adata, "call", inplace=False)
    assert res.notna().all()
    assert np.mean(res.to_numpy() == truth) > 0.99
