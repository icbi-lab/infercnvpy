"""Vectorised / bounded-memory scores must agree with straightforward reference implementations."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import infercnvpy as cnv
from infercnvpy.pl._chromosome_heatmap import chromosome_heatmap_summary

da = pytest.importorskip("dask.array")


def _adata(n=300, p=40, n_groups=5, fmt="csr", seed=0):
    rng = np.random.default_rng(seed)
    X_cnv = rng.normal(size=(n, p)).astype(np.float32)
    X_cnv[rng.random((n, p)) < 0.6] = 0
    X_expr = rng.poisson(2.0, size=(n, p)).astype(np.float32)
    conv = {"csr": sp.csr_matrix, "csc": sp.csc_matrix, "dense": np.asarray}[fmt]
    obs = pd.DataFrame({"g": rng.choice([f"g{i}" for i in range(n_groups)], n)}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=conv(X_expr), obs=obs)
    a.obsm["X_cnv"] = conv(X_cnv)
    a.uns["cnv"] = {"chr_pos": {"chr1": 0, "chr2": p // 2}}
    return a


def _naive_cnv_score(a, groupby="g"):
    X = a.obsm["X_cnv"]
    return {c: np.mean(np.abs(X[a.obs[groupby].values == c, :])) for c in a.obs[groupby].unique()}


def _naive_ith(X, labels):
    out = {}
    for g in pd.unique(labels):
        rows = X[labels == g]
        rows = rows.toarray() if sp.issparse(rows) else rows
        if rows.shape[0] <= 1:
            continue
        q75, q25 = np.percentile(np.corrcoef(rows), [75, 25])
        out[g] = q75 - q25
    return out


@pytest.mark.parametrize("fmt", ["csr", "csc", "dense"])
def test_cnv_score_matches_naive(fmt):
    a = _adata(fmt=fmt)
    got = cnv.tl.cnv_score(a, "g", inplace=False)
    want = _naive_cnv_score(a)
    assert list(got) == list(want)  # same group order
    np.testing.assert_allclose(list(got.values()), list(want.values()), rtol=1e-5)
    cnv.tl.cnv_score(a, "g")
    np.testing.assert_allclose(a.obs["cnv_score"].to_numpy(), [want[g] for g in a.obs["g"]], rtol=1e-5)


@pytest.mark.parametrize("fmt", ["csr", "dense"])
def test_ith_matches_naive(fmt):
    a = _adata(fmt=fmt)
    labels = a.obs["g"].to_numpy()
    X_expr = a.X
    for fn, X in (
        (cnv.tl.ithgex, X_expr),
        (lambda *args, **kw: cnv.tl.ithcna(*args, use_rep="X_cnv", **kw), a.obsm["X_cnv"]),
    ):
        got = fn(a, "g", inplace=False)
        want = _naive_ith(X, labels)
        assert got.keys() == want.keys()
        np.testing.assert_allclose([got[k] for k in want], list(want.values()), rtol=1e-6, atol=1e-7)


def test_ithgex_dask_matches_in_memory():
    a = _adata()
    lazy = a.copy()
    lazy.X = da.from_array(a.X.toarray(), chunks=(70, a.shape[1]))
    got, want = cnv.tl.ithgex(lazy, "g", inplace=False), cnv.tl.ithgex(a, "g", inplace=False)
    np.testing.assert_allclose([got[k] for k in want], list(want.values()), rtol=1e-6)


def test_groups_with_a_single_cell_get_nan_instead_of_failing():
    a = _adata(n=30, n_groups=3)
    a.obs["g"] = ["big"] * 28 + ["lonely", "alone"]
    for fn in (lambda x: cnv.tl.ithgex(x, "g"), lambda x: cnv.tl.ithcna(x, "g"), lambda x: cnv.tl.cnv_score(x, "g")):
        fn(a)  # used to raise a KeyError
    assert a.obs["ithgex"].iloc[:28].notna().all() and a.obs["ithgex"].iloc[28:].isna().all()
    assert a.obs["cnv_score"].notna().all()


def test_cells_without_a_group_are_nan():
    a = _adata(n=40)
    a.obs["g"] = pd.Categorical([None if i % 10 == 0 else f"g{i % 3}" for i in range(40)])
    cnv.tl.cnv_score(a, "g")
    cnv.tl.ithgex(a, "g")
    missing = a.obs["g"].isna().to_numpy()
    assert a.obs.loc[missing, ["cnv_score", "ithgex"]].isna().all().all()
    assert a.obs.loc[~missing, ["cnv_score", "ithgex"]].notna().all().all()


def test_max_cells_subsampling_is_deterministic_and_close():
    rng = np.random.default_rng(3)
    base = rng.normal(size=(1, 30))
    X = (base + 0.5 * rng.normal(size=(1500, 30))).astype(np.float32)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"g": "all"}, index=[f"c{i}" for i in range(1500)]))
    exact = cnv.tl.ithgex(a, "g", inplace=False, max_cells=None)["all"]
    sub1 = cnv.tl.ithgex(a, "g", inplace=False, max_cells=300)["all"]
    sub2 = cnv.tl.ithgex(a, "g", inplace=False, max_cells=300)["all"]
    assert sub1 == sub2 and sub1 != exact
    assert abs(sub1 - exact) < 0.1 * exact


def test_heatmap_downsampling_keeps_color_scale_and_group_balance(monkeypatch):
    import infercnvpy.pl._chromosome_heatmap as hm

    a = _adata(n=600, n_groups=3)
    a.obs["g"] = a.obs["g"].astype("category")
    seen = {}
    real = hm.sc.pl.heatmap

    def spy(tmp_adata, **kw):
        seen["n"], seen["norm"], seen["groups"] = tmp_adata.shape[0], kw["norm"], tmp_adata.obs["g"].value_counts()
        return real(tmp_adata, **kw)

    monkeypatch.setattr(hm.sc.pl, "heatmap", spy)
    cnv.pl.chromosome_heatmap(a, groupby="g", max_cells_per_group=50, show=False)
    assert seen["n"] == 150 and (seen["groups"] <= 50).all()
    norm_small = (seen["norm"].vmin, seen["norm"].vmax)
    cnv.pl.chromosome_heatmap(a, groupby="g", max_cells_per_group=None, show=False)
    assert seen["n"] == 600
    assert norm_small == (seen["norm"].vmin, seen["norm"].vmax), "color scale must not depend on the subsample"


def test_heatmap_summary_matches_naive_group_means(monkeypatch):
    import infercnvpy.pl._chromosome_heatmap as hm

    a = _adata(n=200, n_groups=4)
    a.obs["g"] = a.obs["g"].astype("category")
    seen = {}
    real = hm.sc.pl.heatmap
    monkeypatch.setattr(
        hm.sc.pl, "heatmap", lambda t, **kw: seen.setdefault("X", np.asarray(t.X)) is not None and real(t, **kw)
    )
    chromosome_heatmap_summary(a, groupby="g", show=False)
    groups = a.obs["g"].unique()
    want = np.vstack(
        [
            np.repeat(np.asarray(a.obsm["X_cnv"][a.obs["g"].values == g].mean(axis=0)).reshape(1, -1), 10, axis=0)
            for g in groups
        ]
    )
    np.testing.assert_allclose(seen["X"], want, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize("fmt", ["csr", "csc", "dense"])
def test_group_means_matches_naive(fmt):
    from infercnvpy._util import _group_codes, _group_means

    a = _adata(n=120, p=25, n_groups=6, fmt=fmt)
    a.obs["g"] = pd.Categorical([None if i % 17 == 0 else g for i, g in enumerate(a.obs["g"])])  # incl. missing
    codes, groups = _group_codes(a.obs["g"])
    got = _group_means(a.obsm["X_cnv"], codes, len(groups))
    dense = a.obsm["X_cnv"].toarray() if sp.issparse(a.obsm["X_cnv"]) else a.obsm["X_cnv"]
    want = np.vstack([dense[(a.obs["g"] == g).to_numpy()].mean(axis=0) for g in groups])
    np.testing.assert_allclose(got, want, rtol=1e-5, atol=1e-7)
