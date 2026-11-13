"""The arguments replaced by scanpy-style `use`/`out` keep working, with a FutureWarning."""

import numpy as np
import pytest
import scanpy as sc

import infercnvpy as cnv


@pytest.fixture
def adata():
    rng = np.random.default_rng(0)
    adata = sc.AnnData(
        rng.normal(size=(60, 4)),
        obs={"group": list("AB") * 30},
        obsm={"X_cnv": rng.normal(size=(60, 8)), "X_other": rng.normal(size=(60, 8))},
        uns={"cnv": {"chr_pos": {"chr1": 0, "chr2": 4}}},
    )
    adata.layers["counts"] = adata.X + 1
    return adata


@pytest.mark.parametrize(
    ("old", "new"),
    [
        (
            lambda a: cnv.tl.pca(a, use_rep="other", inplace=False),
            lambda a: cnv.tl.pca(a, use="obsm.X_other", out=True),
        ),
        (
            lambda a: cnv.tl.cnv_score(a, obs_key="group", use_rep="other", inplace=False),
            lambda a: cnv.tl.cnv_score(a, "group", use="obsm.X_other", out=True),
        ),
        (
            lambda a: cnv.tl.ithcna(a, "group", use_rep="X_other", inplace=False),
            lambda a: cnv.tl.ithcna(a, "group", use="obsm.X_other", out=True),
        ),
        (
            lambda a: cnv.tl.ithgex(a, "group", layer="counts", inplace=False),
            lambda a: cnv.tl.ithgex(a, "group", use="layers.counts", out=True),
        ),
    ],
    ids=["pca", "cnv_score", "ithcna", "ithgex"],
)
def test_returned(adata, old, new):
    with pytest.warns(FutureWarning, match="is deprecated"):
        res_old = old(adata)
    np.testing.assert_equal(res_old, new(adata))


def test_key_added(adata):
    with pytest.warns(FutureWarning, match="argument key_added is deprecated"):
        cnv.tl.pca(adata, key_added="old_pca")
    with pytest.warns(FutureWarning, match="argument key_added is deprecated"):
        cnv.tl.cnv_score(adata, "group", key_added="old_score")
    assert "X_old_pca" in adata.obsm
    assert "old_score" in adata.obs


def test_graph(adata):
    cnv.tl.pca(adata)
    with pytest.warns(FutureWarning, match="argument use_rep is deprecated"):
        cnv.pp.neighbors(adata, use_rep="cnv_pca", key_added="old")
    assert "old" in adata.obsp
    with pytest.warns(FutureWarning, match="argument inplace is deprecated"):
        res = cnv.tl.leiden(adata, inplace=False, use="obsp.old")
    assert "cnv_leiden" in res.obs
    assert "cnv_leiden" not in adata.obs
    with pytest.warns(FutureWarning, match="argument neighbors_key is deprecated"):
        cnv.tl.umap(adata, neighbors_key="old")
    assert "X_cnv_umap" in adata.obsm


def test_both_old_and_new(adata):
    with pytest.raises(TypeError, match="not both"), pytest.warns(FutureWarning):
        cnv.tl.pca(adata, use="obsm.X_cnv", use_rep="cnv")


def test_heatmap_use_rep(adata):
    with pytest.warns(FutureWarning, match="argument use_rep is deprecated"):
        cnv.pl.chromosome_heatmap(adata, groupby="group", use_rep="cnv", show=False)
