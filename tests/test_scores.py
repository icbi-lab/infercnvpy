import numpy as np
import pytest

import infercnvpy as cnv


def test_ithgex(adata_ithgex):
    res = cnv.tl.ithgex(adata_ithgex, "group", inplace=False)
    assert res["A"] == 0
    assert res["B"] == pytest.approx(1.2628, abs=0.001)


def test_ithgex_use_raw(adata_ithgex):
    """`use_raw=None` (the default) uses `.X` even when `.raw` is set; `use_raw=True` uses `.raw`."""
    expected_x = cnv.tl.ithgex(adata_ithgex, "group", inplace=False)
    adata_ithgex.raw = adata_ithgex.copy()
    adata_ithgex.X = np.random.default_rng(0).normal(size=adata_ithgex.shape)
    assert cnv.tl.ithgex(adata_ithgex, "group", inplace=False) != expected_x
    assert cnv.tl.ithgex(adata_ithgex, "group", use_raw=True, inplace=False) == expected_x


def test_ithcna(adata_ithgex):
    res = cnv.tl.ithcna(adata_ithgex, "group", inplace=False)
    assert res["A"] == pytest.approx(1.053, abs=0.001)
    assert res["B"] == 0


def test_cnv_score(adata_ithgex):
    res = cnv.tl.cnv_score(adata_ithgex, "group", inplace=False)
    assert res["A"] == pytest.approx(2.25, abs=0.001)
    assert res["B"] == pytest.approx(2.5, abs=0.001)


def test_cnv_score_records_params(adata_ithgex):
    cnv.tl.cnv_score(adata_ithgex, "group", key_added="score")
    assert adata_ithgex.uns["score"] == {"params": {"groupby": "group", "use_rep": "cnv"}}


@pytest.mark.parametrize("func", [cnv.tl.ithgex, cnv.tl.ithcna])
def test_ith_single_cell_group(adata_ithgex, func):
    """A group with one cell has no cell pairs: its score is NaN, the other groups are unaffected."""
    adata_ithgex.obs["group"] = list("AAAAABBC")
    res = func(adata_ithgex, "group", inplace=False)
    assert np.isnan(res["C"])
    assert not np.isnan(res["A"])

    func(adata_ithgex, "group", key_added="ith")
    ith = adata_ithgex.obs["ith"]
    assert np.isnan(ith["c8"])
    assert not ith[["c1", "c6"]].isna().any()
