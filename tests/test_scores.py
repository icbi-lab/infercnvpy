import numpy as np
import pytest

import infercnvpy as cnv


def test_ithgex(adata_ithgex):
    res = cnv.tl.ithgex(adata_ithgex, "group", inplace=False)
    assert res["A"] == 0
    assert res["B"] == pytest.approx(1.2628, abs=0.001)


def test_ithcna(adata_ithgex):
    res = cnv.tl.ithcna(adata_ithgex, "group", inplace=False)
    assert res["A"] == pytest.approx(1.053, abs=0.001)
    assert res["B"] == 0


def test_cnv_score(adata_ithgex):
    res = cnv.tl.cnv_score(adata_ithgex, "group", inplace=False)
    assert res["A"] == pytest.approx(2.25, abs=0.001)
    assert res["B"] == pytest.approx(2.5, abs=0.001)


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
