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
