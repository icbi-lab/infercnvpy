import pytest

import infercnvpy as cnv
from infercnvpy.datasets import _DATASETS, DATASET_ENV_VAR


def test_oligodendroglioma():
    adata = cnv.datasets.oligodendroglioma()
    assert adata.shape == (184, 9450)
    assert {"chromosome", "start", "end"} <= set(adata.var.columns)


@pytest.mark.parametrize("name", _DATASETS)
def test_registry_entry(name):
    (file,) = _DATASETS[name].files
    assert file.sha256 is not None
    # files come only from exampledata.scverse.org, without fallback downloads
    assert not file.fallback_urls


@pytest.mark.internet
def test_maynard2020_3k(tmp_path, monkeypatch):
    monkeypatch.setenv(DATASET_ENV_VAR, str(tmp_path))
    adata = cnv.datasets.maynard2020_3k()
    assert adata.shape == tuple(_DATASETS["maynard2020_3k"].metadata["shape"])
    assert (tmp_path / "maynard2020_3k.h5ad").is_file()
