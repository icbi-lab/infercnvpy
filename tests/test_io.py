import logging

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import scanpy as sc
import scanpy.queries

import infercnvpy.io._scevan
from infercnvpy.io import genomic_position_from_biomart, genomic_position_from_gtf, read_scevan
from infercnvpy.tl import infercnv


@pytest.mark.parametrize(
    "gtf,kwargs,expected_genes",
    [
        # weird gtfparse issue
        # ["chr1_ensembl.gtf", {}, 42],
        ["chr21_gencode.gtf", {}, 56],
        ["chr21_gencode.gtf", {"adata_gene_id": "gene_ids", "gtf_gene_id": "gene_id"}, 116],
    ],
)
def test_get_genomic_position_from_gtf(adata_pbmc3k, gtf, kwargs, testdata, expected_genes):
    adata = adata_pbmc3k
    genomic_position_from_gtf(testdata / gtf, adata, **kwargs)
    # those entries that are not null should start with "chr"
    assert all(adata.var["chromosome"].dropna().str.startswith("chr"))
    # start and end are equally populatedj
    npt.assert_array_equal(adata.var["start"].isnull().values, adata.var["end"].isnull().values)
    # most genes are covered. Matching ENSG we get almost all, using gene symbols only ~17k on this dataset
    assert np.sum(~adata.var["start"].isnull()) == expected_genes
    # can run infercnv on the result
    infercnv(adata)


@pytest.mark.parametrize("gtf", ["chr1_ensembl.gtf", "chr21_gencode.gtf"])
def test_get_genomic_position_from_gtf_without_adata(gtf, testdata):
    gene_pos = genomic_position_from_gtf(testdata / gtf)
    assert {"chromosome", "start", "end", "gene_id", "gene_name"} <= set(gene_pos.columns)
    # the Ensembl GTF has no "chr" prefix, so it gets added
    assert gene_pos["chromosome"].str.startswith("chr").all()
    assert not gene_pos["gene_id"].str.contains(r"\.\d+$").any()
    assert not gene_pos["gene_name"].duplicated().any()


def _genepos_message(caplog, func):
    [msg] = [r.getMessage() for r in caplog.records if r.getMessage().startswith(f"{func}:")]
    return msg


@pytest.mark.parametrize(
    ("gtf_gene_id", "n_duplicated"),
    [
        # gene ids are unique: nothing is dropped, nothing is reported
        ("gene_id", 0),
        # two gene ids share the name "DUP": both rows are dropped
        ("gene_name", 1),
    ],
)
def test_genomic_position_from_gtf_duplicates(tmp_path, caplog, gtf_gene_id, n_duplicated):
    genes = [("ENSG1", "DUP", 100), ("ENSG2", "DUP", 300), ("ENSG3", "UNIQUE", 500)]
    gtf = tmp_path / "dup.gtf"
    gtf.write_text(
        "".join(
            f'chr1\tTEST\tgene\t{start}\t{start + 99}\t.\t+\t.\tgene_id "{gene_id}"; gene_name "{gene_name}";\n'
            for gene_id, gene_name, start in genes
        )
    )
    var_names = ["ENSG1", "ENSG2", "ENSG3"] if gtf_gene_id == "gene_id" else ["DUP", "UNIQUE"]
    adata = sc.AnnData(np.zeros((2, len(var_names))), var=pd.DataFrame(index=var_names))

    with caplog.at_level(logging.INFO, logger="scverse.infercnvpy"):
        genomic_position_from_gtf(gtf, adata, gtf_gene_id=gtf_gene_id)

    msg = _genepos_message(caplog, "genomic_position_from_gtf")
    if n_duplicated:
        assert f"{n_duplicated} with more than one entry for their gene_name" in msg
        assert pd.isna(adata.var.loc["DUP", "chromosome"])
    else:
        assert "more than one entry" not in msg
        assert adata.var["chromosome"].notna().all()


def test_genomic_position_from_biomart_duplicate_message(monkeypatch, caplog):
    """The duplicate message names BioMart as the source, not a GTF file."""
    annot = pd.DataFrame(
        {
            "ensembl_gene_id": ["ENSG1", "ENSG1", "ENSG2"],
            "start_position": [1, 1, 10],
            "end_position": [5, 5, 20],
            "chromosome_name": ["1", "1", "2"],
        }
    )
    monkeypatch.setattr(scanpy.queries, "biomart_annotations", lambda *args, **kwargs: annot.copy())
    adata = sc.AnnData(np.zeros((2, 3)), var=pd.DataFrame(index=["ENSG1", "ENSG2", "ENSG3"]))

    with caplog.at_level(logging.INFO, logger="scverse.infercnvpy"):
        genomic_position_from_biomart(adata)

    msg = _genepos_message(caplog, "genomic_position_from_biomart")
    assert "1 with no entry in BioMart" in msg
    assert "1 with more than one entry for their ensembl_gene_id" in msg
    assert "GTF" not in msg


@pytest.mark.internet
@pytest.mark.usefixtures("biomart_cache")
@pytest.mark.parametrize("kwargs", [{"adata_gene_id": "gene_ids"}, {"biomart_gene_id": "hgnc_symbol"}])
def test_get_genomic_position_from_biomart(adata_pbmc3k, kwargs):
    adata = adata_pbmc3k
    genomic_position_from_biomart(adata, use_cache=True, **kwargs)
    # those entries that are not null should start with "chr"
    assert all(adata.var["chromosome"].dropna().str.startswith("chr"))
    # start and end are equally populatedj
    npt.assert_array_equal(adata.var["start"].isnull().values, adata.var["end"].isnull().values)
    # most genes are covered. Matching ENSG we get almost all, using gene symbols only ~17k on this dataset
    assert np.sum(~adata.var["start"].isnull()) > 15000
    # can run infercnv on the result
    infercnv(adata)


@pytest.fixture
def scevan_dir(tmp_path, monkeypatch):
    """A SCEVAN result directory whose `.RData` files are read by a stand-in for `rdata` (4 of 5 cells kept), and a
    result table."""
    for name in ("s_CNAmtx.RData", "s_count_mtx_annot.RData"):
        (tmp_path / name).touch()
    frames = {
        "CNA_mtx_relat": pd.DataFrame(np.arange(8.0).reshape(2, 4), columns=["c0", "c1", "c2", "c4"]),
        "count_mtx_annot": pd.DataFrame({"seqnames": [1, 2]}),
    }
    monkeypatch.setattr(infercnvpy.io._scevan, "_read_rdata_frame", lambda path, name: frames[name])
    table = pd.DataFrame(
        {"class": ["tumor", "normal"] * 2, "confidentNormal": ["no"] * 4}, index=frames["CNA_mtx_relat"].columns
    )
    table.to_csv(tmp_path / "table.csv")
    return tmp_path


@pytest.mark.parametrize("subset", [True, False])
def test_read_scevan_reports_dropped_cells(scevan_dir, caplog, subset):
    adata = sc.AnnData(obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]))
    with caplog.at_level(logging.INFO, logger="scverse.infercnvpy"):
        read_scevan(adata, scevan_dir, subset=subset)
    msgs = [r.getMessage() for r in caplog.records if r.name == "scverse.infercnvpy"]
    if subset:
        assert msgs == ["read_scevan: kept 4 of 5 cells, dropped 1 that SCEVAN filtered out"]
        assert list(adata.obs_names) == ["c0", "c1", "c2", "c4"]
    else:
        assert msgs == []
        assert np.isnan(adata.obsm["X_scevan"][3]).all()
    assert adata.uns["scevan"]["chr_pos"] == {"chr1": 0, "chr2": 1}


def test_read_scevan_imports_rdata_lazily():
    import subprocess
    import sys

    code = "import sys, infercnvpy; assert 'rdata' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)


def _scevan_adata():
    return sc.AnnData(obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]))


def test_read_scevan_copy(scevan_dir):
    adata = _scevan_adata()
    res = read_scevan(adata, scevan_dir, scevan_res_table=scevan_dir / "table.csv", copy=True)
    assert "X_scevan" in res.obsm
    assert list(res.obs["scevan_class"]) == ["tumor", "normal"] * 2
    assert "X_scevan" not in adata.obsm
    assert read_scevan(adata, scevan_dir) is None
    assert "X_scevan" in adata.obsm


def test_read_scevan_deprecations(scevan_dir):
    adata = _scevan_adata()
    with pytest.warns(FutureWarning, match="Passing `scevan_res_table` to `read_scevan` positionally"):
        read_scevan(adata, scevan_dir, scevan_dir / "table.csv")
    assert "scevan_class" in adata.obs
    with pytest.warns(FutureWarning, match=r"inplace=False\)` returns a modified copy"):
        res = read_scevan(_scevan_adata(), scevan_dir, inplace=False)
    assert isinstance(res, sc.AnnData)
    assert "X_scevan" in res.obsm
