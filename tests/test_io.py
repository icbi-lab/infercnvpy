import numpy as np
import numpy.testing as npt
import pytest

from infercnvpy.io import genomic_position_from_biomart, genomic_position_from_gtf
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
