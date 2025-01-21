import numpy as np
import numpy.testing as npt
import pytest
import scanpy.datasets

from infercnvpy.io import genomic_position_from_biomart, genomic_position_from_gtf
from infercnvpy.tl import infercnv


@pytest.mark.parametrize(
    "adata,gtf,kwargs,expected_genes",
    [
        # weird gtfparse issue
        # [scanpy.datasets.pbmc3k(), "chr1_ensembl.gtf", {}, 42],
        [scanpy.datasets.pbmc3k(), "chr21_gencode.gtf", {}, 56],
        [scanpy.datasets.pbmc3k(), "chr21_gencode.gtf", {"adata_gene_id": "gene_ids", "gtf_gene_id": "gene_id"}, 116],
    ],
)
def test_get_genomic_position_from_gtf(adata, gtf, kwargs, testdata, expected_genes):
    genomic_position_from_gtf(testdata / gtf, adata, **kwargs)
    # those entries that are not null should start with "chr"
    assert all(adata.var["chromosome"].dropna().str.startswith("chr"))
    # start and end are equally populatedj
    npt.assert_array_equal(adata.var["start"].isnull().values, adata.var["end"].isnull().values)
    # most genes are covered. Matching ENSG we get almost all, using gene symbols only ~17k on this dataset
    assert np.sum(~adata.var["start"].isnull()) == expected_genes
    # can run infercnv on the result
    infercnv(adata)


@pytest.mark.parametrize(
    "adata,kwargs",
    [
        [scanpy.datasets.pbmc3k(), {"adata_gene_id": "gene_ids"}],
        [scanpy.datasets.pbmc3k(), {"biomart_gene_id": "hgnc_symbol"}],
    ],
)
def test_get_genomic_position_from_biomart(adata, kwargs):
    genomic_position_from_biomart(adata, use_cache=True, **kwargs)
    # those entries that are not null should start with "chr"
    assert all(adata.var["chromosome"].dropna().str.startswith("chr"))
    # start and end are equally populatedj
    npt.assert_array_equal(adata.var["start"].isnull().values, adata.var["end"].isnull().values)
    # most genes are covered. Matching ENSG we get almost all, using gene symbols only ~17k on this dataset
    assert np.sum(~adata.var["start"].isnull()) > 15000
    # can run infercnv on the result
    infercnv(adata)
