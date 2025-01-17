import pytest
import scanpy.datasets
from infercnvpy.io import genomic_position_from_gtf, genomic_position_from_biomart


def test_get_genomic_position_from_gtf():
    assert False


@pytest.mark.parametrize("adata,kwargs", [[scanpy.datasets.pbmc3k(), {}]])
def test_get_genomic_poisition_from_biomart(adata, kwargs):
    genomic_position_from_biomart(adata, **kwargs)
    assert all(adata.var["chromosome"].str.startswith("chr"))
    assert all(~adata.var["start"].isnull())
    assert all(~adata.var["end"].isnull())
