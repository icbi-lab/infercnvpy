import pytest
from . import TESTDATA
import scanpy as sc
import scipy.sparse


@pytest.fixture(params=["array", "csr", "csc"])
def adata_oligodendroma(request):
    adata = sc.read_h5ad(TESTDATA / "adata_oligodendroma.h5ad")

    if request.param == "array":
        adata.X = adata.X.toarray()
    elif request.param == "csr":
        adata.X = adata.X.tocsr()
    elif request.param == "csc":
        adata.X = adata.X.tocsc()
    else:
        raise AssertionError

    return adata
