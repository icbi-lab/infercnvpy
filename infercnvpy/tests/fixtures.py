import pytest
import infercnvpy as cnv


@pytest.fixture(params=["array", "csr", "csc"])
def adata_oligodendroma(request):
    adata = cnv.datasets.oligodendroglioma()

    if request.param == "array":
        adata.X = adata.X.toarray()
    elif request.param == "csr":
        adata.X = adata.X.tocsr()
    elif request.param == "csc":
        adata.X = adata.X.tocsc()
    else:
        raise AssertionError

    return adata
