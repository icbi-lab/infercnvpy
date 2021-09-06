import pytest
import infercnvpy as cnv
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import pandas as pd


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_oligodendroma(request):
    """Adata with raw counts in .X parametrized to be either sparse or dense."""
    adata = cnv.datasets.oligodendroglioma()

    adata.X = request.param(adata.X.toarray())

    return adata


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_infercnv(request):
    """Adata with infercnv computed and results stored in `.obsm["X_cnv"]`.
    The matrix in obsm is parametrized to be either sparse or dense."""
    adata = cnv.datasets.oligodendroglioma()
    cnv.tl.infercnv(adata)
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)

    adata.obsm["X_cnv"] = request.param(adata.obsm["X_cnv"].toarray())

    return adata


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_mock(request):
    obs = pd.DataFrame().assign(cat=["foo", "foo", "bar", "baz", "bar"])
    X = request.param(
        np.array(
            [
                [1, 1, 1, 2],
                [2, 1, 2, 2],
                [5, 5, 5, 5],
                [7, 5, 5, 7],
                [9, 9, 9, 9],
            ]
        )
    )
    adata = sc.AnnData(X=X, obs=obs)

    return adata
