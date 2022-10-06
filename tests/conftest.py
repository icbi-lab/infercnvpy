import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

import infercnvpy as cnv


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


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_ithgex(request):
    return sc.AnnData(
        X=request.param(
            np.array(
                [
                    [1, 1, 1, 1, 1, 1, 2, 3],
                    [2, 2, 2, 2, 2, 2, 8, 0],
                    [3, 3, 3, 3, 3, 10, 3, 7],
                ]
            ).T
        ),
        obsm={
            "X_cnv": request.param(
                np.array(
                    [
                        [1, 1, 1, 2, 2, 1, 1, 1],
                        [2, 2, 2, 1, 1, 2, 2, 2],
                        [4, 4, 4, 2, 2, 3, 3, 3],
                        [2, 2, 2, 4, 4, 4, 4, 4],
                    ]
                ).T
            )
        },
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]).assign(
            group=list("AAAAABBB"),
        ),
        var=pd.DataFrame(index=["x", "y", "z"]),
    )
