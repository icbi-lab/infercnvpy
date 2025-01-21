from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

import infercnvpy as cnv


@pytest.fixture()
def testdata():
    return Path(__file__).parent / "data"


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


@pytest.fixture
def adata_full_mock():
    np.random.seed(0)  # Set the seed to a fixed value
    obs = pd.DataFrame().assign(sample=["sample1", "sample2", "sample3", "sample4"])  # Added "sample4"
    var = pd.DataFrame().assign(
        gene=["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"],
        start=[100, 200, 300, 400, 500, 0, 100, 200, 300, 400],
        end=[199, 299, 399, 499, 599, 99, 199, 299, 399, 499],
        chromosome=["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr2"],
    )
    var.index = var["gene"]
    X = sp.csr_matrix(np.random.randint(low=0, high=50, size=(4, 10)))  # Changed size to (4, 10)
    adata = sc.AnnData(X=X, obs=obs, var=var)

    return adata


@pytest.fixture
def gene_res_actual():
    df = pd.DataFrame(
        {
            "gene1": [0.75, -1.00, 0.00, 0.00],
            "gene2": [0.00, 0.00, 0.75, 0.00],
            "gene3": [0.000000, 0.000000, 0.91666667, 0.000000],
            "gene4": [0.00, 0.00, 1.25, 0.00],
            "gene5": [-0.75, 0.00, 1.25, 0.00],
            "gene6": [0.000000, 0.000000, 0.000000, 0.921875],
            "gene7": [0.000000, 0.000000, 0.000000, 0.703125],
            "gene8": [0.0, 0.0, 0.0, 0.0],
            "gene9": [0.0, 0.0, 0.0, 0.0],
            "gene10": [0.75, 0.00, 0.00, 0.00],
        }
    )
    df.index = df.index.astype(str)
    return df


@pytest.fixture
def x_res_actual():
    arr = np.array(
        [
            [1.00, 0.00, 0.00, 0.00, 0.00, 1.00],
            [-1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
            [0.00, 1.25, 1.25, 0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00, 0.875, 0.00, 0.00],
        ]
    )
    return arr


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
