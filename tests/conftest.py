import os
from collections.abc import Callable, Generator, Iterable
from pathlib import Path

import matplotlib

matplotlib.use("agg")

import numpy as np
import pandas as pd
import pooch
import pytest
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData, read_h5ad
from matplotlib import pyplot as plt
from matplotlib.testing.compare import compare_images
from packaging.version import Version
from scverse_misc.logging import config as log_config

import infercnvpy as cnv
from infercnvpy.datasets import DATASET_ENV_VAR

IMAGES = Path(__file__).parent / "_images"
MATRIX_TYPES = [np.array, sp.csr_matrix, sp.csc_matrix]


@pytest.fixture(autouse=True)
def _logging_state():
    """Restore the shared scverse verbosity after each test."""
    verbosity_orig = log_config.verbosity
    yield
    log_config.verbosity = verbosity_orig


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help="Run tests that retrieve stuff from the internet. This increases test time.",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: Iterable[pytest.Item]) -> None:
    if config.getoption("--internet-tests"):
        return
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        if "internet" in item.keywords:
            item.add_marker(skip_internet)


@pytest.fixture
def testdata():
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def _oligodendroglioma() -> AnnData:
    return cnv.datasets.oligodendroglioma()


@pytest.fixture
def adata_oligodendroglioma(_oligodendroglioma) -> AnnData:
    """A fresh copy of the oligodendroglioma dataset with raw counts in .X (csr)."""
    return _oligodendroglioma.copy()


@pytest.fixture(params=MATRIX_TYPES)
def adata_oligodendroma(request, _oligodendroglioma) -> AnnData:
    """Adata with raw counts in .X parametrized to be either sparse or dense."""
    adata = _oligodendroglioma.copy()
    adata.X = request.param(adata.X.toarray())
    return adata


@pytest.fixture(scope="session")
def _adata_infercnv(_oligodendroglioma) -> AnnData:
    adata = _oligodendroglioma.copy()
    cnv.tl.infercnv(adata)
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)
    return adata


@pytest.fixture(params=MATRIX_TYPES)
def adata_infercnv(request, _adata_infercnv) -> AnnData:
    """Adata with infercnv computed and results stored in `.obsm["X_cnv"]`.

    infercnv runs once per session; the matrix in obsm is parametrized to be either sparse or dense.
    """
    adata = _adata_infercnv.copy()
    adata.obsm["X_cnv"] = request.param(adata.obsm["X_cnv"].toarray())
    return adata


@pytest.fixture(scope="session")
def _pbmc3k() -> AnnData:
    """scanpy's pbmc3k (symbols and Ensembl ids in var), fetched from scanpy's copy on the scverse CDN."""
    path = pooch.retrieve(
        "https://exampledata.scverse.org/scanpy/pbmc3k_raw.h5ad",
        known_hash="sha256:89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1",
        fname="pbmc3k_raw.h5ad",
        path=os.environ.get(DATASET_ENV_VAR) or pooch.os_cache("infercnvpy"),
    )
    return read_h5ad(path)


@pytest.fixture
def adata_pbmc3k(_pbmc3k) -> AnnData:
    return _pbmc3k.copy()


@pytest.fixture
def biomart_cache(tmp_path_factory, monkeypatch) -> None:
    """pybiomart writes `.pybiomart.sqlite` into the working directory; keep it in a session tmp dir."""
    path = tmp_path_factory.getbasetemp() / "biomart"
    path.mkdir(exist_ok=True)
    monkeypatch.chdir(path)


@pytest.fixture
def image_comparer(tmp_path) -> Generator[Callable[..., None]]:
    """Save the current figure and compare it to `tests/_images/<name>.png` (scanpy's pattern).

    To add or update a baseline, copy the `actual` image named in the failure message to `tests/_images`.
    """
    # baselines are drawn with matplotlib 3.11; 3.10 lays out the heatmap slightly differently
    if Version(matplotlib.__version__) < Version("3.11"):
        pytest.skip("image baselines need matplotlib>=3.11")
    plt.rcdefaults()

    def save_and_compare(name: str, *, tol: float = 15) -> None:
        __tracebackhide__ = True
        expected = IMAGES / f"{name}.png"
        actual = tmp_path / f"{name}.png"
        plt.savefig(actual, dpi=40)
        plt.close("all")
        if not expected.is_file():
            pytest.fail(f"No baseline at {expected}. The actual image is at {actual}.")
        if (msg := compare_images(str(expected), str(actual), tol=tol)) is not None:
            pytest.fail(str(msg))

    yield save_and_compare
    plt.close("all")


@pytest.fixture(params=MATRIX_TYPES)
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
    rng = np.random.RandomState(0)  # legacy generator: the expected values below were computed with it
    obs = pd.DataFrame().assign(sample=["sample1", "sample2", "sample3", "sample4"])  # Added "sample4"
    var = pd.DataFrame().assign(
        gene=["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"],
        start=[100, 200, 300, 400, 500, 0, 100, 200, 300, 400],
        end=[199, 299, 399, 499, 599, 99, 199, 299, 399, 499],
        chromosome=["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr2"],
    )
    var.index = var["gene"]
    X = sp.csr_matrix(rng.randint(low=0, high=50, size=(4, 10)))  # Changed size to (4, 10)
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


@pytest.fixture(params=MATRIX_TYPES)
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
