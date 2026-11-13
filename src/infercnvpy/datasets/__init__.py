"""Example datasets."""

import importlib.resources as pkg_resources
from os import PathLike
from textwrap import dedent, indent
from typing import cast

import pooch
from anndata import AnnData, read_h5ad

from . import data

DATASET_ENV_VAR = "INFERCNVPY_DATA_DIR"

_EXAMPLEDATA = pooch.create(
    path=pooch.os_cache("infercnvpy"),
    base_url="https://github.com/scverse/infercnvpy/releases/download/d0.1.0/",
    env=DATASET_ENV_VAR,
    registry={
        "maynard2020_3k.h5ad": "md5:8a97c55dd72fb92d8e5ab72a52e23679",
    },
)

_POOCH_INFO = dedent(
    f"""\
    .. note::
        infercnvpy datasets are managed through `Pooch <https://github.com/fatiando/pooch>`_.

        By default, the dataset will be downloaded into your operating system's default
        cache directory (See :func:`pooch.os_cache` for more details). If it has already been
        downloaded, it will be retrieved from the cache.

        You can override the default cache dir by setting the `{DATASET_ENV_VAR}` environment variable
        to a path of your preference.
    """
)


def oligodendroglioma() -> AnnData:
    """The original inferCNV example dataset.

    Derived from :cite:`Tirosh2016`.
    """
    with pkg_resources.path(data, "oligodendroglioma.h5ad") as p:
        return read_h5ad(p)


def maynard2020_3k() -> AnnData:
    """Return the dataset from :cite:`Maynard2020` as AnnData object, downsampled to 3000 cells.

    In brief, this data set was processed as follows:
        * raw data downloaded from ENA
        * gene expression quantified using Salmon and the nf-core/rnaseq pipeline.
        * basic quality control (min_counts=20k, max_counts=5M, min_genes=1k, max_mitochondrial_fraction=0.2)
        * filtered to 6000 HVG using `sc.pp.highly_variable_genes(..., flavor="seurat_v3")`
        * raw counts processed using scVI, providing sample information as batch key.
        * cell types manually annotated based on marker genes and leiden clustering and subclustering.
        * downsampled to 3000 cells.

    `adata.X` contains the `log1p` transformed, cpm-normalized raw counts.
    The `scVI` latent representation is stored in `adata.obsm["X_scVI"]`.
    A UMAP for the 3000 cells is precomputed.

    {pooch_info}
    """
    fname = cast(PathLike, _EXAMPLEDATA.fetch("maynard2020_3k.h5ad", progressbar=True))
    return read_h5ad(fname)


maynard2020_3k.__doc__ = maynard2020_3k.__doc__.format(pooch_info=indent(_POOCH_INFO, " " * 4).lstrip())
