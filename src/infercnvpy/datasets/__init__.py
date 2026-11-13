"""Example datasets."""

import os
from pathlib import Path
from textwrap import dedent, indent

import pooch
from anndata import AnnData, read_h5ad
from scverse_misc.datasets import DatasetEntry, DownloadCB, fetch, parse_registry, register_loader

DATASET_ENV_VAR = "INFERCNVPY_DATA_DIR"

_BASE_URL, _DATASETS = parse_registry(Path(__file__).parent / "registry.yaml")


# scverse-misc keys loaders by type name across all packages in the process, so ours uses the package name.
@register_loader("infercnvpy")
def _load_h5ad(entry: DatasetEntry, target: Path, download: DownloadCB, /, **kwargs: object) -> AnnData:
    # Download into the cache dir itself rather than `target` (`<cache dir>/<type>`),
    # so files cached by earlier infercnvpy versions are reused.
    return read_h5ad(download(entry.file(suffix=".h5ad"), dest=target.parent))


def _fetch(name: str) -> AnnData:
    cache_dir = os.environ.get(DATASET_ENV_VAR) or pooch.os_cache("infercnvpy")
    return fetch(_DATASETS[name], cache_dir, base_url=_BASE_URL)


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

    {pooch_info}
    """
    return _fetch("oligodendroglioma")


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
    return _fetch("maynard2020_3k")


for _func in (oligodendroglioma, maynard2020_3k):
    _func.__doc__ = (_func.__doc__ or "").format(pooch_info=indent(_POOCH_INFO, " " * 4).lstrip())
