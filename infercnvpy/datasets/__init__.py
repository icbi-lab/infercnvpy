"""Example datasets"""
import importlib.resources as pkg_resources
from . import data
import scanpy as sc
from anndata import AnnData
from scanpy import settings
from scanpy.readwrite import read


def oligodendroglioma() -> AnnData:
    """The original inferCNV example dataset.

    Derived from :cite:`Tirosh2016`.
    """
    with pkg_resources.path(data, "oligodendroglioma.h5ad") as p:
        return sc.read_h5ad(p)


def maynard2020_3k() -> AnnData:
    """\
    Return the dataset from :cite:`Maynard2020` as AnnData object, downsampled
    to 3000 cells.

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
    """
    url = "https://github.com/icbi-lab/infercnvpy/releases/download/d0.1.0/maynard2020_3k.h5ad"
    filename = settings.datasetdir / "maynard2020_3k.h5ad"
    adata = read(filename, backup_url=url)
    return adata
