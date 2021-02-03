import scanpy as sc
from anndata import AnnData


def umap(adata, **kwargs):
    """Plot the CNV UMAP.

    Thin wrapper around :func:`scanpy.pl.embedding`.
    """
    return sc.pl.embedding(adata, "cnv_umap", **kwargs)


def tsne(adata, **kwargs):
    """Plot the CNV t-SNE.

    Thin wrapper around :func:`scanpy.pl.embedding`.
    """
    return sc.pl.embedding(adata, "cnv_tsne", **kwargs)


def chromosome_heatmap(adata):
    pass
