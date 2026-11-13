import scanpy as sc

from ._chromosome_heatmap import chromosome_heatmap, chromosome_heatmap_summary

__all__ = ["chromosome_heatmap", "chromosome_heatmap_summary", "tsne", "umap"]


def umap(adata, **kwargs):
    """Plot the CNV UMAP.

    Thin wrapper around :func:`scanpy.pl.embedding`.

    Parameters
    ----------
    adata
        annotated data matrix with the result of :func:`infercnvpy.tl.umap` in `adata.obsm["X_cnv_umap"]`
    **kwargs
        Additional arguments passed to :func:`scanpy.pl.embedding`, e.g. `color`.

    Returns
    -------
    Whatever :func:`scanpy.pl.embedding` returns: `None`, or the axes if `show=False`.
    """
    return sc.pl.embedding(adata, "cnv_umap", **kwargs)


def tsne(adata, **kwargs):
    """Plot the CNV t-SNE.

    Thin wrapper around :func:`scanpy.pl.embedding`.

    Parameters
    ----------
    adata
        annotated data matrix with the result of :func:`infercnvpy.tl.tsne` in `adata.obsm["X_cnv_tsne"]`
    **kwargs
        Additional arguments passed to :func:`scanpy.pl.embedding`, e.g. `color`.

    Returns
    -------
    Whatever :func:`scanpy.pl.embedding` returns: `None`, or the axes if `show=False`.
    """
    return sc.pl.embedding(adata, "cnv_tsne", **kwargs)
