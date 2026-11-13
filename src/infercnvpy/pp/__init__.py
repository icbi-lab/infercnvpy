import scanpy as sc
from anndata import AnnData

from infercnvpy import tl
from infercnvpy._logging import log

__all__ = ["neighbors"]


def neighbors(
    adata: AnnData,
    use_rep: str = "cnv_pca",
    key_added: str = "cnv_neighbors",
    inplace: bool = True,
    **kwargs,
):
    """Compute the neighborhood graph based on the result from :func:`infercnvpy.tl.infercnv`.

    Parameters
    ----------
    adata
        AnnData object
    use_rep
        Key under which the PCA of the results of :func:`infercnvpy.tl.infercnv`
        are stored in anndata. If not present, attempts to run :func:`infercnvpy.tl.pca`
        with default parameters.
    key_added
        Distances are stored in .obsp[key_added+’_distances’] and connectivities in
        .obsp[key_added+’_connectivities’].
    inplace
        If `True`, store the neighborhood graph in adata, otherwise return
        a copy of adata that holds the graph.
    **kwargs
        Arguments passed to :func:`scanpy.pp.neighbors`.

    Returns
    -------
    None | anndata.AnnData
        If `inplace=True`, returns `None` and sets the fields below on adata.
        Otherwise returns a copy of adata with these fields.

    Sets the following fields
    -------------------------
    `adata.obsp["{key_added}_distances"]` : :class:`scipy.sparse.csr_matrix`
        Distances between each cell and its nearest neighbors.
    `adata.obsp["{key_added}_connectivities"]` : :class:`scipy.sparse.csr_matrix`
        Weighted adjacency matrix of the neighborhood graph.
    `adata.uns[key_added]` : :class:`dict`
        Keys of the two matrices and the parameters, as recorded by
        :func:`scanpy.pp.neighbors`.
    `adata.obsm["X_cnv_pca"]` : :class:`numpy.ndarray`
        Only if `use_rep="cnv_pca"` and it is missing: the PCA computed with
        :func:`infercnvpy.tl.pca` at default settings. This is set on adata even
        if `inplace=False`.
    """
    if f"X_{use_rep}" not in adata.obsm and use_rep == "cnv_pca":
        tl.pca(adata)
        log.info("neighbors: computed X_cnv_pca with default parameters, as it was not in adata.obsm")

    return sc.pp.neighbors(adata, use_rep=f"X_{use_rep}", key_added=key_added, copy=not inplace, **kwargs)
