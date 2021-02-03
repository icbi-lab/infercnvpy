from anndata import AnnData
import scanpy as sc
from scanpy import logging
from .. import tl


def neighbors(
    adata: AnnData,
    use_rep: str = "cnv_pca",
    key_added: str = "cnv_neighbors",
    inplace: bool = True,
    **kwargs,
):
    """Compute the neighborhood graph based on the result from
    :func:`infercnvpy.tl.infercnv`.

    Parameters
    ----------
    use_rep
        Key under which the PCA of the results of :func:`infercnvpy.tl.infercnv`
        are stored in anndata. If not present, attempts to run :func:`infercnvpy.tl.pca`
        with default parameters.
    key_added
        Distances are stored in .obsp[key_added+’_distances’] and connectivities in
        .obsp[key_added+’_connectivities’].
    inplace
        If `True`, store the neighborhood graph in adata, otherwise return
        the distance and connectivity matrices.
    **kwargs
        Arguments passed to :func:`scanpy.pp.neighbors`.

    Returns
    -------
    Depending on the value of inplace, updates anndata or returns the distance
    and connectivity matrices.
    """
    if f"X_{use_rep}" not in adata.obsm and use_rep == "cnv_pca":
        logging.warning(
            "X_cnv_pca not found in adata.obsm. Computing PCA with default parameters"
        )  # type: ignore
        tl.pca(adata)

    return sc.pp.neighbors(
        adata, use_rep=f"X_{use_rep}", key_added=key_added, copy=not inplace, **kwargs
    )
