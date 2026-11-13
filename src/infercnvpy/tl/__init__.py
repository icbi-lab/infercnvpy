import numpy as np
import scanpy as sc
from anndata import AnnData

from infercnvpy._logging import log

from ._copykat import copykat
from ._infercnv import infercnv
from ._scores import cnv_score, ithcna, ithgex

__all__ = ["cnv_score", "copykat", "infercnv", "ithcna", "ithgex", "leiden", "pca", "tsne", "umap"]


def leiden(
    adata: AnnData,
    neighbors_key: str = "cnv_neighbors",
    key_added: str = "cnv_leiden",
    inplace: bool = True,
    **kwargs,
):
    """Perform leiden clustering on the CNV neighborhood graph.

    Thin wrapper around :func:`scanpy.tl.leiden`.

    Parameters
    ----------
    adata
        annotated data matrix
    neighbors_key
        Key under which the result of :func:`infercnvpy.pp.neighbors` is stored
        in adata
    key_added
        Key under which the cluster labels will be stored in `adata.obs`.
    inplace
        If True, store the result in adata, otherwise return a copy of adata with the result.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.leiden`.

    Returns
    -------
    Depending on the value of `inplace`, either returns `None` or a copy of `adata`
    with the cluster labels in `adata.obs[key_added]`.
    """
    return sc.tl.leiden(
        adata,
        neighbors_key=neighbors_key,
        key_added=key_added,
        copy=not inplace,
        **kwargs,
    )


def pca(
    adata: AnnData,
    svd_solver: str = "arpack",
    zero_center: bool = False,
    inplace: bool = True,
    use_rep: str = "cnv",
    key_added: str = "cnv_pca",
    **kwargs,
) -> np.ndarray | None:
    """Compute the PCA on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.pp.pca`.

    Parameters
    ----------
    adata
        annotated data matrix
    svd_solver
        See :func:`scanpy.pp.pca`.
    zero_center
        See :func:`scanpy.pp.pca`.
    inplace
        If True, store the result in adata.obsm. Otherwise return the PCA matrix.
    use_rep
        Key under which the result of infercnv is stored in adata
    key_added
        Key under which the result will be stored in adata.obsm if `inplace=True`.
    **kwargs
        Additional arguments passed to :func:`scanpy.pp.pca`.
    """
    if f"X_{use_rep}" not in adata.obsm:
        raise KeyError(f"X_{use_rep} is not in adata.obsm. Did you run `tl.infercnv`?")

    pca_res = sc.tl.pca(
        adata.obsm[f"X_{use_rep}"],
        svd_solver=svd_solver,
        zero_center=zero_center,
        **kwargs,
    )
    if inplace:
        adata.obsm[f"X_{key_added}"] = pca_res
        return None
    return pca_res


def umap(
    adata: AnnData,
    neighbors_key: str = "cnv_neighbors",
    key_added: str = "cnv_umap",
    inplace: bool = True,
    **kwargs,
):
    """Compute the UMAP on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.tl.umap`

    Parameters
    ----------
    adata
        annotated data matrix
    neighbors_key
        Key under which the result of :func:`infercnvpy.pp.neighbors` is stored
        in adata
    key_added
        Key under which the result UMAP will be stored in adata.obsm
    inplace
        If True, store the result in adata.obsm, otherwise return the result of UMAP.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.umap`.
    """
    tmp_adata = sc.tl.umap(adata, neighbors_key=neighbors_key, copy=True, **kwargs)

    if inplace:
        adata.obsm[f"X_{key_added}"] = tmp_adata.obsm["X_umap"]
        return None
    return tmp_adata.obsm["X_umap"]


def tsne(
    adata: AnnData,
    use_rep: str = "cnv_pca",
    key_added: str = "cnv_tsne",
    inplace: bool = True,
    **kwargs,
):
    """Compute the t-SNE on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.tl.tsne`

    Parameters
    ----------
    adata
        annotated data matrix
    use_rep
        Key under which the result of :func:`infercnvpy.tl.pca` is stored
        in adata
    key_added
        Key under which the result of t-SNE will be stored in adata.obsm
    inplace
        If True, store the result in adata.obsm, otherwise return the result of t-SNE.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.tsne`.
    """
    if f"X_{use_rep}" not in adata.obsm and use_rep == "cnv_pca":
        pca(adata)
        log.info("tsne: computed X_cnv_pca with default parameters, as it was not in adata.obsm")
    tmp_adata = sc.tl.tsne(adata, use_rep=f"X_{use_rep}", copy=True, **kwargs)

    if inplace:
        adata.obsm[f"X_{key_added}"] = tmp_adata.obsm["X_tsne"]
        return None
    return tmp_adata.obsm["X_tsne"]
