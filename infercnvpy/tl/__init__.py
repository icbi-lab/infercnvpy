from typing import Union
from ._infercnv import infercnv, cnv_score
import numpy as np
from anndata import AnnData
import scanpy as sc
from scanpy import logging


def leiden(
    adata: AnnData,
    neighbors_key: str = "cnv_neighbors",
    key_added: str = "cnv_leiden",
    inplace: bool = True,
    **kwargs,
):
    """Perform leiden clustering on the CNV neighborhood graph.

    Thin wrapper around :func:`scanpy.tl.leiden`.
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
) -> Union[np.ndarray, None]:
    """Compute the PCA on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.tl.pca`.

    Parameters
    ----------
    adata
        annotated data matrix
    svd_solver
        See :func:`scanpy.tl.pca`.
    zero_center
        See :func:`scanpy.tl.pca`.
    inplace
        If True, store the result in adata.obsm. Otherwise return the PCA matrix.
    use_rep
        Key under which the result of infercnv is stored in adata
    key_added
        Key under which the result will be stored in adata.obsm if `inplace=True`.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.pca`.
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
    else:
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
    else:
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
        logging.warning(
            "X_cnv_pca not found in adata.obsm. Computing PCA with default parameters"
        )  # type: ignore
        pca(adata)
    tmp_adata = sc.tl.tsne(adata, use_rep=f"X_cnv_pca", copy=True, **kwargs)

    if inplace:
        adata.obsm[f"X_{key_added}"] = tmp_adata.obsm["X_tsne"]
    else:
        return tmp_adata.obsm["X_tsne"]
