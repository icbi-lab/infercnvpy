from typing import Union
from ._infercnv import infercnv
import numpy as np
from anndata import AnnData
import scanpy as sc


def leiden(adata, *args, **kwargs):
    pass


def pca(
    adata: AnnData,
    svd_solver: str = "arpack",
    zero_center: bool = False,
    inplace: bool = False,
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


def umap(adata):
    pass
