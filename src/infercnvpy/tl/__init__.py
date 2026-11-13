from typing import Literal

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
    *,
    neighbors_key: str = "cnv_neighbors",
    key_added: str = "cnv_leiden",
    inplace: bool = True,
    flavor: Literal["igraph", "leidenalg"] | None = None,
    n_iterations: int | None = None,
    directed: bool | None = None,
    **kwargs,
):
    """Perform leiden clustering on the CNV neighborhood graph.

    Thin wrapper around :func:`scanpy.tl.leiden`. Unlike scanpy, it uses the
    `igraph` implementation by default, which is faster and does not raise scanpy's
    flavor warning. Pass `flavor="leidenalg"` for the previous default.

    .. versionchanged:: 0.7.0
        The default implementation changed from `leidenalg` to `igraph`, with
        `n_iterations=2` and `directed=False`, ahead of the same change in scanpy.
        Cluster labels differ from earlier versions.

    Parameters
    ----------
    adata
        Annotated data matrix.
    neighbors_key
        Key under which the result of :func:`infercnvpy.pp.neighbors` is stored in adata.
    key_added
        Key under which the cluster labels are stored in `adata.obs`.
    inplace
        If `True`, store the result in adata, otherwise return a copy of adata with the result.
    flavor
        Leiden implementation, see :func:`scanpy.tl.leiden`. Default: `"igraph"`, or
        `"leidenalg"` when `partition_type` is passed, which only leidenalg supports.
    n_iterations
        Number of iterations of the Leiden algorithm. Default: `2` for `flavor="igraph"`,
        scanpy's default (`-1`, until convergence) for `flavor="leidenalg"`.
    directed
        Whether to treat the graph as directed. Default: `False` for `flavor="igraph"`
        (the only value igraph supports), scanpy's default for `flavor="leidenalg"`.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.leiden`, e.g. `resolution`
        or `random_state`.

    Returns
    -------
    If `inplace=True`, returns `None` and sets the fields below on adata.
    Otherwise returns a copy of adata with these fields.
    Sets the following fields:

    `adata.obs[key_added]` : :class:`pandas.Categorical`
        Cluster label of each cell.
    `adata.uns[key_added]["params"]` : :class:`dict`
        Parameters of the clustering, as recorded by :func:`scanpy.tl.leiden`.
    """
    if flavor is None:
        flavor = "leidenalg" if "partition_type" in kwargs else "igraph"
    if flavor == "igraph":
        n_iterations = 2 if n_iterations is None else n_iterations
        directed = False if directed is None else directed
    if n_iterations is not None:
        kwargs["n_iterations"] = n_iterations
    if directed is not None:
        kwargs["directed"] = directed
    return sc.tl.leiden(
        adata,
        neighbors_key=neighbors_key,
        key_added=key_added,
        copy=not inplace,
        flavor=flavor,
        **kwargs,
    )


def pca(
    adata: AnnData,
    *,
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

    Returns
    -------
    If `inplace=True`, returns `None` and sets the field below.
    Otherwise returns the PCA coordinates.
    Sets the following fields:

    `adata.obsm["X_{key_added}"]` : :class:`numpy.ndarray`
        PCA coordinates of the CNV matrix, cells x components.
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
    *,
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

    Returns
    -------
    If `inplace=True`, returns `None` and sets the field below.
    Otherwise returns the UMAP coordinates.
    Sets the following fields:

    `adata.obsm["X_{key_added}"]` : :class:`numpy.ndarray`
        UMAP coordinates, cells x `n_components` (default 2).
    """
    tmp_adata = sc.tl.umap(adata, neighbors_key=neighbors_key, copy=True, **kwargs)

    if inplace:
        adata.obsm[f"X_{key_added}"] = tmp_adata.obsm["X_umap"]
        return None
    return tmp_adata.obsm["X_umap"]


def tsne(
    adata: AnnData,
    *,
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

    Returns
    -------
    If `inplace=True`, returns `None` and sets the field below.
    Otherwise returns the t-SNE coordinates.
    Sets the following fields:

    `adata.obsm["X_{key_added}"]` : :class:`numpy.ndarray`
        t-SNE coordinates, cells x 2.
    `adata.obsm["X_cnv_pca"]` : :class:`numpy.ndarray`
        Only if `use_rep="cnv_pca"` and it is missing: the PCA computed with
        :func:`infercnvpy.tl.pca` at default settings. This is set even if
        `inplace=False`.
    """
    if f"X_{use_rep}" not in adata.obsm and use_rep == "cnv_pca":
        pca(adata)
        log.info("tsne: computed X_cnv_pca with default parameters, as it was not in adata.obsm")
    tmp_adata = sc.tl.tsne(adata, use_rep=f"X_{use_rep}", copy=True, **kwargs)

    if inplace:
        adata.obsm[f"X_{key_added}"] = tmp_adata.obsm["X_tsne"]
        return None
    return tmp_adata.obsm["X_tsne"]
