from typing import Any, Literal

import numpy as np
import scanpy as sc
from anndata import AnnData
from anndata.acc import A
from scverse_misc import deprecated_arg

from infercnvpy._acc import (
    DEPR_INPLACE,
    DEPR_INPLACE_COPY,
    DEPR_KEY_ADDED,
    DEPR_NEIGHBORS_KEY,
    DEPR_USE_REP,
    _dest,
    _get,
    _key,
    _legacy,
    _obsm_key,
    _write_obsm,
)
from infercnvpy._logging import log

from ._copykat import copykat
from ._infercnv import infercnv
from ._scores import cnv_score, ithcna, ithgex

__all__ = ["cnv_score", "copykat", "infercnv", "ithcna", "ithgex", "leiden", "pca", "tsne", "umap"]

DEFAULT_GRAPH = A.obsp["cnv_neighbors"]


def _legacy_graph(adata: AnnData, use: Any, neighbors_key: str | None) -> Any:
    """Translate the deprecated `neighbors_key` to the graph accessor scanpy's `use` takes."""
    if neighbors_key is None:
        return use
    return _legacy(use, neighbors_key, A.obsp[adata.uns[neighbors_key]["connectivities_key"]], name="neighbors_key")


def _ensure_pca(adata: AnnData, use: Any, what: str) -> None:
    """Run :func:`pca` with default parameters if `use` is the missing default PCA."""
    if _dest(use) == ("obsm", key := _obsm_key("cnv_pca")) and key not in adata.obsm:
        pca(adata)
        log.info("%s: computed %s with default parameters, as it was not in adata.obsm", what, key)


@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE_COPY)
def leiden(
    adata: AnnData,
    resolution: float = 1,
    *,
    use: Any = None,
    out: Any = False,
    flavor: Literal["igraph", "leidenalg"] | None = None,
    n_iterations: int | None = None,
    directed: bool | None = None,
    neighbors_key: str | None = None,
    key_added: str | None = None,
    inplace: bool | None = None,
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
    resolution
        Resolution of the clustering, see :func:`scanpy.tl.leiden`.
    use
        The neighborhood graph to cluster, e.g. `A.obsp["cnv_neighbors"]` (the default,
        written by :func:`infercnvpy.pp.neighbors`). Strings are resolved by
        :meth:`anndata.acc.AdAcc.resolve`, e.g. `"obsp.cnv_neighbors"`.
    out
        Where to write the cluster labels, e.g. `A.obs["cnv_leiden"]` (the default).
        If `True`, return the labels instead of writing them.
    flavor
        Leiden implementation, see :func:`scanpy.tl.leiden`. Default: `"igraph"`, or
        `"leidenalg"` when `partition_type` is passed, which only leidenalg supports.
    n_iterations
        Number of iterations of the Leiden algorithm. Default: `2` for `flavor="igraph"`,
        scanpy's default (`-1`, until convergence) for `flavor="leidenalg"`.
    directed
        Whether to treat the graph as directed. Default: `False` for `flavor="igraph"`
        (the only value igraph supports), scanpy's default for `flavor="leidenalg"`.
    neighbors_key
        Deprecated, use `use`.
    key_added
        Deprecated, use `out`.
    inplace
        Deprecated, copy `adata` before calling this function instead.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.leiden`, e.g. `rng`.

    Returns
    -------
    The cluster labels if `out=True`, else `None`. By default sets the following fields:

    `adata.obs["cnv_leiden"]` : :class:`pandas.Categorical`
        Cluster label of each cell.
    `adata.uns["cnv_leiden"]["params"]` : :class:`dict`
        Parameters of the clustering, as recorded by :func:`scanpy.tl.leiden`.
    """
    use = _legacy_graph(adata, use, neighbors_key)
    out = _legacy(out, key_added, A.obs[key_added] if key_added else None, name="key_added", default=False)
    if inplace is False:
        adata = adata.copy()
    if flavor is None:
        flavor = "leidenalg" if "partition_type" in kwargs else "igraph"
    if flavor == "igraph":
        n_iterations = 2 if n_iterations is None else n_iterations
        directed = False if directed is None else directed
    if n_iterations is not None:
        kwargs["n_iterations"] = n_iterations
    if directed is not None:
        kwargs["directed"] = directed
    res = sc.tl.leiden(
        adata,
        resolution,
        use=DEFAULT_GRAPH if use is None else use,
        out=A.obs["cnv_leiden"] if out is False else out,
        flavor=flavor,
        **kwargs,
    )
    return adata if inplace is False else res


@deprecated_arg("use_rep", DEPR_USE_REP)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE)
def pca(
    adata: AnnData,
    *,
    use: Any = None,
    out: Any = False,
    svd_solver: str = "arpack",
    zero_center: bool = False,
    use_rep: str | None = None,
    key_added: str | None = None,
    inplace: bool | None = None,
    **kwargs,
) -> np.ndarray | None:
    """Compute the PCA on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.pp.pca`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    use
        The CNV matrix, e.g. `A.obsm["X_cnv"]` (the default). Strings are resolved by
        :meth:`anndata.acc.AdAcc.resolve`, e.g. `"obsm.X_cnv"`.
    out
        Where to write the PCA coordinates, e.g. `A.obsm["X_cnv_pca"]` (the default).
        If `True`, return them instead of writing them.
    svd_solver
        See :func:`scanpy.pp.pca`.
    zero_center
        See :func:`scanpy.pp.pca`.
    use_rep
        Deprecated, use `use`.
    key_added
        Deprecated, use `out`.
    inplace
        Deprecated, use `out=True` to return the result.
    **kwargs
        Additional arguments passed to :func:`scanpy.pp.pca`.

    Returns
    -------
    The PCA coordinates if `out=True`, else `None`. By default sets the following field:

    `adata.obsm["X_cnv_pca"]` : :class:`numpy.ndarray`
        PCA coordinates of the CNV matrix, cells x components.
    """
    use = _legacy(use, use_rep, A.obsm[f"X_{use_rep}"], name="use_rep")
    out = _legacy(out, key_added, A.obsm[f"X_{key_added}"], name="key_added", default=False)
    out = True if inplace is False else out
    use = A.obsm[_obsm_key("cnv")] if use is None else use
    if (key := _key(use, "obsm")) not in adata.obsm:
        raise KeyError(f"{key} is not in adata.obsm. Did you run `tl.infercnv`?")
    res = sc.pp.pca(_get(adata, use), svd_solver=svd_solver, zero_center=zero_center, **kwargs)
    return _write_obsm(adata, res, out, "cnv_pca")


@deprecated_arg("neighbors_key", DEPR_NEIGHBORS_KEY)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE)
def umap(
    adata: AnnData,
    *,
    use: Any = None,
    out: Any = False,
    neighbors_key: str | None = None,
    key_added: str | None = None,
    inplace: bool | None = None,
    **kwargs,
) -> np.ndarray | None:
    """Compute the UMAP on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.tl.umap`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    use
        The neighborhood graph to lay out, e.g. `A.obsp["cnv_neighbors"]` (the default,
        written by :func:`infercnvpy.pp.neighbors`). Strings are resolved by
        :meth:`anndata.acc.AdAcc.resolve`.
    out
        Where to write the embedding, e.g. `A.obsm["X_cnv_umap"]` (the default).
        If `True`, return it instead of writing it.
    neighbors_key
        Deprecated, use `use`.
    key_added
        Deprecated, use `out`.
    inplace
        Deprecated, use `out=True` to return the result.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.umap`.

    Returns
    -------
    The UMAP coordinates if `out=True`, else `None`. By default sets the following field:

    `adata.obsm["X_cnv_umap"]` : :class:`numpy.ndarray`
        UMAP coordinates, cells x `n_components` (default 2).
    """
    use = _legacy_graph(adata, use, neighbors_key)
    out = _legacy(out, key_added, A.obsm[f"X_{key_added}"], name="key_added", default=False)
    out = True if inplace is False else out
    res = sc.tl.umap(adata, use=DEFAULT_GRAPH if use is None else use, out=True, **kwargs)
    return _write_obsm(adata, res, out, "cnv_umap")


@deprecated_arg("use_rep", DEPR_USE_REP)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE)
def tsne(
    adata: AnnData,
    *,
    use: Any = None,
    out: Any = False,
    use_rep: str | None = None,
    key_added: str | None = None,
    inplace: bool | None = None,
    **kwargs,
) -> np.ndarray | None:
    """Compute the t-SNE on the result of :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.tl.tsne`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    use
        The representation to embed, e.g. `A.obsm["X_cnv_pca"]` (the default; computed
        with :func:`infercnvpy.tl.pca` at default settings if missing). Strings are
        resolved by :meth:`anndata.acc.AdAcc.resolve`.
    out
        Where to write the embedding, e.g. `A.obsm["X_cnv_tsne"]` (the default).
        If `True`, return it instead of writing it.
    use_rep
        Deprecated, use `use`.
    key_added
        Deprecated, use `out`.
    inplace
        Deprecated, use `out=True` to return the result.
    **kwargs
        Additional arguments passed to :func:`scanpy.tl.tsne`.

    Returns
    -------
    The t-SNE coordinates if `out=True`, else `None`. By default sets the following fields:

    `adata.obsm["X_cnv_tsne"]` : :class:`numpy.ndarray`
        t-SNE coordinates, cells x 2.
    `adata.obsm["X_cnv_pca"]` : :class:`numpy.ndarray`
        Only if `use` is the default and it is missing: the PCA computed with
        :func:`infercnvpy.tl.pca` at default settings. This is set even if `out=True`.
    """
    use = _legacy(use, use_rep, A.obsm[f"X_{use_rep}"], name="use_rep")
    out = _legacy(out, key_added, A.obsm[f"X_{key_added}"], name="key_added", default=False)
    out = True if inplace is False else out
    use = A.obsm[_obsm_key("cnv_pca")] if use is None else use
    _ensure_pca(adata, use, "tsne")
    res = sc.tl.tsne(adata, use=use, out=True, **kwargs)
    return _write_obsm(adata, res, out, "cnv_tsne")
