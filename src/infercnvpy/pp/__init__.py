from typing import Any

import scanpy as sc
from anndata import AnnData
from anndata.acc import A
from scverse_misc import deprecated_arg

from infercnvpy import tl
from infercnvpy._acc import DEPR_INPLACE_COPY, DEPR_KEY_ADDED, DEPR_USE_REP, _legacy, _obsm_key

__all__ = ["neighbors"]


@deprecated_arg("use_rep", DEPR_USE_REP)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("inplace", DEPR_INPLACE_COPY)
def neighbors(
    adata: AnnData,
    *,
    use: Any = None,
    out: Any = False,
    use_rep: str | None = None,
    key_added: str | None = None,
    inplace: bool | None = None,
    **kwargs,
):
    """Compute the neighborhood graph based on the result from :func:`infercnvpy.tl.infercnv`.

    Thin wrapper around :func:`scanpy.pp.neighbors`.

    .. versionchanged:: 0.7.0
        The connectivities are stored in `adata.obsp["cnv_neighbors"]` instead of
        `adata.obsp["cnv_neighbors_connectivities"]`, as in scanpy>=1.13.

    Parameters
    ----------
    adata
        Annotated data matrix.
    use
        The representation to compute neighbors on, e.g. `A.obsm["X_cnv_pca"]` (the default;
        computed with :func:`infercnvpy.tl.pca` at default settings if missing).
        Strings are resolved by :meth:`anndata.acc.AdAcc.resolve`, e.g. `"obsm.X_cnv_pca"`.
    out
        Where to write the graph, e.g. `A.obsp["cnv_neighbors"]` (the default), see
        :func:`scanpy.pp.neighbors`. If `True`, return it instead of writing it.
    use_rep
        Deprecated, use `use`.
    key_added
        Deprecated, use `out`.
    inplace
        Deprecated, copy `adata` before calling this function instead.
    **kwargs
        Arguments passed to :func:`scanpy.pp.neighbors`.

    Returns
    -------
    The graph if `out=True`, else `None`. By default sets the following fields:

    `adata.obsp["cnv_neighbors"]` : :class:`scipy.sparse.csr_matrix`
        Weighted adjacency matrix of the neighborhood graph.
    `adata.obsp["cnv_neighbors_distances"]` : :class:`scipy.sparse.csr_matrix`
        Distances between each cell and its nearest neighbors.
    `adata.uns["cnv_neighbors"]` : :class:`dict`
        Keys of the two matrices and the parameters, as recorded by :func:`scanpy.pp.neighbors`.
    `adata.obsm["X_cnv_pca"]` : :class:`numpy.ndarray`
        Only if `use` is the default and it is missing: the PCA computed with
        :func:`infercnvpy.tl.pca` at default settings. This is set on adata even if
        `inplace=False`.
    """
    use = _legacy(use, use_rep, A.obsm[f"X_{use_rep}"], name="use_rep")
    out = _legacy(out, key_added, A.obsp[key_added] if key_added else None, name="key_added", default=False)
    use = A.obsm[_obsm_key("cnv_pca")] if use is None else use
    tl._ensure_pca(adata, use, "neighbors")
    if inplace is False:
        adata = adata.copy()
    res = sc.pp.neighbors(adata, use=use, out=A.obsp["cnv_neighbors"] if out is False else out, **kwargs)
    return adata if inplace is False else res
