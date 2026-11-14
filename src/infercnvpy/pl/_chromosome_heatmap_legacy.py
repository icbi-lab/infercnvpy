"""`mode="legacy"` of the chromosome heatmaps: the implementation before 0.7.0, through :func:`scanpy.pl.heatmap`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib.colors import TwoSlopeNorm

from infercnvpy._util import _as_df, _as_matrix

try:
    from scanpy.plotting._utils import savefig_or_show
except ImportError:  # scanpy >= 1.13
    from scanpy.plotting.legacy._utils import savefig_or_show

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.colors import Colormap


def legacy_heatmap(
    adata: AnnData,
    *,
    groupby: str,
    use_rep: str,
    cmap: str | Colormap,
    figsize: tuple[float, float],
    show: bool | None,
    save: str | bool | None,
    **kwargs,
) -> dict[str, Axes] | None:
    """Every cell as a row of :func:`scanpy.pl.heatmap`, in the original order within each group."""
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`")
    tmp_adata = AnnData(X=X, obs=_as_df(adata.obs), uns=adata.uns)
    tmp_adata = tmp_adata[_stable_group_order(tmp_adata, groupby, kwargs.get("dendrogram", False))]

    # re-sort, as saving & loading anndata destroys the order
    chr_pos_dict = dict(sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    # center color map at 0
    X = _as_matrix(tmp_adata.X)
    tmp_data = X if isinstance(X, np.ndarray) else X.data
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    if vmin is None:
        vmin = np.nanmin(tmp_data)
    if vmax is None:
        vmax = np.nanmax(tmp_data)
    kwargs["norm"] = TwoSlopeNorm(0, vmin=vmin, vmax=vmax)

    # add chromosome annotations
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]], strict=False))

    return_ax_dic = sc.pl.heatmap(
        tmp_adata,
        var_names=tmp_adata.var.index.values,
        groupby=groupby,
        figsize=figsize,
        cmap=cmap,
        show_gene_labels=False,
        var_group_positions=var_group_positions,
        var_group_labels=list(chr_pos_dict.keys()),
        show=False,
        **kwargs,
    )

    return_ax_dic["heatmap_ax"].vlines(chr_pos[1:], lw=0.6, ymin=0, ymax=tmp_adata.shape[0])

    savefig_or_show("heatmap", show=show, save=save)
    show = sc.settings.autoshow if show is None else show
    if not show:
        return return_ax_dic
    return None


def _stable_group_order(adata: AnnData, groupby: str, dendrogram: bool | str) -> np.ndarray:
    """The cells sorted by group with a stable sort, in the group order `scanpy.pl.heatmap` will use.

    `scanpy.pl.heatmap` sorts the cells with pandas' default (unstable) sort, so the order of the cells within a
    group depends on the platform. pandas does not reorder an index that is already sorted, so cells passed in this
    order keep their original order within each group. With a dendrogram, the groups are in its order; it is
    computed here, as `scanpy.pl.heatmap` would, if `adata.uns` does not have it yet.
    """
    categories = _as_df(adata.obs)[groupby].astype("category").cat.categories
    if dendrogram:
        key = dendrogram if isinstance(dendrogram, str) else f"dendrogram_{groupby}"
        if key not in adata.uns:
            sc.tl.dendrogram(adata, groupby, key_added=key)
        categories = pd.Index(adata.uns[key]["categories_ordered"])
    rank = pd.Categorical(_as_df(adata.obs)[groupby], categories=categories).codes.astype(np.int64)
    rank[rank < 0] = len(categories)  # cells without a group last, as pandas sorts NaN
    return np.argsort(rank, kind="stable")


def legacy_summary(
    adata: AnnData,
    *,
    groupby: str,
    use_rep: str,
    cmap: str | Colormap,
    figsize: tuple[float, float],
    show: bool | None,
    save: str | bool | None,
    **kwargs,
) -> dict[str, Axes] | None:
    """Every group mean as 10 rows of :func:`scanpy.pl.heatmap`."""
    # this dirty hack repeats each row 10 times, since scanpy
    # heatmap cannot really handle it if there's just one observation
    # per row. Scanpy matrixplot is not an option, since it plots each
    # gene individually.
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`")
    groups = adata.obs[groupby].unique()
    tmp_obs = pd.DataFrame()
    tmp_obs[groupby] = np.hstack([np.repeat(x, 10) for x in groups])
    tmp_obs.index = tmp_obs.index.astype(str)

    def _get_group_mean(group):
        group_mean = np.mean(X[adata.obs[groupby].values == group, :], axis=0)
        if len(group_mean.shape) == 1:
            # derived from an array instead of sparse matrix -> 1 dim instead of 2
            group_mean = group_mean[np.newaxis, :]
        return group_mean

    tmp_adata = sc.AnnData(
        X=np.vstack([np.repeat(_get_group_mean(group), 10, axis=0) for group in groups]), obs=tmp_obs, uns=adata.uns
    )

    chr_pos_dict = dict(sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    # center color map at 0
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    if vmin is None:
        vmin = np.min(tmp_adata.X)
    if vmax is None:
        vmax = np.max(tmp_adata.X)
    kwargs["norm"] = TwoSlopeNorm(0, vmin=vmin, vmax=vmax)

    # add chromosome annotations
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]], strict=False))

    return_ax_dic = sc.pl.heatmap(
        tmp_adata,
        var_names=tmp_adata.var.index.values,
        groupby=groupby,
        figsize=figsize,
        cmap=cmap,
        show_gene_labels=False,
        var_group_positions=var_group_positions,
        var_group_labels=list(chr_pos_dict.keys()),
        show=False,
        **kwargs,
    )

    return_ax_dic["heatmap_ax"].vlines(chr_pos[1:], lw=0.6, ymin=-1, ymax=tmp_adata.shape[0])

    savefig_or_show("heatmap", show=show, save=save)
    show = sc.settings.autoshow if show is None else show
    if not show:
        return return_ax_dic
    return None
