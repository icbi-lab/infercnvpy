from typing import Any

import matplotlib.axes
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata.acc import A
from matplotlib.colors import Colormap, TwoSlopeNorm
from scanpy.plotting.legacy._utils import savefig_or_show  # TODO should not use private scanpy functions
from scverse_misc import deprecated_arg

from infercnvpy._acc import DEPR_USE_REP, _key, _legacy, _obsm_key, _uns_entry
from infercnvpy._util import _as_df, _as_matrix


def _cnv(adata: AnnData, use: Any, use_rep: str | None) -> tuple[Any, dict[str, int]]:
    """The CNV matrix and its chromosome start positions, sorted (saving and loading anndata loses the order)."""
    use = _legacy(use, use_rep, A.obsm[f"X_{use_rep}"], name="use_rep")
    key = _obsm_key("cnv") if use is None else _key(use, "obsm")
    chr_pos = dict(sorted(_uns_entry(adata, key)["chr_pos"].items(), key=lambda x: x[1]))
    return adata.obsm[key], chr_pos


@deprecated_arg("use_rep", DEPR_USE_REP)
def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use: Any = None,
    use_rep: str | None = None,
    cmap: str | Colormap = "bwr",
    figsize: tuple[int, int] = (16, 10),
    show: bool | None = None,
    save: str | bool | None = None,
    **kwargs,
) -> dict[str, matplotlib.axes.Axes] | None:
    """Plot a heatmap of smoothed gene expression by chromosome.

    Wrapper around :func:`scanpy.pl.heatmap`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use
        The CNV matrix, e.g. `A.obsm["X_cnv"]` (the default). Its chromosome positions are
        read from `adata.uns` under the same name without the `X_` prefix.
    use_rep
        Deprecated, use `use`.
    cmap
        colormap to use
    figsize
        (width, height) tuple in inches
    show
        Whether to show the figure or to return the axes objects.
    save
        If `True` or a `str`, save the figure. A string is appended to the default filename.
        Infer the filetype if ending on `{'.pdf', '.png', '.svg'}`.
    **kwargs
        Arguments passed on to :func:`scanpy.pl.heatmap`.

    Returns
    -------
    If `show` is False, a dictionary of axes.

    """
    if groupby == "cnv_leiden" and "cnv_leiden" not in adata.obs.columns:
        raise ValueError("'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?")
    X, chr_pos_dict = _cnv(adata, use, use_rep)
    tmp_adata = AnnData(X=_as_matrix(X), obs=_as_df(adata.obs), uns=adata.uns)
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


@deprecated_arg("use_rep", DEPR_USE_REP)
def chromosome_heatmap_summary(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use: Any = None,
    use_rep: str | None = None,
    cmap: str | Colormap = "bwr",
    figsize: tuple[int, int] = (16, 10),
    show: bool | None = None,
    save: str | bool | None = None,
    **kwargs,
) -> dict[str, matplotlib.axes.Axes] | None:
    """Plot a heatmap of average of the smoothed gene expression by chromosome per category in groupby.

    Wrapper around :func:`scanpy.pl.heatmap`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use
        The CNV matrix, e.g. `A.obsm["X_cnv"]` (the default). Its chromosome positions are
        read from `adata.uns` under the same name without the `X_` prefix.
    use_rep
        Deprecated, use `use`.
    cmap
        colormap to use
    figsize
        (width, height) tuple in inches
    show
        Whether to show the figure or to return the axes objects.
    save
        If `True` or a `str`, save the figure. A string is appended to the default filename.
        Infer the filetype if ending on `{'.pdf', '.png', '.svg'}`.
    **kwargs
        Arguments passed on to :func:`scanpy.pl.heatmap`.

    Returns
    -------
    If `show` is False, a dictionary of axes.

    """
    if groupby == "cnv_leiden" and "cnv_leiden" not in adata.obs.columns:
        raise ValueError("'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?")

    # TODO this dirty hack repeats each row 10 times, since scanpy
    # heatmap cannot really handle it if there's just one observation
    # per row. Scanpy matrixplot is not an option, since it plots each
    # gene individually.
    X, chr_pos_dict = _cnv(adata, use, use_rep)
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
