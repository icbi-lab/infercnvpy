import scanpy as sc
from anndata import AnnData
import matplotlib.axes
from matplotlib.colors import Colormap, TwoSlopeNorm
from typing import Dict, Union, Tuple, Optional
import numpy as np
from scanpy.plotting._utils import savefig_or_show
import pandas as pd
from scipy.sparse import issparse


def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: Union[str, Colormap] = "bwr",
    figsize: Tuple[int, int] = (16, 10),
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    **kwargs,
) -> Optional[Dict[str, matplotlib.axes.Axes]]:
    """
    Plot a heatmap of smoothed gene expression by chromosome.

    Wrapper around :func:`scanpy.pl.heatmap`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use_rep
        Key under which the result from :func:`infercnvpy.tl.infercnv` are stored.
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
        raise ValueError(
            "'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?"
        )
    tmp_adata = AnnData(X=adata.obsm["X_cnv"], obs=adata.obs)

    # re-sort, as saving & loading anndata destroys the order
    chr_pos_dict = dict(
        sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1])
    )
    chr_pos = list(chr_pos_dict.values())

    # center color map at 0
    norm = TwoSlopeNorm(0, vmin=np.min(tmp_adata.X), vmax=np.max(tmp_adata.X))

    # add chromosome annotations
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]]))

    return_ax_dic = sc.pl.heatmap(
        tmp_adata,
        var_names=tmp_adata.var.index.values,
        groupby=groupby,
        figsize=figsize,
        cmap=cmap,
        show_gene_labels=False,
        var_group_positions=var_group_positions,
        var_group_labels=list(chr_pos_dict.keys()),
        norm=norm,
        show=False,
        **kwargs,
    )

    return_ax_dic["heatmap_ax"].vlines(
        chr_pos[1:], lw=0.6, ymin=0, ymax=tmp_adata.shape[0]
    )

    savefig_or_show("heatmap", show=show, save=save)
    show = sc.settings.autoshow if show is None else show
    if not show:
        return return_ax_dic


def chromosome_heatmap_summary(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: Union[str, Colormap] = "bwr",
    figsize: Tuple[int, int] = (16, 10),
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    **kwargs,
) -> Optional[Dict[str, matplotlib.axes.Axes]]:
    """
    Plot a heatmap of average of the smoothed gene expression by chromosome per
    category in groupby.

    Wrapper around :func:`scanpy.pl.heatmap`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use_rep
        Key under which the result from :func:`infercnvpy.tl.infercnv` are stored.
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
        raise ValueError(
            "'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?"
        )

    # TODO this dirty hack repeats each row 10 times, since scanpy
    # heatmap cannot really handle it if there's just one observation
    # per row. Scanpy matrixplot is not an option, since it plots each
    # gene individually.
    groups = adata.obs[groupby].unique()
    tmp_obs = pd.DataFrame()
    tmp_obs[groupby] = np.hstack([np.repeat(x, 10) for x in groups])

    def _get_group_mean(group):
        group_mean = np.mean(
            adata.obsm[f"X_{use_rep}"][adata.obs[groupby] == group, :], axis=0
        )
        if len(group_mean.shape) == 1:
            # derived from an array instead of sparse matrix -> 1 dim instead of 2
            group_mean = group_mean[np.newaxis, :]
        return group_mean

    tmp_adata = sc.AnnData(
        X=np.vstack(
            [np.repeat(_get_group_mean(group), 10, axis=0) for group in groups]
        ),
        obs=tmp_obs,
    )

    chr_pos_dict = dict(
        sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1])
    )
    chr_pos = list(chr_pos_dict.values())

    # center color map at 0
    norm = TwoSlopeNorm(0, vmin=np.min(tmp_adata.X), vmax=np.max(tmp_adata.X))

    # add chromosome annotations
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]]))

    return_ax_dic = sc.pl.heatmap(
        tmp_adata,
        var_names=tmp_adata.var.index.values,
        groupby=groupby,
        figsize=figsize,
        cmap=cmap,
        show_gene_labels=False,
        var_group_positions=var_group_positions,
        var_group_labels=list(chr_pos_dict.keys()),
        norm=norm,
        show=False,
        **kwargs,
    )

    return_ax_dic["heatmap_ax"].vlines(
        chr_pos[1:], lw=0.6, ymin=-1, ymax=tmp_adata.shape[0]
    )

    savefig_or_show("heatmap", show=show, save=save)
    show = sc.settings.autoshow if show is None else show
    if not show:
        return return_ax_dic
