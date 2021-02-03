import scanpy as sc
from anndata import AnnData
import matplotlib.axes
from matplotlib.colors import Colormap, TwoSlopeNorm
from typing import Dict, Union, Tuple, Optional
import numpy as np
from scanpy.plotting._utils import savefig_or_show


def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: Union[str, Colormap] = "vlag",
    figsize: Tuple[int, int] = (16, 10),
    show: Optional[bool] = None,
    save: Union[str, bool, None] = None,
    **kwargs
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
    chr_pos = list(adata.uns[use_rep]["chr_pos"].values())

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
        var_group_labels=list(adata.uns["cnv"]["chr_pos"].keys()),
        norm=norm,
        show=False,
        **kwargs
    )

    return_ax_dic["heatmap_ax"].vlines(
        chr_pos[1:], lw=0.6, ymin=0, ymax=tmp_adata.shape[0]
    )

    savefig_or_show("heatmap", show=show, save=save)
    show = sc.settings.autoshow if show is None else show
    if not show:
        return return_ax_dic
