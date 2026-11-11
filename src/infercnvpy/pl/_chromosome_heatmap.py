import matplotlib.axes
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib.colors import Colormap, TwoSlopeNorm
from scanpy.plotting._utils import savefig_or_show
from scipy.sparse import issparse

from infercnvpy._util import _group_codes, _group_means, _group_members, _subsample


def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: str | Colormap = "bwr",
    figsize: tuple[int, int] = (16, 10),
    max_cells_per_group: int | None = 2000,
    random_state: int | None = 0,
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
    use_rep
        Key under which the result from :func:`infercnvpy.tl.infercnv` are stored.
    cmap
        colormap to use
    figsize
        (width, height) tuple in inches
    max_cells_per_group
        Plot at most this many randomly selected cells of every group. A heatmap of hundreds of thousands
        of cells is slow, needs a lot of memory and cannot be resolved by the eye anyway. The color scale is
        always determined from all cells. Set to `None` to plot all cells.
    random_state
        Seed for the selection of cells if a group has more than `max_cells_per_group` cells.
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
    X_cnv = adata.obsm[f"X_{use_rep}"]
    obs = adata.obs.to_memory() if hasattr(adata.obs, "to_memory") else adata.obs  # lazy obs (`read_lazy`)
    rng = np.random.default_rng(random_state)
    codes, groups = _group_codes(obs[groupby])
    keep = np.concatenate(
        [_subsample(members, max_cells_per_group, rng) for members in _group_members(codes, len(groups))]
        + [np.flatnonzero(codes < 0)]
    )
    keep.sort()
    if len(keep) < adata.shape[0]:
        tmp_adata = AnnData(X=X_cnv[keep], obs=obs.iloc[keep], uns=adata.uns)
    else:
        tmp_adata = AnnData(X=X_cnv, obs=obs, uns=adata.uns)

    # re-sort, as saving & loading anndata destroys the order
    chr_pos_dict = dict(sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    # center color map at 0 (using all cells, also if only some of them are plotted)
    tmp_data = X_cnv.data if issparse(X_cnv) else X_cnv
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


def chromosome_heatmap_summary(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
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
        raise ValueError("'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?")

    # TODO this dirty hack repeats each row 10 times, since scanpy
    # heatmap cannot really handle it if there's just one observation
    # per row. Scanpy matrixplot is not an option, since it plots each
    # gene individually.
    codes, groups = _group_codes(adata.obs[groupby])
    tmp_obs = pd.DataFrame()
    tmp_obs[groupby] = np.hstack([np.repeat(x, 10) for x in groups])
    tmp_obs.index = tmp_obs.index.astype(str)

    # mean of each group with a single pass over the data
    group_means = _group_means(adata.obsm[f"X_{use_rep}"], codes, len(groups))
    tmp_adata = sc.AnnData(X=np.repeat(group_means, 10, axis=0), obs=tmp_obs, uns=adata.uns)

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
