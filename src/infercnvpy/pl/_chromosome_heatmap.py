from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm, Colormap, ListedColormap, Normalize, TwoSlopeNorm
from matplotlib.path import Path as MplPath
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_size import Fixed

from infercnvpy._logging import log
from infercnvpy._numba import n_kernel_chunks, numba_threads
from infercnvpy._util import _as_df, _as_matrix
from infercnvpy.pl._heatmap_kernels import binned_sums, pixel_rows, sampled_rows

if TYPE_CHECKING:
    from collections.abc import Callable

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

# Sizes in inches, as `scanpy.pl.heatmap` lays them out at figsize=(16, 10).
_GROUP_W = 0.155  # colour blocks of the groups
_CBAR_W = 0.155
_DENDRO_W = 0.775
_CHR_H = 0.115  # chromosome brackets
_PAD_W = 0.03
_PAD_H = 0.05
_CBAR_MAX_H = 4.0

#: :func:`scanpy.pl.heatmap` arguments that `mode="fast"` implements; any other one switches to `mode="legacy"`.
FAST_KWARGS = frozenset(
    {"swap_axes", "standard_scale", "var_group_rotation", "dendrogram", "vmin", "vmax", "vcenter", "norm"}
    | {"show_gene_labels", "log"}
)


def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: str | Colormap = "bwr",
    figsize: tuple[float, float] = (16, 10),
    mode: Literal["fast", "legacy"] = "fast",
    aggregate: Literal["sample", "mean"] = "sample",
    ax: Axes | None = None,
    n_jobs: int | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    return_fig: bool | None = None,
    **kwargs,
) -> Figure | dict[str, Axes] | None:
    """Plot a heatmap of smoothed gene expression by chromosome.

    Cells are sorted by `groupby` and each group gets a share of the heatmap's height proportional to its size.

    .. versionchanged:: 0.7.0
       Rendered by infercnvpy instead of :func:`scanpy.pl.heatmap` (`mode="fast"`): each pixel row of the figure
       shows one cell of the cells that fall into it, or their mean, so large data plots in seconds. The previous
       implementation is available as `mode="legacy"`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use_rep
        Key under which the result from :func:`infercnvpy.tl.infercnv` are stored. With `mode="fast"`,
        `adata.obsm[f"X_{use_rep}"]` can also be a dask array or a backed matrix, which is read block by block.
    cmap
        colormap to use
    figsize
        (width, height) tuple in inches. Ignored if `ax` is given.
    mode
        `"fast"` draws the heatmap with infercnvpy at the resolution of the figure. `"legacy"` draws every cell
        with :func:`scanpy.pl.heatmap`, which needs the CNV matrix as a dense table and several times its memory.
        `mode="fast"` implements the :func:`scanpy.pl.heatmap` arguments `swap_axes`, `standard_scale`,
        `var_group_rotation`, `dendrogram`, `vmin`, `vmax`, `vcenter`, `norm`, `show_gene_labels` (ignored) and
        `log`; with any other argument (or a numeric, not categorical, `groupby` column), it logs an info message
        and uses `mode="legacy"`.
    aggregate
        What a pixel row shows when there are more cells than pixel rows (`mode="fast"` only). `"sample"` shows the
        cell in the middle of the cells that fall into the row, as the full-resolution plot downscaled with
        nearest-neighbour interpolation, with the colour scale from the values of all cells. `"mean"` shows their
        mean, with the colour scale from the range of the drawn (averaged) values unless `vmin` and `vmax` are
        given.
    ax
        Draw the heatmap into this axes (`mode="fast"` only). The group colours, chromosome labels, dendrogram and
        colour bar go into axes split off from it.
    n_jobs
        Number of threads for reading the CNV matrix (`mode="fast"`). `None` uses all cores.
    show
        Whether to show the figure. Defaults to :attr:`scanpy.settings.autoshow`.
    save
        If `True` or a `str`, save the figure to :attr:`scanpy.settings.figdir`. A string is appended to the default
        filename `heatmap`; it sets the file type if it ends in `.pdf`, `.png` or `.svg`, otherwise
        :attr:`scanpy.settings.file_format_figs` does.
    return_fig
        Return the :class:`~matplotlib.figure.Figure`.
    **kwargs
        Arguments of :func:`scanpy.pl.heatmap`, see `mode`. By default, the colour scale is a
        :class:`~matplotlib.colors.TwoSlopeNorm` centred at `vcenter=0`. `dendrogram=True` orders the groups by
        hierarchical clustering of their mean CNV profiles (Pearson correlation, complete linkage, as
        :func:`scanpy.tl.dendrogram`) and draws the dendrogram next to the heatmap.

    Returns
    -------
    The figure if `return_fig` is true. Otherwise, if `show` is false, a dictionary of axes with the keys
    `heatmap_ax`, `groupby_ax` (group colours), `gene_groups_ax` (chromosome labels), `colorbar_ax` (not with
    `mode="legacy"`) and, with `dendrogram`, `dendrogram_ax`.
    """
    if aggregate not in ("sample", "mean"):
        raise ValueError(f"`aggregate` must be 'sample' or 'mean', got {aggregate!r}.")
    if _use_legacy("chromosome_heatmap", adata, groupby, mode, kwargs):
        return _legacy(_legacy_heatmap, adata, groupby, use_rep, cmap, figsize, ax, show, save, return_fig, kwargs)
    opts = _Options.pop(kwargs)
    fig, ax = _figure(ax, figsize)
    cat = _categories(adata, groupby)
    palette = _group_colors(adata, groupby, len(cat.categories))
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`", lazy=True)
    dendro = None
    with numba_threads(n_jobs):
        if opts.dendrogram:
            cat, palette, dendro = _order_by_dendrogram(
                adata, groupby, cat, palette, lambda: _group_means(X, cat.codes.astype(np.int64), len(cat.categories))
            )
        codes = cat.codes.astype(np.int64)
        order, bounds, rows = pixel_rows(codes, len(cat.categories), _n_pixel_rows(fig, ax, swap=opts.swap_axes))
        if aggregate == "sample":
            values, lo, hi = sampled_rows(X, order[(bounds[:-1] + bounds[1:]) // 2])
        else:
            values = binned_sums(X, order, bounds) / np.diff(bounds)[:, None]
    values = opts.transform(values)
    if aggregate == "mean":
        lo, hi = np.nanmin(values), np.nanmax(values)
        center = 0 if opts.vcenter is None else opts.vcenter
        lo, hi = min(lo, center - 1e-9), max(hi, center + 1e-9)  # TwoSlopeNorm needs vmin < vcenter < vmax
    norm = opts.color_norm(lo, hi)
    axes = _draw(fig, ax, adata, groupby, cat, palette, values, rows, dendro, use_rep, norm, cmap, opts)
    return _finish(fig, axes, show=show, save=save, return_fig=return_fig)


def chromosome_heatmap_summary(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: str | Colormap = "bwr",
    figsize: tuple[float, float] = (16, 10),
    mode: Literal["fast", "legacy"] = "fast",
    ax: Axes | None = None,
    n_jobs: int | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    return_fig: bool | None = None,
    **kwargs,
) -> Figure | dict[str, Axes] | None:
    """Plot a heatmap of average of the smoothed gene expression by chromosome per category in groupby.

    Each group is one row of equal height.

    .. versionchanged:: 0.7.0
       Rendered by infercnvpy instead of :func:`scanpy.pl.heatmap` (`mode="fast"`). The previous implementation is
       available as `mode="legacy"`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        group the cells by a categorical variable from adata.obs. It usually makes
        sense to either group by unsupervised clustering obtained from
        :func:`infercnvpy.tl.leiden` (the default) or a cell-type label.
    use_rep
        Key under which the result from :func:`infercnvpy.tl.infercnv` are stored. With `mode="fast"`,
        `adata.obsm[f"X_{use_rep}"]` can also be a dask array or a backed matrix, which is read block by block.
    cmap
        colormap to use
    figsize
        (width, height) tuple in inches. Ignored if `ax` is given.
    mode
        `"fast"` computes the group means in one pass over the CNV matrix and draws them with infercnvpy.
        `"legacy"` uses :func:`scanpy.pl.heatmap`. `mode="fast"` implements the :func:`scanpy.pl.heatmap` arguments
        `swap_axes`, `standard_scale`, `var_group_rotation`, `dendrogram`, `vmin`, `vmax`, `vcenter`, `norm`,
        `show_gene_labels` (ignored) and `log`; with any other argument (or a numeric, not categorical, `groupby`
        column), it logs an info message and uses `mode="legacy"`.
    ax
        Draw the heatmap into this axes (`mode="fast"` only). The group colours, chromosome labels, dendrogram and
        colour bar go into axes split off from it.
    n_jobs
        Number of threads for reading the CNV matrix (`mode="fast"`). `None` uses all cores.
    show
        Whether to show the figure. Defaults to :attr:`scanpy.settings.autoshow`.
    save
        If `True` or a `str`, save the figure to :attr:`scanpy.settings.figdir`. A string is appended to the default
        filename `heatmap`; it sets the file type if it ends in `.pdf`, `.png` or `.svg`, otherwise
        :attr:`scanpy.settings.file_format_figs` does.
    return_fig
        Return the :class:`~matplotlib.figure.Figure`.
    **kwargs
        Arguments of :func:`scanpy.pl.heatmap`, see `mode`. By default, the colour scale is a
        :class:`~matplotlib.colors.TwoSlopeNorm` centred at `vcenter=0` from the smallest to the largest group mean.

    Returns
    -------
    The figure if `return_fig` is true. Otherwise, if `show` is false, a dictionary of axes with the keys
    `heatmap_ax`, `groupby_ax` (group colours), `gene_groups_ax` (chromosome labels), `colorbar_ax` (not with
    `mode="legacy"`) and, with `dendrogram`, `dendrogram_ax`.
    """
    if _use_legacy("chromosome_heatmap_summary", adata, groupby, mode, kwargs):
        return _legacy(_legacy_summary, adata, groupby, use_rep, cmap, figsize, ax, show, save, return_fig, kwargs)
    opts = _Options.pop(kwargs)
    fig, ax = _figure(ax, figsize)
    cat = _categories(adata, groupby)
    palette = _group_colors(adata, groupby, len(cat.categories))
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`", lazy=True)
    with numba_threads(n_jobs):
        means = _group_means(X, cat.codes.astype(np.int64), len(cat.categories))
    dendro = None
    if opts.dendrogram:
        present = np.flatnonzero(np.bincount(cat.codes[cat.codes >= 0], minlength=len(cat.categories)))
        before = list(cat.categories[present])
        cat, palette, dendro = _order_by_dendrogram(adata, groupby, cat, palette, lambda: means)
        means = means[[before.index(c) for c in cat.categories if c in before]]
    rows = (np.bincount(cat.codes[cat.codes >= 0], minlength=len(cat.categories)) > 0).astype(np.int64)
    norm = opts.color_norm(float(np.min(means)), float(np.max(means)))
    axes = _draw(fig, ax, adata, groupby, cat, palette, opts.transform(means), rows, dendro, use_rep, norm, cmap, opts)
    return _finish(fig, axes, show=show, save=save, return_fig=return_fig)


# --- options of mode="fast" -------------------------------------------------------------------------------------------


@dataclass
class _Options:
    """The :func:`scanpy.pl.heatmap` arguments that `mode="fast"` implements."""

    swap_axes: bool = False
    standard_scale: Literal["var", "obs"] | None = None
    var_group_rotation: float | None = None
    dendrogram: bool = False
    vmin: float | None = None
    vmax: float | None = None
    vcenter: float | None = None
    norm: Normalize | None = None
    log: bool = False

    @classmethod
    def pop(cls, kwargs: dict[str, Any]) -> _Options:
        kwargs.pop("show_gene_labels", None)  # there are too many windows to label
        opts = cls(**kwargs)
        kwargs.clear()
        return opts

    def transform(self, values: np.ndarray) -> np.ndarray:
        """`log` and `standard_scale` of the drawn rows, as :func:`scanpy.pl.heatmap` does for all cells."""
        if self.log:
            values = np.log1p(values)
        if self.standard_scale is not None:
            axis = 0 if self.standard_scale == "var" else 1
            values = values - np.nanmin(values, axis=axis, keepdims=True)
            with np.errstate(invalid="ignore", divide="ignore"):
                values = np.nan_to_num(values / np.nanmax(values, axis=axis, keepdims=True), nan=0)
        return values

    def color_norm(self, lo: float, hi: float) -> Normalize:
        """The user's `norm`, or a `TwoSlopeNorm` at `vcenter` (0) from `vmin` to `vmax` (default `lo` and `hi`)."""
        if self.norm is not None:
            return self.norm
        return TwoSlopeNorm(
            0 if self.vcenter is None else self.vcenter,
            vmin=lo if self.vmin is None else self.vmin,
            vmax=hi if self.vmax is None else self.vmax,
        )


# --- mode selection -------------------------------------------------------------------------------------------------


def _use_legacy(func: str, adata: AnnData, groupby: str, mode: str, kwargs: dict[str, Any]) -> bool:
    """Whether to draw with `mode="legacy"`: asked for, or needed for an argument `mode="fast"` does not implement."""
    if mode not in ("fast", "legacy"):
        raise ValueError(f"`mode` must be 'fast' or 'legacy', got {mode!r}.")
    if groupby not in adata.obs.columns:
        if groupby == "cnv_leiden":
            raise ValueError("'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?")
        raise ValueError(f"{groupby!r} is not in `adata.obs`.")
    if mode == "legacy":
        return True
    reasons = sorted(set(kwargs) - FAST_KWARGS)
    if not isinstance(kwargs.get("dendrogram", False), bool | np.bool_):
        reasons.append("dendrogram=<key>")
    if kwargs.get("standard_scale") not in (None, "var", "obs"):
        reasons.append(f"standard_scale={kwargs['standard_scale']!r}")
    col = _as_df(adata.obs)[groupby]
    if pd.api.types.is_numeric_dtype(col) and not isinstance(col.dtype, pd.CategoricalDtype):
        reasons.append("a numeric `groupby` column")
    if reasons:
        log.info("%s: mode='fast' does not support %s, using mode='legacy'", func, ", ".join(reasons))
        return True
    return False


def _legacy(
    fn: Callable[..., dict[str, Axes] | None],
    adata: AnnData,
    groupby: str,
    use_rep: str,
    cmap: str | Colormap,
    figsize: tuple[float, float],
    ax: Axes | None,
    show: bool | None,
    save: str | bool | None,
    return_fig: bool | None,
    kwargs: dict[str, Any],
) -> Figure | dict[str, Axes] | None:
    if ax is not None:
        raise ValueError("`ax` needs mode='fast' (and only arguments that it supports, see `mode`).")
    common = {"groupby": groupby, "use_rep": use_rep, "cmap": cmap, "figsize": figsize, "save": save}
    if not return_fig:
        return fn(adata, show=show, **common, **kwargs)
    axes = fn(adata, show=False, **common, **kwargs)
    assert axes is not None
    fig = axes["heatmap_ax"].get_figure(root=True)
    assert fig is not None
    return fig


def _legacy_heatmap(adata: AnnData, **kwargs) -> dict[str, Axes] | None:
    from infercnvpy.pl._chromosome_heatmap_legacy import legacy_heatmap

    return legacy_heatmap(adata, **kwargs)


def _legacy_summary(adata: AnnData, **kwargs) -> dict[str, Axes] | None:
    from infercnvpy.pl._chromosome_heatmap_legacy import legacy_summary

    return legacy_summary(adata, **kwargs)


# --- data -----------------------------------------------------------------------------------------------------------


def _figure(ax: Axes | None, figsize: tuple[float, float]) -> tuple[Figure, Axes]:
    if ax is None:
        ax = plt.figure(figsize=figsize).add_subplot()
    fig = ax.get_figure(root=True)
    assert fig is not None
    return fig, ax


def _n_pixel_rows(fig: Figure, ax: Axes, *, swap: bool) -> int:
    """Extent of `ax` along the cells, in pixels at the larger of the figure and savefig resolution."""
    dpi = fig.dpi
    if not isinstance(save_dpi := plt.rcParams["savefig.dpi"], str):
        dpi = max(dpi, save_dpi)
    pos, (width, height) = ax.get_position(), fig.get_size_inches()
    return max(1, int(np.ceil((pos.width * width if swap else pos.height * height) * dpi)))


def _categories(adata: AnnData, groupby: str) -> pd.Categorical:
    col = _as_df(adata.obs)[groupby].astype("category")  # keeps the order of categories that are already categorical
    return pd.Categorical(col, categories=col.cat.categories)


def _group_colors(adata: AnnData, groupby: str, n: int) -> list:
    """The colours scanpy stores for `groupby`, or the default palette scanpy would pick."""
    colors = adata.uns.get(f"{groupby}_colors")
    if colors is not None and len(colors) >= n:
        return list(colors[:n])
    if n <= len(cycle := plt.rcParams["axes.prop_cycle"].by_key()["color"]):
        return cycle[:n]
    for palette in (sc.pl.palettes.default_20, sc.pl.palettes.default_28, sc.pl.palettes.default_102):
        if n <= len(palette):
            return list(palette[:n])
    return ["grey"] * n


def _reduce_groups(sums: np.ndarray, counts: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Means of the non-empty groups from the sums and cell counts of their `rows` runs of cells."""
    starts = np.concatenate([[0], np.cumsum(rows[rows > 0])[:-1]])
    return np.add.reduceat(sums, starts, axis=0) / np.add.reduceat(counts, starts)[:, None]


def _group_means(X: Any, codes: np.ndarray, n_groups: int) -> np.ndarray:
    """Means of the non-empty groups in one pass, with several runs per group so that all threads take part."""
    order, bounds, rows = pixel_rows(codes, n_groups, max(n_groups, n_kernel_chunks(len(codes))))
    return _reduce_groups(binned_sums(X, order, bounds), np.diff(bounds), rows)


def _order_by_dendrogram(
    adata: AnnData, groupby: str, cat: pd.Categorical, palette: list, means: Callable[[], np.ndarray]
) -> tuple[pd.Categorical, list, dict[str, Any] | None]:
    """Reorder the categories (and their colours) as the leaves of the dendrogram of the non-empty groups.

    The dendrogram is `adata.uns[f"dendrogram_{groupby}"]` from :func:`scanpy.tl.dendrogram` if that has exactly
    these groups, as for `mode="legacy"`, otherwise it is computed from the groups' mean CNV profiles (`means()`).
    Unused categories go last. With fewer than two groups, there is no dendrogram.
    """
    present = np.flatnonzero(np.bincount(cat.codes[cat.codes >= 0], minlength=len(cat.categories)))
    if len(present) < 2:
        return cat, palette, None
    labels = [str(c) for c in cat.categories[present]]
    stored = adata.uns.get(f"dendrogram_{groupby}", {}).get("dendrogram_info")
    if stored is not None and sorted(map(str, stored["ivl"])) == sorted(labels):
        dendro = {**stored, "leaves": [labels.index(str(c)) for c in stored["ivl"]]}
    else:
        dendro = _dendrogram(means(), labels)
    new = [*present[dendro["leaves"]], *np.setdiff1d(np.arange(len(cat.categories)), present)]
    return cat.reorder_categories(cat.categories[new]), [palette[i] for i in new], dendro


def _dendrogram(means: np.ndarray, labels: list[str]) -> dict[str, Any]:
    """`scipy.cluster.hierarchy.dendrogram` of the groups, computed like :func:`scanpy.tl.dendrogram`."""
    import scipy.cluster.hierarchy as sch
    from scipy.spatial import distance

    corr = np.clip(pd.DataFrame(means.T).corr().to_numpy(), -1, 1)
    z = sch.linkage(distance.squareform(1 - corr, checks=False), method="complete")
    return sch.dendrogram(z, labels=labels, no_plot=True)


# --- drawing --------------------------------------------------------------------------------------------------------


def _draw(
    fig: Figure,
    ax: Axes,
    adata: AnnData,
    groupby: str,
    cat: pd.Categorical,
    palette: list,
    values: np.ndarray,
    rows: np.ndarray,
    dendro: dict[str, Any] | None,
    use_rep: str,
    norm: Normalize,
    cmap: str | Colormap,
    opts: _Options,
) -> dict[str, Axes]:
    """Draw `values` (pixel rows x windows, `rows` rows per category of `cat`) with its annotations.

    Cells run down the heatmap, or along it with `swap_axes`. `palette` has a colour per category, `dendro` the
    dendrogram of the non-empty groups in their order.
    """
    keep = np.flatnonzero(rows)
    labels = [str(c) for c in cat.categories[keep]]
    colors = [palette[g] for g in keep]
    rows = rows[keep]
    n_rows, n_cols = values.shape
    ends = np.cumsum(rows)
    centers = ends - rows / 2 - 0.5
    swap = opts.swap_axes

    chr_pos_dict = dict(sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    im = ax.imshow(values.T if swap else values, aspect="auto", interpolation="nearest", norm=norm, cmap=cmap)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    ax.grid(visible=False)
    line: dict[str, Any] = {"lw": 1, "color": "black", "zorder": 10, "clip_on": False}
    if swap:
        ax.set_xlim(-0.5, n_rows - 0.5)
        ax.set_ylim(n_cols - 0.5, -0.5)
        ax.vlines(ends[:-1] - 0.5, -0.5, n_cols - 0.5, **line)
        ax.hlines(chr_pos[1:], lw=0.6, xmin=-0.5, xmax=n_rows - 0.5)
    else:
        ax.set_ylim(n_rows - 0.5, -0.5)
        ax.set_xlim(-0.5, n_cols - 0.5)
        ax.hlines(ends[:-1] - 0.5, -0.5, n_cols - 0.5, **line)
        ax.vlines(chr_pos[1:], lw=0.6, ymin=-0.5, ymax=n_rows - 0.5)

    divider = make_axes_locatable(ax)
    axes = {"heatmap_ax": ax}
    along_cells = {"sharex": ax} if swap else {"sharey": ax}
    along_windows = {"sharey": ax} if swap else {"sharex": ax}

    chr_ax = axes["gene_groups_ax"] = divider.append_axes(
        "right" if swap else "top", size=Fixed(_CHR_H), pad=Fixed(_PAD_W if swap else _PAD_H), **along_windows
    )
    _chromosome_brackets(chr_ax, list(chr_pos_dict), chr_pos, n_cols, rotation=opts.var_group_rotation, right=swap)

    group_ax = axes["groupby_ax"] = divider.append_axes(
        "bottom" if swap else "left", size=Fixed(_GROUP_W), pad=Fixed(_PAD_H if swap else _PAD_W), **along_cells
    )
    blocks = np.repeat(np.arange(len(rows)), rows)[:, None]
    group_ax.imshow(
        blocks.T if swap else blocks,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(colors),
        norm=BoundaryNorm(np.arange(len(rows) + 1) - 0.5, len(rows)),
    )
    group_ax.tick_params(left=False, bottom=False, labelleft=not swap, labelbottom=swap, labelsize="small")
    if swap:
        group_ax.set_xticks(centers, labels, rotation=0 if max(len(x) for x in labels) < 3 else 90)
        group_ax.set_xlabel(groupby)
    else:
        group_ax.set_yticks(centers, labels)
        group_ax.set_ylabel(groupby)
    group_ax.set_frame_on(False)
    group_ax.grid(visible=False)

    if dendro is not None:
        dendro_ax = axes["dendrogram_ax"] = divider.append_axes(
            "top" if swap else "right", size=Fixed(_DENDRO_W), pad=Fixed(_PAD_H if swap else _PAD_W), **along_cells
        )
        leaf_pos = np.arange(5, 10 * len(rows), 10, dtype=float)  # where scipy puts the leaves
        for xs, ys in zip(dendro["icoord"], dendro["dcoord"], strict=True):
            pos = np.interp(xs, leaf_pos, centers)
            dendro_ax.plot(*((pos, ys) if swap else (ys, pos)), color="#555555")
        dendro_ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        dendro_ax.set_frame_on(False)
        dendro_ax.grid(visible=False)

    cbar_holder = divider.append_axes("right", size=Fixed(_CBAR_W), pad=Fixed(_PAD_W))
    cbar_holder.axis("off")
    cbar_ax = axes["colorbar_ax"] = cbar_holder.inset_axes((0, 0, 1, min(1, _CBAR_MAX_H / fig.get_size_inches()[1])))
    fig.colorbar(im, cax=cbar_ax)
    return axes


def _chromosome_brackets(
    ax: Axes, labels: list[str], starts: list[int], n_cols: int, *, rotation: float | None, right: bool
) -> None:
    """Brackets over (or, with `right`, next to) the windows of each chromosome with its label, as scanpy draws."""
    verts, codes = [], []
    if rotation is None:
        rotation = 90 if max(len(x) for x in labels) > 4 else 0
    for label, start, end in zip(labels, starts, [*starts[1:], n_cols], strict=True):
        lo, hi = start - 0.3, end + 0.3
        if right:
            verts += [(0, lo), (0.15, lo), (0.15, hi), (0, hi)]
            text = label if (hi - lo) * 2 >= len(label) else f"{label[: int((hi - lo) * 2)]}."
            ax.text(0.6, (lo + hi) / 2, text, ha="right", va="center", rotation=270, fontsize="small")
        else:
            verts += [(lo, 0), (lo, 0.6), (hi, 0.6), (hi, 0)]
            ax.text((lo + hi) / 2, 1.1, label, ha="center", va="bottom", rotation=rotation)
        codes += [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO]
    ax.add_patch(patches.PathPatch(MplPath(verts, codes), facecolor="none", lw=1.5))
    if right:
        ax.set_xlim(0, 1)
    else:
        ax.set_ylim(0, 1)
    ax.axis("off")


def _finish(
    fig: Figure, axes: dict[str, Axes], *, show: bool | None, save: str | bool | None, return_fig: bool | None
) -> Figure | dict[str, Axes] | None:
    if sc.settings.autosave if save is None else save:
        _save(fig, save)
    if return_fig:
        return fig
    if sc.settings.autoshow if show is None else show:
        plt.show()
        return None
    return axes


def _save(fig: Figure, save: str | bool | None) -> None:
    """Save as `<figdir>/heatmap<save><plot_suffix>.<ext>`, the file name `scanpy.pl.heatmap` used."""
    name, ext = "heatmap", sc.settings.file_format_figs
    if isinstance(save, str):
        if (suffix := Path(save).suffix) in {".svg", ".pdf", ".png"}:
            save, ext = save.removesuffix(suffix), suffix[1:]
        name += save
    figdir = Path(sc.settings.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / f"{name}{sc.settings.plot_suffix}.{ext}", bbox_inches="tight")
