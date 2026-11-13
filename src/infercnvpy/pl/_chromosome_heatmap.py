from __future__ import annotations

import re
import textwrap
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm, Colormap, ListedColormap, TwoSlopeNorm
from matplotlib.path import Path as MplPath
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_size import Fixed

from infercnvpy._numba import n_kernel_chunks, numba_threads
from infercnvpy._util import _as_df, _as_matrix
from infercnvpy.pl._heatmap_kernels import _pixel_rows, binned_sums

if TYPE_CHECKING:
    from collections.abc import Callable

    from anndata import AnnData
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

F = TypeVar("F", bound="Callable[..., Any]")

# Sizes in inches, as `scanpy.pl.heatmap` lays them out at figsize=(16, 10).
_GROUP_W = 0.155  # colour blocks of the groups
_CBAR_W = 0.155
_DENDRO_W = 0.775
_CHR_H = 0.115  # chromosome brackets
_PAD_W = 0.03
_PAD_H = 0.05
_CBAR_MAX_H = 4.0

_PARAMS_DOC = """\
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
        (width, height) tuple in inches. Ignored if `ax` is given.
    ax
        Draw the heatmap into this axes. The group colours, chromosome labels, dendrogram and colour bar go into
        axes split off from it.
    n_jobs
        Number of threads for averaging the cells. `None` uses all cores.
    show
        Whether to show the figure. Defaults to :attr:`scanpy.settings.autoshow`.
    save
        If `True` or a `str`, save the figure to :attr:`scanpy.settings.figdir`. A string is appended to the default
        filename `heatmap`; it sets the file type if it ends in `.pdf`, `.png` or `.svg`, otherwise
        :attr:`scanpy.settings.file_format_figs` does.
    return_fig
        Return the :class:`~matplotlib.figure.Figure`.
    **kwargs
        `vmin` and `vmax` set the limits of the colour scale (centred at 0), `dendrogram=True` orders the groups by
        hierarchical clustering of their mean CNV profiles (Pearson correlation, complete linkage, as
        :func:`scanpy.tl.dendrogram`) and draws the dendrogram on the right. Other arguments are passed on to
        :meth:`matplotlib.axes.Axes.imshow`.

    Returns
    -------
    The figure if `return_fig` is true. Otherwise, if `show` is false, a dictionary of axes with the keys
    `heatmap_ax`, `groupby_ax` (group colours), `gene_groups_ax` (chromosome labels), `colorbar_ax` and, with
    `dendrogram=True`, `dendrogram_ax`.
"""


def _with_params_doc(fn: F) -> F:
    """Fill in `{params}`, indented like that line (Python >= 3.13 strips the indentation of docstrings)."""
    if fn.__doc__ and (m := re.search(r"^( *)\{params\}$", fn.__doc__, re.MULTILINE)):  # no doc with `python -OO`
        params = textwrap.indent(textwrap.dedent(_PARAMS_DOC).strip(), m[1])
        fn.__doc__ = fn.__doc__.replace(m[0], params)
    return fn


@_with_params_doc
def chromosome_heatmap(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: str | Colormap = "bwr",
    figsize: tuple[float, float] = (16, 10),
    ax: Axes | None = None,
    n_jobs: int | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    return_fig: bool | None = None,
    **kwargs,
) -> Figure | dict[str, Axes] | None:
    """Plot a heatmap of smoothed gene expression by chromosome.

    Cells are sorted by `groupby` and each group gets a share of the heatmap's height proportional to its size.
    Each pixel row shows the mean of the cells that fall into it.

    .. versionchanged:: 0.7.0
       Rendered by infercnvpy instead of :func:`scanpy.pl.heatmap`; cells are averaged into pixel rows, so large
       data plots in seconds.

    Parameters
    ----------
    {params}
    """
    fig, ax = _figure(ax, figsize)
    cat = _categories(adata, groupby)
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`")
    order, bounds, rows = _pixel_rows(cat.codes.astype(np.int64), len(cat.categories), _n_pixel_rows(fig, ax))
    with numba_threads(n_jobs):
        sums, lo, hi = binned_sums(X, order, bounds)
    values = sums / np.diff(bounds)[:, None]
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    norm = TwoSlopeNorm(0, vmin=lo if vmin is None else vmin, vmax=hi if vmax is None else vmax)
    axes = _draw(
        fig, ax, adata, groupby, cat, values, rows, _group_means(sums, bounds, rows), use_rep, norm, cmap, kwargs
    )
    return _finish(fig, axes, show=show, save=save, return_fig=return_fig)


@_with_params_doc
def chromosome_heatmap_summary(
    adata: AnnData,
    *,
    groupby: str = "cnv_leiden",
    use_rep: str = "cnv",
    cmap: str | Colormap = "bwr",
    figsize: tuple[float, float] = (16, 10),
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
       Rendered by infercnvpy instead of :func:`scanpy.pl.heatmap`.

    Parameters
    ----------
    {params}
    """
    fig, ax = _figure(ax, figsize)
    cat = _categories(adata, groupby)
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`adata.obsm['X_{use_rep}']`")
    codes = cat.codes.astype(np.int64)
    # more runs than groups, so that big groups are summed by several threads
    order, bounds, rows = _pixel_rows(codes, len(cat.categories), max(len(cat.categories), n_kernel_chunks(len(codes))))
    with numba_threads(n_jobs):
        sums, _, _ = binned_sums(X, order, bounds)
    means = _group_means(sums, bounds, rows)
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    norm = TwoSlopeNorm(0, vmin=means.min() if vmin is None else vmin, vmax=means.max() if vmax is None else vmax)
    axes = _draw(fig, ax, adata, groupby, cat, means, (rows > 0).astype(np.int64), means, use_rep, norm, cmap, kwargs)
    return _finish(fig, axes, show=show, save=save, return_fig=return_fig)


def _figure(ax: Axes | None, figsize: tuple[float, float]) -> tuple[Figure, Axes]:
    if ax is None:
        ax = plt.figure(figsize=figsize).add_subplot()
    fig = ax.get_figure(root=True)
    assert fig is not None
    return fig, ax


def _n_pixel_rows(fig: Figure, ax: Axes) -> int:
    """Height of `ax` in pixels at the larger of the figure and savefig resolution."""
    dpi = fig.dpi
    if not isinstance(save_dpi := plt.rcParams["savefig.dpi"], str):
        dpi = max(dpi, save_dpi)
    return max(1, int(np.ceil(ax.get_position().height * fig.get_size_inches()[1] * dpi)))


def _categories(adata: AnnData, groupby: str) -> pd.Categorical:
    if groupby not in adata.obs.columns:
        if groupby == "cnv_leiden":
            raise ValueError("'cnv_leiden' is not in `adata.obs`. Did you run `tl.leiden()`?")
        raise ValueError(f"{groupby!r} is not in `adata.obs`.")
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


def _group_means(sums: np.ndarray, bounds: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Means of the non-empty groups from the sums of their runs of cells."""
    starts = np.concatenate([[0], np.cumsum(rows[rows > 0])[:-1]])
    return np.add.reduceat(sums, starts, axis=0) / np.add.reduceat(np.diff(bounds), starts)[:, None]


def _dendrogram(means: np.ndarray, labels: list[str]) -> dict[str, Any]:
    """`scipy.cluster.hierarchy.dendrogram` of the groups, computed like :func:`scanpy.tl.dendrogram`."""
    import scipy.cluster.hierarchy as sch
    from scipy.spatial import distance

    corr = np.clip(pd.DataFrame(means.T).corr().to_numpy(), -1, 1)
    z = sch.linkage(distance.squareform(1 - corr, checks=False), method="complete")
    return sch.dendrogram(z, labels=labels, no_plot=True)


def _draw(
    fig: Figure,
    ax: Axes,
    adata: AnnData,
    groupby: str,
    cat: pd.Categorical,
    values: np.ndarray,
    rows: np.ndarray,
    group_means: np.ndarray,
    use_rep: str,
    norm: TwoSlopeNorm,
    cmap: str | Colormap,
    kwargs: dict[str, Any],
) -> dict[str, Axes]:
    """Draw `values` (pixel rows x windows, grouped by `rows` rows per category of `cat`) with its annotations."""
    keep = np.flatnonzero(rows)
    labels = [str(c) for c in cat.categories[keep]]
    palette = _group_colors(adata, groupby, len(cat.categories))
    colors = [palette[g] for g in keep]
    rows = rows[keep]
    dendro = None
    if kwargs.pop("dendrogram", False) and len(keep) > 1:
        dendro = _dendrogram(group_means, labels)
        leaves = dendro["leaves"]
        ends = np.cumsum(rows)
        values = np.concatenate([values[ends[g] - rows[g] : ends[g]] for g in leaves])
        labels, colors, rows = [labels[g] for g in leaves], [colors[g] for g in leaves], rows[leaves]
    n_rows, n_cols = values.shape

    chr_pos_dict = dict(sorted(adata.uns[use_rep]["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    kwargs.setdefault("interpolation", "nearest")
    im = ax.imshow(values, aspect="auto", norm=norm, cmap=cmap, **kwargs)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    ax.grid(visible=False)
    ends = np.cumsum(rows)
    ax.hlines(ends[:-1] - 0.5, -0.5, n_cols - 0.5, lw=1, color="black", zorder=10, clip_on=False)
    ax.vlines(chr_pos[1:], lw=0.6, ymin=-0.5, ymax=n_rows - 0.5)

    divider = make_axes_locatable(ax)
    axes = {"heatmap_ax": ax}

    chr_ax = axes["gene_groups_ax"] = divider.append_axes("top", size=Fixed(_CHR_H), pad=Fixed(_PAD_H), sharex=ax)
    _chromosome_brackets(chr_ax, list(chr_pos_dict), chr_pos, n_cols)

    group_ax = axes["groupby_ax"] = divider.append_axes("left", size=Fixed(_GROUP_W), pad=Fixed(_PAD_W), sharey=ax)
    group_ax.imshow(
        np.repeat(np.arange(len(rows)), rows)[:, None],
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(colors),
        norm=BoundaryNorm(np.arange(len(rows) + 1) - 0.5, len(rows)),
    )
    group_ax.set_yticks(ends - rows / 2 - 0.5, labels)
    group_ax.tick_params(axis="y", left=False, labelsize="small")
    group_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    group_ax.set_ylabel(groupby)
    group_ax.set_frame_on(False)
    group_ax.grid(visible=False)

    if dendro is not None:
        dendro_ax = axes["dendrogram_ax"] = divider.append_axes(
            "right", size=Fixed(_DENDRO_W), pad=Fixed(_PAD_W), sharey=ax
        )
        leaf_x = np.arange(5, 10 * len(rows), 10, dtype=float)  # where scipy puts the leaves
        for xs, ys in zip(dendro["icoord"], dendro["dcoord"], strict=True):
            dendro_ax.plot(ys, np.interp(xs, leaf_x, ends - rows / 2 - 0.5), color="#555555")
        dendro_ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        dendro_ax.set_frame_on(False)
        dendro_ax.grid(visible=False)

    cbar_holder = divider.append_axes("right", size=Fixed(_CBAR_W), pad=Fixed(_PAD_W))
    cbar_holder.axis("off")
    cbar_ax = axes["colorbar_ax"] = cbar_holder.inset_axes((0, 0, 1, min(1, _CBAR_MAX_H / fig.get_size_inches()[1])))
    fig.colorbar(im, cax=cbar_ax)
    return axes


def _chromosome_brackets(ax: Axes, labels: list[str], starts: list[int], n_cols: int) -> None:
    """Brackets over the windows of each chromosome with its label on top (as `scanpy.pl.heatmap`)."""
    verts, codes = [], []
    rotation = 90 if max(len(x) for x in labels) > 4 else 0
    for label, start, end in zip(labels, starts, [*starts[1:], n_cols], strict=True):
        left, right = start - 0.3, end + 0.3
        verts += [(left, 0), (left, 0.6), (right, 0.6), (right, 0)]
        codes += [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO]
        ax.text((left + right) / 2, 1.1, label, ha="center", va="bottom", rotation=rotation)
    ax.add_patch(patches.PathPatch(MplPath(verts, codes), facecolor="none", lw=1.5))
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
