from __future__ import annotations

import re
from collections.abc import Sequence

import numpy as np
import scanpy as sc
from anndata import AnnData

from infercnvpy._util import _as_df, _as_matrix

_REGION = re.compile(r"^(?P<chrom>[^:]+)(?::(?P<start>[\d_,]+)-(?P<end>[\d_,]+))?$")


def _region_columns(adata: AnnData, use_rep: str, region: str) -> np.ndarray:
    """Columns of `X_{use_rep}` in `region` ("chr3", or "chr3:90000000-200000000" by window overlap)."""
    if (m := _REGION.match(region)) is None:
        raise ValueError(f'`region` must look like "chr3" or "chr3:90000000-200000000", got {region!r}.')
    uns = adata.uns[use_rep]
    chr_pos = dict(sorted(uns["chr_pos"].items(), key=lambda x: x[1]))
    if (chrom := m["chrom"]) not in chr_pos:
        raise ValueError(f"Chromosome {chrom!r} is not in `adata.uns[{use_rep!r}]['chr_pos']`.")
    bounds = [*chr_pos.values(), adata.obsm[f"X_{use_rep}"].shape[1]]
    i = list(chr_pos).index(chrom)
    cols = np.arange(bounds[i], bounds[i + 1])
    if m["start"] is not None:
        if "window_start" not in uns:
            raise ValueError(
                f"`adata.uns[{use_rep!r}]` has no window positions. Rerun `tl.infercnv` to select a genomic range."
            )
        start, end = (int(m[k].replace("_", "").replace(",", "")) for k in ("start", "end"))
        overlap = (np.asarray(uns["window_end"])[cols] >= start) & (np.asarray(uns["window_start"])[cols] <= end)
        cols = cols[overlap]
    if not len(cols):
        raise ValueError(f"No windows of `X_{use_rep}` in {region!r}.")
    return cols


def spatial_cnv(
    adata: AnnData,
    *,
    region: str | Sequence[str] | None = None,
    color: str | Sequence[str] | None = None,
    use_rep: str = "cnv",
    basis: str = "spatial",
    cmap: str = "bwr",
    **kwargs,
):
    """Plot the mean CNV of genomic regions, or columns of `adata.obs`, on the tissue.

    Thin wrapper around :func:`scanpy.pl.embedding` with `basis="spatial"`. Each region
    becomes a panel coloured by the mean of `adata.obsm["X_{use_rep}"]` over its windows,
    with a diverging colour map centred at 0 (pass `vcenter=None` to turn that off).

    Parameters
    ----------
    adata
        Annotated data matrix with coordinates in `adata.obsm["spatial"]`.
    region
        A chromosome (`"chr3"`) or a genomic range (`"chr3:90_000_000-200_000_000"`, all
        windows that overlap it), or a list of them. Ranges need the window positions that
        :func:`infercnvpy.tl.infercnv` stores in `adata.uns[use_rep]`.
    color
        Columns of `adata.obs` to plot as additional panels, e.g. labels from
        :func:`infercnvpy.tl.spatial_regions`.
    use_rep
        Key of the CNV matrix, `adata.obsm["X_{use_rep}"]`, e.g. `"cnv_smoothed"` from
        :func:`infercnvpy.tl.smooth_cnv`.
    basis
        Key of the coordinates in `adata.obsm` (`"spatial"` refers to `obsm["spatial"]`).
    cmap
        Colour map of the regions (and other continuous panels).
    **kwargs
        Passed to :func:`scanpy.pl.embedding`, e.g. `size`, `ncols`, `vmin`/`vmax` or `show`.

    Returns
    -------
    Whatever :func:`scanpy.pl.embedding` returns: `None`, or the axes if `show=False`.
    """
    regions = [region] if isinstance(region, str) else list(region or [])
    colors = [color] if isinstance(color, str) else list(color or [])
    if not regions and not colors:
        raise ValueError("Pass at least one `region` or `color`.")
    obs = _as_df(adata.obs).copy()
    if regions:
        X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`obsm['X_{use_rep}']`")
        for r in regions:
            obs[f"{use_rep} {r}"] = np.asarray(X[:, _region_columns(adata, use_rep, r)].mean(axis=1)).ravel()
    tmp = AnnData(obs=obs, uns=adata.uns)
    tmp.obsm[basis] = np.asarray(adata.obsm[basis])
    if regions:
        kwargs.setdefault("vcenter", 0)
    return sc.pl.embedding(tmp, basis, color=[f"{use_rep} {r}" for r in regions] + colors, cmap=cmap, **kwargs)
