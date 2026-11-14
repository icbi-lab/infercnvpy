"""Scores to summarize and assess copy number variation."""

import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from infercnvpy._logging import log, log_result
from infercnvpy._numba import numba_threads
from infercnvpy._util import _as_df, _as_matrix, _choose_mtx_rep, _is_backed, _is_dask, _params_for_uns

from ._ith import corr_iqr


def cnv_score(
    adata: AnnData,
    groupby: str = "cnv_leiden",
    *,
    use_rep: str = "cnv",
    key_added: str = "cnv_score",
    inplace: bool = True,
    obs_key=None,
) -> Mapping[Any, np.number] | None:
    """Assign each cnv cluster a CNV score.

    Clusters with a high score are likely affected by copy number aberrations.
    Based on this score, cells can be divided into tumor/normal cells.

    The score is currently simply defined as the mean of the absolute values of the
    CNV scores in each cluster.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        Key under which the clustering is stored in adata.obs. Usually
        the result of :func:`infercnvpy.tl.leiden`, but could also be
        other grouping information, e.g. sample or patient information.
    use_rep
        Key under which the result of :func:`infercnvpy.tl.infercnv` is stored
        in adata. `adata.obsm[f"X_{use_rep}"]` can also be a dask array or a backed
        matrix, which is read block by block.
    key_added
        Key under which the score will be stored in `adata.obs`.
    inplace
        If True, store the result in adata, otherwise return it.
    obs_key
        Deprecated alias for `groupby`.

    Returns
    -------
    If `inplace=True`, returns `None` and sets the fields below.
    Otherwise returns a dictionary with the score per group.
    Sets the following fields:

    `adata.obs[key_added]` : float
        Score of the group each cell belongs to.
    `adata.uns[key_added]["params"]` : :class:`dict`
        Parameters of this call (`groupby`, `use_rep`).
    """
    if obs_key is not None:
        warnings.warn(
            "The obs_key argument has been renamed to `groupby` for consistency with "
            "other functions and will be removed in the future. ",
            category=FutureWarning,
            stacklevel=2,
        )
        groupby = obs_key

    if groupby not in adata.obs.columns and groupby == "cnv_leiden":
        raise ValueError("`cnv_leiden` not found in `adata.obs`. Did you run `tl.leiden`?")
    X = _as_matrix(adata.obsm[f"X_{use_rep}"], f"`obsm['X_{use_rep}']`", lazy=True)
    groups = _as_df(adata.obs)[groupby]
    if _is_dask(X) or _is_backed(X):
        cluster_score = _abs_group_means(X, groups)
    else:
        cluster_score = {cluster: np.mean(np.abs(X[groups.values == cluster, :])) for cluster in groups.unique()}
    log_result("cnv_score", "scored %d groups of %r", len(cluster_score), groupby)

    if inplace:
        score_array = np.array([cluster_score[c] for c in adata.obs[groupby]])
        adata.obs[key_added] = score_array
        adata.uns[key_added] = {"params": _params_for_uns(groupby=groupby, use_rep=use_rep)}
        return None
    return cluster_score


def _abs_group_means(X, groups: pd.Series) -> dict[Any, np.number]:
    """Mean absolute value of the rows of the lazy matrix `X` in each group, in one pass over its row blocks."""
    from infercnvpy.tl._infercnv import _Rows

    uniq = groups.unique()
    n = len(uniq)
    codes = pd.Index(uniq).get_indexer(groups)
    codes[pd.isna(groups.to_numpy())] = n  # cells without a group (a NaN group gets no cells, as in memory)
    sums = np.zeros(n + 1)
    source = _Rows(X, 10_000)
    blocks = list(range(source.n_blocks))
    for i, (Xb, row0) in zip(blocks, source.read_ahead(blocks), strict=True):
        s, e = int(source.bounds[i]), int(source.bounds[i + 1])
        row_sums = np.asarray(abs(Xb[row0 : row0 + e - s]).sum(axis=1), dtype=np.float64).ravel()
        sums += np.bincount(codes[s:e], weights=row_sums, minlength=n + 1)
    counts = np.bincount(codes, minlength=n + 1)[:n] * X.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums[:n] / counts
    return dict(zip(uniq, means, strict=True))


def ithgex(
    adata: AnnData,
    groupby: str,
    *,
    use_raw: bool | None = None,
    layer: str | None = None,
    inplace: bool = True,
    key_added: str = "ithgex",
    n_jobs: int | None = None,
    n_pairs: int | None = None,
    random_state: int | None = 0,
    exclude_self_pairs: bool = False,
) -> Mapping[str, float] | None:
    """Compute the ITHGEX diversity score based on gene expression :cite:`Wu2021`.

    A high score indicates a high diversity of gene expression profiles of cells
    within a group.

    The score is defined as follows:

        Intratumoral heterogeneity scores based on CNAs and gene expressions
        The calculations of intratumoral heterogeneity scores were inspired by a
        previous study and modified as follows. First, to calculate ITHCNA, we used
        the relative expression value matrix generated by inferCNV and calculated the
        pairwise cell–cell distances using Pearson's correlation coefficients for each
        patient. ITHCNA was defined as interquartile range (IQR) of the distribution for
        all malignant cell pairs' Pearson's correlation coefficients. **Similarly, we also
        used gene expression profiles of cancer cells of each patient to construct the
        distribution of the intratumoral distances. ITHGEX was assigned as the IQR of the
        distribution.**

        (from :cite:`Wu2021`)

    The correlations are taken over all ordered pairs of cells in a group, including each
    cell with itself (self-pairs, see `exclude_self_pairs`), as in earlier versions. Cells
    with a constant profile, e.g. all-zero rows, have no defined correlation and are left out.

    .. versionchanged:: 0.7.0
        The correlations are computed block by block with a numba kernel instead of as one
        cells x cells matrix, so large groups no longer run out of memory. The quantiles match
        numpy to about 1e-6. Cells with a constant profile (e.g. no CNV signal at all) have no
        defined correlation and are left out; before, any such cell made the group's score `NaN`.
        The number of cells left out is logged per group.

    .. versionadded:: 0.7.0
        `n_pairs` and `exclude_self_pairs`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        calculate diversity for each group defined in this category.
    use_raw
        If `True`, use gene expression from `adata.raw.X`. Default: use `adata.X`
        (or `layer`, if given), even when `adata.raw` is present.
    layer
        Use gene expression from `adata.layers[layer]`
    inplace
        If True, store the result in adata, otherwise return it.
    key_added
        Key under which the score will be stored in `adata.obs`.
    n_jobs
        Number of threads. Default: all cores (numba's setting, see :func:`numba.set_num_threads`).
    n_pairs
        Estimate the score of large groups from this many randomly drawn cell pairs instead of
        all pairs, e.g. `1_000_000`. Pairs are drawn uniformly from the full correlation matrix,
        so the estimate is consistent (it converges to the exact score as `n_pairs` grows); its
        error shrinks like `1 / sqrt(n_pairs)` and the cost no longer grows with the group size.
        Groups with fewer pairs are computed exactly. Default: all pairs (exact).
    random_state
        Seed for drawing the pairs with `n_pairs`.
    exclude_self_pairs
        Leave out the correlation of each cell with itself. These self-pairs are always 1;
        including them (the default, as in earlier versions) raises the upper quartile for small
        groups, while for large groups the difference is negligible.

    Returns
    -------
    If `inplace=True`, returns `None` and sets the field below.
    Otherwise returns a dictionary with one value per group. Groups with a
    single cell have no cell pairs to correlate and get a score of `NaN`.
    Sets the following fields:

    `adata.obs[key_added]` : float
        Score of the group each cell belongs to, `NaN` for single-cell groups.
    """
    X = _as_matrix(_choose_mtx_rep(adata, use_raw, layer))
    return _ith(adata, X, groupby, key_added, inplace, n_jobs, n_pairs, random_state, exclude_self_pairs)


def ithcna(
    adata: AnnData,
    groupby: str,
    *,
    use_rep: str = "X_cnv",
    key_added: str = "ithcna",
    inplace: bool = True,
    n_jobs: int | None = None,
    n_pairs: int | None = None,
    random_state: int | None = 0,
    exclude_self_pairs: bool = False,
) -> Mapping[str, float] | None:
    """Compute the ITHCNA diversity score based on copy number variation :cite:`Wu2021`.

    A high score indicates a high diversity of CNV profiles of cells
    within a group.

    The score is defined as follows:

        Intratumoral heterogeneity scores based on CNAs and gene expressions
        The calculations of intratumoral heterogeneity scores were inspired by a
        previous study and modified as follows. First, to calculate ITHCNA, we used
        the relative expression value matrix generated by inferCNV and calculated the
        pairwise cell–cell distances using Pearson's correlation coefficients for each
        patient. ITHCNA was defined as interquartile range (IQR) of the distribution for
        all malignant cell pairs' Pearson's correlation coefficients.

        (from :cite:`Wu2021`)

    The correlations are taken over all ordered pairs of cells in a group, including each
    cell with itself (self-pairs, see `exclude_self_pairs`), as in earlier versions. Cells
    with a constant profile, e.g. the all-zero rows of a CNV matrix where every value lies
    in the noise band of :func:`infercnvpy.tl.infercnv`, have no defined correlation and are
    left out.

    .. versionchanged:: 0.7.0
        The correlations are computed block by block with a numba kernel instead of as one
        cells x cells matrix, so large groups no longer run out of memory. The quantiles match
        numpy to about 1e-6. Cells with a constant profile (e.g. no CNV signal at all) have no
        defined correlation and are left out; before, any such cell made the group's score `NaN`.
        The number of cells left out is logged per group.

    .. versionadded:: 0.7.0
        `n_pairs` and `exclude_self_pairs`.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        calculate diversity for each group defined in this category.
    use_rep
        Key under which the result of :func:`infercnvpy.tl.infercnv` is stored
        in adata.
    key_added
        Key under which the score will be stored in `adata.obs`.
    inplace
        If True, store the result in adata, otherwise return it.
    n_jobs
        Number of threads. Default: all cores (numba's setting, see :func:`numba.set_num_threads`).
    n_pairs
        Estimate the score of large groups from this many randomly drawn cell pairs instead of
        all pairs, e.g. `1_000_000`. Pairs are drawn uniformly from the full correlation matrix,
        so the estimate is consistent (it converges to the exact score as `n_pairs` grows); its
        error shrinks like `1 / sqrt(n_pairs)` and the cost no longer grows with the group size.
        Groups with fewer pairs are computed exactly. Default: all pairs (exact).
    random_state
        Seed for drawing the pairs with `n_pairs`.
    exclude_self_pairs
        Leave out the correlation of each cell with itself. These self-pairs are always 1;
        including them (the default, as in earlier versions) raises the upper quartile for small
        groups, while for large groups the difference is negligible.

    Returns
    -------
    If `inplace=True`, returns `None` and sets the field below.
    Otherwise returns a dictionary with one value per group. Groups with a
    single cell have no cell pairs to correlate and get a score of `NaN`.
    Sets the following fields:

    `adata.obs[key_added]` : float
        Score of the group each cell belongs to, `NaN` for single-cell groups.
    """
    X = _as_matrix(adata.obsm[use_rep], f"`obsm[{use_rep!r}]`")
    return _ith(adata, X, groupby, key_added, inplace, n_jobs, n_pairs, random_state, exclude_self_pairs)


def _ith(
    adata: AnnData,
    X,
    groupby: str,
    key_added: str,
    inplace: bool,
    n_jobs: int | None,
    n_pairs: int | None,
    random_state: int | None,
    exclude_self_pairs: bool,
):
    """IQR of the cell-cell correlations of the rows of `X` within each group of `adata.obs[groupby]`."""
    codes = adata.obs[groupby].astype("category").cat
    with numba_threads(n_jobs):
        rng = np.random.default_rng(random_state)
        sizes = np.bincount(codes.codes[codes.codes >= 0], minlength=len(codes.categories))
        if n_pairs is not None and (
            sampled := [str(g) for g, n in zip(codes.categories, sizes, strict=True) if n * n > n_pairs]
        ):
            log.info("%s: estimated from %d sampled cell pairs for groups %s", key_added, n_pairs, ", ".join(sampled))
        scores, heavy = {}, []
        for i, group in enumerate(codes.categories):
            scores[group], n_dropped = corr_iqr(
                X[np.flatnonzero(codes.codes == i)], n_pairs=n_pairs, rng=rng, exclude_self_pairs=exclude_self_pairs
            )
            if n_dropped:
                log.info(
                    "%s: left out %d of %d cells of group %s with a constant profile (no defined correlation)",
                    *(key_added, n_dropped, sizes[i], group),
                )
                if n_dropped > 0.2 * sizes[i]:
                    heavy.append(f"{group} ({n_dropped} of {sizes[i]})")
    if heavy:
        log.warning(
            "%s: more than 20%% of the cells of groups %s have a constant profile and were left out",
            *(key_added, ", ".join(heavy)),
            deep="in a noise-filtered CNV matrix these are cells without any CNV signal",
        )
    if not inplace:
        return scores
    adata.obs[key_added] = np.array([*scores.values(), np.nan])[codes.codes]  # code -1 (no group) -> NaN
    return None
