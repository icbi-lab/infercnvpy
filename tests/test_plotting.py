import logging

import dask
import dask.array as da
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

import infercnvpy as cnv
from infercnvpy.pl._chromosome_heatmap import _group_means
from infercnvpy.pl._heatmap_kernels import _rows_per_group, binned_sums, pixel_rows, sampled_rows

AX_KEYS = {"heatmap_ax", "groupby_ax", "gene_groups_ax", "colorbar_ax"}
LOGGER = "scverse.infercnvpy"
both_plots = pytest.mark.parametrize("plot", [cnv.pl.chromosome_heatmap, cnv.pl.chromosome_heatmap_summary])


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def adata_small() -> AnnData:
    """60 cells in uneven groups (one category unused, one cell without a group) x 12 windows on 3 chromosomes."""
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 12)).astype(np.float32)
    X[np.abs(X) < 0.8] = 0
    groups = rng.choice(["b", "a", "c"], size=60, p=[0.6, 0.3, 0.1]).astype(object)
    groups[7] = None
    obs = pd.DataFrame({"group": pd.Categorical(groups, categories=["c", "b", "unused", "a"])})
    obs.index = obs.index.astype(str)
    return AnnData(
        obs=obs,
        obsm={"X_cnv": sp.csr_matrix(X)},
        uns={"cnv": {"chr_pos": {"chr2": 5, "chr1": 0, "chr3": 9}}, "group_colors": ["red", "green", "blue", "black"]},
    )


def _image(axes) -> np.ndarray:
    return np.asarray(axes["heatmap_ax"].images[0].get_array())


def _dask(X, block_rows=7):
    """A dask array with row blocks of `block_rows` (CSR blocks for sparse `X`)."""
    if not sp.issparse(X):
        return da.from_array(X, chunks=(block_rows, -1))
    X = sp.csr_matrix(X)
    meta = sp.csr_matrix((0, X.shape[1]), dtype=X.dtype)
    blocks = [
        da.from_delayed(
            dask.delayed(X[s : s + block_rows]), shape=X[s : s + block_rows].shape, dtype=X.dtype, meta=meta
        )
        for s in range(0, X.shape[0], block_rows)
    ]
    return da.concatenate(blocks, axis=0)


# --- data ---------------------------------------------------------------------------------------------------------


def test_rows_per_group():
    counts = np.array([500, 3, 0, 97, 1])
    np.testing.assert_array_equal(_rows_per_group(counts, 100), [83, 1, 0, 16, 1])  # proportional, >= 1 if not empty
    np.testing.assert_array_equal(_rows_per_group(counts, counts.sum()), counts)  # one row per cell


@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix, sp.csc_matrix, _dask, lambda X: _dask(sp.csr_matrix(X))])
@pytest.mark.parametrize("n_rows", [7, 20, 1000])
def test_pixel_rows_match_numpy(fmt, n_rows):
    rng = np.random.default_rng(1)
    X = rng.normal(size=(120, 9))
    X[np.abs(X) < 1] = 0
    codes = rng.choice([0, 1, 3, -1], size=120, p=[0.7, 0.2, 0.05, 0.05])  # group 2 empty, -1 no group
    order, bounds, rows = pixel_rows(codes, 4, n_rows)

    np.testing.assert_array_equal(np.sort(order), np.flatnonzero(codes >= 0))
    assert np.all(np.diff(codes[order]) >= 0)  # sorted by group
    assert np.all(np.diff(bounds) > 0)  # no empty pixel row
    assert rows[2] == 0
    assert np.all(rows[[0, 1, 3]] >= 1)
    assert rows.sum() == len(bounds) - 1
    starts = np.concatenate([[0], np.cumsum(np.bincount(codes[order], minlength=4))])
    assert set(bounds[:-1]) >= set(starts[:-1][rows > 0])  # each pixel row holds cells of one group

    runs = list(zip(bounds[:-1], bounds[1:], strict=True))
    sums = binned_sums(fmt(X), order, bounds)
    np.testing.assert_allclose(sums, np.stack([X[order[a:b]].sum(axis=0) for a, b in runs]))

    picks = order[(bounds[:-1] + bounds[1:]) // 2]
    picked, lo, hi = sampled_rows(fmt(X), picks)
    np.testing.assert_array_equal(picked, X[picks])
    stored = X if fmt is np.asarray or fmt is _dask else X[X != 0]  # sparse: stored values only, as before 0.7.0
    assert (lo, hi) == (stored.min(), stored.max())

    expected = pd.DataFrame(X[codes >= 0]).groupby(codes[codes >= 0]).mean().to_numpy()
    np.testing.assert_allclose(_group_means(fmt(X), codes, 4), expected)


# --- fast vs legacy -----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix])
@pytest.mark.parametrize("kwargs", [{}, {"standard_scale": "var"}, {"standard_scale": "obs", "log": True}])
@both_plots
def test_fast_draws_what_legacy_draws(adata_small, fmt, kwargs, plot):
    """With a row per cell, `mode="fast"` draws the cells of `mode="legacy"` in the same order and colours."""
    adata_small = adata_small[adata_small.obs["group"].notna()].copy()  # legacy cannot draw cells without group
    adata_small.obsm["X_cnv"] = fmt(adata_small.obsm["X_cnv"].toarray())
    legacy = plot(adata_small, groupby="group", show=False, mode="legacy", **kwargs)
    fast = plot(adata_small, groupby="group", show=False, **kwargs)
    assert set(fast) == AX_KEYS
    legacy_values = _image(legacy)
    if plot is cnv.pl.chromosome_heatmap_summary:
        # legacy repeats each group 10 times, and sorts the groups by name instead of the order of the categories
        legacy_labels = [t.get_text() for t in legacy["groupby_ax"].get_yticklabels()]
        legacy_values = legacy_values[::10][[legacy_labels.index(c) for c in ["c", "b", "a"]]]
    np.testing.assert_allclose(_image(fast), legacy_values, rtol=1e-6, atol=1e-6)
    for key in ("vmin", "vcenter", "vmax"):
        assert getattr(fast["heatmap_ax"].images[0].norm, key) == pytest.approx(
            getattr(legacy["heatmap_ax"].images[0].norm, key)
        )
    assert [t.get_text() for t in fast["groupby_ax"].get_yticklabels()] == ["c", "b", "a"]


def test_dendrogram_reuses_stored_one(adata_small):
    adata_small = adata_small[adata_small.obs["group"].notna()].copy()
    adata_small.obsm["X_cnv"] = adata_small.obsm["X_cnv"].toarray()  # scanpy's dendrogram needs dense data here
    legacy = cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False, mode="legacy", dendrogram=True)
    assert "dendrogram_group" in adata_small.uns  # scanpy stores it
    fast = cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False, dendrogram=True)
    assert set(fast) == AX_KEYS | {"dendrogram_ax"}
    np.testing.assert_allclose(_image(fast), _image(legacy), rtol=1e-6)


def test_dendrogram_from_mean_profiles(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", dendrogram=True, show=False)
    tmp = AnnData(adata_small.obsm["X_cnv"].toarray(), obs=adata_small.obs)
    tmp = tmp[tmp.obs["group"].notna()].copy()
    tmp.obs["group"] = tmp.obs["group"].cat.remove_unused_categories()
    expected = sc.tl.dendrogram(tmp, "group", var_names=tmp.var_names, inplace=False)["categories_ordered"]
    assert [t.get_text() for t in axes["groupby_ax"].get_yticklabels()] == expected
    summary = cnv.pl.chromosome_heatmap_summary(adata_small, groupby="group", dendrogram=True, show=False)
    assert [t.get_text() for t in summary["groupby_ax"].get_yticklabels()] == expected
    X, g = adata_small.obsm["X_cnv"].toarray(), adata_small.obs["group"]
    np.testing.assert_allclose(_image(summary), [X[(g == c).to_numpy()].mean(axis=0) for c in expected], rtol=1e-6)


@both_plots
def test_unsupported_argument_falls_back_to_legacy(adata_small, caplog, plot):
    adata_small = adata_small[adata_small.obs["group"].notna()].copy()
    with caplog.at_level(logging.INFO, logger=LOGGER):
        axes = plot(adata_small, groupby="group", show=False, interpolation="none")
    assert "colorbar_ax" not in axes  # drawn by scanpy
    [msg] = [r.getMessage() for r in caplog.records if r.name == LOGGER]
    assert msg.endswith("mode='fast' does not support interpolation, using mode='legacy'")


def test_legacy_rejects_ax(adata_small):
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="`ax` needs mode='fast'"):
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", ax=ax, mode="legacy")


def test_legacy_return_fig(adata_small):
    adata_small = adata_small[adata_small.obs["group"].notna()].copy()
    assert isinstance(cnv.pl.chromosome_heatmap(adata_small, groupby="group", mode="legacy", return_fig=True), Figure)


# --- fast mode ----------------------------------------------------------------------------------------------------


def test_sample_picks_middle_cell(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(4, 0.2), show=False)
    values = _image(axes)
    assert 3 <= len(values) < 59
    X = adata_small.obsm["X_cnv"].toarray()
    assert all(any(np.array_equal(row, cell) for cell in X) for row in values)  # whole cells, not means
    norm = axes["heatmap_ax"].images[0].norm
    assert (norm.vmin, norm.vmax) == (X[X != 0].min(), X[X != 0].max())  # colour scale of all cells


def test_mean(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(4, 0.2), aggregate="mean", show=False)
    values = _image(axes)
    codes = adata_small.obs["group"].cat.codes.to_numpy().astype(np.int64)
    order, bounds, _ = pixel_rows(codes, 4, len(values))
    X = adata_small.obsm["X_cnv"].toarray()
    expected = [X[order[a:b]].mean(axis=0) for a, b in zip(bounds[:-1], bounds[1:], strict=True)]
    np.testing.assert_allclose(values, expected, rtol=1e-6)
    norm = axes["heatmap_ax"].images[0].norm
    assert (norm.vmin, norm.vmax) == pytest.approx((np.min(expected), np.max(expected)))
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", aggregate="mean", vmin=-3, vmax=3, show=False)
    assert (axes["heatmap_ax"].images[0].norm.vmin, axes["heatmap_ax"].images[0].norm.vmax) == (-3, 3)


@pytest.mark.parametrize("aggregate", ["sample", "mean"])
@pytest.mark.parametrize("sparse", [False, True])
def test_dask_matches_in_memory(adata_small, aggregate, sparse):
    X = adata_small.obsm["X_cnv"]
    expected = _image(
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(4, 0.2), aggregate=aggregate, show=False)
    )
    adata_small.obsm["X_cnv"] = _dask(X if sparse else X.toarray())
    actual = _image(
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(4, 0.2), aggregate=aggregate, show=False)
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-6)
    summary = _image(cnv.pl.chromosome_heatmap_summary(adata_small, groupby="group", show=False))
    adata_small.obsm["X_cnv"] = X
    np.testing.assert_allclose(
        summary, _image(cnv.pl.chromosome_heatmap_summary(adata_small, groupby="group", show=False))
    )


def test_swap_axes(adata_small):
    plain = cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False)
    swapped = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(16, 10), swap_axes=True, show=False)
    assert set(swapped) == AX_KEYS
    np.testing.assert_array_equal(_image(swapped), _image(plain).T)  # 59 cells fit both ways
    assert [t.get_text() for t in swapped["groupby_ax"].get_xticklabels()] == ["c", "b", "a"]


def test_norm_options(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", vmin=-5, vmax=2, vcenter=0.5, show=False)
    norm = axes["heatmap_ax"].images[0].norm
    assert (norm.vmin, norm.vcenter, norm.vmax) == (-5, 0.5, 2)
    custom = Normalize(-1, 1)
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", norm=custom, show_gene_labels=True, show=False)
    assert axes["heatmap_ax"].images[0].norm is custom


def test_var_group_rotation(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", var_group_rotation=45, show=False)
    assert [t.get_rotation() for t in axes["gene_groups_ax"].texts] == [45] * 3
    assert [t.get_text() for t in axes["gene_groups_ax"].texts] == ["chr1", "chr2", "chr3"]


def test_group_colors(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False)
    blocks = axes["groupby_ax"].images[0]
    assert [tuple(blocks.cmap(i)) for i in range(3)] == [
        (1.0, 0.0, 0.0, 1.0),
        (0.0, 0.5019607843137255, 0.0, 1.0),
        (0.0, 0.0, 0.0, 1.0),
    ]


@both_plots
def test_ax_and_return_fig(adata_small, plot):
    fig, (ax, other) = plt.subplots(1, 2)
    assert plot(adata_small, groupby="group", ax=ax, return_fig=True) is fig
    axes = plot(adata_small, groupby="group", ax=other, show=False)
    assert axes["heatmap_ax"] is other
    assert all(a.get_figure(root=True) is fig for a in axes.values())


def test_return_fig_without_ax(adata_small):
    fig = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(5, 3), return_fig=True)
    assert isinstance(fig, Figure)
    assert tuple(fig.get_size_inches()) == (5, 3)


def test_save(adata_small, tmp_path):
    figdir = sc.settings.figdir
    sc.settings.figdir = tmp_path  # a property without deleter, so no monkeypatch
    try:
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False, save=True)
        cnv.pl.chromosome_heatmap_summary(adata_small, groupby="group", show=False, save="_summary.png")
    finally:
        sc.settings.figdir = figdir
    ext = sc.settings.file_format_figs
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted([f"heatmap.{ext}", "heatmap_summary.png"])


def test_bad_arguments(adata_small):
    with pytest.raises(ValueError, match="tl.leiden"):
        cnv.pl.chromosome_heatmap(adata_small)
    with pytest.raises(ValueError, match="`mode`"):
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", mode="slow")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="`aggregate`"):
        cnv.pl.chromosome_heatmap(adata_small, groupby="group", aggregate="max")  # type: ignore[arg-type]


# --- on infercnv results ------------------------------------------------------------------------------------------


def test_plot_chromosome_heatmap(adata_infercnv):
    cnv.pl.chromosome_heatmap(adata_infercnv, show=False)


def test_plot_chromosome_heatmap_summary(adata_infercnv):
    cnv.pl.chromosome_heatmap_summary(adata_infercnv, show=False)


# Grouped by the annotated cell types, not leiden clusters, so the baseline does not depend on the igraph version.
@both_plots
def test_chromosome_heatmap_image(_adata_infercnv, image_comparer, plot):
    plot(_adata_infercnv.copy(), groupby="cell_type", figsize=(8, 5), show=False)
    image_comparer(plot.__name__)
