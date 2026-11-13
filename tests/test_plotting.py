import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

import infercnvpy as cnv
from infercnvpy.pl._chromosome_heatmap import _group_means
from infercnvpy.pl._heatmap_kernels import _pixel_rows, _rows_per_group, binned_sums

AX_KEYS = {"heatmap_ax", "groupby_ax", "gene_groups_ax", "colorbar_ax"}


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


def test_rows_per_group():
    counts = np.array([500, 3, 0, 97, 1])
    rows = _rows_per_group(counts, 100)
    np.testing.assert_array_equal(rows, [83, 1, 0, 16, 1])  # proportional, >= 1 per non-empty group
    np.testing.assert_array_equal(_rows_per_group(counts, counts.sum()), counts)  # one row per cell


@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix, sp.csc_matrix])
@pytest.mark.parametrize("n_rows", [7, 20, 1000])
def test_binned_sums_match_numpy(fmt, n_rows):
    rng = np.random.default_rng(1)
    X = rng.normal(size=(120, 9))
    X[np.abs(X) < 1] = 0
    codes = rng.choice([0, 1, 3, -1], size=120, p=[0.7, 0.2, 0.05, 0.05])  # group 2 empty, -1 no group
    order, bounds, rows = _pixel_rows(codes, 4, n_rows)

    np.testing.assert_array_equal(np.sort(order), np.flatnonzero(codes >= 0))
    assert np.all(np.diff(codes[order]) >= 0)  # sorted by group
    assert np.all(np.diff(bounds) > 0)  # no empty pixel row
    assert rows[2] == 0
    assert np.all(rows[[0, 1, 3]] >= 1)
    assert rows.sum() == len(bounds) - 1
    # each pixel row holds cells of a single group
    starts = np.concatenate([[0], np.cumsum(np.bincount(codes[order], minlength=4))])
    assert set(bounds[:-1]) >= set(starts[:-1][rows > 0])

    sums, lo, hi = binned_sums(fmt(X), order, bounds)
    expected = np.stack([X[order[a:b]].sum(axis=0) for a, b in zip(bounds[:-1], bounds[1:], strict=True)])
    np.testing.assert_allclose(sums, expected)
    kept = X[order]
    if fmt is np.asarray:
        assert (lo, hi) == (kept.min(), kept.max())
    else:  # stored values only, as before 0.7.0
        assert (lo, hi) == (kept[kept != 0].min(), kept[kept != 0].max())

    means = _group_means(sums, bounds, rows)
    expected = pd.DataFrame(X[codes >= 0]).groupby(codes[codes >= 0]).mean().to_numpy()
    np.testing.assert_allclose(means, expected)


def test_summary_shows_group_means(adata_small):
    axes = cnv.pl.chromosome_heatmap_summary(adata_small, groupby="group", show=False)
    assert set(axes) == AX_KEYS
    X = adata_small.obsm["X_cnv"].toarray()
    g = adata_small.obs["group"]
    expected = np.stack([X[(g == c).to_numpy()].mean(axis=0) for c in ["c", "b", "a"]])
    np.testing.assert_allclose(axes["heatmap_ax"].images[0].get_array(), expected, rtol=1e-6)
    assert [t.get_text() for t in axes["groupby_ax"].get_yticklabels()] == ["c", "b", "a"]
    plt.close("all")


def test_heatmap_rows(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", show=False)
    assert set(axes) == AX_KEYS
    values = axes["heatmap_ax"].images[0].get_array()
    assert values.shape == (59, 12)  # tall figure: one row per cell with a group
    X = adata_small.obsm["X_cnv"].toarray()
    order = np.argsort(adata_small.obs["group"].cat.codes.to_numpy(), kind="stable")[1:]
    np.testing.assert_allclose(values, X[order], rtol=1e-6)
    group_colors = axes["groupby_ax"].images[0]
    assert [tuple(group_colors.cmap(i)) for i in range(3)] == [
        (1.0, 0.0, 0.0, 1.0),
        (0.0, 0.5019607843137255, 0.0, 1.0),
        (0.0, 0.0, 0.0, 1.0),
    ]
    assert [t.get_text() for t in axes["gene_groups_ax"].texts] == ["chr1", "chr2", "chr3"]
    plt.close("all")


def test_heatmap_averages_into_pixel_rows(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(4, 0.2), show=False)
    n_rows = axes["heatmap_ax"].images[0].get_array().shape[0]
    assert 3 <= n_rows < 59
    plt.close("all")


def test_vmin_vmax(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", vmin=-5, vmax=2, show=False)
    norm = axes["heatmap_ax"].images[0].norm
    assert (norm.vmin, norm.vcenter, norm.vmax) == (-5, 0, 2)
    plt.close("all")


@pytest.mark.parametrize("plot", [cnv.pl.chromosome_heatmap, cnv.pl.chromosome_heatmap_summary])
def test_ax_and_return_fig(adata_small, plot):
    fig, (ax, other) = plt.subplots(1, 2)
    assert plot(adata_small, groupby="group", ax=ax, return_fig=True) is fig
    axes = plot(adata_small, groupby="group", ax=other, show=False)
    assert axes["heatmap_ax"] is other
    assert all(a.get_figure(root=True) is fig for a in axes.values())
    plt.close("all")


def test_return_fig_without_ax(adata_small):
    fig = cnv.pl.chromosome_heatmap(adata_small, groupby="group", figsize=(5, 3), return_fig=True)
    assert isinstance(fig, Figure)
    assert tuple(fig.get_size_inches()) == (5, 3)
    plt.close("all")


def test_dendrogram_matches_scanpy(adata_small):
    axes = cnv.pl.chromosome_heatmap(adata_small, groupby="group", dendrogram=True, show=False)
    assert set(axes) == AX_KEYS | {"dendrogram_ax"}
    tmp = AnnData(adata_small.obsm["X_cnv"].toarray(), obs=adata_small.obs)
    tmp = tmp[tmp.obs["group"].notna()].copy()
    tmp.obs["group"] = tmp.obs["group"].cat.remove_unused_categories()
    expected = sc.tl.dendrogram(tmp, "group", var_names=tmp.var_names, inplace=False)["categories_ordered"]
    assert [t.get_text() for t in axes["groupby_ax"].get_yticklabels()] == expected
    plt.close("all")


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
    plt.close("all")


def test_missing_groupby(adata_small):
    with pytest.raises(ValueError, match="tl.leiden"):
        cnv.pl.chromosome_heatmap(adata_small)


def test_plot_chromosome_heatmap(adata_infercnv):
    cnv.pl.chromosome_heatmap(adata_infercnv, show=False)


def test_plot_chromosome_heatmap_summary(adata_infercnv):
    cnv.pl.chromosome_heatmap_summary(adata_infercnv, show=False)


# Grouped by the annotated cell types, not leiden clusters, so the baseline does not depend on the igraph version.
@pytest.mark.parametrize("plot", [cnv.pl.chromosome_heatmap, cnv.pl.chromosome_heatmap_summary])
def test_chromosome_heatmap_image(_adata_infercnv, image_comparer, plot):
    plot(_adata_infercnv.copy(), groupby="cell_type", figsize=(8, 5), show=False)
    image_comparer(plot.__name__)
