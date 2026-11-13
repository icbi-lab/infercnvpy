import pytest

import infercnvpy as cnv


def test_plot_chromosome_heatmap(adata_infercnv):
    cnv.pl.chromosome_heatmap(adata_infercnv, show=False)


def test_plot_chromosome_heatmap_summary(adata_infercnv):
    cnv.pl.chromosome_heatmap_summary(adata_infercnv, show=False)


# Grouped by the annotated cell types, not leiden clusters, so the baseline does not depend on the igraph version.
@pytest.mark.parametrize("plot", [cnv.pl.chromosome_heatmap, cnv.pl.chromosome_heatmap_summary])
def test_chromosome_heatmap_image(_adata_infercnv, image_comparer, plot):
    plot(_adata_infercnv.copy(), groupby="cell_type", figsize=(8, 5), show=False)
    image_comparer(plot.__name__)
