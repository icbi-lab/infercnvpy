from .fixtures import adata_infercnv
import infercnvpy as cnv


def test_plot_chromosome_heatmap(adata_infercnv):
    cnv.pl.chromosome_heatmap(adata_infercnv, show=False)


def test_plot_chromosome_heatmap_summary(adata_infercnv):
    cnv.pl.chromosome_heatmap_summary(adata_infercnv, show=False)
