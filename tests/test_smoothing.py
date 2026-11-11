"""Tests of the smoothing weight matrices (pyramid weights, `step`, per-gene averages) that `infercnv` uses."""

import numpy as np
import pandas as pd

import infercnvpy.tl._infercnv as m


def _smooth(x, n, step, genes):
    """Apply the layout matrices like `_infercnv_chunk_kernel` does (without centering / clipping)."""
    var = pd.DataFrame({"chromosome": "chr1", "start": np.arange(len(genes))}, index=genes)
    layout = m._build_layout(var, n, step, with_gene_map=True)
    (b,) = layout.blocks  # one chromosome
    smoothed = x @ b.W
    gene_values = pd.DataFrame(smoothed @ b.G, columns=var.index[layout.order[layout.covered_genes]])
    return smoothed, gene_values


def test_running_mean_pyramid_and_step():
    """Pyramid weights (1,2,3,2,1)/9 for window 5, and only every 2nd window is kept."""
    x = np.arange(1, 9, dtype=float)[None, :] ** 2
    result, _ = _smooth(x, 5, 2, np.array([f"g{i}" for i in range(8)]))
    expected = np.convolve(x[0], [1, 2, 3, 2, 1], mode="valid")[::2] / 9
    np.testing.assert_allclose(result[0], expected)


def test_gene_averages_from_layout():
    """Gene values are the mean of all windows that contain the gene (windows of size 3 over 5 genes)."""
    var = pd.DataFrame({"chromosome": "chr1", "start": np.arange(5)}, index=[f"gene{i}" for i in range(1, 6)])
    layout = m._build_layout(var, 3, 1, with_gene_map=True)
    smoothed_x = np.array([[2, 3, 4], [4, 4, 6], [6, 2, 1]])
    got = pd.DataFrame(smoothed_x @ layout.blocks[0].G, columns=var.index[layout.order[layout.covered_genes]])
    expected = pd.DataFrame(
        {
            "gene1": [2.0, 4.0, 6.0],
            "gene2": [2.5, 4.0, 4.0],
            "gene3": [3.000000, 4.666667, 3.000000],
            "gene4": [3.5, 5.0, 1.5],
            "gene5": [4.0, 6.0, 1.0],
        }
    )
    pd.testing.assert_frame_equal(got, expected, check_exact=False, atol=1e-6)
