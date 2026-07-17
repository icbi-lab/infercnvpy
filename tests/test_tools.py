import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import scanpy as sc

import infercnvpy as cnv
from infercnvpy.tl._infercnv import _get_reference


def test_get_reference_key_and_cat(adata_mock):
    """Test that reference is correctly calculated given a reference key and category"""
    actual = _get_reference(adata_mock, "cat", ["foo", "baz"], None, layer=None)
    npt.assert_almost_equal(
        actual,
        np.array(
            [
                [1.5, 1, 1.5, 2],
                [7, 5, 5, 7],
            ]
        ),
    )


def test_get_reference_no_reference(adata_mock):
    """If no reference is specified, the mean of the entire adata object is taken"""
    actual = _get_reference(adata_mock, None, None, None, layer=None)
    npt.assert_almost_equal(actual, np.array([[4.8, 4.2, 4.4, 5]]), decimal=5)


def test_get_reference_given_reference(adata_mock):
    """Predefined reference takes precendence over reference_key and reference_cat"""
    reference = np.array([1, 2, 3, 4])
    actual = _get_reference(adata_mock, "foo", "bar", reference, layer=None)
    npt.assert_equal(reference, actual[0, :])

    with pytest.raises(ValueError):
        reference = np.array([1, 2, 3])
        actual = _get_reference(adata_mock, "foo", "bar", reference, layer=None)


@pytest.mark.parametrize(
    "reference_key,reference_cat",
    [
        (None, None),
        ("cell_type", ["Microglia/Macrophage", "Oligodendrocytes (non-malignant)"]),
    ],
)
def test_infercnv(adata_oligodendroma, reference_key, reference_cat):
    cnv.tl.infercnv(adata_oligodendroma, reference_key=reference_key, reference_cat=reference_cat)
    assert "X_cnv" in adata_oligodendroma.obsm_keys(), "cnv not in adata.obsm"
    assert "cnv" in adata_oligodendroma.uns_keys(), "cnv not in adata.uns"
    ## Test that the gene values are not calculated by default
    assert "gene_values_cnv" not in adata_oligodendroma.layers.keys(), "gene_values_cnv in .layers"


def test_infercnv_gene_values(adata_oligodendroma):
    cnv.tl.infercnv(adata_oligodendroma, calculate_gene_values=True)
    assert "X_cnv" in adata_oligodendroma.obsm_keys(), "cnv not in adata.obsm"
    assert "cnv" in adata_oligodendroma.uns_keys(), "cnv not in adata.uns"
    assert "gene_values_cnv" in adata_oligodendroma.layers.keys(), "gene_values_cnv not in .layers"


def test_running_mean_n_less_than_genes():
    # Create a 2D numpy array
    x = np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
    n = 3
    step = 1
    gene_list = np.array(["gene1", "gene2", "gene3", "gene4", "gene5"])

    # Call the function with the test parameters
    result, convolved_gene_values = cnv.tl._infercnv._running_mean(x, n, step, gene_list, calculate_gene_values=True)

    # Define the expected output
    expected_result = np.array([[2, 3, 4], [7, 8, 9]])
    expected_gene_values = pd.DataFrame(
        {
            "gene1": [2.0, 7.0],
            "gene2": [2.5, 7.5],
            "gene3": [3.0, 8.0],
            "gene4": [3.5, 8.5],
            "gene5": [4.0, 9.0],
        }
    )

    # Assert that the function output is as expected
    np.testing.assert_array_equal(result, expected_result)
    # Assert that the gene dfs are equal
    pd.testing.assert_frame_equal(convolved_gene_values, expected_gene_values)


def test_running_mean_n_greater_than_genes():
    # Create a 2D numpy array
    x = np.array([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
    n = 7
    step = 1
    gene_list = np.array(["gene1", "gene2", "gene3", "gene4", "gene5"])

    # Call the function with the test parameters
    result, convolved_gene_values = cnv.tl._infercnv._running_mean(x, n, step, gene_list, calculate_gene_values=True)

    # Define the expected output
    expected_result = np.array([[3], [8]])
    expected_gene_values = pd.DataFrame(
        {
            "gene1": [3.0, 8.0],
            "gene2": [3.0, 8.0],
            "gene3": [3.0, 8.0],
            "gene4": [3.0, 8.0],
            "gene5": [3.0, 8.0],
        }
    )

    # Assert that the function output is as expected
    np.testing.assert_array_equal(result, expected_result)
    # Assert that the gene dfs are equal
    pd.testing.assert_frame_equal(convolved_gene_values, expected_gene_values)


def test_calculate_gene_averages():
    convolved_gene_names = np.array(
        [["gene1", "gene2", "gene3"], ["gene2", "gene3", "gene4"], ["gene3", "gene4", "gene5"]]
    )

    smoothed_x = np.array([[2, 3, 4], [4, 4, 6], [6, 2, 1]])

    convolved_gene_values = cnv.tl._infercnv._calculate_gene_averages(convolved_gene_names, smoothed_x)

    pd.testing.assert_frame_equal(
        convolved_gene_values,
        pd.DataFrame(
            {
                "gene1": [2.0, 4.0, 6.0],
                "gene2": [2.5, 4.0, 4.0],
                "gene3": [3.000000, 4.666667, 3.000000],
                "gene4": [3.5, 5.0, 1.5],
                "gene5": [4.0, 6.0, 1.0],
            }
        ),
    )


def test_infercnv_chunk_with_gene_values(adata_full_mock, gene_res_actual, x_res_actual):
    reference = _get_reference(adata_full_mock, reference_key=None, reference_cat=None, reference=None, layer=None)
    var = adata_full_mock.var.loc[:, ["chromosome", "start", "end"]]
    tmp_x = adata_full_mock.X

    chr_pos, x_res, gene_res = cnv.tl._infercnv._infercnv_chunk(
        tmp_x, var, reference, lfc_cap=1, window_size=3, step=1, dynamic_threshold=1, calculate_gene_values=True
    )

    gene_res.index = gene_res.index.astype(str)
    pd.testing.assert_frame_equal(gene_res, gene_res_actual)
    np.testing.assert_array_equal(x_res.toarray(), x_res_actual)
    assert chr_pos == {"chr1": 0, "chr2": 3}, "chr_pos is not as expected"


def test_infercnv_chunk_default(adata_full_mock, x_res_actual):
    reference = _get_reference(adata_full_mock, reference_key=None, reference_cat=None, reference=None, layer=None)
    var = adata_full_mock.var.loc[:, ["chromosome", "start", "end"]]
    tmp_x = adata_full_mock.X

    chr_pos, x_res, gene_res = cnv.tl._infercnv._infercnv_chunk(
        tmp_x, var, reference, lfc_cap=1, window_size=3, step=1, dynamic_threshold=1
    )

    assert gene_res is None
    np.testing.assert_array_equal(x_res.toarray(), x_res_actual)
    assert chr_pos == {"chr1": 0, "chr2": 3}, "chr_pos is not as expected"


def test_infercnv_more_than_2_chunks(adata_full_mock, x_res_actual):
    chr_pos, res, per_gene_mtx = cnv.tl.infercnv(
        adata_full_mock,
        reference_key=None,
        reference_cat=None,
        reference=None,
        chunksize=2,
        lfc_clip=1,
        window_size=3,
        step=1,
        dynamic_threshold=1,
        calculate_gene_values=True,
        inplace=False,
    )

    ## each chunk will contain 2 samples, this simulatanously tests the chunking and the merging
    np.testing.assert_array_equal(per_gene_mtx[0], np.array([0.75, 0.0, 0.0, 0.0, -0.75, 0.0, 0.0, 0.0, 0.0, 0.75]))
    np.testing.assert_array_equal(per_gene_mtx[3], np.array([0, 0, 0, 0, 0, 0.921875, 0.703125, 0, 0, 0]))
    np.testing.assert_array_equal(res.toarray(), x_res_actual)
    assert chr_pos == {"chr1": 0, "chr2": 3}, "chr_pos is not as expected"


def test_infercnv_manual_reference(adata_oligodendroma):
    cnv.tl.infercnv(adata_oligodendroma, reference=np.ones(adata_oligodendroma.shape[1]))


@pytest.mark.skip(
    reason="rpy2 segfaults on the CI. I don't know why and don't have the time for a painful debugging session."
)
def test_copykat(adata_oligodendroma):
    sc.pp.subsample(adata_oligodendroma, n_obs=50)
    cnv.tl.copykat(adata_oligodendroma)


def test_workflow(adata_oligodendroma):
    cnv.tl.infercnv(adata_oligodendroma)
    cnv.tl.pca(adata_oligodendroma)
    cnv.pp.neighbors(adata_oligodendroma)
    cnv.tl.tsne(adata_oligodendroma)
    cnv.tl.umap(adata_oligodendroma)
    cnv.tl.leiden(adata_oligodendroma)
    cnv.tl.cnv_score(adata_oligodendroma)

    cnv.pl.umap(adata_oligodendroma, color=["cnv_score", "cnv_leiden"], show=False)
    cnv.pl.tsne(adata_oligodendroma, color=["cnv_score", "cnv_leiden"], show=False)
    cnv.pl.chromosome_heatmap(adata_oligodendroma, show=False)
    cnv.pl.chromosome_heatmap_summary(adata_oligodendroma, show=False)


def test_layer_parameter(adata_oligodendroma):
    adata_oligodendroglioma = cnv.datasets.oligodendroglioma()

    # create lognorm layer but leave X unchanged
    adata_oligodendroglioma.layers["LogNormalize"] = adata_oligodendroglioma.X.copy()
    sc.pp.normalize_total(adata_oligodendroglioma, layer="LogNormalize", target_sum=1e4)
    sc.pp.log1p(adata_oligodendroglioma, layer="LogNormalize")

    # copy to test if results change with layer option
    adata_oligodendroglioma2 = adata_oligodendroglioma.copy()
    adata_oligodendroglioma2.X = adata_oligodendroglioma.layers["LogNormalize"]

    cnv.tl.infercnv(adata_oligodendroglioma, layer="LogNormalize")
    cnv.tl.infercnv(adata_oligodendroglioma2, layer=None)

    X_cnv = adata_oligodendroglioma.obsm["X_cnv"].toarray()
    X_cnv2 = adata_oligodendroglioma2.obsm["X_cnv"].toarray()

    assert np.all(X_cnv == X_cnv2), "Different results found with infercnv layer parameter"


def test_calculate_gene_values_speed(benchmark, adata_oligodendroma):
    # Benchmark with calculate_gene_values=True
    benchmark(
        cnv.tl.infercnv,
        adata_oligodendroma,
        calculate_gene_values=True,
    )


def test_default_infercnv(benchmark, adata_oligodendroma):
    benchmark(
        cnv.tl.infercnv,
        adata_oligodendroma,
        calculate_gene_values=False,
    )
