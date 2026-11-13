import warnings

import anndata as ad
import leidenalg
import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
import scanpy as sc
from anndata.acc import A

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
    with pytest.warns(UserWarning, match="mean of all cells as reference") as record:
        actual = _get_reference(adata_mock, None, None, None, layer=None)
    assert record[0].filename == __file__  # points at the caller, not at infercnvpy
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
    assert "X_cnv" in adata_oligodendroma.obsm, "cnv not in adata.obsm"
    assert "cnv" in adata_oligodendroma.uns, "cnv not in adata.uns"
    ## Test that the gene values are not calculated by default
    assert "gene_values_cnv" not in adata_oligodendroma.layers.keys(), "gene_values_cnv in .layers"


def test_infercnv_gene_values(adata_oligodendroma):
    cnv.tl.infercnv(adata_oligodendroma, calculate_gene_values=True)
    assert "X_cnv" in adata_oligodendroma.obsm, "cnv not in adata.obsm"
    assert "cnv" in adata_oligodendroma.uns, "cnv not in adata.uns"
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


def test_infercnv_records_params(adata_full_mock, tmp_path):
    cnv.tl.infercnv(adata_full_mock, reference_key="sample", reference_cat="sample1", window_size=3, step=1)
    assert set(adata_full_mock.uns["cnv"]) == {"chr_pos", "params"}
    assert adata_full_mock.uns["cnv"]["params"] == {
        "reference_key": "sample",
        "reference_cat": "sample1",
        "lfc_clip": 3,
        "window_size": 3,
        "step": 1,
        "dynamic_threshold": 1.5,
        "exclude_chromosomes": ["chrX", "chrY"],
        "chunksize": 5000,
        "layer": None,
        "calculate_gene_values": False,
    }
    # a disabled step is recorded as None, and params round-trip through h5ad
    cnv.tl.infercnv(adata_full_mock, reference_key="sample", reference_cat="sample1", dynamic_threshold=None)
    adata_full_mock.write_h5ad(tmp_path / "cnv.h5ad")
    assert ad.read_h5ad(tmp_path / "cnv.h5ad").uns["cnv"]["params"]["dynamic_threshold"] is None


def test_infercnv_no_chr_prefix(adata_full_mock):
    """Chromosome names without a "chr" prefix (Ensembl style) give a clear error."""
    adata_full_mock.var["chromosome"] = adata_full_mock.var["chromosome"].str.removeprefix("chr")
    with pytest.raises(ValueError, match='start with "chr"'):
        cnv.tl.infercnv(adata_full_mock, window_size=3, step=1)


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


def test_tsne_use():
    """`tl.tsne` embeds the representation passed as `use`, not always `X_cnv_pca`."""
    rng = np.random.default_rng(0)
    adata = sc.AnnData(np.zeros((60, 1)), obsm={"X_other": rng.normal(size=(60, 5))})
    res = cnv.tl.tsne(adata, use="obsm.X_other", out=True)
    assert res.shape == (60, 2)
    assert "X_cnv_pca" not in adata.obsm


@pytest.fixture
def adata_cnv_neighbors():
    rng = np.random.default_rng(0)
    adata = sc.AnnData(np.zeros((100, 1)), obsm={"X_cnv_pca": rng.normal(size=(100, 5))})
    cnv.pp.neighbors(adata)
    return adata


def test_neighbors_layout(adata_cnv_neighbors):
    """The graph uses scanpy>=1.13's layout, and `out` picks its name."""
    assert {"cnv_neighbors", "cnv_neighbors_distances"} <= set(adata_cnv_neighbors.obsp)
    assert adata_cnv_neighbors.uns["cnv_neighbors"]["connectivities_key"] == "cnv_neighbors"
    cnv.pp.neighbors(adata_cnv_neighbors, out=A.obsp["other"])
    assert {"other", "other_distances"} <= set(adata_cnv_neighbors.obsp)
    cnv.tl.leiden(adata_cnv_neighbors, use="obsp.other", out="obs.other_leiden")
    assert adata_cnv_neighbors.obs["other_leiden"].nunique() > 1


def test_leiden_defaults_to_igraph(adata_cnv_neighbors):
    cnv.tl.leiden(adata_cnv_neighbors)
    assert adata_cnv_neighbors.uns["cnv_leiden"]["params"]["n_iterations"] == 2
    cnv.tl.leiden(adata_cnv_neighbors, n_iterations=5)
    assert adata_cnv_neighbors.uns["cnv_leiden"]["params"]["n_iterations"] == 5


@pytest.mark.parametrize(
    "kwargs",
    [
        {"flavor": "leidenalg"},
        # only leidenalg supports partition_type, so it is used without an explicit flavor
        {"partition_type": leidenalg.RBConfigurationVertexPartition},
    ],
)
def test_leiden_leidenalg(adata_cnv_neighbors, kwargs):
    # scanpy adds a UserWarning recommending igraph
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        cnv.tl.leiden(adata_cnv_neighbors, **kwargs)
    # leidenalg keeps scanpy's own default of iterating until convergence
    assert adata_cnv_neighbors.uns["cnv_leiden"]["params"]["n_iterations"] == -1
    assert adata_cnv_neighbors.obs["cnv_leiden"].nunique() > 1


def test_leiden_no_flavor_warning(adata_cnv_neighbors):
    """At default settings, scanpy's leiden flavor FutureWarning is not raised."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        cnv.tl.leiden(adata_cnv_neighbors)
    assert not [w for w in record if issubclass(w.category, FutureWarning) or "leiden" in str(w.message)]
    assert adata_cnv_neighbors.obs["cnv_leiden"].nunique() > 1


def test_layer_parameter(adata_oligodendroglioma):

    # create lognorm layer but leave X unchanged
    adata_oligodendroglioma.layers["LogNormalize"] = adata_oligodendroglioma.X.copy()
    lognorm = A.layers["LogNormalize"]
    sc.pp.normalize_total(adata_oligodendroglioma, use=lognorm, out=lognorm, target_sum=1e4)
    sc.pp.log1p(adata_oligodendroglioma, use=lognorm, out=lognorm)

    # copy to test if results change with layer option
    adata_oligodendroglioma2 = adata_oligodendroglioma.copy()
    adata_oligodendroglioma2.X = adata_oligodendroglioma.layers["LogNormalize"]

    cnv.tl.infercnv(adata_oligodendroglioma, layer="LogNormalize")
    cnv.tl.infercnv(adata_oligodendroglioma2, layer=None)

    X_cnv = adata_oligodendroglioma.obsm["X_cnv"].toarray()
    X_cnv2 = adata_oligodendroglioma2.obsm["X_cnv"].toarray()

    assert np.all(X_cnv == X_cnv2), "Different results found with infercnv layer parameter"


# The benchmarks are skipped by default (`--benchmark-skip` in pyproject.toml); run them with `--benchmark-only`.
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


def test_default_keys_follow_scanpy_preset(monkeypatch):
    """Under scanpy's 2.0 preview preset, default `.obsm` keys have no `X_` prefix."""
    # setting the real preset requires the `scanpy[scanpy2]` extras
    monkeypatch.setattr("infercnvpy._acc._v2_preset", lambda: True)
    rng = np.random.default_rng(0)
    adata = sc.AnnData(np.zeros((60, 1)), obsm={"cnv": rng.normal(size=(60, 8))})
    cnv.pp.neighbors(adata)
    assert "cnv_pca" in adata.obsm
    assert "X_cnv_pca" not in adata.obsm
