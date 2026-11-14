import warnings

import leidenalg
import numpy as np
import numpy.testing as npt
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


@pytest.mark.parametrize("chunksize", [1, 3, 5000])
def test_infercnv_mock_values(adata_full_mock, gene_res_actual, x_res_actual, chunksize):
    """Values of infercnvpy <= 0.6 with all cells in one chunk, whose threshold is the std of all cells."""
    chr_pos, res, per_gene_mtx = cnv.tl.infercnv(
        adata_full_mock,
        lfc_clip=1,
        window_size=3,
        step=1,
        dynamic_threshold=1,
        noise_reference="all",
        chunksize=chunksize,
        calculate_gene_values=True,
        inplace=False,
    )
    np.testing.assert_array_equal(res.toarray(), x_res_actual)
    np.testing.assert_allclose(per_gene_mtx, gene_res_actual.to_numpy(), rtol=1e-6)
    assert res.dtype == per_gene_mtx.dtype == np.float32
    assert chr_pos == {"chr1": 0, "chr2": 3}, "chr_pos is not as expected"


def test_infercnv_records_params(adata_full_mock, tmp_path):
    cnv.tl.infercnv(adata_full_mock, reference_key="sample", reference_cat="sample1", window_size=3, step=1)
    assert set(adata_full_mock.uns["cnv"]) == {"chr_pos", "window_start", "window_end", "params"}
    assert adata_full_mock.uns["cnv"]["params"] == {
        "reference_key": "sample",
        "reference_cat": "sample1",
        "lfc_clip": 3,
        "window_size": 3,
        "step": 1,
        "dynamic_threshold": 1.5,
        "noise_reference": "reference",
        "exclude_chromosomes": ["chrX", "chrY"],
        "chunksize": 5000,
        "layer": None,
        "calculate_gene_values": False,
        "min_mean": None,
    }
    # a disabled step is recorded as None, and params round-trip through h5ad
    cnv.tl.infercnv(adata_full_mock, reference_key="sample", reference_cat="sample1", dynamic_threshold=None)
    adata_full_mock.write_h5ad(tmp_path / "cnv.h5ad")
    assert sc.read_h5ad(tmp_path / "cnv.h5ad").uns["cnv"]["params"]["dynamic_threshold"] is None


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


def test_tsne_use_rep():
    """`tl.tsne` embeds the representation passed as `use_rep`, not always `X_cnv_pca`."""
    rng = np.random.default_rng(0)
    adata = sc.AnnData(np.zeros((60, 1)), obsm={"X_other": rng.normal(size=(60, 5))})
    res = cnv.tl.tsne(adata, use_rep="other", inplace=False)
    assert res.shape == (60, 2)
    assert "X_cnv_pca" not in adata.obsm


@pytest.fixture
def adata_cnv_neighbors():
    rng = np.random.default_rng(0)
    adata = sc.AnnData(np.zeros((100, 1)), obsm={"X_cnv_pca": rng.normal(size=(100, 5))})
    cnv.pp.neighbors(adata)
    return adata


def test_neighbors_copy(adata_cnv_neighbors):
    """`copy=True` returns a copy of adata with the graph and leaves adata unchanged."""
    res = cnv.pp.neighbors(adata_cnv_neighbors, key_added="other", copy=True)
    assert isinstance(res, sc.AnnData)
    assert "other_connectivities" in res.obsp
    assert "other_connectivities" not in adata_cnv_neighbors.obsp


@pytest.mark.parametrize(
    ("func", "kwargs", "added"),
    [
        (cnv.pp.neighbors, {"key_added": "other"}, lambda a: a.obsp["other_connectivities"]),
        (cnv.tl.leiden, {"key_added": "other"}, lambda a: a.obs["other"]),
    ],
)
def test_inplace_false_still_returns_a_copy(adata_cnv_neighbors, func, kwargs, added):
    """`inplace=False` keeps returning a modified copy for one release, with a FutureWarning."""
    with pytest.warns(FutureWarning, match=r"inplace=False\)` returns a modified copy.*copy=True"):
        res = func(adata_cnv_neighbors, inplace=False, **kwargs)
    assert isinstance(res, sc.AnnData)
    added(res)
    with pytest.raises(KeyError):
        added(adata_cnv_neighbors)
    with pytest.raises(ValueError, match="not both"):
        func(adata_cnv_neighbors, inplace=False, copy=True, **kwargs)


def test_leiden_copy(adata_cnv_neighbors):
    res = cnv.tl.leiden(adata_cnv_neighbors, copy=True)
    assert "cnv_leiden" in res.obs
    assert "cnv_leiden" not in adata_cnv_neighbors.obs


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
    # scanpy>=1.12 adds a UserWarning recommending igraph
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


def test_positional_arguments_from_0_6_still_work(adata_cnv_neighbors):
    """The arguments that became keyword-only in 0.7.0 still work positionally, with a FutureWarning."""
    adata = adata_cnv_neighbors
    adata.obsm["X_cnv"] = np.random.default_rng(1).normal(size=(adata.n_obs, 20))
    deprecated = pytest.warns(FutureWarning, match="positionally is deprecated")
    with deprecated:
        cnv.pp.neighbors(adata, "cnv_pca", "nb")
    assert "nb_connectivities" in adata.obsp
    with deprecated:
        cnv.tl.leiden(adata, "nb", "cl")
    assert "cl" in adata.obs
    with deprecated:
        res = cnv.tl.pca(adata, "arpack", False, False, "cnv", n_comps=3)
    assert res.shape == (adata.n_obs, 3)
    with deprecated:
        cnv.tl.umap(adata, "nb", "um")
    assert "X_um" in adata.obsm
    with deprecated:
        cnv.tl.tsne(adata, "cnv_pca", "ts")
    assert "X_ts" in adata.obsm
    with pytest.raises(TypeError, match="positional arguments"):
        cnv.pp.neighbors(adata, "cnv_pca", "nb", True, "extra")


def test_copykat_accepts_0_6_positional_arguments():
    """copykat needs rpy2 and R (not installed for the tests): check that the positional arguments get through."""
    import inspect

    import infercnvpy.tl._copykat as ck

    names = list(inspect.signature(ck.copykat).parameters)[1:]
    assert names[:3] == ["gene_ids", "organism", "segmentation_cut"]
    import importlib.util

    if importlib.util.find_spec("rpy2") is not None:
        pytest.skip("with rpy2 installed, copykat would run R")
    with pytest.warns(FutureWarning, match="`gene_ids`, `organism`"), pytest.raises(ImportError, match="rpy2"):
        cnv.tl.copykat(sc.AnnData(np.ones((2, 2))), "S", "human")
