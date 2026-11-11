"""The sparse (no-densify) and dense smoothing paths must agree with each other and with a naive reference."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import infercnvpy.tl._infercnv as m


def _naive(x, reference, var, lfc_cap, window_size, step):
    """Straightforward (slow) re-implementation of steps 1-3 using np.convolve."""
    if reference.shape[0] == 1:
        x = x - reference[0]
    else:
        lo, hi = reference.min(0), reference.max(0)
        x = np.where(x > hi, x - hi, 0) + np.where(x < lo, x - lo, 0)
    x = np.clip(x, -lfc_cap, lfc_cap)
    out = []
    for chr_ in m._natural_sort([c for c in var["chromosome"].unique() if c.startswith("chr") and c != "chrM"]):
        genes = np.flatnonzero((var["chromosome"] == chr_).to_numpy())
        genes = genes[np.argsort(var["start"].to_numpy()[genes], kind="stable")]
        xs = x[:, genes]
        if window_size < xs.shape[1]:
            r = np.arange(1, window_size + 1)
            pyr = np.minimum(r, r[::-1])
            sm = np.array([np.convolve(row, pyr, mode="valid") for row in xs]) / pyr.sum()
            sm = sm[:, ::step]
        else:
            sm = xs.mean(axis=1, keepdims=True)
        out.append(sm)
    return np.hstack(out)


@pytest.fixture(scope="module")
def data():
    rng = np.random.default_rng(1)
    n_cells, n_genes = 30, 400
    var = pd.DataFrame(
        {
            "chromosome": rng.choice(["chr1", "chr2", "chr10", "chrM", "chr3"], n_genes, p=[0.4, 0.3, 0.2, 0.02, 0.08]),
            "start": rng.integers(0, 10**6, n_genes),
        },
        index=[f"g{i}" for i in range(n_genes)],
    )
    dense = rng.poisson(0.3, size=(n_cells, n_genes)).astype(np.float32) * rng.random(
        (n_cells, n_genes), dtype=np.float32
    )
    return var, dense


@pytest.mark.parametrize("n_ref", [1, 3])
@pytest.mark.parametrize("window_size,step", [(30, 5), (7, 1), (1000, 10)])
@pytest.mark.parametrize("path", ["sparse", "dense"])
def test_smoothing_matches_naive(data, n_ref, window_size, step, path):
    var, dense = data
    rng = np.random.default_rng(2)
    reference = (rng.random((n_ref, dense.shape[1])) * 0.5).astype(np.float32)
    layout = m._build_layout(var, window_size, step)
    ordered = dense[:, layout.order]
    ref = reference[:, layout.order]
    expected = _naive(dense, reference, var, 1.0, window_size, step)

    if path == "sparse":
        got = m._smooth_sparse(sp.csr_matrix(ordered), ref, layout, 1.0, np.dtype(np.float32))
    else:
        got = m._smooth_dense(ordered, ref, layout, 1.0, np.dtype(np.float32))
    np.testing.assert_allclose(got, expected, rtol=1e-4, atol=1e-5)


def test_sparse_path_handles_explicit_zeros_and_ints(data):
    var, dense = data
    layout = m._build_layout(var, 20, 5)
    ref = np.full((1, dense.shape[1]), 0.5, dtype=np.float32)[:, layout.order]
    xi = np.round(dense[:, layout.order] * 3).astype(np.int64)
    csr = sp.csr_matrix(xi)
    csr.data[::7] = 0  # explicitly stored zeros
    a = m._infercnv_chunk_kernel(csr, ref, layout, 1.0, None)[0].toarray()
    b = m._infercnv_chunk_kernel(csr.toarray(), ref, layout, 1.0, None)[0].toarray()
    np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)


def test_path_choice_does_not_change_result(data, monkeypatch):
    var, dense = data
    layout = m._build_layout(var, 20, 5)
    ref = np.zeros((2, dense.shape[1]), dtype=np.float32)
    ref[1] = 0.4
    x = sp.csr_matrix(dense[:, layout.order])
    monkeypatch.setattr(m, "_SPARSE_PATH_MAX_DENSITY", 2.0)  # force sparse
    a = m._infercnv_chunk_kernel(x, ref[:, layout.order], layout, 1.0, 1.0)[0].toarray()
    monkeypatch.setattr(m, "_SPARSE_PATH_MAX_DENSITY", 0.0)  # force dense
    b = m._infercnv_chunk_kernel(x, ref[:, layout.order], layout, 1.0, 1.0)[0].toarray()
    np.testing.assert_allclose(a, b, rtol=1e-4, atol=1e-5)
