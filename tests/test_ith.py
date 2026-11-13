import numpy as np
import pytest
import scipy.sparse as sp

from infercnvpy.tl._ith import corr_iqr


def _numpy_iqr(X):
    X = X.toarray() if sp.issparse(X) else X
    q75, q25 = np.percentile(np.corrcoef(X), [75, 25])
    return q75 - q25


def _cnv_like(n, rng):
    """Two subclones with sparse gains/losses plus noise, thresholded like X_cnv."""
    base = rng.normal(size=(2, 300)) * (rng.random((2, 300)) < 0.3)
    X = base[rng.integers(0, 2, n)] + rng.normal(scale=0.5, size=(n, 300))
    X[np.abs(X) < 0.6] = 0
    return X


@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix])
@pytest.mark.parametrize("tile", [64, 10_000])  # histogram kernel over several blocks / exact single block
def test_matches_numpy(fmt, tile):
    X = _cnv_like(500, np.random.default_rng(0))
    assert corr_iqr(fmt(X), tile=tile) == pytest.approx(_numpy_iqr(X), abs=1e-5)


def test_constant_cells_are_left_out():
    X = _cnv_like(300, np.random.default_rng(1))
    with_zero_rows = np.vstack([X, np.zeros((5, X.shape[1]))])
    assert np.isnan(_numpy_iqr(with_zero_rows))  # numpy propagates the undefined correlation
    for tile in (64, 10_000):
        assert corr_iqr(with_zero_rows, tile=tile) == pytest.approx(_numpy_iqr(X), abs=1e-5)


@pytest.mark.parametrize("n", [0, 1])
def test_fewer_than_two_cells(n):
    assert np.isnan(corr_iqr(np.ones((n, 3)) * np.arange(3)))


def test_sampled_pairs_estimate_the_exact_score():
    X = _cnv_like(3000, np.random.default_rng(2))
    exact = corr_iqr(X)
    estimates = [corr_iqr(X, n_pairs=200_000, rng=np.random.default_rng(s)) for s in range(5)]
    assert np.mean(estimates) == pytest.approx(exact, abs=2e-3)
    assert np.std(estimates) < 5e-3
    # small groups (fewer matrix entries than n_pairs) stay exact
    assert corr_iqr(X[:100], n_pairs=200_000) == corr_iqr(X[:100])
