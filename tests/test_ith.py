import numpy as np
import pytest
import scipy.sparse as sp

from infercnvpy.tl._ith import corr_iqr


def _numpy_iqr(X, *, exclude_self_pairs=False):
    X = X.toarray() if sp.issparse(X) else X
    C = np.corrcoef(X)
    if exclude_self_pairs:
        C = C[~np.eye(len(C), dtype=bool)]
    q75, q25 = np.percentile(C, [75, 25])
    return q75 - q25


def _iqr(X, **kwargs):
    return corr_iqr(X, **kwargs)[0]


def _cnv_like(n, rng):
    """Two subclones with sparse gains/losses plus noise, thresholded like X_cnv."""
    base = rng.normal(size=(2, 300)) * (rng.random((2, 300)) < 0.3)
    X = base[rng.integers(0, 2, n)] + rng.normal(scale=0.5, size=(n, 300))
    X[np.abs(X) < 0.6] = 0
    return X


@pytest.mark.parametrize("exclude_self_pairs", [False, True])
@pytest.mark.parametrize("fmt", [np.asarray, sp.csr_matrix])
@pytest.mark.parametrize("tile", [64, 10_000])  # histogram kernel over several blocks / exact single block
def test_matches_numpy(fmt, tile, exclude_self_pairs):
    X = _cnv_like(500, np.random.default_rng(0))
    expected = _numpy_iqr(X, exclude_self_pairs=exclude_self_pairs)
    assert _iqr(fmt(X), tile=tile, exclude_self_pairs=exclude_self_pairs) == pytest.approx(expected, abs=1e-5)


def test_self_pairs_matter_for_small_groups():
    """Self-pairs (always 1) raise the upper quartile of a small group; the default keeps them, as before 0.7.0."""
    X = _cnv_like(12, np.random.default_rng(3))
    assert _iqr(X) == _numpy_iqr(X)
    assert _iqr(X, exclude_self_pairs=True) == _numpy_iqr(X, exclude_self_pairs=True)
    assert _iqr(X, exclude_self_pairs=True) != pytest.approx(_iqr(X), abs=1e-3)


def test_constant_cells_are_left_out():
    X = _cnv_like(300, np.random.default_rng(1))
    with_zero_rows = np.vstack([X, np.zeros((5, X.shape[1]))])
    assert np.isnan(_numpy_iqr(with_zero_rows))  # numpy propagates the undefined correlation
    for tile in (64, 10_000):
        iqr, n_dropped = corr_iqr(with_zero_rows, tile=tile)
        assert iqr == pytest.approx(_numpy_iqr(X), abs=1e-5)
        assert n_dropped == 5


@pytest.mark.parametrize("n", [0, 1])
def test_fewer_than_two_cells(n):
    assert np.isnan(_iqr(np.ones((n, 3)) * np.arange(3)))


@pytest.mark.parametrize("exclude_self_pairs", [False, True])
def test_sampled_pairs_estimate_the_exact_score(exclude_self_pairs):
    X = _cnv_like(3000, np.random.default_rng(2))
    kw = {"exclude_self_pairs": exclude_self_pairs}
    exact = _iqr(X, **kw)
    estimates = [_iqr(X, n_pairs=200_000, rng=np.random.default_rng(s), **kw) for s in range(5)]
    assert np.mean(estimates) == pytest.approx(exact, abs=2e-3)
    assert np.std(estimates) < 5e-3
    # small groups (fewer matrix entries than n_pairs) stay exact
    assert _iqr(X[:100], n_pairs=200_000, **kw) == _iqr(X[:100], **kw)


def test_sampled_pairs_without_self_pairs_estimate_the_off_diagonal_quantiles():
    """With few cells the diagonal is a large share of the matrix, so the sample must leave it out."""
    X = _cnv_like(40, np.random.default_rng(4))
    exact = _numpy_iqr(X, exclude_self_pairs=True)
    estimates = [_iqr(X, n_pairs=1500, rng=np.random.default_rng(s), exclude_self_pairs=True) for s in range(40)]
    assert np.mean(estimates) == pytest.approx(exact, abs=0.01)
