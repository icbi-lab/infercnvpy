"""Interquartile range of the cell-cell Pearson correlations within a group, without the cells x cells matrix.

`np.percentile(np.corrcoef(X), [75, 25])` takes quantiles over all n^2 entries of the correlation matrix: the n
diagonal ones and every pair twice. Here the rows are standardized once (float32, centred, unit norm), the
correlations of each block of rows with all later rows come from one BLAS product, and a parallel numba kernel
adds that block's upper triangle to thread-local histograms with weight 2. The quantiles are read off the
cumulative counts (linear interpolation inside a bin, then numpy's linear interpolation between ranks), which
matches numpy to about 1e-6. Memory is the standardized matrix (cells x features, float32) plus one block,
instead of cells^2 float64.

With `n_pairs`, the quantiles are estimated from that many entries of the correlation matrix drawn uniformly at
random (row and column independently, so the diagonal and both triangles appear as often as in the full matrix).
The sample quantiles are unbiased estimates of the exact ones, their error shrinks like 1/sqrt(n_pairs), and the
cost no longer depends on the group size.
"""

from __future__ import annotations

from typing import Any

import numba
import numpy as np
import scipy.sparse as sp

from infercnvpy._numba import _on_main_thread, njit_parallel, prange

N_BINS = 1 << 20  # bin width 2 / 2^20 ~ 2e-6
TILE = 2048  # rows per block: a block of correlations is TILE x cells float32


@njit_parallel
def _standardize_dense(X, out):  # pragma: no cover
    """Centre each row and scale it to unit norm; a constant row becomes NaN."""
    n, g = X.shape
    for i in prange(n):
        m = 0.0
        for c in range(g):
            m += X[i, c]
        m /= g
        ss = 0.0
        for c in range(g):
            v = X[i, c] - m
            out[i, c] = v
            ss += v * v
        inv = 1.0 / np.sqrt(ss) if ss > 0 else np.nan
        for c in range(g):
            out[i, c] *= inv


@njit_parallel
def _standardize_csr(indptr, indices, data, n_cols, out):  # pragma: no cover
    """Like `_standardize_dense`, straight from CSR into a dense float32 row (no float64 densification)."""
    n = len(indptr) - 1
    for i in prange(n):
        s = 0.0
        q = 0.0
        for k in range(indptr[i], indptr[i + 1]):
            s += data[k]
            q += data[k] * data[k]
        m = s / n_cols
        ss = q - s * m  # sum of squared deviations, implicit zeros included
        inv = 1.0 / np.sqrt(ss) if ss > 0 else np.nan
        fill = -m * inv
        for c in range(n_cols):
            out[i, c] = fill
        for k in range(indptr[i], indptr[i + 1]):
            out[i, indices[k]] = (data[k] - m) * inv


@njit_parallel
def _bin_tile(C, row0, hists):  # pragma: no cover
    """Add the entries of correlation block `C` (rows `row0 + a`, columns `row0 + b`) above the diagonal."""
    rows, cols = C.shape
    n_threads, n_bins = hists.shape
    scale = 0.5 * n_bins
    for t in prange(n_threads):
        h = hists[t]
        for a in range(t, rows, n_threads):  # strided rows balance the triangle
            for b in range(a + 1, cols):
                r = C[a, b]
                if r == r:  # pairs with a constant cell are NaN and skipped
                    k = int((r + 1.0) * scale)
                    h[min(max(k, 0), n_bins - 1)] += 2


@njit_parallel
def _sampled_corr(Z, rows, cols, out):  # pragma: no cover
    """Correlations of the sampled matrix entries (rows[k], cols[k]) of standardized rows `Z`."""
    g = Z.shape[1]
    for k in prange(len(rows)):
        i = rows[k]
        j = cols[k]
        if i == j:
            out[k] = 1.0
            continue
        r = 0.0
        for c in range(g):
            r += Z[i, c] * Z[j, c]
        out[k] = r


def _quantiles(hist: np.ndarray, qs: tuple[float, ...]) -> list[float]:
    """Numpy's default (linear) quantiles of the values summarised by `hist` over [-1, 1]."""
    n_bins = len(hist)
    cum = np.cumsum(hist)
    total = int(cum[-1])
    width = 2.0 / n_bins

    def value(k: int) -> float:
        b = int(np.searchsorted(cum, k, side="right"))
        before = int(cum[b - 1]) if b else 0
        return -1.0 + width * (b + (k - before + 0.5) / hist[b])  # ranks spread evenly inside the bin

    out = []
    for q in qs:
        h = q * (total - 1)
        lo = int(np.floor(h))
        out.append(value(lo) + (h - lo) * (value(min(lo + 1, total - 1)) - value(lo)))
    return out


def corr_iqr(X: Any, *, tile: int = TILE, n_pairs: int | None = None, rng: np.random.Generator | None = None) -> float:
    """IQR of the Pearson correlations between all rows of `X` (numpy array or scipy sparse matrix).

    Rows with constant values have no defined correlation and are left out. Returns NaN with fewer than two rows.
    With `n_pairs`, groups with more than `n_pairs` matrix entries are estimated from a uniform sample of entries.
    """
    Z = np.empty(X.shape, dtype=np.float32)
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        _standardize_csr(X.indptr, X.indices, X.data.astype(np.float32, copy=False), X.shape[1], Z)
    else:
        _standardize_dense(np.ascontiguousarray(X, dtype=np.float32), Z)
    keep = ~np.isnan(Z[:, 0]) if X.shape[1] else np.zeros(X.shape[0], dtype=bool)
    Z = Z[keep]
    n = Z.shape[0]
    if n <= 1:
        return np.nan
    if n_pairs is not None and n * n > n_pairs:
        rng = np.random.default_rng() if rng is None else rng
        rows = rng.integers(0, n, n_pairs)
        cols = rng.integers(0, n, n_pairs)
        r = np.empty(n_pairs, dtype=np.float32)
        _sampled_corr(Z, rows, cols, r)
        q75, q25 = np.percentile(r, [75, 25])
        return float(q75 - q25)
    if n <= tile:
        # one block: the exact float64 quantiles cost at most tile^2 values, as before 0.7.0
        rows = X[np.flatnonzero(keep)]
        q75, q25 = np.percentile(np.corrcoef(rows.toarray() if sp.issparse(rows) else rows), [75, 25])
        return q75 - q25
    hists = np.zeros((numba.get_num_threads() if _on_main_thread() else 1, N_BINS), dtype=np.int64)
    for r0 in range(0, n, tile):
        _bin_tile(Z[r0 : r0 + tile] @ Z[r0:].T, r0, hists)
    hist = hists.sum(axis=0)
    hist[-1] += n  # the diagonal: every cell correlates 1 with itself
    q75, q25 = _quantiles(hist, (0.75, 0.25))
    return q75 - q25
