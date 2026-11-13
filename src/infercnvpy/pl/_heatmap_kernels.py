"""numba kernel of :func:`infercnvpy.pl.chromosome_heatmap`: sums of the cells that fall into each pixel row.

The cells are sorted by group and every group is split into contiguous runs of cells, one per pixel row
(`_pixel_rows`). `_bin_sums_csr` / `_bin_sums_dense` add up the rows of each run straight from the CNV matrix (CSR
without densifying), one pixel row per `prange` iteration, and record the smallest and largest value they see for the
colour scale. A figure is ~1000 pixels tall, so the result is a small dense matrix whatever the number of cells.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.sparse as sp

from infercnvpy._numba import njit_parallel, prange


@njit_parallel
def _bin_sums_csr(indptr, indices, data, order, bounds, sums, mins, maxs):  # pragma: no cover
    """Add rows `order[bounds[b]:bounds[b + 1]]` into `sums[b]` (zeroed); stored values only for `mins`/`maxs`."""
    for b in prange(len(bounds) - 1):
        lo = np.inf
        hi = -np.inf
        row = sums[b]
        for k in range(bounds[b], bounds[b + 1]):
            i = order[k]
            for p in range(indptr[i], indptr[i + 1]):
                v = data[p]
                row[indices[p]] += v
                if v < lo:
                    lo = v
                if v > hi:
                    hi = v
        mins[b] = lo
        maxs[b] = hi


@njit_parallel
def _bin_sums_dense(X, order, bounds, sums, mins, maxs):  # pragma: no cover
    """Like `_bin_sums_csr` for a dense matrix."""
    n_cols = X.shape[1]
    for b in prange(len(bounds) - 1):
        lo = np.inf
        hi = -np.inf
        row = sums[b]
        for k in range(bounds[b], bounds[b + 1]):
            x = X[order[k]]
            for c in range(n_cols):
                v = x[c]
                row[c] += v
                if v < lo:
                    lo = v
                if v > hi:
                    hi = v
        mins[b] = lo
        maxs[b] = hi


def _rows_per_group(counts: np.ndarray, n_rows: int) -> np.ndarray:
    """Split `n_rows` pixel rows over groups in proportion to `counts` (largest remainder).

    A non-empty group gets at least one row, and no group more rows than cells.
    """
    exact = counts * (n_rows / counts.sum())
    rows = np.floor(exact).astype(np.int64)
    rest = n_rows - rows.sum()
    rows[np.argsort(rows - exact, kind="stable")[:rest]] += 1
    return np.clip(rows, counts > 0, counts)


def _pixel_rows(codes: np.ndarray, n_groups: int, n_rows: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Order the cells by group (`codes`, -1 is left out) and cut each group into runs of cells, one per pixel row.

    Returns the cell order, the start of each run in that order (plus the end), and the number of rows per group.
    """
    order = np.argsort(codes, kind="stable")
    order = order[codes[order] >= 0]
    counts = np.bincount(codes[order], minlength=n_groups)
    rows = _rows_per_group(counts, min(n_rows, len(order)))
    starts = np.concatenate([[0], np.cumsum(counts)])
    bounds = [starts[g] + np.arange(rows[g]) * counts[g] // rows[g] for g in range(n_groups) if rows[g]]
    return order, np.concatenate([*bounds, [len(order)]]).astype(np.int64), rows


def binned_sums(X: Any, order: np.ndarray, bounds: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Sum of the rows of `X` in each run `order[bounds[b]:bounds[b + 1]]`, and the min and max over those rows.

    `X` is a numpy array or scipy sparse matrix (for sparse input, the min and max of the stored values).
    """
    n_bins = len(bounds) - 1
    sums = np.zeros((n_bins, X.shape[1]), dtype=np.float64)
    mins = np.empty(n_bins, dtype=np.float64)
    maxs = np.empty(n_bins, dtype=np.float64)
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        _bin_sums_csr(X.indptr, X.indices, X.data, order, bounds, sums, mins, maxs)
    else:
        _bin_sums_dense(np.asarray(X), order, bounds, sums, mins, maxs)
    return sums, float(mins.min(initial=np.inf)), float(maxs.max(initial=-np.inf))
