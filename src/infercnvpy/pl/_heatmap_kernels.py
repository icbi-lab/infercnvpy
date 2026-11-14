"""Reading the CNV matrix for :func:`infercnvpy.pl.chromosome_heatmap`, one pixel row per run of cells.

The cells are sorted by group and every group is split into contiguous runs of cells, one per pixel row
(`pixel_rows`). A pixel row then shows either one cell of its run (`sampled_rows`: the cell in the middle, which is
what matplotlib's nearest-neighbour downscaling shows when every cell is drawn) or the mean of the run (`binned_sums`:
`_bin_sums_csr` / `_bin_sums_dense` add up the rows of each run straight from the matrix, CSR without densifying).
A figure is ~1000 pixels tall, so the result is a small dense matrix whatever the number of cells.

The matrix is read in row blocks (`infercnvpy.tl._infercnv._Rows`): an in-memory matrix as one block with the parallel
kernel, a dask array block by block and a backed matrix `CHUNKSIZE` rows at a time, with the serial kernel.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.sparse as sp

from infercnvpy._numba import njit_parallel, prange
from infercnvpy.tl._infercnv import _Rows

CHUNKSIZE = 10_000  # rows per block of a backed matrix


@njit_parallel
def _bin_sums_csr(indptr, indices, data, order, bounds, bins, sums):  # pragma: no cover
    """Add rows `order[bounds[j]:bounds[j + 1]]` to `sums[bins[j]]`."""
    for j in prange(len(bounds) - 1):
        row = sums[bins[j]]
        for k in range(bounds[j], bounds[j + 1]):
            i = order[k]
            for p in range(indptr[i], indptr[i + 1]):
                row[indices[p]] += data[p]


@njit_parallel
def _bin_sums_dense(X, order, bounds, bins, sums):  # pragma: no cover
    """Like `_bin_sums_csr` for a dense matrix."""
    n_cols = X.shape[1]
    for j in prange(len(bounds) - 1):
        row = sums[bins[j]]
        for k in range(bounds[j], bounds[j + 1]):
            x = X[order[k]]
            for c in range(n_cols):
                row[c] += x[c]


def _rows_per_group(counts: np.ndarray, n_rows: int) -> np.ndarray:
    """Split `n_rows` pixel rows over groups in proportion to `counts` (largest remainder).

    A non-empty group gets at least one row, and no group more rows than cells.
    """
    exact = counts * (n_rows / counts.sum())
    rows = np.floor(exact).astype(np.int64)
    rest = n_rows - rows.sum()
    rows[np.argsort(rows - exact, kind="stable")[:rest]] += 1
    return np.clip(rows, counts > 0, counts)


def pixel_rows(codes: np.ndarray, n_groups: int, n_rows: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Order the cells by group (`codes`, -1 is left out) and cut each group into runs of cells, one per pixel row.

    Within a group, the cells keep their order in `codes`. (:func:`scanpy.pl.heatmap` sorts with pandas' unstable
    default sort, so its order within a group is arbitrary.) Returns the cell order, the start of each run in that order (plus the end), and the
    number of rows per group.
    """
    keep = np.flatnonzero(codes >= 0)
    order = keep[codes[keep].argsort(kind="stable")]
    counts = np.bincount(codes[order], minlength=n_groups)
    rows = _rows_per_group(counts, min(n_rows, len(order)))
    starts = np.concatenate([[0], np.cumsum(counts)])
    bounds = [starts[g] + np.arange(rows[g]) * counts[g] // rows[g] for g in range(n_groups) if rows[g]]
    return order, np.concatenate([*bounds, [len(order)]]).astype(np.int64), rows


def _blocks(X: Any) -> _Rows:
    return _Rows(X, max(X.shape[0], 1) if isinstance(X, np.ndarray) or sp.issparse(X) else CHUNKSIZE)


def _each_block(source: _Rows, fn) -> None:
    """Call `fn(X, row0, n_rows, start)` for every block of `source`, in order (backed and dask blocks read ahead)."""
    blocks = list(range(source.n_blocks))
    for i, (Xb, row0) in zip(blocks, source.read_ahead(blocks), strict=True):
        s, e = int(source.bounds[i]), int(source.bounds[i + 1])
        fn(Xb, row0, e - s, s)


def binned_sums(X: Any, order: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    """Sum of the rows of `X` in each run `order[bounds[b]:bounds[b + 1]]` (float64, runs x columns).

    `X` is a numpy array, scipy sparse matrix, dask array or backed matrix.
    """
    n_bins = len(bounds) - 1
    sums = np.zeros((n_bins, X.shape[1]), dtype=np.float64)
    row_bin = np.full(X.shape[0], -1, dtype=np.int64)
    row_bin[order] = np.repeat(np.arange(n_bins), np.diff(bounds))
    source = _blocks(X)
    serial = source.kind != "in memory"

    def add(Xb, row0, n, start):
        local = np.flatnonzero(row_bin[start : start + n] >= 0)
        b = row_bin[start + local]
        by_bin = np.argsort(b, kind="stable")
        bins, first = np.unique(b[by_bin], return_index=True)
        o = (row0 + local[by_bin]).astype(np.int64)
        bnd = np.append(first, len(o)).astype(np.int64)
        if sp.issparse(Xb):
            _bin_sums_csr(Xb.indptr, Xb.indices, Xb.data, o, bnd, bins, sums, serial=serial)
        else:
            _bin_sums_dense(Xb, o, bnd, bins, sums, serial=serial)

    _each_block(source, add)
    return sums


def sampled_rows(X: Any, picks: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Rows `picks` of `X` as a dense float64 array, and the min and max of all values of `X`.

    For sparse input, the min and max are those of the stored values (as before 0.7.0). NaNs are ignored.
    """
    out = np.empty((len(picks), X.shape[1]), dtype=np.float64)
    sorted_picks = np.argsort(picks, kind="stable")
    ranges = []

    def take(Xb, row0, n, start):
        values = Xb.data if sp.issparse(Xb) else Xb[row0 : row0 + n]
        if values.size:
            ranges.append((np.nanmin(values), np.nanmax(values)))
        lo, hi = np.searchsorted(picks[sorted_picks], [start, start + n])
        dest = sorted_picks[lo:hi]
        rows = Xb[row0 + picks[dest] - start]
        out[dest] = rows.toarray() if sp.issparse(rows) else rows

    source = _blocks(X)
    _each_block(source, take)
    if not ranges:
        return out, np.nan, np.nan
    return out, float(min(r[0] for r in ranges)), float(max(r[1] for r in ranges))
