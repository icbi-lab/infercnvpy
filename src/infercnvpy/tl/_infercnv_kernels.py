"""numba kernels of :func:`infercnvpy.tl.infercnv`, one row block of cells at a time.

Per cell, `_cnv_csr` / `_cnv_dense` compute the pyramid-weighted running mean of every kept window, per chromosome
and in genomic gene order, of the clipped log fold changes. The input column of a value is mapped to its gene position
with `col_map` (-1 for dropped genes), so the input is never densified or copied. Two ways give the same sums:

* a sparse cell (at most `max_nnz` values that differ from the baseline, the clipped log fold change of an expression
  of 0): start from the weighted window sums of the baseline, computed once, and add each value's difference to the
  baseline to the windows that contain its gene, with the gene's pyramid weights. O(values x windows per gene);
* any other cell: fill a scratch row with the baseline and the cell's values. A pyramid of `n` genes is the
  convolution of two boxes of `a = ceil(n / 2)` and `b = n + 1 - a` genes, so a window is the difference of two
  entries of a prefix sum of a prefix sum: O(genes), whatever `window_size` is.

Then the median of the cell is subtracted and the windows are written as float32, with the float64 sum and sum of
squares of the written values for the noise threshold and the number of values outside the noise band.

Scratch buffers are allocated once per call, one row per chunk of cells, outside the `prange` loop: allocating inside
the loop body serialises the threads on numba's allocator lock. `np.median` would allocate too, hence `_median`.
"""

from __future__ import annotations

import numba
import numpy as np

from infercnvpy._numba import njit_parallel, prange


@numba.njit(cache=True, nogil=True)
def _lfc(v, lo, hi, cap):
    """Clipped log fold change of `v` against the reference bounds `[lo, hi]` (`lo == hi` for one reference)."""
    if v > hi:
        d = v - hi
    elif v < lo:
        d = v - lo
    else:
        return np.float32(0)
    if d > cap:
        return cap
    if d < -cap:
        return -cap
    return d


@numba.njit(cache=True, nogil=True)
def _select(a, n, k):
    """Partially sort `a[:n]` in place so that `a[k]` is its k-th smallest value, and return that (quickselect)."""
    lo = 0
    hi = n - 1
    while hi > lo:
        mid = (lo + hi) // 2
        if a[mid] < a[lo]:
            a[mid], a[lo] = a[lo], a[mid]
        if a[hi] < a[lo]:
            a[hi], a[lo] = a[lo], a[hi]
        if a[hi] < a[mid]:
            a[hi], a[mid] = a[mid], a[hi]
        pivot = a[mid]
        i = lo
        j = hi
        while i <= j:
            while a[i] < pivot:
                i += 1
            while a[j] > pivot:
                j -= 1
            if i <= j:
                a[i], a[j] = a[j], a[i]
                i += 1
                j -= 1
        if k <= j:
            hi = j
        elif k >= i:
            lo = i
        else:
            break
    return a[k]


@numba.njit(cache=True, nogil=True)
def _median(values, buf):
    """Median of `values` (as `np.median`), using `buf` (same size) as scratch."""
    n = values.size
    for i in range(n):
        buf[i] = values[i]
    k = n // 2
    upper = _select(buf, n, k)
    if n % 2:
        return upper
    lower = buf[0]
    for i in range(1, k):  # after the selection, buf[:k] holds the k smallest values
        if buf[i] > lower:
            lower = buf[i]
    return (lower + upper) / 2


@numba.njit(cache=True, nogil=True)
def _smooth(x, p1, p2, win, g0s, gs, As, Bs, Ks, w0s, step):
    """Running means of the kept windows of one cell (gene-ordered `x`) into `win`, from prefix sums."""
    for c in range(g0s.size):
        g0 = g0s[c]
        g = gs[c]
        a = As[c]
        b = Bs[c]
        w0 = w0s[c]
        p1[0] = 0.0
        for i in range(g):
            p1[i + 1] = p1[i] + x[g0 + i]
        p2[0] = 0.0
        for j in range(g - a + 1):
            p2[j + 1] = p2[j] + (p1[j + a] - p1[j])
        weight = a * b
        for w in range(Ks[c]):
            s = w * step
            win[w0 + w] = (p2[s + b] - p2[s]) / weight


@numba.njit(cache=True, nogil=True)
def _scatter(win, j, d, first, weights):
    """Add `d` times the pyramid weights of gene `j` to the weighted window sums `win`, from window `first[j]` on.

    `weights[j]` is 0 past the windows that contain gene `j`, and `win` has room for them past its last window.
    """
    w = first[j]
    wj = weights[j]
    for m in range(wj.size):
        win[w + m] += d * wj[m]


@numba.njit(cache=True, nogil=True)
def _normalise(win, n_windows, window_weight):
    """Divide the weighted window sums by the sums of their weights: the running means."""
    for w in range(n_windows):
        win[w] /= window_weight[w]


@numba.njit(cache=True, nogil=True)
def _center_emit(win, buf, thr, centre, out, r, rowstats, counts):
    """Subtract the median from the running means `win`, write them to `out[r]` as float32.

    Also `rowstats[r]`: sum and sum of squares of the written values, and `counts[r]`: the number of them outside the
    noise band `centre ± thr`. `out` may have no rows: then only these are computed.
    """
    med = _median(win, buf)
    total = 0.0
    squares = 0.0
    kept = 0
    write = out.shape[0] > 0
    for w in range(win.size):
        v = np.float32(win[w] - med)
        if write:
            out[r, w] = v
        total += v
        squares += np.float64(v) * v
        if v != 0 and abs(np.float64(v) - centre) >= thr:
            kept += 1
    rowstats[r, 0] = total
    rowstats[r, 1] = squares
    counts[r] = kept


@njit_parallel
def _cnv_csr(
    data, indices, indptr, rows, col_map, lo, hi, base, cap, g0s, gs, As, Bs, Ks, w0s, step,
    first, weights, window_weight, base_sums, max_nnz, thr, centre, out, rowstats, counts, n_chunks,
):  # fmt: skip
    """Centred windows of the CSR rows `rows` into `out`, see `_center_emit`.

    A row with at most `max_nnz` stored values starts from `base_sums`, the weighted window sums of the baseline `base`,
    and adds the difference to the baseline of each stored value to the windows that contain it (`_scatter`). Other
    rows are smoothed with prefix sums. `n_chunks` is the number of row chunks the threads share (see
    `n_kernel_chunks`).
    """
    n_rows = rows.size
    n_windows = out.shape[1]
    max_g = gs.max() + 1
    x = np.empty((n_chunks, base.size), np.float32)
    p1 = np.empty((n_chunks, max_g))
    p2 = np.empty((n_chunks, max_g))
    win = np.empty((n_chunks, n_windows + weights.shape[1]))
    buf = np.empty((n_chunks, n_windows))
    for t in prange(n_chunks):
        xt = x[t]
        wt = win[t]
        for r in range(t * n_rows // n_chunks, (t + 1) * n_rows // n_chunks):
            k0 = indptr[rows[r]]
            k1 = indptr[rows[r] + 1]
            if k1 - k0 <= max_nnz:
                wt[:n_windows] = base_sums
                for k in range(k0, k1):
                    j = col_map[indices[k]]
                    if j >= 0:
                        d = np.float64(_lfc(np.float32(data[k]), lo[j], hi[j], cap)) - base[j]
                        if d != 0:
                            _scatter(wt, j, d, first, weights)
                _normalise(wt, n_windows, window_weight)
            else:
                xt[:] = base
                for k in range(k0, k1):
                    j = col_map[indices[k]]
                    if j >= 0:
                        xt[j] = _lfc(np.float32(data[k]), lo[j], hi[j], cap)
                _smooth(xt, p1[t], p2[t], wt, g0s, gs, As, Bs, Ks, w0s, step)
            _center_emit(wt[:n_windows], buf[t], thr, centre, out, r, rowstats, counts)


@njit_parallel
def _cnv_dense(
    X, rows, col_map, lo, hi, base, cap, g0s, gs, As, Bs, Ks, w0s, step,
    first, weights, window_weight, base_sums, max_nnz, thr, centre, out, rowstats, counts, n_chunks,
):  # fmt: skip
    """Same as `_cnv_csr` for the rows `rows` of the dense matrix `X`; values that differ from the baseline count."""
    n_rows = rows.size
    n_windows = out.shape[1]
    max_g = gs.max() + 1
    x = np.empty((n_chunks, base.size), np.float32)
    p1 = np.empty((n_chunks, max_g))
    p2 = np.empty((n_chunks, max_g))
    win = np.empty((n_chunks, n_windows + weights.shape[1]))
    buf = np.empty((n_chunks, n_windows))
    pos = np.empty((n_chunks, max_nnz + 1), np.int64)
    diff = np.empty((n_chunks, max_nnz + 1))
    n_cols = X.shape[1]
    for t in prange(n_chunks):
        xt = x[t]
        wt = win[t]
        pt = pos[t]
        dt = diff[t]
        for r in range(t * n_rows // n_chunks, (t + 1) * n_rows // n_chunks):
            m = 0
            for i in range(n_cols):
                j = col_map[i]
                if j >= 0:
                    v = _lfc(np.float32(X[rows[r], i]), lo[j], hi[j], cap)
                    xt[j] = v
                    d = np.float64(v) - base[j]
                    if d != 0 and m <= max_nnz:
                        pt[m] = j
                        dt[m] = d
                        m += 1
            if m <= max_nnz:
                wt[:n_windows] = base_sums
                for q in range(m):
                    _scatter(wt, pt[q], dt[q], first, weights)
                _normalise(wt, n_windows, window_weight)
            else:
                _smooth(xt, p1[t], p2[t], wt, g0s, gs, As, Bs, Ks, w0s, step)
            _center_emit(wt[:n_windows], buf[t], thr, centre, out, r, rowstats, counts)


@njit_parallel
def _count_kept(out, thr, centre, counts):
    """`counts[r]`: entries of row `r` that are non-zero and outside the noise band `centre ± thr`."""
    n_rows, n_windows = out.shape
    for r in prange(n_rows):
        c = 0
        for w in range(n_windows):
            v = out[r, w]
            if v != 0 and abs(np.float64(v) - centre) >= thr:
                c += 1
        counts[r] = c


@njit_parallel
def _fill_kept(out, thr, centre, indptr, data, indices):
    """Fill the CSR arrays (`indptr` from `_count_kept`) with the entries outside the noise band."""
    n_rows, n_windows = out.shape
    for r in prange(n_rows):
        k = indptr[r]
        for w in range(n_windows):
            v = out[r, w]
            if v != 0 and abs(np.float64(v) - centre) >= thr:
                data[k] = v
                indices[k] = w
                k += 1


@njit_parallel
def _gene_values(out, gs, ns, Ks, w0s, g0s, step, order, n_covered, thr, centre, gv, n_chunks):
    """Per gene: mean of the kept windows that contain it, median-centred over all covered genes, denoised.

    `out` holds the centred windows of a block of cells; `ns` the window length per chromosome (the number of genes
    for a chromosome that is a single window). Values are written into `gv` (cells x all genes of `adata`, prefilled
    with NaN) at the column `order[gene]`. A window prefix sum makes this O(1) per gene.
    """
    n_rows = out.shape[0]
    pw = np.empty((n_chunks, Ks.max() + 1))
    vals = np.empty((n_chunks, n_covered))
    buf = np.empty((n_chunks, n_covered))
    for t in prange(n_chunks):
        pwt = pw[t]
        vt = vals[t]
        for r in range(t * n_rows // n_chunks, (t + 1) * n_rows // n_chunks):
            q = 0
            for c in range(gs.size):
                K = Ks[c]
                n = ns[c]
                pwt[0] = 0.0
                for w in range(K):
                    pwt[w + 1] = pwt[w] + out[r, w0s[c] + w]
                for i in range(min(gs[c], (K - 1) * step + n)):
                    first = 0 if i < n else (i - n) // step + 1  # first window that contains gene i
                    last = min(K - 1, i // step)
                    vt[q] = (pwt[last + 1] - pwt[first]) / (last - first + 1)
                    q += 1
            med = _median(vt, buf[t])
            q = 0
            for c in range(gs.size):
                for i in range(min(gs[c], (Ks[c] - 1) * step + ns[c])):
                    v = vt[q] - med
                    q += 1
                    gv[r, order[g0s[c] + i]] = 0 if abs(v - centre) < thr else v


@njit_parallel
def _group_sums_csr(data, indices, indptr, rows, groups, sums):
    """`sums[t, groups[i]] += X[rows[i]]` for the CSR matrix `X`, where `t` is the chunk of `i` (one per `sums[t]`)."""
    n_chunks = sums.shape[0]
    n = rows.size
    for t in prange(n_chunks):
        st = sums[t]
        for i in range(t * n // n_chunks, (t + 1) * n // n_chunks):
            s = st[groups[i]]
            r = rows[i]
            for k in range(indptr[r], indptr[r + 1]):
                s[indices[k]] += data[k]


@njit_parallel
def _group_sums_dense(X, rows, groups, sums):
    """Same as `_group_sums_csr` for the dense matrix `X`."""
    n_chunks = sums.shape[0]
    n = rows.size
    for t in prange(n_chunks):
        st = sums[t]
        for i in range(t * n // n_chunks, (t + 1) * n // n_chunks):
            s = st[groups[i]]
            r = rows[i]
            for j in range(X.shape[1]):
                s[j] += X[r, j]


@njit_parallel
def _n_unsorted_rows(indices, indptr):
    """Number of CSR rows whose column indices are not strictly increasing (unsorted or duplicate)."""
    n = 0
    for r in prange(indptr.size - 1):
        for k in range(indptr[r] + 1, indptr[r + 1]):
            if indices[k] <= indices[k - 1]:
                n += 1
                break
    return n
