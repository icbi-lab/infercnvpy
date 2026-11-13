"""numba kernels of :func:`infercnvpy.tl.infercnv`, one row block of cells at a time.

Per cell, `_cnv_csr` / `_cnv_dense`:

1. fill a scratch row in genomic gene order with the clipped log fold change of an expression of 0, then scatter the
   cell's values into it through `col_map` (input column -> gene position, -1 for dropped genes), so the input is
   never densified or copied;
2. compute the pyramid-weighted running mean of every kept window per chromosome. A pyramid of `n` genes is the
   convolution of two boxes of `a = ceil(n / 2)` and `b = n + 1 - a` genes, so a window is the difference of two
   entries of a prefix sum of a prefix sum: O(1) per kept window, whatever `window_size` is;
3. subtract the median of the cell and write the windows as float32, with the float64 sum and sum of squares of the
   written values for the noise threshold.

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
def _smooth_center_emit(x, p1, p2, win, buf, g0s, gs, As, Bs, Ks, w0s, step, out, r, rowstats):
    """Running means of the kept windows of one cell (gene-ordered `x`), median-centred, written to `out[r]`."""
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
    med = _median(win, buf)
    total = 0.0
    squares = 0.0
    for w in range(win.size):
        v = np.float32(win[w] - med)
        out[r, w] = v
        total += v
        squares += np.float64(v) * v
    rowstats[r, 0] = total
    rowstats[r, 1] = squares


@njit_parallel
def _cnv_csr(
    data, indices, indptr, row0, col_map, lo, hi, base, cap, g0s, gs, As, Bs, Ks, w0s, step, out, rowstats, n_chunks
):
    """Centred windows of the CSR rows `row0 : row0 + out.shape[0]` into `out`, their sum/sum of squares into `rowstats`.

    `n_chunks` is the number of row chunks the threads share (see `n_kernel_chunks`).
    """
    n_rows, n_windows = out.shape
    max_g = gs.max() + 1
    x = np.empty((n_chunks, base.size), np.float32)
    p1 = np.empty((n_chunks, max_g))
    p2 = np.empty((n_chunks, max_g))
    win = np.empty((n_chunks, n_windows))
    buf = np.empty((n_chunks, n_windows))
    for t in prange(n_chunks):
        xt = x[t]
        for r in range(t * n_rows // n_chunks, (t + 1) * n_rows // n_chunks):
            xt[:] = base
            for k in range(indptr[row0 + r], indptr[row0 + r + 1]):
                j = col_map[indices[k]]
                if j >= 0:
                    xt[j] = _lfc(np.float32(data[k]), lo[j], hi[j], cap)
            _smooth_center_emit(xt, p1[t], p2[t], win[t], buf[t], g0s, gs, As, Bs, Ks, w0s, step, out, r, rowstats)


@njit_parallel
def _cnv_dense(X, row0, col_map, lo, hi, cap, g0s, gs, As, Bs, Ks, w0s, step, out, rowstats, n_chunks):
    """Same as `_cnv_csr` for the rows `row0 : row0 + out.shape[0]` of the dense matrix `X`."""
    n_rows, n_windows = out.shape
    max_g = gs.max() + 1
    x = np.empty((n_chunks, lo.size), np.float32)
    p1 = np.empty((n_chunks, max_g))
    p2 = np.empty((n_chunks, max_g))
    win = np.empty((n_chunks, n_windows))
    buf = np.empty((n_chunks, n_windows))
    n_cols = X.shape[1]
    for t in prange(n_chunks):
        xt = x[t]
        for r in range(t * n_rows // n_chunks, (t + 1) * n_rows // n_chunks):
            for i in range(n_cols):
                j = col_map[i]
                if j >= 0:
                    xt[j] = _lfc(np.float32(X[row0 + r, i]), lo[j], hi[j], cap)
            _smooth_center_emit(xt, p1[t], p2[t], win[t], buf[t], g0s, gs, As, Bs, Ks, w0s, step, out, r, rowstats)


@njit_parallel
def _count_kept(out, thr, centre, counts):
    """`counts[r + 1]`: entries of row `r` that are non-zero and outside the noise band `centre ± thr`."""
    n_rows, n_windows = out.shape
    for r in prange(n_rows):
        c = 0
        for w in range(n_windows):
            v = out[r, w]
            if v != 0 and abs(np.float64(v) - centre) >= thr:
                c += 1
        counts[r + 1] = c


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
