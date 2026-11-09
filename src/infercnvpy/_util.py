import numpy as np


def _ensure_array(a):
    """If a is a matrix, turn it into an array."""
    if isinstance(a, np.matrix):
        return a.A
    else:
        return a


def _choose_mtx_rep(adata, use_raw=False, layer=None):
    """Get gene expression from anndata depending on use_raw and layer."""
    is_layer = layer is not None
    if use_raw and is_layer:
        raise ValueError(
            f"Cannot use expression from both layer and raw. You provided:'use_raw={use_raw}' and 'layer={layer}'"
        )
    if is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    else:
        return adata.X


def _fetch_rows(X, start: int, stop: int):
    """Materialise rows `start:stop` of `X` as an in-memory `ndarray` or scipy sparse matrix.

    Works for in-memory arrays, scipy sparse matrices, dask arrays (also with sparse blocks, e.g. built with
    `da.from_delayed`, as in lazily loaded spatialdata-io tables) and backed AnnData matrices
    (`h5py.Dataset`, `anndata` sparse datasets). Only the requested rows are read.
    """
    import scipy.sparse as sp

    rows = X[start:stop]
    if hasattr(rows, "compute"):  # dask
        rows = rows.compute()
    if hasattr(rows, "to_scipy_sparse"):  # pydata/sparse COO
        rows = rows.to_scipy_sparse()
    if sp.issparse(rows):
        return rows.tocsr()
    return np.asarray(rows)


def _is_lazy(X) -> bool:
    """True for matrices that are not in memory: dask arrays and backed (on-disk) matrices."""
    import scipy.sparse as sp

    return not isinstance(X, np.ndarray) and not sp.issparse(X)


def _read_units(X, chunksize: int) -> list[tuple[int, int]]:
    """Consecutive row ranges `(start, stop)` of `X` that can be read efficiently, one read each.

    For in-memory and backed matrices these are `chunksize` rows. For dask arrays a read always loads whole
    blocks, so units are made of whole blocks (several small ones are grouped up to `chunksize` rows; a block larger
    than `chunksize` is a unit of its own): reading `chunksize` rows out of a larger block would otherwise load, and
    discard, that entire block again for every chunk.
    """
    n = X.shape[0]
    if hasattr(X, "compute") and getattr(X, "chunks", None) is not None:
        sizes = [int(c) for c in X.chunks[0]]
        bounds = np.concatenate([[0], np.cumsum(sizes)])
        units, i = [], 0
        while i < len(sizes):
            j = i + 1
            while j < len(sizes) and bounds[j + 1] - bounds[i] <= chunksize:
                j += 1
            if bounds[j] > bounds[i]:
                units.append((int(bounds[i]), int(bounds[j])))
            i = j
        return units
    return [(s, min(s + chunksize, n)) for s in range(0, n, chunksize)]


def _group_codes(values) -> tuple[np.ndarray, list]:
    """Integer code per observation and the groups in order of appearance. Missing values get the code -1."""
    import pandas as pd

    if hasattr(values, "to_pandas"):  # lazy (xarray) column of `anndata.experimental.read_lazy`
        values = values.to_pandas()
    codes, uniques = pd.factorize(values)
    return codes.astype(np.int64, copy=False), list(uniques)


def _group_members(codes: np.ndarray, n_groups: int) -> list[np.ndarray]:
    """Sorted row indices of every group, in O(n) instead of one boolean mask per group."""
    valid = np.flatnonzero(codes >= 0)
    order = valid[np.argsort(codes[valid], kind="stable")]
    bounds = np.concatenate([[0], np.cumsum(np.bincount(codes[valid], minlength=n_groups))])
    return [order[bounds[i] : bounds[i + 1]] for i in range(n_groups)]


def _group_means(X, codes: np.ndarray, n_groups: int) -> np.ndarray:
    """Mean over the rows of `X` of each group (`n_groups x n_cols`).

    The rows of every group are selected with precomputed indices (instead of a boolean mask over all rows per
    group), which makes this independent of the number of groups for sparse matrices; dense matrices use
    a sparse indicator matrix product.
    """
    import scipy.sparse as sp

    counts = np.bincount(codes[codes >= 0], minlength=n_groups)
    if sp.issparse(X):
        X = X.tocsr()  # row selection is only cheap for CSR
        sums = np.zeros((n_groups, X.shape[1]))
        for i, members in enumerate(_group_members(codes, n_groups)):
            if len(members):
                sums[i] = np.asarray(X[members].sum(axis=0)).ravel()
    else:
        valid = np.flatnonzero(codes >= 0)
        indicator = sp.csr_matrix((np.ones(len(valid)), (codes[valid], valid)), shape=(n_groups, X.shape[0]))
        sums = np.asarray(indicator @ np.asarray(X))
    with np.errstate(invalid="ignore", divide="ignore"):
        return sums / counts[:, None]


def _subsample(idx: np.ndarray, max_n: int | None, rng: np.random.Generator) -> np.ndarray:
    """At most `max_n` of the (sorted) indices `idx`, drawn without replacement and returned sorted."""
    if max_n is None or len(idx) <= max_n:
        return idx
    return np.sort(rng.choice(idx, size=max_n, replace=False))


def _take_rows(X, idx: np.ndarray):
    """Rows `idx` (sorted) of `X` as in-memory array/sparse matrix; also for dask and backed matrices."""
    import scipy.sparse as sp

    rows = X[idx]
    if hasattr(rows, "compute"):
        rows = rows.compute()
    if hasattr(rows, "to_scipy_sparse"):
        rows = rows.to_scipy_sparse()
    return rows.tocsr() if sp.issparse(rows) else np.asarray(rows)
