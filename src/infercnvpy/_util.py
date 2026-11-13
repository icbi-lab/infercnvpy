import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _as_matrix(x, what: str = "the expression matrix", *, lazy: bool = False):
    """Narrow an AnnData array to an in-memory numpy array or scipy sparse matrix.

    With `lazy`, also accept a dask array or a backed (on-disk) matrix, see `_is_dask` and `_is_backed`.
    """
    if isinstance(x, np.ndarray) or sp.issparse(x) or (lazy and (_is_dask(x) or _is_backed(x))):
        return x
    kinds = "a numpy array or scipy sparse matrix in memory"
    if lazy:
        kinds = "a numpy array, scipy sparse matrix, dask array or backed matrix"
    raise TypeError(f"infercnvpy needs {what} as {kinds}, got {type(x).__name__}.")


def _is_dask(x) -> bool:
    """Whether `x` is a dask array (without importing dask)."""
    da = sys.modules.get("dask.array")
    return da is not None and isinstance(x, da.Array)


def _is_backed(x) -> bool:
    """Whether `x` is the matrix of a backed AnnData (an h5py dataset or an anndata sparse dataset)."""
    import h5py
    from anndata.abc import CSCDataset, CSRDataset

    return isinstance(x, h5py.Dataset | CSRDataset | CSCDataset)


def _as_df(df) -> pd.DataFrame:
    """Narrow `obs`/`var` to a pandas DataFrame (AnnData also allows lazy `Dataset2D`)."""
    if isinstance(df, pd.DataFrame):
        return df
    raise TypeError(f"infercnvpy needs `obs` and `var` as pandas DataFrames, got {type(df).__name__}.")


def _params_for_uns(**params) -> dict:
    """Prepare call parameters for `adata.uns[...]["params"]`: store sequences as lists."""
    return {k: list(v) if isinstance(v, tuple | list) else v for k, v in params.items()}
