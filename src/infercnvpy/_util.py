import numpy as np
import pandas as pd
import scipy.sparse as sp


def _ensure_array(a):
    """If a is a matrix, turn it into an array."""
    if isinstance(a, np.matrix):
        return a.A
    else:
        return a


def _as_matrix(x, what: str = "the expression matrix") -> np.ndarray | sp.spmatrix | sp.sparray:
    """Narrow an AnnData array to an in-memory numpy array or scipy sparse matrix."""
    if isinstance(x, np.ndarray) or sp.issparse(x):
        return x
    raise TypeError(
        f"infercnvpy needs {what} as a numpy array or scipy sparse matrix in memory, got {type(x).__name__}."
    )


def _as_df(df) -> pd.DataFrame:
    """Narrow `obs`/`var` to a pandas DataFrame (AnnData also allows lazy `Dataset2D`)."""
    if isinstance(df, pd.DataFrame):
        return df
    raise TypeError(f"infercnvpy needs `obs` and `var` as pandas DataFrames, got {type(df).__name__}.")


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


def _params_for_uns(**params) -> dict:
    """Prepare call parameters for `adata.uns[...]["params"]`: drop `None`, store sequences as lists."""
    return {k: list(v) if isinstance(v, tuple | list) else v for k, v in params.items() if v is not None}
