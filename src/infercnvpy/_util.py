import inspect
import sys
import warnings
from collections.abc import Callable
from functools import wraps
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from infercnvpy._logging import SKIP_OWN_FRAMES


def _as_matrix(x, what: str = "the expression matrix", *, lazy: bool = False):
    """Narrow an AnnData array to an in-memory numpy array or scipy sparse matrix.

    With `lazy`, also accept a dask array or a backed (on-disk) matrix, see `_is_dask` and `_is_backed`. Without it,
    such a matrix raises an error that tells the user how to load it.
    """
    if isinstance(x, np.ndarray) or sp.issparse(x) or (lazy and (_is_dask(x) or _is_backed(x))):
        return x
    if _is_dask(x) or _is_backed(x):
        kind = "a dask array" if _is_dask(x) else "a backed (on-disk) matrix"
        raise TypeError(
            f"This function needs {what} in memory, got {kind}. Load it first: call `.compute()` on a dask array "
            "and store the result back in adata (e.g. `adata.obsm['X_cnv'] = adata.obsm['X_cnv'].compute()`), or "
            "use `adata = adata.to_memory()` for a backed AnnData."
        )
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
    """Prepare call parameters for `adata.uns[...]["params"]`: store sequences as lists."""
    return {k: list(v) if isinstance(v, tuple | list) else v for k, v in params.items()}


def _deprecated_positional[F: Callable[..., Any]](*names: str, until: str = "0.8.0") -> Callable[[F], F]:
    """Accept the keyword-only arguments `names` positionally, in this order, with a FutureWarning.

    For arguments that became keyword-only: calls that pass them positionally (after the positional
    parameters of the new signature) keep working until `until`.
    """

    def decorate(fn: F) -> F:
        n_positional = sum(
            p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD) for p in inspect.signature(fn).parameters.values()
        )

        @wraps(fn)
        def wrapper(*args, **kwargs):
            if len(args) > n_positional:
                extra = args[n_positional:]
                if len(extra) > len(names):
                    raise TypeError(
                        f"{fn.__name__}() takes {n_positional} positional arguments (up to {n_positional + len(names)} "
                        f"with deprecated ones), but {len(args)} were given."
                    )
                given = names[: len(extra)]
                if twice := sorted(set(given) & kwargs.keys()):
                    raise TypeError(f"{fn.__name__}() got multiple values for {', '.join(twice)}.")
                warnings.warn(
                    f"Passing {', '.join(f'`{n}`' for n in given)} to `{fn.__name__}` positionally is deprecated and "
                    f"will raise an error from infercnvpy {until}. Pass "
                    f"{', '.join(f'`{n}=...`' for n in given)} as keyword argument{'s' * (len(given) > 1)} instead.",
                    FutureWarning,
                    skip_file_prefixes=SKIP_OWN_FRAMES,
                )
                kwargs.update(zip(given, extra, strict=True))
                args = args[:n_positional]
            return fn(*args, **kwargs)

        return wrapper  # type: ignore[return-value]

    return decorate


def _copy_for_inplace_false(func: str, inplace: bool, copy: bool, *, future: str) -> bool:
    """Whether to work on a copy of adata: `copy=True`, or the deprecated `inplace=False` (with a FutureWarning).

    `inplace=False` returned a modified copy of adata before 0.8.0; from 0.8.0 it returns the results (`future`), as in
    the other functions.
    """
    if copy and not inplace:
        raise ValueError(f"{func}: pass either `copy=True` or `inplace=False`, not both.")
    if not inplace:
        warnings.warn(
            f"`{func}(..., inplace=False)` returns a modified copy of adata. From infercnvpy 0.8.0 it will return "
            f"{future} instead, as the other functions do. Pass `copy=True` to keep getting a modified copy.",
            FutureWarning,
            skip_file_prefixes=SKIP_OWN_FRAMES,
        )
        return True
    return copy
