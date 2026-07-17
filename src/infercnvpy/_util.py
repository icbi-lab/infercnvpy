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
