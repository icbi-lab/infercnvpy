import numpy as np

# Determine the right tqdm version, see https://github.com/tqdm/tqdm/issues/1082
try:
    import ipywidgets  # type: ignore # NOQA
    from tqdm.auto import tqdm
except ModuleNotFoundError:
    from tqdm import tqdm  # NOQA


def _ensure_array(a):
    """If a is a matrix, turn it into an array."""
    if isinstance(a, np.matrix):
        return a.A
    else:
        return a
