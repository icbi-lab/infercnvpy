import numpy as np
import scipy.ndimage
import scipy.sparse


def _running_mean(
    x: scipy.sparse.spmatrix, n: int = 50, chunksize: int = 5000, step: int = 10
) -> np.ndarray:
    """
    Compute the running mean along rows of an array.

    Densifies the matrix. Use `step` and `chunksize` to save memory.

    Parameters
    ----------
    n
        Length of the running window
    chunksize
        Process matrix in chunks of `chunksize` rows. Use this if your matrix
        doesn't fit into memory while dense.
    step
        only compute running windows ever `step` columns, e.g. if step is 10
        0:100, 10:110, 20:120 etc. Saves memory.
    """
    x = x.tocsr()
    return np.vstack(
        [
            scipy.ndimage.uniform_filter1d(
                x[i : i + chunksize, :].toarray(), size=n, axis=1, mode="nearest"
            )[:, np.arange(0, x.shape[1], step)]
            for i in range(0, x.shape[0], chunksize)
        ]
    )
