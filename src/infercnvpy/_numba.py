"""numba helpers: a parallel kernel with a serial twin, and `n_jobs` as the number of numba threads.

A ``parallel=True`` kernel only runs on the main thread. Called from any other thread (a dask worker, a thread pool),
the serial twin runs instead: numba's default ``workqueue`` threading layer (the only one in the macOS wheels) aborts
the process on concurrent use, and a caller that brings its own threads should not get numba threads on top.
"""

from __future__ import annotations

import threading
from contextlib import contextmanager
from functools import update_wrapper
from types import FunctionType
from typing import TYPE_CHECKING, Any

import numba

if TYPE_CHECKING:
    from builtins import range as prange  # numba's prange iterates like range; its own type says nothing about that
    from collections.abc import Callable, Generator
else:
    from numba import prange

__all__ = ["n_kernel_chunks", "njit_parallel", "numba_threads", "prange"]


class njit_parallel:
    """Compile a function with numba twice, ``parallel=True`` and serial, both ``cache=True`` and ``nogil=True``.

    Calling it runs the parallel version on the main thread and the serial one on any other thread, or always the
    serial one with ``serial=True`` (for code that runs inside dask tasks). The two versions are distinct function
    objects with distinct names, so they don't share a cache entry.
    """

    def __init__(self, fn: Callable[..., Any], /):
        assert isinstance(fn, FunctionType)
        self.parallel = numba.njit(_renamed(fn, "parallel"), parallel=True, cache=True, nogil=True)
        self.serial = numba.njit(_renamed(fn, "serial"), cache=True, nogil=True)
        update_wrapper(self, fn)

    def __call__(self, *args, serial: bool = False):
        return (self.serial if serial or not _on_main_thread() else self.parallel)(*args)


def _on_main_thread() -> bool:
    return threading.current_thread() is threading.main_thread()


def _renamed(fn: FunctionType, suffix: str) -> FunctionType:
    new = FunctionType(fn.__code__, fn.__globals__, fn.__name__, fn.__defaults__, fn.__closure__)
    update_wrapper(new, fn)
    new.__qualname__ = f"{fn.__qualname__}-{suffix}"
    return new


def n_kernel_chunks(n_rows: int, *, serial: bool = False, per_thread: int = 4) -> int:
    """Number of row chunks for a `prange` kernel: `per_thread` per numba thread, or 1 if it runs serially."""
    if serial or not _on_main_thread():
        return 1
    return max(1, min(n_rows, per_thread * numba.get_num_threads()))


@contextmanager
def numba_threads(n_jobs: int | None) -> Generator[None]:
    """Run the block with `n_jobs` numba threads and restore the previous number afterwards.

    `None` keeps numba's setting (all cores unless ``NUMBA_NUM_THREADS`` says otherwise), `-1` means all
    ``NUMBA_NUM_THREADS``, and larger values are capped at that.
    """
    if n_jobs is None:
        yield
        return
    if n_jobs == 0 or n_jobs < -1:
        raise ValueError(f"`n_jobs` must be a positive number of threads, -1 (all) or None, got {n_jobs}.")
    n_max = numba.config.NUMBA_NUM_THREADS  # type: ignore[attr-defined]
    previous = numba.get_num_threads()
    numba.set_num_threads(n_max if n_jobs == -1 else min(n_jobs, n_max))
    try:
        yield
    finally:
        numba.set_num_threads(previous)
