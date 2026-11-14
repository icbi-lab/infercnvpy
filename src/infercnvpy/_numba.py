"""numba helpers: a parallel kernel with a serial twin, and `n_jobs` as the number of numba threads.

A ``parallel=True`` kernel only runs on the main thread. Called from any other thread (a dask worker, a thread pool),
the serial twin runs instead: numba's default ``workqueue`` threading layer (the only one in the macOS wheels) aborts
the process on concurrent use, and a caller that brings its own threads should not get numba threads on top.

All kernels use numba's disk cache, with data files named after their cache key (see `_KeyedCacheFile`).
"""

from __future__ import annotations

import hashlib
import threading
from contextlib import contextmanager
from functools import update_wrapper
from types import FunctionType
from typing import TYPE_CHECKING, Any

import numba
from numba.core.caching import FunctionCache, IndexDataCacheFile

if TYPE_CHECKING:
    from builtins import range as prange  # numba's prange iterates like range; its own type says nothing about that
    from collections.abc import Callable, Generator
else:
    from numba import prange

__all__ = ["n_kernel_chunks", "njit", "njit_parallel", "numba_threads", "prange"]


class _KeyedCacheFile(IndexDataCacheFile):
    """numba's cache index and data files, with each data file named after its key (signature, CPU, bytecode).

    numba numbers the data files of a function (``.1.nbc``, ``.2.nbc``, ...) by the free slots in the index as the
    process read it. Processes that compile different signatures of one function at the same time (pytest-xdist
    workers, a pool of jobs on one environment) read the same index and pick the same number. The index can then end
    up mapping one signature to another signature's machine code, and every later call with it fails with "can't
    unbox array from PyObject" until the cache is deleted. With a name per key, a data file only ever holds the code
    for its own key. Concurrent index writes can still drop an entry, which only costs one recompilation.
    """

    def save(self, key, data):
        digest = hashlib.sha256(repr(key).encode()).hexdigest()[:16]
        data_name = self._data_name_pattern.replace("{number:d}", digest)
        self._save_data(data_name, data)  # before the index, so that no index entry points to a missing file
        overloads = self._load_index()
        if overloads.get(key) != data_name:
            overloads[key] = data_name
            self._save_index(overloads)


class _KeyedFunctionCache(FunctionCache):
    def __init__(self, py_func):
        super().__init__(py_func)
        self._cache_file = _KeyedCacheFile(
            self._cache_path, self._impl.filename_base, self._impl.locator.get_source_stamp()
        )


def njit(fn: Callable[..., Any], /, **options: Any) -> Any:
    """`numba.njit` with ``nogil=True`` and the disk cache of `_KeyedCacheFile`."""
    dispatcher = numba.njit(fn, nogil=True, **options)
    dispatcher._cache = _KeyedFunctionCache(dispatcher.py_func)  # type: ignore[attr-defined]  # what `cache=True` sets
    return dispatcher


class njit_parallel:
    """Compile a function with `njit` twice, ``parallel=True`` and serial.

    Calling it runs the parallel version on the main thread and the serial one on any other thread, or always the
    serial one with ``serial=True`` (for code that runs inside dask tasks). The two versions are distinct function
    objects with distinct names, so they don't share a cache entry.
    """

    def __init__(self, fn: Callable[..., Any], /):
        assert isinstance(fn, FunctionType)
        self.parallel = njit(_renamed(fn, "parallel"), parallel=True)
        self.serial = njit(_renamed(fn, "serial"))
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
