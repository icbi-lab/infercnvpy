"""The numba disk cache stays correct when processes compile different signatures of one kernel at the same time."""

import os
import subprocess
import sys
import textwrap

KERNEL = """
from infercnvpy._numba import njit_parallel

@njit_parallel
def total(x):
    s = 0.0
    for v in x:
        s += v
    return s
"""

# Two processes A and B compile `total` for float32 and float64 at the same time, in this order: A writes the index,
# B (which read the index before A wrote it) writes the index and its data, then A's data write lands.
# One process plays both, with two copies of the kernel module and numba's file access patched to that order.
RACE = """
import importlib.util
import sys

import numpy as np
from numba.core.caching import IndexDataCacheFile


def kernel():
    spec = importlib.util.spec_from_file_location("kernels", sys.argv[1])
    module = sys.modules["kernels"] = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.total


a, b = kernel(), kernel()
save_data, load_index = IndexDataCacheFile._save_data, IndexDataCacheFile._load_index
late = []
IndexDataCacheFile._save_data = lambda self, name, data: late.append((self, name, data))
a(np.ones(3, dtype=np.float32))
IndexDataCacheFile._save_data = save_data
IndexDataCacheFile._load_index = lambda self: {}
b(np.ones(3))
IndexDataCacheFile._load_index = load_index
for args in late:
    save_data(*args)
"""

CALL = """
import numpy as np

import kernels as module

assert module.total(np.ones(3, dtype=np.float32)) == 3
assert module.total(np.ones(3)) == 3
"""


def test_concurrent_compiles_share_the_cache(tmp_path):
    kernels = tmp_path / "kernels.py"
    kernels.write_text(KERNEL)
    env = {**os.environ, "NUMBA_CACHE_DIR": str(tmp_path / "cache")}
    for script in (RACE, CALL):
        res = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(script), str(kernels)],
            env=env,
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0, res.stderr
