"""infercnvpy's logger, the ``scverse.infercnvpy`` child of the shared scverse logger.

One setting controls infercnvpy and every other package on the shared logger::

    from scverse_misc.logging import config

    config.verbosity = "info"  # default "warning"
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from scverse_misc.logging import get_logger

if TYPE_CHECKING:
    from datetime import datetime

log = get_logger("infercnvpy", timed=True)

#: ``warnings.warn(..., skip_file_prefixes=SKIP_OWN_FRAMES)`` points at the caller's line, not ours.
SKIP_OWN_FRAMES = (os.path.dirname(__file__) + os.sep,)


def log_result(
    func: str,
    msg: str,
    *args: object,
    dropped: int = 0,
    total: int = 0,
    remedy: str | None = None,
    time: datetime | None = None,
) -> None:
    """Log the one info line of `func` with its result, as a warning when half or more of `total` were dropped.

    `msg` is a %-format string for `args`; it should mention the `dropped` count when there is one.
    The `remedy` goes to ``deep=``, so it shows only when the verbosity is more detailed than the line's level.
    """
    heavy = dropped > 0 and 2 * dropped >= total
    emit = log.warning if heavy else log.info
    emit("%s: " + msg, func, *args, deep=remedy if dropped else None, time=time)
