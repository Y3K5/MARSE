"""Deprecated alias of :mod:`marse.analysis.calibration`, kept for one release."""

import warnings

from marse.analysis.calibration import *  # noqa: F403
from marse.analysis.calibration import __all__  # noqa: F401

warnings.warn(
    "marse.calibration has moved to marse.analysis.calibration; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
