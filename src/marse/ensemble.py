"""Deprecated alias of :mod:`marse.analysis.ensemble`, kept for one release."""

import warnings

from marse.analysis.ensemble import *  # noqa: F403
from marse.analysis.ensemble import __all__  # noqa: F401

warnings.warn(
    "marse.ensemble has moved to marse.analysis.ensemble; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
