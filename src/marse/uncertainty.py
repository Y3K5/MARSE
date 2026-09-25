"""Deprecated alias of :mod:`marse.analysis.uncertainty`, kept for one release."""

import warnings

from marse.analysis.uncertainty import *  # noqa: F403
from marse.analysis.uncertainty import __all__  # noqa: F401

warnings.warn(
    "marse.uncertainty has moved to marse.analysis.uncertainty; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
