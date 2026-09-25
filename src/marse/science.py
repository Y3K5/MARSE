"""Deprecated alias of :mod:`marse.evidence.science`, kept for one release."""

import warnings

from marse.evidence.science import *  # noqa: F403
from marse.evidence.science import __all__  # noqa: F401

warnings.warn(
    "marse.science has moved to marse.evidence.science; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
