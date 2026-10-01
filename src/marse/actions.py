"""Deprecated alias of :mod:`marse.experimental.host.actions`, kept for one release."""

import warnings

from marse.experimental.host.actions import *  # noqa: F403
from marse.experimental.host.actions import __all__  # noqa: F401

warnings.warn(
    "marse.actions has moved to marse.experimental.host.actions; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
