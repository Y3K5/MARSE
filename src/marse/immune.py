"""Deprecated alias of :mod:`marse.experimental.host.immune`, kept for one release."""

import warnings

from marse.experimental.host.immune import *  # noqa: F403
from marse.experimental.host.immune import __all__  # noqa: F401

warnings.warn(
    "marse.immune has moved to marse.experimental.host.immune; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
