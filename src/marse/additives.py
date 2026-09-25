"""Deprecated alias of :mod:`marse.microbes.additives`, kept for one release."""

import warnings

from marse.microbes.additives import *  # noqa: F403
from marse.microbes.additives import __all__  # noqa: F401

warnings.warn(
    "marse.additives has moved to marse.microbes.additives; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
