"""Deprecated alias of :mod:`marse.microbes.niche`, kept for one release."""

import warnings

from marse.microbes.niche import *  # noqa: F403
from marse.microbes.niche import __all__  # noqa: F401

warnings.warn(
    "marse.niche has moved to marse.microbes.niche; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
