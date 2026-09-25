"""Deprecated alias of :mod:`marse.microbes.genotype`, kept for one release."""

import warnings

from marse.microbes.genotype import *  # noqa: F403
from marse.microbes.genotype import __all__  # noqa: F401

warnings.warn(
    "marse.genotype has moved to marse.microbes.genotype; "
    "the old import path will be removed in the next release",
    DeprecationWarning,
    stacklevel=2,
)
