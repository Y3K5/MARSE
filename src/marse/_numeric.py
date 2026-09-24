"""Helpers shared by the numerical modules."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

type FloatOrArray = float | NDArray[np.float64]


def as_array(value: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(value, dtype=float)


def unwrap(array: NDArray[np.float64]) -> FloatOrArray:
    """Return a NumPy scalar for 0-d results and the array itself otherwise."""
    return array[()]


def require(condition: bool | np.bool_, message: str) -> None:
    if not condition:
        raise ValueError(message)
