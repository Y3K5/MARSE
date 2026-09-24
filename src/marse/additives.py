"""Continuous small-molecule fields and phenomenological dose responses."""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from numpy.typing import ArrayLike, NDArray

__all__ = ["AdditiveEffect", "hill_response", "apply_effect"]


def hill_response(concentration: ArrayLike, half_effect: float, coefficient: float) -> NDArray:
    """Return a bounded Hill occupancy for a non-negative concentration."""
    if half_effect <= 0 or coefficient <= 0:
        raise ValueError("half_effect and coefficient must be positive")
    value = np.maximum(np.asarray(concentration, dtype=float), 0.0)
    return value**coefficient / (half_effect**coefficient + value**coefficient)


@dataclass(frozen=True, slots=True)
class AdditiveEffect:
    """Species response to one additive, from inhibition to stimulation."""

    additive: str
    minimum_multiplier: float
    maximum_multiplier: float
    half_effect: float
    coefficient: float = 1.0
    direction: str = "increasing"

    def __post_init__(self) -> None:
        if not self.additive.strip():
            raise ValueError("additive must not be empty")
        if self.minimum_multiplier < 0 or self.maximum_multiplier < 0:
            raise ValueError("effect multipliers must be non-negative")
        if self.half_effect <= 0 or self.coefficient <= 0:
            raise ValueError("half_effect and coefficient must be positive")
        if self.direction not in {"increasing", "decreasing"}:
            raise ValueError("direction must be 'increasing' or 'decreasing'")


def apply_effect(effect: AdditiveEffect, concentration: ArrayLike) -> NDArray:
    """Evaluate an additive effect as a multiplicative phenotype factor."""
    occupancy = hill_response(concentration, effect.half_effect, effect.coefficient)
    if effect.direction == "decreasing":
        occupancy = 1.0 - occupancy
    return effect.minimum_multiplier + (
        effect.maximum_multiplier - effect.minimum_multiplier
    ) * occupancy
