"""Bounded immune and molecular interaction primitives.

These functions operate on explicit concentration fields and biomass arrays.
They are intentionally phenomenological: they do not infer receptor biology,
cell trajectories, or clinical outcomes from sparse parameters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from marse.additives import hill_response

__all__ = [
    "ImmuneError",
    "ImmuneInteraction",
    "ImmunePressureResult",
    "MolecularNeutralizer",
    "apply_immune_pressure",
]


class ImmuneError(ValueError):
    """An immune or molecular interaction is invalid."""


@dataclass(frozen=True, slots=True)
class ImmuneInteraction:
    """Effector-mediated loss for one named species."""

    species: str
    effector: str
    maximum_kill_per_h: float
    half_effect: float
    hill_coefficient: float = 1.0
    susceptibility: float = 1.0

    def __post_init__(self) -> None:
        if not self.species.strip() or not self.effector.strip():
            raise ImmuneError("species and effector names must not be empty")
        if self.maximum_kill_per_h < 0 or self.half_effect <= 0:
            raise ImmuneError("kill rate must be non-negative and half_effect positive")
        if self.hill_coefficient <= 0 or self.susceptibility < 0:
            raise ImmuneError("hill coefficient must be positive and susceptibility non-negative")


@dataclass(frozen=True, slots=True)
class MolecularNeutralizer:
    """Reduce an effector field by a competing molecular concentration."""

    effector: str
    molecule: str
    neutralization_fraction: float
    half_effect: float
    hill_coefficient: float = 1.0

    def __post_init__(self) -> None:
        if not self.effector.strip() or not self.molecule.strip():
            raise ImmuneError("effector and molecule names must not be empty")
        if not 0 <= self.neutralization_fraction <= 1:
            raise ImmuneError("neutralization_fraction must be in [0, 1]")
        if self.half_effect <= 0 or self.hill_coefficient <= 0:
            raise ImmuneError("neutralizer half_effect and hill_coefficient must be positive")


@dataclass(frozen=True, slots=True)
class ImmunePressureResult:
    """Biomass after pressure and the removed biomass field."""

    biomass: NDArray[np.float64]
    removed: NDArray[np.float64]


def apply_immune_pressure(
    biomass: ArrayLike,
    *,
    effectors: dict[str, ArrayLike],
    molecules: dict[str, ArrayLike] | None = None,
    interaction: ImmuneInteraction,
    neutralizers: tuple[MolecularNeutralizer, ...] = (),
    dt: float,
) -> ImmunePressureResult:
    """Apply deterministic, bounded effector killing to one biomass field."""
    if dt < 0:
        raise ImmuneError("dt must be non-negative")
    target = np.maximum(np.asarray(biomass, dtype=float), 0.0)
    if interaction.effector not in effectors:
        raise ImmuneError(f"missing effector field '{interaction.effector}'")
    effector = np.maximum(np.asarray(effectors[interaction.effector], dtype=float), 0.0)
    if effector.shape != target.shape:
        raise ImmuneError("effector and biomass fields must have matching shapes")
    effective_effector = effector.copy()
    for neutralizer in neutralizers:
        if neutralizer.effector != interaction.effector:
            continue
        if molecules is None or neutralizer.molecule not in molecules:
            raise ImmuneError(f"missing neutralizer molecule '{neutralizer.molecule}'")
        molecule = np.maximum(np.asarray(molecules[neutralizer.molecule], dtype=float), 0.0)
        if molecule.shape != target.shape:
            raise ImmuneError("molecule and biomass fields must have matching shapes")
        occupancy = hill_response(molecule, neutralizer.half_effect, neutralizer.hill_coefficient)
        effective_effector *= 1.0 - neutralizer.neutralization_fraction * occupancy
    kill_rate = (
        interaction.maximum_kill_per_h
        * interaction.susceptibility
        * hill_response(effective_effector, interaction.half_effect, interaction.hill_coefficient)
    )
    updated = target * np.exp(-dt * kill_rate)
    return ImmunePressureResult(updated, target - updated)
