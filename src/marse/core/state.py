"""Canonical simulation state.

The state is everything needed to resume a run: the clock, the biomass of each
organism, the substrate concentration, and the cumulative substrate consumed.
That last field is not redundant — it is what makes conservation checkable at
runtime rather than only in tests (docs/theory.md, section 9.5).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

__all__ = ["SimulationState"]


@dataclass(frozen=True, slots=True)
class SimulationState:
    """State at one instant. Arrays are ordered as the configuration lists organisms."""

    time_h: float
    step: int
    biomass_g_per_l: NDArray[np.float64]
    substrate_mm: float
    substrate_consumed_mm: float

    @classmethod
    def initial(cls, biomass: NDArray[np.float64], substrate_mm: float) -> SimulationState:
        return cls(
            time_h=0.0,
            step=0,
            biomass_g_per_l=np.asarray(biomass, dtype=float).copy(),
            substrate_mm=float(substrate_mm),
            substrate_consumed_mm=0.0,
        )

    def advanced(
        self,
        *,
        time_h: float,
        step: int,
        biomass: NDArray[np.float64],
        substrate_mm: float,
        substrate_consumed_mm: float,
    ) -> SimulationState:
        return replace(
            self,
            time_h=time_h,
            step=step,
            biomass_g_per_l=np.asarray(biomass, dtype=float),
            substrate_mm=float(substrate_mm),
            substrate_consumed_mm=float(substrate_consumed_mm),
        )

    @property
    def total_biomass_g_per_l(self) -> float:
        return float(np.sum(self.biomass_g_per_l))

    def to_dict(self, organism_names: tuple[str, ...]) -> dict[str, Any]:
        """Plain data for checkpoints and summaries."""
        return {
            "time_h": self.time_h,
            "step": self.step,
            "substrate_mm": self.substrate_mm,
            "substrate_consumed_mm": self.substrate_consumed_mm,
            "biomass_g_per_l": {
                name: float(value)
                for name, value in zip(organism_names, self.biomass_g_per_l, strict=True)
            },
            "total_biomass_g_per_l": self.total_biomass_g_per_l,
        }
