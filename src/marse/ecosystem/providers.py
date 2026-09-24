"""Versioned provider contracts for the ecosystem update pipeline.

Providers are deliberately small and numerical. They do not own configuration
or orchestration; the ecosystem model supplies fields, units, and boundaries.
This makes alternative transport or mechanics implementations replaceable
without changing the canonical state representation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "EcosystemProviders",
    "ExplicitTransportProvider",
    "FieldBoundary",
    "NoBoundary",
    "TransportProvider",
]


class FieldBoundary(Protocol):
    boundary_value: float | None
    boundary_edges: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class NoBoundary:
    """Boundary contract for biomass spreading, which has no fixed edges."""

    boundary_value: float | None = None
    boundary_edges: tuple[str, ...] = ()


class TransportProvider(Protocol):
    """Advance one non-negative scalar field by explicit transport."""

    version: str

    def advance(
        self,
        field: NDArray[np.float64],
        *,
        diffusivity: float,
        dt: float,
        cell_size_um: float,
        decay_per_h: float = 0.0,
        boundary: FieldBoundary,
    ) -> NDArray[np.float64]:
        """Return the transported field after one timestep."""


def _laplacian(field: NDArray[np.float64]) -> NDArray[np.float64]:
    padded = np.pad(field, 1, mode="edge")
    return padded[1:-1, :-2] + padded[1:-1, 2:] + padded[:-2, 1:-1] + padded[2:, 1:-1] - 4.0 * field


def _apply_boundary(field: NDArray[np.float64], boundary: FieldBoundary) -> None:
    if boundary.boundary_value is None:
        return
    value = boundary.boundary_value
    if "top" in boundary.boundary_edges:
        field[0, :] = value
    if "bottom" in boundary.boundary_edges:
        field[-1, :] = value
    if "left" in boundary.boundary_edges:
        field[:, 0] = value
    if "right" in boundary.boundary_edges:
        field[:, -1] = value


@dataclass(frozen=True, slots=True)
class ExplicitTransportProvider:
    """Reference finite-difference transport provider used by MARSE v1."""

    version: str = "explicit_transport_2d_v1"

    def advance(
        self,
        field: NDArray[np.float64],
        *,
        diffusivity: float,
        dt: float,
        cell_size_um: float,
        decay_per_h: float = 0.0,
        boundary: FieldBoundary,
    ) -> NDArray[np.float64]:
        updated = field + dt * diffusivity / cell_size_um**2 * _laplacian(field)
        if decay_per_h:
            updated *= np.exp(-dt * decay_per_h)
        _apply_boundary(updated, boundary)
        return updated


@dataclass(frozen=True, slots=True)
class EcosystemProviders:
    """The declared numerical providers for one ecosystem run."""

    transport: TransportProvider = ExplicitTransportProvider()
    biomass_transport: TransportProvider = ExplicitTransportProvider(
        version="explicit_biomass_spreading_2d_v1"
    )
    reaction_version: str = "local_reactions_2d_v1"
    biomass_version: str = "logistic_biomass_2d_v1"
    diagnostics_version: str = "frame_diagnostics_v1"

    @property
    def versions(self) -> dict[str, str]:
        return {
            "transport": self.transport.version,
            "biomass_transport": self.biomass_transport.version,
            "reactions": self.reaction_version,
            "biomass": self.biomass_version,
            "diagnostics": self.diagnostics_version,
        }
