"""Spatial domains.

Phase 2 covers the one-dimensional case: depth through a flat biofilm or
sediment, with a well-mixed bulk liquid above it and an impermeable substratum
below. That geometry is enough for the questions MARSE answers first — how far
a solute penetrates, and what fraction of the population sits beyond it — and
it is the geometry the published benchmark problems use.

Two and three dimensions arrive with biofilm structure in Phase 4; the
interfaces here are deliberately narrow so that adding them does not disturb
what already works.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse._numeric import require

__all__ = ["Grid1D"]


@dataclass(frozen=True, slots=True)
class Grid1D:
    """A uniform one-dimensional grid through the depth of a slab.

    Depth runs from 0 at the surface in contact with the bulk liquid to
    ``thickness`` at the substratum. Solved nodes sit at
    ``dx, 2 dx, ..., thickness`` with uniform spacing ``dx``. The surface at
    ``z = 0`` is a boundary value rather than an unknown, so it is not among
    them; the last node lies on the impermeable base, where the no-flux
    condition applies.

    Units are the caller's, but must be consistent with the diffusivity used
    with them. MARSE works in micrometres and seconds internally.
    """

    thickness: float
    cells: int

    def __post_init__(self) -> None:
        require(self.thickness > 0, "thickness must be positive")
        require(isinstance(self.cells, int) and self.cells >= 2, "cells must be an integer >= 2")

    @property
    def dx(self) -> float:
        """Spacing between nodes."""
        return self.thickness / self.cells

    @property
    def depths(self) -> NDArray[np.float64]:
        """Depth of each node, from ``dx`` to ``thickness`` inclusive."""
        return np.linspace(self.dx, self.thickness, self.cells)

    def integrate(self, surface_value: float, node_values: NDArray[np.float64]) -> float:
        """Integrate a quantity over the full depth, from the surface to the base.

        The trapezoidal rule over ``z = 0, dx, ..., thickness``. The surface
        value has to be supplied because it is a boundary condition rather than
        a solved node, and omitting it would bias the integral by half a
        spacing — enough to break a flux balance that should hold exactly.
        """
        values = np.asarray(node_values, dtype=float)
        require(values.shape == (self.cells,), f"expected {self.cells} node values")
        interior = float(np.sum(values[:-1]))
        return float((0.5 * surface_value + interior + 0.5 * values[-1]) * self.dx)
