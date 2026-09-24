"""Steady-state reaction and diffusion of a solute through a slab.

Solves

    D d2C/dz2 = R(C),      C(0) = surface,      dC/dz = 0 at the base,

with Monod uptake ``R(C) = k_max C / (K + C)``. The equation is nonlinear in C,
so it is solved by Newton iteration; each step is a tridiagonal system, solved
directly by the Thomas algorithm rather than a general matrix routine.

**Why steady state rather than time stepping.** Solutes equilibrate two to
three orders of magnitude faster than biomass changes (docs/theory.md §5.3),
and an explicit time step for diffusion on a biofilm-resolution grid is limited
to milliseconds (§9.1). Marching the solute field would therefore spend
millions of steps reproducing an equilibrium that can be solved for directly.
This is the approach taken by the established biofilm models.

The solver is checked against three independent analytical limits — no uptake,
first-order uptake and zero-order uptake — plus a flux balance and a grid
convergence study (validation case V3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from marse._numeric import require
from marse.spatial.domain import Grid1D

__all__ = ["SoluteProfile", "solve_steady_state"]


class ConvergenceError(RuntimeError):
    """Newton iteration did not converge; the result would be meaningless."""


@dataclass(frozen=True, slots=True)
class SoluteProfile:
    """The solved concentration field and the quantities derived from it."""

    grid: Grid1D
    concentration: NDArray[np.float64]
    surface: float
    iterations: int

    @property
    def depths(self) -> NDArray[np.float64]:
        return self.grid.depths

    def uptake(self, max_uptake: ArrayLike, half_saturation: float) -> NDArray[np.float64]:
        """Volumetric uptake rate at each node."""
        c = self.concentration
        return np.asarray(max_uptake, dtype=float) * c / (half_saturation + c)

    def total_uptake(self, max_uptake: ArrayLike, half_saturation: float) -> float:
        """Uptake integrated over the full depth, per unit surface area.

        With per-node capacity the surface rate uses the capacity of the
        topmost node, since the boundary itself holds no biomass of its own.
        """
        capacity = np.asarray(max_uptake, dtype=float)
        topmost = float(capacity if capacity.ndim == 0 else capacity[0])
        surface_rate = topmost * self.surface / (half_saturation + self.surface)
        return self.grid.integrate(surface_rate, self.uptake(capacity, half_saturation))

    def surface_flux(self, diffusivity: float) -> float:
        """Flux into the slab, from a second-order one-sided derivative at the surface.

        At steady state this equals :meth:`total_uptake`; comparing them is the
        strongest available check that the solution is physically consistent,
        because it ties the boundary condition to the interior chemistry.
        """
        dx, c = self.grid.dx, self.concentration
        # Three-point one-sided derivative at z = 0 over the uniformly spaced
        # points surface, c[0], c[1]; the sign makes a downward flux positive.
        return float(diffusivity * (3.0 * self.surface - 4.0 * c[0] + c[1]) / (2.0 * dx))

    def penetration_depth(self, threshold_fraction: float = 0.01) -> float:
        """Depth at which the solute first falls below a fraction of its surface value.

        Returns the slab thickness when it never does, meaning the solute
        reaches the base. This is a measured property of the numerical profile,
        not the analytical estimate in
        :func:`marse.validation.analytical.zero_order_penetration_depth`;
        comparing the two is the point of validation case V3.
        """
        require(0.0 < threshold_fraction < 1.0, "threshold_fraction must lie in (0, 1)")
        below = np.flatnonzero(self.concentration < threshold_fraction * self.surface)
        if below.size == 0:
            return self.grid.thickness
        return float(self.depths[below[0]])

    @property
    def anoxic_fraction(self) -> float:
        """Fraction of the slab holding less than 1% of the surface concentration."""
        return 1.0 - self.penetration_depth() / self.grid.thickness


def _solve_tridiagonal(
    lower: NDArray[np.float64],
    diagonal: NDArray[np.float64],
    upper: NDArray[np.float64],
    rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Thomas algorithm. ``lower[0]`` and ``upper[-1]`` are unused."""
    n = diagonal.size
    c, d = np.zeros(n), np.zeros(n)
    beta = diagonal[0]
    require(beta != 0.0, "singular system")
    c[0], d[0] = upper[0] / beta, rhs[0] / beta
    for i in range(1, n):
        beta = diagonal[i] - lower[i] * c[i - 1]
        if beta == 0.0:
            raise ConvergenceError("tridiagonal system became singular")
        c[i] = upper[i] / beta if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / beta
    x = np.empty(n)
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def solve_steady_state(
    grid: Grid1D,
    *,
    diffusivity: float,
    surface: float,
    max_uptake: ArrayLike,
    half_saturation: float,
    tolerance: float = 1e-12,
    max_iterations: int = 200,
) -> SoluteProfile:
    """Solve for the steady concentration profile through ``grid``.

    Parameters use the caller's units and must be mutually consistent: with
    depth in micrometres, ``diffusivity`` is um^2/s, ``surface`` and
    ``half_saturation`` are mM, and ``max_uptake`` is mM/s.

    ``max_uptake`` may be a single value, for uniform demand, or one value per
    node. The per-node form is what couples this solver to a biofilm whose
    biomass varies with depth: uptake capacity is proportional to the local
    biomass, so a sparsely colonised region consumes less and lets the solute
    travel further. Zero is permitted and gives the pure-diffusion case.
    """
    require(diffusivity > 0, "diffusivity must be positive")
    require(surface >= 0, "surface concentration must be non-negative")
    require(half_saturation > 0, "half_saturation must be positive")
    require(tolerance > 0, "tolerance must be positive")

    capacity = np.asarray(max_uptake, dtype=float)
    require(bool(np.all(capacity >= 0)), "max_uptake must be non-negative")
    require(
        capacity.ndim == 0 or capacity.shape == (grid.cells,),
        f"max_uptake must be a single value or one per node; "
        f"got shape {capacity.shape} for {grid.cells} nodes",
    )
    max_uptake = capacity

    n, dx = grid.cells, grid.dx
    conduction = diffusivity / dx**2

    # The no-flux base is imposed by mirroring the last interior node, which
    # doubles its contribution to the final row of the Jacobian.
    lower = np.full(n, conduction)
    lower[-1] = 2.0 * conduction
    upper = np.full(n, conduction)

    concentration = np.full(n, surface, dtype=float)
    for iteration in range(1, max_iterations + 1):
        left = np.concatenate(([surface], concentration[:-1]))
        right = np.concatenate((concentration[1:], [concentration[-2]]))
        saturation = half_saturation + concentration
        uptake = max_uptake * concentration / saturation
        residual = conduction * (left - 2.0 * concentration + right) - uptake
        slope = max_uptake * half_saturation / saturation**2  # dR/dC

        step = _solve_tridiagonal(lower, -2.0 * conduction - slope, upper, -residual)
        # Concentrations cannot go negative; clamping keeps the uptake term
        # well defined while Newton is still far from the solution.
        concentration = np.maximum(concentration + step, 0.0)
        if not np.all(np.isfinite(concentration)):
            raise ConvergenceError(f"the profile became non-finite at iteration {iteration}")
        if np.max(np.abs(step)) <= tolerance * max(surface, 1.0):
            return SoluteProfile(grid, concentration, float(surface), iteration)

    raise ConvergenceError(
        f"Newton iteration did not converge in {max_iterations} steps; "
        "the last correction was "
        f"{float(np.max(np.abs(step))):.3e}"
    )
