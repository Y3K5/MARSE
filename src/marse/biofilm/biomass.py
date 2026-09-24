"""Biomass distributed through the depth of a biofilm, and where it actually grows.

This is the join between the two halves MARSE has built so far: the spatial
solute field (:mod:`marse.spatial.diffusion`) and the growth kinetics
(:mod:`marse.microbes.growth`). Until now growth responded to a single
well-mixed concentration. Here each depth sees its own concentration, and
therefore grows at its own rate.

The coupling runs one way per solve, and that is deliberate. Uptake capacity
is set by the biomass present, the solute field is solved to steady state
against that capacity, and the resulting concentrations give the local growth
rates. Biomass is *not* advanced here: at fixed biomass this is a
well-defined steady problem, which is exactly the form the published
monospecies benchmark uses, and it avoids the biomass-spreading question that
docs/modeling-landscape.md argues should stay a declared, swappable choice.

The result that falls out is stratification: in an active biofilm thicker than
the penetration depth, only a surface layer grows at all, and the interior is
starved however much biomass it holds (docs/theory.md, section 5.1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse._numeric import require
from marse.microbes.growth import monod
from marse.spatial.diffusion import SoluteProfile, solve_steady_state
from marse.spatial.domain import Grid1D

__all__ = ["BiofilmState", "GrowthProfile", "Population", "solve_growth_profile"]


@dataclass(frozen=True, slots=True)
class Population:
    """One species, with the parameters that govern how it uses the substrate.

    ``yield_per_substrate`` is biomass produced per unit substrate consumed, so
    its reciprocal converts a growth rate into an uptake rate (Pirt, 1965;
    docs/theory.md section 3.2).
    """

    name: str
    mu_max: float
    half_saturation: float
    yield_per_substrate: float

    def __post_init__(self) -> None:
        require(bool(self.name.strip()), "population name must not be empty")
        require(self.mu_max > 0, f"{self.name}: mu_max must be positive")
        require(self.half_saturation > 0, f"{self.name}: half_saturation must be positive")
        require(self.yield_per_substrate > 0, f"{self.name}: yield must be positive")


@dataclass(frozen=True, slots=True)
class BiofilmState:
    """Biomass density of each population at each depth.

    ``density`` has shape ``(populations, nodes)`` and uses the same node
    ordering as the grid, so row ``i`` column ``j`` is the density of
    population ``i`` at depth ``grid.depths[j]``.
    """

    grid: Grid1D
    populations: tuple[Population, ...]
    density: NDArray[np.float64]

    def __post_init__(self) -> None:
        require(len(self.populations) > 0, "at least one population is required")
        names = [p.name for p in self.populations]
        require(len(names) == len(set(names)), f"population names must be unique, got {names}")
        expected = (len(self.populations), self.grid.cells)
        require(
            self.density.shape == expected,
            f"density must have shape {expected}, got {self.density.shape}",
        )
        require(bool(np.all(self.density >= 0)), "density must be non-negative")

    @classmethod
    def uniform(
        cls, grid: Grid1D, populations: tuple[Population, ...], densities: NDArray | list[float]
    ) -> BiofilmState:
        """A biofilm of constant composition through its depth."""
        column = np.asarray(densities, dtype=float).reshape(-1, 1)
        require(
            column.shape[0] == len(populations),
            f"expected one density per population, got {column.shape[0]}",
        )
        return cls(grid, populations, np.repeat(column, grid.cells, axis=1))

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(p.name for p in self.populations)

    @property
    def total_density(self) -> NDArray[np.float64]:
        """Combined biomass density at each depth."""
        return np.sum(self.density, axis=0)

    def areal_density(self) -> float:
        """Biomass per unit surface area, integrated over the depth."""
        total = self.total_density
        return self.grid.integrate(float(total[0]), total)


@dataclass(frozen=True, slots=True)
class GrowthProfile:
    """Where the biofilm grows, once the solute field has been solved."""

    state: BiofilmState
    solute: SoluteProfile
    growth_rate: NDArray[np.float64]  # (populations, nodes), per unit time

    @property
    def depths(self) -> NDArray[np.float64]:
        return self.state.grid.depths

    @property
    def production(self) -> NDArray[np.float64]:
        """Biomass produced per unit volume per unit time, at each depth."""
        return self.growth_rate * self.state.density

    def active_zone(self, fraction_of_max: float = 0.5) -> float:
        """Depth beyond which growth has fallen below ``fraction_of_max`` of its peak.

        The measurable counterpart of the stratification seen in real biofilms,
        where protein synthesis is confined to a band tens of micrometres wide
        next to the oxic surface.
        """
        require(0.0 < fraction_of_max < 1.0, "fraction_of_max must lie in (0, 1)")
        combined = np.sum(self.production, axis=0)
        peak = float(np.max(combined))
        if peak <= 0.0:
            return 0.0
        active = np.flatnonzero(combined >= fraction_of_max * peak)
        return float(self.depths[active[-1]])

    @property
    def active_fraction(self) -> float:
        """Share of total biomass production happening above the active-zone depth."""
        combined = np.sum(self.production, axis=0)
        total = self.state.grid.integrate(float(combined[0]), combined)
        if total <= 0.0:
            return 0.0
        within = self.depths <= self.active_zone()
        return float(np.sum(combined[within]) * self.state.grid.dx / total)

    def mean_growth_rate(self) -> NDArray[np.float64]:
        """Biomass-weighted mean growth rate of each population.

        This is what a measurement of the whole biofilm would report, and it is
        well below the surface rate whenever the film is stratified — the
        difference between the two is the practical consequence of gradients.
        """
        weights = self.state.density
        totals = np.sum(weights, axis=1)
        weighted = np.sum(self.growth_rate * weights, axis=1)
        return np.divide(weighted, totals, out=np.zeros_like(weighted), where=totals > 0)


def solve_growth_profile(
    state: BiofilmState,
    *,
    diffusivity: float,
    surface: float,
    tolerance: float = 1e-12,
    max_iterations: int = 200,
) -> GrowthProfile:
    """Solve the solute field for this biomass, then the growth it supports.

    **The time units must match, and this is the easiest thing to get wrong
    here.** Published diffusivities are quoted in um^2/s while growth rates are
    quoted per hour, and mixing the two silently overstates demand by a factor
    of 3600, which shortens the penetration depth sixtyfold — from something
    like 130 um to 2 um, which still looks like a plausible number. MARSE works
    in hours internally (docs/theory.md section 0.1), so a diffusivity taken
    from a table in um^2/s must be multiplied by 3600 before it is passed here.

    All populations must share a half-saturation constant, because the solver
    represents uptake with a single Monod term; their differing capacities are
    carried by the per-node uptake capacity. Mixed affinities need a reaction
    term that sums several Monod expressions, which is a later extension.
    """
    affinities = {p.half_saturation for p in state.populations}
    require(
        len(affinities) == 1,
        f"all populations must currently share a half_saturation; got {sorted(affinities)}",
    )
    half_saturation = affinities.pop()

    # Uptake capacity at each depth: what the biomass there would consume if
    # the substrate were saturating.
    capacity = np.zeros(state.grid.cells)
    for population, density in zip(state.populations, state.density, strict=True):
        capacity += population.mu_max / population.yield_per_substrate * density

    solute = solve_steady_state(
        state.grid,
        diffusivity=diffusivity,
        surface=surface,
        max_uptake=capacity,
        half_saturation=half_saturation,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )

    growth = np.vstack(
        [
            np.asarray(monod(solute.concentration, p.mu_max, p.half_saturation), dtype=float)
            for p in state.populations
        ]
    )
    return GrowthProfile(state=state, solute=solute, growth_rate=growth)
