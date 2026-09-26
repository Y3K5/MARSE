"""Space: the box of voxels a network runs in, and what starts where.

A version 2 document with a ``domain`` runs in space (docs/networks.md,
"Running a network in space"). The domain is a box of cubic voxels over a flat
substratum, height its last axis, with the bulk liquid held above it
(:mod:`marse.spatial.grid`). It states:

- ``voxels`` and ``voxel_um``: the grid, in one, two or three dimensions;
- ``bulk_mol_per_m3``: what the liquid above holds, dissolved components only;
- ``diffusivity_m2_per_s``: one value for every dissolved component;
  particulate components (biomass) do not diffuse, and grow in place until
  Stage 2d adds spreading;
- ``colonies``: hemispheres of a particulate component on the substratum, each
  placed at stated coordinates;
- ``random_colonies``: a number of such hemispheres placed at random, from the
  run's seed. The placement is recorded as a random stream, so a replay puts
  them back in the same voxels.

The run fields' ``initial_mol_per_m3`` fill every voxel, and the colonies then
set their component inside their hemispheres. A colony that covers no voxel
is refused, because it would silently vanish.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from marse.core.config import ConfigError
from marse.core.seeds import SeedRegistry
from marse.schemas._reading import plain, read_object
from marse.schemas.network import (
    COLONY_FIELDS,
    DOMAIN_FIELDS,
    RANDOM_COLONY_FIELDS,
    Network,
    _known,
)
from marse.spatial.grid import Grid

__all__ = ["M2_PER_S_TO_UM2_PER_H", "Colony", "Domain", "RandomColonies", "read_domain"]

M2_PER_S_TO_UM2_PER_H = 1e12 * 3600.0
MAX_VOXELS = 2**24
COLONY_STREAM = "colonies"


@dataclass(frozen=True, slots=True)
class Colony:
    """A hemisphere of one particulate component, standing on the substratum."""

    component: str
    center_um: tuple[float, ...]
    radius_um: float
    concentration_mol_per_m3: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "component": self.component,
            "center_um": list(self.center_um),
            "radius_um": self.radius_um,
            "concentration_mol_per_m3": self.concentration_mol_per_m3,
        }


@dataclass(frozen=True, slots=True)
class RandomColonies:
    """``count`` hemispheres of one particulate component, centred at random on the surface."""

    component: str
    count: int
    radius_um: float
    concentration_mol_per_m3: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "component": self.component,
            "count": self.count,
            "radius_um": self.radius_um,
            "concentration_mol_per_m3": self.concentration_mol_per_m3,
        }


@dataclass(frozen=True, slots=True)
class Domain:
    """The grid, the liquid above it, how each component moves, and where the colonies are."""

    voxels: tuple[int, ...]
    voxel_um: float
    bulk_mol_per_m3: dict[str, float]
    diffusivity_m2_per_s: dict[str, float]
    colonies: tuple[Colony, ...] = ()
    random_colonies: tuple[RandomColonies, ...] = ()

    @property
    def grid(self) -> Grid:
        return Grid(self.voxels, self.voxel_um)

    def diffusivities_um2_per_h(self, names: tuple[str, ...]) -> NDArray[np.float64]:
        return np.array(
            [self.diffusivity_m2_per_s.get(n, 0.0) * M2_PER_S_TO_UM2_PER_H for n in names]
        )

    def bulk(self, names: tuple[str, ...]) -> NDArray[np.float64]:
        return np.array([self.bulk_mol_per_m3.get(n, 0.0) for n in names])

    @property
    def random_streams(self) -> tuple[str, ...]:
        return (COLONY_STREAM,) if self.random_colonies else ()

    def placed_colonies(self, seed: int) -> tuple[Colony, ...]:
        """Every colony with its centre: the stated ones, then the random ones from the seed."""
        placed = list(self.colonies)
        if self.random_colonies:
            rng = SeedRegistry(seed).stream(COLONY_STREAM)
            grid = self.grid
            sizes = [grid.size_um[a] for a in grid.lateral_axes]
            for group in self.random_colonies:
                centres = rng.uniform(0.0, 1.0, size=(group.count, len(sizes))) * np.array(sizes)
                placed.extend(
                    Colony(
                        group.component,
                        tuple(float(x) for x in centre),
                        group.radius_um,
                        group.concentration_mol_per_m3,
                    )
                    for centre in centres
                )
        return tuple(placed)

    def initial_state(
        self, names: tuple[str, ...], uniform: dict[str, float], seed: int
    ) -> NDArray[np.float64]:
        """The starting concentrations, shape (components, *voxels)."""
        grid = self.grid
        state = np.empty((len(names), *grid.shape))
        for j, name in enumerate(names):
            state[j] = uniform.get(name, 0.0)
        index = {name: j for j, name in enumerate(names)}
        for colony in self.placed_colonies(seed):
            state[index[colony.component]][_hemisphere(grid, colony)] = (
                colony.concentration_mol_per_m3
            )
        return state

    def to_dict(self) -> dict[str, Any]:
        return {
            "voxels": list(self.voxels),
            "voxel_um": self.voxel_um,
            "bulk_mol_per_m3": dict(self.bulk_mol_per_m3),
            "diffusivity_m2_per_s": dict(self.diffusivity_m2_per_s),
            "colonies": [c.to_dict() for c in self.colonies],
            "random_colonies": [c.to_dict() for c in self.random_colonies],
        }


def _hemisphere(grid: Grid, colony: Colony) -> NDArray[np.bool_]:
    """Voxels whose centre lies inside the colony; lateral distances wrap around the box."""
    squared = np.zeros(grid.shape)
    dims = grid.dimensions
    for axis, centre in zip(grid.lateral_axes, colony.center_um, strict=True):
        length = grid.size_um[axis]
        offset = np.abs(grid.centers_um(axis) - centre) % length
        offset = np.minimum(offset, length - offset)
        shape = [1] * dims
        shape[axis] = -1
        squared = squared + (offset**2).reshape(shape)
    shape = [1] * dims
    shape[-1] = -1
    squared = squared + (grid.heights_um() ** 2).reshape(shape)
    return squared <= colony.radius_um**2


def _positive(value: Any, where: str) -> float:
    if not value > 0:
        raise ConfigError(f"{where}: must be positive")
    return float(plain(value))


def _colony_common(
    values: dict[str, Any], network: Network, where: str
) -> tuple[str, float, float]:
    components = {c.name: c for c in network.components}
    _known([values["component"]], components, f"{where}.component")
    component = values["component"]
    if components[component].phase != "particulate":
        raise ConfigError(
            f"{where}.component: '{component}' is dissolved; a colony is made of a particulate "
            "component, such as biomass"
        )
    radius = _positive(values["radius_um"], f"{where}.radius_um")
    concentration = _positive(
        values["concentration_mol_per_m3"], f"{where}.concentration_mol_per_m3"
    )
    return component, radius, concentration


def read_domain(raw: Any, network: Network, where: str = "experiment.domain") -> Domain:
    """Read and check a domain against the network it will run."""
    values = read_object(raw, where, DOMAIN_FIELDS)
    try:
        grid = Grid(values["voxels"], 1.0)
    except ValueError as error:
        raise ConfigError(f"{where}.voxels: {error}") from None
    if grid.voxels > MAX_VOXELS:
        raise ConfigError(
            f"{where}.voxels: {grid.voxels:,} voxels is more than the {MAX_VOXELS:,} one run "
            "can hold; use larger voxels or a smaller box"
        )
    voxel_um = _positive(values["voxel_um"], f"{where}.voxel_um")
    grid = Grid(values["voxels"], voxel_um)
    components = {c.name: c for c in network.components}

    bulk = values.get("bulk_mol_per_m3", {})
    _known(list(bulk), components, f"{where}.bulk_mol_per_m3")
    for name, amount in bulk.items():
        if components[name].phase != "dissolved":
            raise ConfigError(
                f"{where}.bulk_mol_per_m3: '{name}' is particulate; the bulk liquid holds "
                "dissolved components"
            )
        if amount < 0:
            raise ConfigError(f"{where}.bulk_mol_per_m3: '{name}' must not be negative")

    diffusivity = values["diffusivity_m2_per_s"]
    _known(list(diffusivity), components, f"{where}.diffusivity_m2_per_s")
    for name, value in diffusivity.items():
        if components[name].phase != "dissolved":
            raise ConfigError(
                f"{where}.diffusivity_m2_per_s: '{name}' is particulate and does not diffuse; "
                "biomass grows in place until colonies can spread (Stage 2d)"
            )
        if value < 0:
            raise ConfigError(f"{where}.diffusivity_m2_per_s: '{name}' must not be negative")
    missing = [
        c.name for c in network.components if c.phase == "dissolved" and c.name not in diffusivity
    ]
    if missing:
        listed = ", ".join(f"'{m}'" for m in missing)
        raise ConfigError(
            f"{where}.diffusivity_m2_per_s: every dissolved component needs a diffusivity to "
            f"run in space; {listed} has none"
        )

    lateral = [grid.size_um[a] for a in grid.lateral_axes]
    colonies = []
    for i, item in enumerate(values.get("colonies", [])):
        here = f"{where}.colonies[{i}]"
        colony_values = read_object(item, here, COLONY_FIELDS)
        component, radius, concentration = _colony_common(colony_values, network, here)
        centre = tuple(float(plain(x)) for x in colony_values["center_um"])
        if len(centre) != len(lateral):
            raise ConfigError(
                f"{here}.center_um: a {grid.dimensions}-D box needs {len(lateral)} lateral "
                f"coordinate(s), got {len(centre)}"
            )
        for x, length in zip(centre, lateral, strict=True):
            if not 0 <= x < length:
                raise ConfigError(f"{here}.center_um: {x:g} lies outside the box, 0 to {length:g}")
        colony = Colony(component, centre, radius, concentration)
        if not _hemisphere(grid, colony).any():
            raise ConfigError(
                f"{here}: covers no voxel; give it a radius of at least one voxel "
                f"({voxel_um:g} um) or centre it on a voxel"
            )
        colonies.append(colony)

    random_colonies = []
    for i, item in enumerate(values.get("random_colonies", [])):
        here = f"{where}.random_colonies[{i}]"
        group_values = read_object(item, here, RANDOM_COLONY_FIELDS)
        component, radius, concentration = _colony_common(group_values, network, here)
        if not lateral:
            raise ConfigError(
                f"{here}: a one-dimensional column has no surface to scatter colonies on; "
                "use colonies instead"
            )
        count = group_values["count"]
        if count < 1:
            raise ConfigError(f"{here}.count: must be at least 1")
        if radius < 0.5 * voxel_um * math.sqrt(grid.dimensions):
            raise ConfigError(
                f"{here}.radius_um: {radius:g} um can miss every voxel centre; use at least "
                f"{0.5 * voxel_um * math.sqrt(grid.dimensions):g} um"
            )
        random_colonies.append(RandomColonies(component, count, radius, concentration))

    return Domain(
        voxels=tuple(values["voxels"]),
        voxel_um=voxel_um,
        bulk_mol_per_m3={
            c.name: float(plain(bulk[c.name])) if c.name in bulk else 0.0
            for c in network.components
            if c.phase == "dissolved"
        },
        diffusivity_m2_per_s={
            c.name: float(plain(diffusivity[c.name]))
            for c in network.components
            if c.phase == "dissolved"
        },
        colonies=tuple(colonies),
        random_colonies=tuple(random_colonies),
    )
