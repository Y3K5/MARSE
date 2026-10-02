"""Space: the box of voxels a network runs in, and what starts where.

A version 2 document with a ``domain`` runs in space (docs/networks.md,
"Running a network in space"). The domain is a box of cubic voxels over a flat
substratum, height its last axis, with the bulk liquid held above it
(:mod:`marse.spatial.grid`). It states:

- ``voxels`` and ``voxel_um``: the grid, in one, two or three dimensions;
- ``bulk_mol_per_m3``: what the liquid above holds, dissolved components only;
- ``diffusivity_m2_per_s``: one value for every dissolved component;
  particulate components (biomass) do not diffuse. They grow in place, unless
  the domain spreads them;
- ``colonies``: hemispheres of a particulate component on the substratum, each
  placed at stated coordinates;
- ``random_colonies``: a number of such hemispheres placed at random, from the
  run's seed. The placement is recorded as a random stream, so a replay puts
  them back in the same voxels;
- ``substratum``, ``liquid``, ``flow``, ``suspension`` and ``adhesion``,
  stated together: cells in the liquid binding to the materials of the
  substratum (docs/environments.md, :mod:`marse.microbes.adhesion`);
- ``film`` and ``mouth``, stated together: the top voxels are a salivary film,
  closed to the air and renewed from the mouth's saliva, which is secreted and
  swallowed as the run goes (docs/environments.md, :mod:`marse.oral`);
- ``diet``, with the mouth: rinses held, drinks sipped and foods that dissolve,
  each for a stated time, and food some of them leave on the teeth;
- ``spreading``: biomass that outgrows its voxel pushes the excess into its
  neighbours, every species sharing the room (:mod:`marse.spatial.spreading`).
  Each species then needs a packing density, and the start must fit. A column
  spreads; boxes in two and three dimensions arrive with Stage 2d, 2d.3.

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
    ADHESION_FIELDS,
    COLONY_FIELDS,
    DOMAIN_FIELDS,
    FILM_FIELDS,
    FLOW_FIELDS,
    INTAKE_FIELDS,
    LIQUID_FIELDS,
    MOUTH_FIELDS,
    PATCH_FIELDS,
    RANDOM_COLONY_FIELDS,
    RETAINED_FIELDS,
    SPREADING_FIELDS,
    SUBSTRATUM_FIELDS,
    SUSPENSION_FIELDS,
    Network,
    _known,
)
from marse.spatial.grid import Grid
from marse.spatial.surface import Patch, Substratum, faces_within

__all__ = [
    "DEFAULT_MIXING_PER_S",
    "M2_PER_S_TO_UM2_PER_H",
    "Adhesion",
    "Colony",
    "Domain",
    "Film",
    "Flow",
    "Intake",
    "Liquid",
    "Mouth",
    "RandomColonies",
    "Retained",
    "Surface",
    "Suspension",
    "read_domain",
]

M2_PER_S_TO_UM2_PER_H = 1e12 * 3600.0
MAX_VOXELS = 2**24
COLONY_STREAM = "colonies"
DEFAULT_MIXING_PER_S = 1.0
"""How fast an intake mixes the film with the mouth's liquid, per second, unless it says."""


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
class Surface:
    """What the substratum is made of: patches of materials under one conditioning film."""

    conditioning_film: str
    patches: tuple[Patch, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "conditioning_film": self.conditioning_film,
            "patches": [
                {"material": p.material, "region_um": list(p.region_um)} for p in self.patches
            ],
        }


@dataclass(frozen=True, slots=True)
class Liquid:
    """The liquid's temperature and viscosity, which set how fast cells diffuse."""

    temperature_c: float
    viscosity_mpa_s: float

    def to_dict(self) -> dict[str, Any]:
        return {"temperature_c": self.temperature_c, "viscosity_mpa_s": self.viscosity_mpa_s}


@dataclass(frozen=True, slots=True)
class Flow:
    """The shear the liquid's flow applies at the substratum, and how far downstream it is."""

    wall_shear_rate_per_s: float
    distance_from_inlet_mm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "wall_shear_rate_per_s": self.wall_shear_rate_per_s,
            "distance_from_inlet_mm": self.distance_from_inlet_mm,
        }


@dataclass(frozen=True, slots=True)
class Suspension:
    """A species suspended in the liquid, and the two components it binds into."""

    reversible: str
    attached: str
    cells_per_ml: float
    cell_diameter_um: float
    carbon_fmol_per_cell: float
    blocked_area_um2: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "reversible": self.reversible,
            "attached": self.attached,
            "cells_per_ml": self.cells_per_ml,
            "cell_diameter_um": self.cell_diameter_um,
            "carbon_fmol_per_cell": self.carbon_fmol_per_cell,
            "blocked_area_um2": self.blocked_area_um2,
        }


@dataclass(frozen=True, slots=True)
class Adhesion:
    """How one species binds to one material: efficiency, detachment and locking."""

    attached: str
    material: str
    efficiency: float
    detachment_per_h: float
    locking_per_h: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "attached": self.attached,
            "material": self.material,
            "efficiency": self.efficiency,
            "detachment_per_h": self.detachment_per_h,
            "locking_per_h": self.locking_per_h,
        }


SCENE = ("substratum", "liquid", "flow", "suspension", "adhesion")
"""The fields that describe cells binding to a surface, stated all together or not at all."""


@dataclass(frozen=True, slots=True)
class Film:
    """The salivary film over the box: its thickness, its mean velocity, and the plaque it crosses.

    The film is the top ``thickness_um`` of the box. It moves over the tooth
    at ``velocity_mm_per_min``, and reaches the site after crossing
    ``plaque_length_mm`` of plaque, so fresh saliva replaces it at the rate
    u(z) / l in each of its voxels (docs/theory.md, section 4.8).
    """

    thickness_um: float
    velocity_mm_per_min: float
    plaque_length_mm: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "thickness_um": self.thickness_um,
            "velocity_mm_per_min": self.velocity_mm_per_min,
            "plaque_length_mm": self.plaque_length_mm,
        }


@dataclass(frozen=True, slots=True)
class Mouth:
    """The mouth's saliva: what the glands secrete, how much the mouth holds, and when it swallows.

    Dawes's (1983) model: the volume grows from ``resting_volume_ml`` at the
    salivary flow until it reaches ``swallow_volume_ml``, when a swallow takes
    it back to the resting volume. The flow is ``unstimulated_flow_ml_per_min``
    plus up to ``stimulated_flow_ml_per_min`` more, half of it at
    ``stimulus_half_mol_per_m3`` of the ``stimulus``. Secreted saliva moves
    from ``saliva_mol_per_m3`` towards ``stimulated_saliva_mol_per_m3`` as the
    flow rises. ``plaque_area_cm2`` is the plaque the box stands for, which
    exchanges with the mouth through its film. Every composition names every
    dissolved component.
    """

    saliva_mol_per_m3: dict[str, float]
    stimulated_saliva_mol_per_m3: dict[str, float] | None
    resting_volume_ml: float
    swallow_volume_ml: float
    unstimulated_flow_ml_per_min: float
    stimulated_flow_ml_per_min: float
    stimulus: str | None
    stimulus_half_mol_per_m3: float | None
    plaque_area_cm2: float
    initial_mol_per_m3: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        written: dict[str, Any] = {"saliva_mol_per_m3": dict(self.saliva_mol_per_m3)}
        if self.stimulated_saliva_mol_per_m3 is not None:
            written["stimulated_saliva_mol_per_m3"] = dict(self.stimulated_saliva_mol_per_m3)
        written |= {
            "resting_volume_ml": self.resting_volume_ml,
            "swallow_volume_ml": self.swallow_volume_ml,
            "unstimulated_flow_ml_per_min": self.unstimulated_flow_ml_per_min,
            "stimulated_flow_ml_per_min": self.stimulated_flow_ml_per_min,
        }
        if self.stimulus is not None:
            written["stimulus"] = self.stimulus
            written["stimulus_half_mol_per_m3"] = self.stimulus_half_mol_per_m3
        written["plaque_area_cm2"] = self.plaque_area_cm2
        written["initial_mol_per_m3"] = dict(self.initial_mol_per_m3)
        return written


@dataclass(frozen=True, slots=True)
class Retained:
    """Food an intake leaves on the teeth: an amount per m² of a particulate component.

    It is placed in the film over ``region_um`` of the substratum, bounds as a
    patch has them, when the intake ends; a process of the network releases
    what dissolves from it.
    """

    component: str
    amount_mol_per_m2: float
    region_um: tuple[float, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "component": self.component,
            "amount_mol_per_m2": self.amount_mol_per_m2,
            "region_um": list(self.region_um),
        }


@dataclass(frozen=True, slots=True)
class Intake:
    """Something taken into the mouth from ``start_h`` for ``duration_min``.

    - a **rinse** adds ``volume_ml`` of liquid of ``composition_mol_per_m3``
      at once; the mouth holds it without swallowing, and expels everything
      above its resting volume at the end;
    - a **drink** of ``volume_ml`` is sipped steadily over the duration, and
      swallowed as the mouth fills;
    - a **food** releases ``released_mmol`` into the saliva steadily over the
      duration, as a sweet sucked slowly does, and adds no liquid.

    While it lasts, the film mixes with the mouth's liquid at ``mixing_per_s``
    in each of its voxels. ``retained`` is food it leaves on the teeth.
    """

    kind: str
    start_h: float
    duration_min: float
    volume_ml: float | None = None
    composition_mol_per_m3: dict[str, float] | None = None
    released_mmol: dict[str, float] | None = None
    mixing_per_s: float = DEFAULT_MIXING_PER_S
    retained: Retained | None = None

    @property
    def end_h(self) -> float:
        return self.start_h + self.duration_min / 60.0

    def to_dict(self) -> dict[str, Any]:
        written: dict[str, Any] = {
            "kind": self.kind,
            "start_h": self.start_h,
            "duration_min": self.duration_min,
        }
        if self.volume_ml is not None:
            written["volume_ml"] = self.volume_ml
        if self.composition_mol_per_m3 is not None:
            written["composition_mol_per_m3"] = dict(self.composition_mol_per_m3)
        if self.released_mmol is not None:
            written["released_mmol"] = dict(self.released_mmol)
        written["mixing_per_s"] = self.mixing_per_s
        if self.retained is not None:
            written["retained"] = self.retained.to_dict()
        return written


DEFAULT_SPREADING_INTERVAL_H = 0.25
"""The default spreading interval, from the measurement of criterion D7 (docs/validation.md)."""


@dataclass(frozen=True, slots=True)
class Spreading:
    """How biomass spreads, and how often.

    Each recording step is cut into equal spans of at most ``interval_h``, and
    the biomass is spread at the end of every span.
    """

    mechanism: str
    interval_h: float

    def to_dict(self) -> dict[str, Any]:
        return {"mechanism": self.mechanism, "interval_h": self.interval_h}


@dataclass(frozen=True, slots=True)
class Domain:
    """The grid, the liquid above it, how each component moves, and where the colonies are."""

    voxels: tuple[int, ...]
    voxel_um: float
    bulk_mol_per_m3: dict[str, float]
    diffusivity_m2_per_s: dict[str, float]
    colonies: tuple[Colony, ...] = ()
    random_colonies: tuple[RandomColonies, ...] = ()
    surface: Surface | None = None
    liquid: Liquid | None = None
    flow: Flow | None = None
    suspension: tuple[Suspension, ...] = ()
    adhesion: tuple[Adhesion, ...] = ()
    film: Film | None = None
    mouth: Mouth | None = None
    diet: tuple[Intake, ...] = ()
    spreading: Spreading | None = None

    @property
    def grid(self) -> Grid:
        return Grid(self.voxels, self.voxel_um)

    @property
    def substratum(self) -> Substratum | None:
        """The material of every face of the substratum, if cells bind to it."""
        return None if self.surface is None else Substratum(self.grid, self.surface.patches)

    def adhesion_of(self, attached: str, material: str) -> Adhesion:
        return next(a for a in self.adhesion if a.attached == attached and a.material == material)

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
        written: dict[str, Any] = {
            "voxels": list(self.voxels),
            "voxel_um": self.voxel_um,
            "bulk_mol_per_m3": dict(self.bulk_mol_per_m3),
            "diffusivity_m2_per_s": dict(self.diffusivity_m2_per_s),
            "colonies": [c.to_dict() for c in self.colonies],
            "random_colonies": [c.to_dict() for c in self.random_colonies],
        }
        # Only a scene with binding writes these, so a domain without one reads
        # back, and hashes, exactly as it did before surfaces existed.
        if self.surface is not None and self.liquid is not None and self.flow is not None:
            written["substratum"] = self.surface.to_dict()
            written["liquid"] = self.liquid.to_dict()
            written["flow"] = self.flow.to_dict()
            written["suspension"] = [s.to_dict() for s in self.suspension]
            written["adhesion"] = [a.to_dict() for a in self.adhesion]
        # Likewise, only a scene in the mouth writes these.
        if self.film is not None and self.mouth is not None:
            written["film"] = self.film.to_dict()
            written["mouth"] = self.mouth.to_dict()
            if self.diet:
                written["diet"] = [intake.to_dict() for intake in self.diet]
        # And only a domain that spreads writes this.
        if self.spreading is not None:
            written["spreading"] = self.spreading.to_dict()
        return written


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
                "biomass grows in place, or moves when the domain spreads it"
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

    stated = [key for key in SCENE if key in values]
    scene: dict[str, Any] = {}
    if stated:
        missing = [key for key in SCENE if key not in values]
        if missing:
            raise ConfigError(
                f"{where}: cells binding to a surface need {', '.join(SCENE)} together; "
                f"{', '.join(missing)} {'is' if len(missing) == 1 else 'are'} missing"
            )
        scene = _read_scene(values, grid, network, where)
    oral = (
        _read_oral(values, grid, network, where)
        if any(key in values for key in ("film", "mouth", "diet"))
        else {}
    )
    if oral and any(amount for amount in bulk.values()):
        raise ConfigError(
            f"{where}.bulk_mol_per_m3: under a film, the liquid comes from the mouth; give its "
            "composition as mouth.saliva_mol_per_m3 and leave the bulk out"
        )
    spreading = (
        _read_spreading(values["spreading"], grid, network, bool(oral), f"{where}.spreading")
        if "spreading" in values
        else None
    )

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
        **scene,
        **oral,
        spreading=spreading,
    )


def _read_spreading(raw: Any, grid: Grid, network: Network, oral: bool, where: str) -> Spreading:
    """The spreading block, and what a network must state to spread."""
    values = read_object(raw, where, SPREADING_FIELDS)
    if grid.dimensions != 1:
        raise ConfigError(
            f"{where}: a column spreads; spreading in two and three dimensions arrives with "
            "Stage 2d, increment 2d.3 (docs/stage-2d-plan.md)"
        )
    if grid.shape[-1] < 2:
        raise ConfigError(
            f"{where}: a column that spreads needs at least two voxels, since biomass must stay "
            "out of the top one"
        )
    if oral:
        raise ConfigError(
            f"{where}: spreading under a salivary film is not supported yet; it arrives when the "
            "oral scenes need it"
        )
    interval = (
        float(plain(values["interval_h"]))
        if "interval_h" in values
        else DEFAULT_SPREADING_INTERVAL_H
    )
    if not interval > 0:
        raise ConfigError(f"{where}.interval_h: must be positive")
    species = {p.growth.biomass: p.name for p in network.processes if p.growth is not None}
    for name, process in species.items():
        if network.component(name).density_mol_per_m3 is None:
            raise ConfigError(
                f"{where}: '{name}' grows in process '{process}', so it takes up room as it "
                "spreads; give it a density_mol_per_m3, its concentration when it alone fills "
                "a voxel"
            )
    return Spreading(values["mechanism"], interval)


def _composition(
    given: dict[str, Any], network: Network, where: str, holder: str = "saliva holds"
) -> dict[str, float]:
    """A liquid's composition, every dissolved component named, the unnamed ones at zero."""
    components = {c.name: c for c in network.components}
    _known(list(given), components, where)
    for name, amount in given.items():
        if components[name].phase != "dissolved":
            raise ConfigError(f"{where}: '{name}' is particulate; {holder} dissolved components")
        if amount < 0:
            raise ConfigError(f"{where}: '{name}' must not be negative")
    return {
        c.name: float(plain(given[c.name])) if c.name in given else 0.0
        for c in network.components
        if c.phase == "dissolved"
    }


def _read_oral(values: dict[str, Any], grid: Grid, network: Network, where: str) -> dict[str, Any]:
    """The film and the mouth, which come together, checked against the box and each other."""
    if "diet" in values and "mouth" not in values:
        raise ConfigError(
            f"{where}.diet: a diet is taken into the mouth; give the domain a film and a mouth"
        )
    missing = [key for key in ("film", "mouth") if key not in values]
    if missing:
        raise ConfigError(
            f"{where}: a salivary film and the mouth that renews it come together; "
            f"{missing[0]} is missing"
        )
    here = f"{where}.film"
    film_values = read_object(values["film"], here, FILM_FIELDS)
    thickness = _positive(film_values["thickness_um"], f"{here}.thickness_um")
    layers = thickness / grid.voxel_um
    if abs(layers - round(layers)) > 1e-9 * layers or round(layers) < 1:
        raise ConfigError(
            f"{here}.thickness_um: {thickness:g} um is not a whole number of voxels of "
            f"{grid.voxel_um:g} um"
        )
    if round(layers) >= grid.shape[-1]:
        raise ConfigError(
            f"{here}.thickness_um: the film fills the box; make the box taller than "
            f"{thickness:g} um, so there is room below it"
        )
    film = Film(
        thickness,
        _positive(film_values["velocity_mm_per_min"], f"{here}.velocity_mm_per_min"),
        _positive(film_values["plaque_length_mm"], f"{here}.plaque_length_mm"),
    )

    here = f"{where}.mouth"
    m = read_object(values["mouth"], here, MOUTH_FIELDS)
    components = {c.name: c for c in network.components}
    saliva = _composition(m["saliva_mol_per_m3"], network, f"{here}.saliva_mol_per_m3")
    stimulated = None
    if "stimulated_saliva_mol_per_m3" in m:
        stimulated = _composition(
            m["stimulated_saliva_mol_per_m3"], network, f"{here}.stimulated_saliva_mol_per_m3"
        )
    initial = saliva
    if "initial_mol_per_m3" in m:
        initial = _composition(m["initial_mol_per_m3"], network, f"{here}.initial_mol_per_m3")
    resting = _positive(m["resting_volume_ml"], f"{here}.resting_volume_ml")
    swallow = _positive(m["swallow_volume_ml"], f"{here}.swallow_volume_ml")
    if not swallow > resting:
        raise ConfigError(f"{here}.swallow_volume_ml: must be more than resting_volume_ml")
    unstimulated = _positive(
        m["unstimulated_flow_ml_per_min"], f"{here}.unstimulated_flow_ml_per_min"
    )
    stimulated_flow = float(plain(m.get("stimulated_flow_ml_per_min", 0)))
    if stimulated_flow < 0:
        raise ConfigError(f"{here}.stimulated_flow_ml_per_min: must not be negative")
    stimulus, half = m.get("stimulus"), None
    if stimulated_flow > 0 and stimulus is None:
        raise ConfigError(
            f"{here}: stimulated_flow_ml_per_min needs the stimulus that raises the flow, "
            "such as sugar"
        )
    if stimulus is not None:
        _known([stimulus], components, f"{here}.stimulus")
        if components[stimulus].phase != "dissolved":
            raise ConfigError(
                f"{here}.stimulus: '{stimulus}' is particulate; the mouth tastes what is dissolved"
            )
        if "stimulus_half_mol_per_m3" not in m:
            raise ConfigError(f"{here}: a stimulus needs stimulus_half_mol_per_m3")
        half = _positive(m["stimulus_half_mol_per_m3"], f"{here}.stimulus_half_mol_per_m3")
    elif "stimulus_half_mol_per_m3" in m:
        raise ConfigError(f"{here}.stimulus_half_mol_per_m3: applies only with a stimulus")
    area = _positive(m["plaque_area_cm2"], f"{here}.plaque_area_cm2")
    film_ml = area * thickness * 1e-4  # cm2 x um, in mL
    if not film_ml < resting:
        raise ConfigError(
            f"{here}.plaque_area_cm2: the film over {area:g} cm2 of plaque holds {film_ml:g} mL, "
            f"no less than the resting volume; the mouth's other surfaces hold saliva too"
        )
    mouth = Mouth(
        saliva_mol_per_m3=saliva,
        stimulated_saliva_mol_per_m3=stimulated,
        resting_volume_ml=resting,
        swallow_volume_ml=swallow,
        unstimulated_flow_ml_per_min=unstimulated,
        stimulated_flow_ml_per_min=stimulated_flow,
        stimulus=stimulus,
        stimulus_half_mol_per_m3=half,
        plaque_area_cm2=area,
        initial_mol_per_m3=initial,
    )
    diet = _read_diet(values.get("diet", []), grid, network, f"{where}.diet")
    return {"film": film, "mouth": mouth, "diet": diet}


def _read_diet(items: list[Any], grid: Grid, network: Network, where: str) -> tuple[Intake, ...]:
    """The intakes, in order and one at a time, each with what it brings for its kind."""
    diet: list[Intake] = []
    for i, item in enumerate(items):
        here = f"{where}[{i}]"
        v = read_object(item, here, INTAKE_FIELDS)
        kind = v["kind"]
        start = float(plain(v["start_h"]))
        if start < 0:
            raise ConfigError(f"{here}.start_h: must not be negative")
        duration = _positive(v["duration_min"], f"{here}.duration_min")
        volume: float | None = None
        composition: dict[str, float] | None = None
        released: dict[str, float] | None = None
        if kind == "food":
            for key in ("volume_ml", "composition_mol_per_m3"):
                if key in v:
                    raise ConfigError(
                        f"{here}.{key}: food adds no liquid; give what it releases into the "
                        "saliva as released_mmol"
                    )
            if "released_mmol" not in v:
                raise ConfigError(f"{here}: food needs released_mmol, what it releases")
            released = _composition(
                v["released_mmol"], network, f"{here}.released_mmol", "food releases"
            )
        else:
            if "released_mmol" in v:
                raise ConfigError(
                    f"{here}.released_mmol: a {kind} brings what it holds in its liquid; give "
                    "it as composition_mol_per_m3"
                )
            if "volume_ml" not in v:
                raise ConfigError(f"{here}: a {kind} needs volume_ml")
            volume = _positive(v["volume_ml"], f"{here}.volume_ml")
            composition = _composition(
                v.get("composition_mol_per_m3", {}),
                network,
                f"{here}.composition_mol_per_m3",
                f"a {kind} holds",
            )
        mixing = DEFAULT_MIXING_PER_S
        if "mixing_per_s" in v:
            mixing = float(plain(v["mixing_per_s"]))
        if mixing < 0:
            raise ConfigError(f"{here}.mixing_per_s: must not be negative")
        retained = None
        if "retained" in v:
            retained = _read_retained(v["retained"], grid, network, f"{here}.retained")
        if diet and start < diet[-1].end_h - 1e-12 * max(1.0, diet[-1].end_h):
            raise ConfigError(
                f"{here}.start_h: starts at {start:g} h, before intake {i - 1} ends at "
                f"{diet[-1].end_h:g} h; list the intakes in order, one at a time"
            )
        diet.append(Intake(kind, start, duration, volume, composition, released, mixing, retained))
    return tuple(diet)


def _read_retained(raw: Any, grid: Grid, network: Network, where: str) -> Retained:
    """Food left on the teeth: particulate, released by a process, over part of the substratum."""
    v = read_object(raw, where, RETAINED_FIELDS)
    components = {c.name: c for c in network.components}
    name = v["component"]
    _known([name], components, f"{where}.component")
    if components[name].phase != "particulate":
        raise ConfigError(
            f"{where}.component: '{name}' is dissolved; food left on the teeth is particulate, "
            "and a process of the network releases what dissolves from it"
        )
    if not any(p.coefficient(name) < 0 for p in network.processes):
        raise ConfigError(
            f"{where}.component: no process consumes '{name}', so nothing would dissolve from "
            "it; declare one, such as a first-order reaction to the sugar it holds"
        )
    amount = _positive(v["amount_mol_per_m2"], f"{where}.amount_mol_per_m2")
    region = tuple(bound for axis in grid.lateral_axes for bound in (0.0, grid.size_um[axis]))
    if "region_um" in v:
        region = tuple(float(plain(x)) for x in v["region_um"])
    try:
        faces_within(grid, region)
    except ValueError as error:
        raise ConfigError(f"{where}.region_um: {error}") from None
    return Retained(name, amount, region)


def _same_formula(a: Any, b: Any) -> bool:
    return sorted(a.formula.counts) == sorted(b.formula.counts) and a.formula.charge == (
        b.formula.charge
    )


def _read_scene(values: dict[str, Any], grid: Grid, network: Network, where: str) -> dict[str, Any]:
    """The substratum, liquid, flow, suspension and adhesion, checked against each other."""
    components = {c.name: c for c in network.components}

    here = f"{where}.substratum"
    surface_values = read_object(values["substratum"], here, SUBSTRATUM_FIELDS)
    patches = []
    for i, item in enumerate(surface_values["patches"]):
        patch_values = read_object(item, f"{here}.patches[{i}]", PATCH_FIELDS)
        region = tuple(float(plain(x)) for x in patch_values["region_um"])
        patches.append(Patch(patch_values["material"], region))
    try:
        substratum = Substratum(grid, patches)
    except ValueError as error:
        raise ConfigError(f"{here}: {error}") from None
    surface = Surface(surface_values["conditioning_film"], tuple(patches))

    here = f"{where}.liquid"
    liquid_values = read_object(values["liquid"], here, LIQUID_FIELDS)
    temperature = float(plain(liquid_values["temperature_c"]))
    if not temperature > -273.15:
        raise ConfigError(f"{here}.temperature_c: must be above absolute zero, -273.15")
    viscosity = _positive(liquid_values["viscosity_mpa_s"], f"{here}.viscosity_mpa_s")
    liquid = Liquid(temperature, viscosity)

    here = f"{where}.flow"
    flow_values = read_object(values["flow"], here, FLOW_FIELDS)
    flow = Flow(
        _positive(flow_values["wall_shear_rate_per_s"], f"{here}.wall_shear_rate_per_s"),
        _positive(flow_values["distance_from_inlet_mm"], f"{here}.distance_from_inlet_mm"),
    )

    suspension: list[Suspension] = []
    used: dict[str, str] = {}
    if not values["suspension"]:
        raise ConfigError(f"{where}.suspension: name at least one species that binds")
    for i, item in enumerate(values["suspension"]):
        here = f"{where}.suspension[{i}]"
        s = read_object(item, here, SUSPENSION_FIELDS)
        reversible, attached = s["reversible"], s["attached"]
        _known([reversible, attached], components, here)
        if reversible == attached:
            raise ConfigError(f"{here}: reversible and attached must be two components")
        for role, name_ in (("reversible", reversible), ("attached", attached)):
            if components[name_].phase != "particulate":
                raise ConfigError(
                    f"{here}.{role}: '{name_}' is dissolved; bound cells are particulate"
                )
            if name_ in used:
                raise ConfigError(f"{here}.{role}: '{name_}' already binds in {used[name_]}")
            used[name_] = here
        if not _same_formula(components[reversible], components[attached]):
            raise ConfigError(
                f"{here}: '{reversible}' and '{attached}' must have the same formula, because "
                "locking turns one into the other and must conserve every element"
            )
        touching = [p.name for p in network.processes if p.coefficient(reversible) != 0]
        if touching:
            raise ConfigError(
                f"{here}.reversible: '{reversible}' takes part in {', '.join(touching)}; "
                "reversibly bound cells only detach or lock, and grow once locked"
            )
        suspension.append(
            Suspension(
                reversible,
                attached,
                _positive(s["cells_per_ml"], f"{here}.cells_per_ml"),
                _positive(s["cell_diameter_um"], f"{here}.cell_diameter_um"),
                _positive(s["carbon_fmol_per_cell"], f"{here}.carbon_fmol_per_cell"),
                _positive(s["blocked_area_um2"], f"{here}.blocked_area_um2"),
            )
        )

    adhesion: list[Adhesion] = []
    species = {s.attached for s in suspension}
    seen: set[tuple[str, str]] = set()
    for i, item in enumerate(values["adhesion"]):
        here = f"{where}.adhesion[{i}]"
        a = read_object(item, here, ADHESION_FIELDS)
        attached, material = a["attached"], a["material"]
        if attached not in species:
            raise ConfigError(
                f"{here}.attached: '{attached}' is not the attached component of any species "
                "in the suspension"
            )
        if material not in substratum.materials:
            raise ConfigError(
                f"{here}.material: '{material}' is not a material of the substratum "
                f"({', '.join(substratum.materials)})"
            )
        if (attached, material) in seen:
            raise ConfigError(f"{here}: '{attached}' on '{material}' is stated twice")
        seen.add((attached, material))
        efficiency = float(plain(a["efficiency"]))
        if not 0 <= efficiency <= 1:
            raise ConfigError(f"{here}.efficiency: a fraction, between 0 and 1")
        detachment = float(plain(a["detachment_per_h"]))
        locking = float(plain(a["locking_per_h"]))
        for field_name, rate in (("detachment_per_h", detachment), ("locking_per_h", locking)):
            if rate < 0:
                raise ConfigError(f"{here}.{field_name}: must not be negative")
        adhesion.append(Adhesion(attached, material, efficiency, detachment, locking))
    missing_pairs = [
        f"'{s}' on '{m}'"
        for s in sorted(species)
        for m in substratum.materials
        if (s, m) not in seen
    ]
    if missing_pairs:
        raise ConfigError(
            f"{where}.adhesion: every species needs its binding on every material; missing "
            f"{', '.join(missing_pairs)}"
        )
    return {
        "surface": surface,
        "liquid": liquid,
        "flow": flow,
        "suspension": tuple(suspension),
        "adhesion": tuple(adhesion),
    }
