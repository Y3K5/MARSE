"""A small, explicit 2D ecosystem model.

The model is intentionally a foundation rather than a claim of whole-organism
realism. Species grow from nutrient fields, consume those fields, and diffuse
through a rectangular domain. A seeded mutation stream can switch a cell's
growth multiplier; this provides a reproducible hook for richer phenotype and
lineage models later.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from marse.core.seeds import SeedRegistry

__all__ = [
    "EcosystemConfig",
    "EcosystemError",
    "EcosystemFrame",
    "EcosystemResult",
    "EcosystemState",
    "NutrientConfig",
    "SeedRegion",
    "SpeciesConfig",
    "load_experiment",
    "run",
]


class EcosystemError(ValueError):
    """The ecosystem configuration or numerical state is invalid."""


def _number(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EcosystemError(f"{where}: expected a number")
    result = float(value)
    if not np.isfinite(result):
        raise EcosystemError(f"{where}: must be finite")
    return result


@dataclass(frozen=True, slots=True)
class SeedRegion:
    """A circular initial colony expressed in normalized domain coordinates."""

    x: float
    y: float
    radius: float
    biomass: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.x <= 1.0 or not 0.0 <= self.y <= 1.0:
            raise EcosystemError("seed region center coordinates must be in [0, 1]")
        if self.radius <= 0.0 or self.radius > 1.0:
            raise EcosystemError("seed region radius must be in (0, 1]")
        if self.biomass < 0.0:
            raise EcosystemError("seed region biomass must be non-negative")


@dataclass(frozen=True, slots=True)
class NutrientConfig:
    """One diffusing nutrient or metabolite field."""

    name: str
    initial: float
    diffusivity: float
    boundary_value: float | None = None
    boundary_edges: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise EcosystemError("nutrient.name: must not be empty")
        if self.initial < 0:
            raise EcosystemError("nutrient.initial: must be non-negative")
        if self.diffusivity < 0:
            raise EcosystemError("nutrient.diffusivity: must be non-negative")
        if self.boundary_value is not None and self.boundary_value < 0:
            raise EcosystemError("nutrient.boundary_value must be non-negative")
        valid_edges = {"top", "bottom", "left", "right"}
        if any(edge not in valid_edges for edge in self.boundary_edges):
            raise EcosystemError("nutrient.boundary_edges contains an unknown edge")
        if self.boundary_value is None and self.boundary_edges:
            raise EcosystemError("nutrient.boundary_edges requires boundary_value")


@dataclass(frozen=True, slots=True)
class SpeciesConfig:
    """Species-level growth and nutrient-use parameters."""

    name: str
    initial_biomass: float
    maximum_growth_per_h: float
    half_saturation: tuple[float, ...]
    yield_per_nutrient: tuple[float, ...]
    mutation_probability: float = 0.0
    mutation_growth_multiplier: float = 1.0
    seed_regions: tuple[SeedRegion, ...] = ()

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise EcosystemError("species.name: must not be empty")
        if self.initial_biomass < 0 or self.maximum_growth_per_h < 0:
            raise EcosystemError(f"species '{self.name}': biomass and growth must be non-negative")
        if len(self.half_saturation) == 0:
            raise EcosystemError(f"species '{self.name}': at least one nutrient is required")
        if len(self.half_saturation) != len(self.yield_per_nutrient):
            raise EcosystemError(f"species '{self.name}': nutrient parameter lengths differ")
        if any(value <= 0 for value in self.half_saturation + self.yield_per_nutrient):
            raise EcosystemError(f"species '{self.name}': nutrient parameters must be positive")
        if not 0 <= self.mutation_probability <= 1:
            raise EcosystemError(f"species '{self.name}': mutation probability must be in [0, 1]")
        if self.mutation_growth_multiplier <= 0:
            raise EcosystemError(f"species '{self.name}': mutation multiplier must be positive")


@dataclass(frozen=True, slots=True)
class EcosystemConfig:
    """Validated rectangular 2D experiment configuration."""

    experiment_id: str
    width: int
    height: int
    cell_size_um: float
    duration_h: float
    timestep_h: float
    seed: int
    nutrients: tuple[NutrientConfig, ...]
    species: tuple[SpeciesConfig, ...]
    mutation_interval_h: float | None = None

    def __post_init__(self) -> None:
        if not self.experiment_id.strip():
            raise EcosystemError("experiment_id: must not be empty")
        if self.width < 2 or self.height < 2:
            raise EcosystemError("width and height must each be at least 2")
        if self.cell_size_um <= 0 or self.duration_h <= 0 or self.timestep_h <= 0:
            raise EcosystemError("cell_size_um, duration_h and timestep_h must be positive")
        if self.timestep_h > self.duration_h:
            raise EcosystemError("timestep_h must not exceed duration_h")
        if not isinstance(self.seed, int) or isinstance(self.seed, bool) or self.seed < 0:
            raise EcosystemError("seed must be a non-negative integer")
        if not self.nutrients or not self.species:
            raise EcosystemError("at least one nutrient and one species are required")
        if len({n.name for n in self.nutrients}) != len(self.nutrients):
            raise EcosystemError("nutrient names must be unique")
        if len({s.name for s in self.species}) != len(self.species):
            raise EcosystemError("species names must be unique")
        nutrient_count = len(self.nutrients)
        if any(len(s.half_saturation) != nutrient_count for s in self.species):
            raise EcosystemError("each species must define every nutrient parameter")
        if self.mutation_interval_h is not None and self.mutation_interval_h <= 0:
            raise EcosystemError("mutation_interval_h must be positive")

    @property
    def steps(self) -> int:
        return int(np.ceil(round(self.duration_h / self.timestep_h, 9)))


@dataclass(frozen=True, slots=True)
class EcosystemState:
    """Canonical state arrays: biomass is species,y,x; nutrients are nutrient,y,x."""

    time_h: float
    step: int
    biomass: NDArray[np.float64]
    nutrients: NDArray[np.float64]
    mutations: NDArray[np.int64]


@dataclass(frozen=True, slots=True)
class EcosystemFrame:
    """One visualization-ready snapshot with JSON-compatible metadata."""

    time_h: float
    biomass: NDArray[np.float64]
    nutrients: NDArray[np.float64]
    mutations: NDArray[np.int64]

    def to_dict(
        self, species_names: tuple[str, ...], nutrient_names: tuple[str, ...]
    ) -> dict[str, Any]:
        return {
            "time_h": self.time_h,
            "species": {name: self.biomass[i].tolist() for i, name in enumerate(species_names)},
            "nutrients": {
                name: self.nutrients[i].tolist() for i, name in enumerate(nutrient_names)
            },
            "mutations": self.mutations.tolist(),
        }


@dataclass(frozen=True, slots=True)
class EcosystemResult:
    config: EcosystemConfig
    frames: tuple[EcosystemFrame, ...]
    final_state: EcosystemState

    def write_frames(self, path: str | Path) -> Path:
        destination = Path(path)
        payload = {
            "experiment_id": self.config.experiment_id,
            "width": self.config.width,
            "height": self.config.height,
            "cell_size_um": self.config.cell_size_um,
            "species": [s.name for s in self.config.species],
            "nutrients": [n.name for n in self.config.nutrients],
            "frames": [
                frame.to_dict(
                    tuple(s.name for s in self.config.species),
                    tuple(n.name for n in self.config.nutrients),
                )
                for frame in self.frames
            ],
        }
        destination.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8")
        return destination


def _tuple_numbers(raw: Any, where: str) -> tuple[float, ...]:
    if not isinstance(raw, list) or not raw:
        raise EcosystemError(f"{where}: expected a non-empty list")
    return tuple(_number(value, f"{where}[{i}]") for i, value in enumerate(raw))


def ecosystem_from_dict(raw: dict[str, Any]) -> EcosystemConfig:
    """Build a validated ecosystem configuration from JSON-compatible data."""
    if not isinstance(raw, dict):
        raise EcosystemError("experiment: expected a JSON object")
    known = set(EcosystemConfig.__dataclass_fields__)
    unknown = set(raw) - known
    if unknown:
        raise EcosystemError(f"experiment: unknown field(s) {sorted(unknown)}")
    for name in (
        "experiment_id",
        "width",
        "height",
        "cell_size_um",
        "duration_h",
        "timestep_h",
        "seed",
        "nutrients",
        "species",
    ):
        if name not in raw:
            raise EcosystemError(f"experiment: missing required field '{name}'")
    nutrients = raw["nutrients"]
    species = raw["species"]
    if not isinstance(nutrients, list) or not isinstance(species, list):
        raise EcosystemError("nutrients and species must be lists")
    parsed_nutrients = tuple(
        NutrientConfig(
            name=str(item.get("name", "")),
            initial=_number(item.get("initial"), f"nutrients[{i}].initial"),
            diffusivity=_number(item.get("diffusivity"), f"nutrients[{i}].diffusivity"),
            boundary_value=(
                _number(item["boundary_value"], f"nutrients[{i}].boundary_value")
                if item.get("boundary_value") is not None
                else None
            ),
            boundary_edges=tuple(str(edge) for edge in item.get("boundary_edges", ())),
        )
        if isinstance(item, dict)
        else (_ for _ in ()).throw(EcosystemError(f"nutrients[{i}]: expected an object"))
        for i, item in enumerate(nutrients)
    )
    nutrient_count = len(parsed_nutrients)
    parsed_species = []
    for i, item in enumerate(species):
        if not isinstance(item, dict):
            raise EcosystemError(f"species[{i}]: expected an object")
        parsed_species.append(
            SpeciesConfig(
                name=str(item.get("name", "")),
                initial_biomass=_number(
                    item.get("initial_biomass"),
                    f"species[{i}].initial_biomass",
                ),
                maximum_growth_per_h=_number(
                    item.get("maximum_growth_per_h"), f"species[{i}].maximum_growth_per_h"
                ),
                half_saturation=_tuple_numbers(
                    item.get("half_saturation"), f"species[{i}].half_saturation"
                ),
                yield_per_nutrient=_tuple_numbers(
                    item.get("yield_per_nutrient"), f"species[{i}].yield_per_nutrient"
                ),
                mutation_probability=_number(
                    item.get("mutation_probability", 0.0), f"species[{i}].mutation_probability"
                ),
                mutation_growth_multiplier=_number(
                    item.get("mutation_growth_multiplier", 1.0),
                    f"species[{i}].mutation_growth_multiplier",
                ),
                seed_regions=tuple(
                    SeedRegion(
                        x=_number(region.get("x"), f"species[{i}].seed_regions[{j}].x"),
                        y=_number(region.get("y"), f"species[{i}].seed_regions[{j}].y"),
                        radius=_number(
                            region.get("radius"),
                            f"species[{i}].seed_regions[{j}].radius",
                        ),
                        biomass=_number(
                            region.get("biomass"),
                            f"species[{i}].seed_regions[{j}].biomass",
                        ),
                    )
                    if isinstance(region, dict)
                    else (_ for _ in ()).throw(
                        EcosystemError(f"species[{i}].seed_regions[{j}]: expected an object")
                    )
                    for j, region in enumerate(item.get("seed_regions", ()))
                ),
            )
        )
        if len(parsed_species[-1].half_saturation) != nutrient_count:
            raise EcosystemError(f"species[{i}]: expected one nutrient parameter per nutrient")
    seed = raw["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise EcosystemError("seed: expected a non-negative integer")
    return EcosystemConfig(
        experiment_id=str(raw["experiment_id"]),
        width=int(raw["width"]),
        height=int(raw["height"]),
        cell_size_um=_number(raw["cell_size_um"], "cell_size_um"),
        duration_h=_number(raw["duration_h"], "duration_h"),
        timestep_h=_number(raw["timestep_h"], "timestep_h"),
        seed=seed,
        nutrients=parsed_nutrients,
        species=tuple(parsed_species),
        mutation_interval_h=(
            _number(raw["mutation_interval_h"], "mutation_interval_h")
            if raw.get("mutation_interval_h") is not None
            else None
        ),
    )


def load_experiment(path: str | Path) -> EcosystemConfig:
    """Read and validate an ecosystem experiment from JSON."""
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise EcosystemError(f"{Path(path).name}: invalid JSON ({error.msg})") from error
    return ecosystem_from_dict(raw)


def _laplacian(field: NDArray[np.float64]) -> NDArray[np.float64]:
    padded = np.pad(field, 1, mode="edge")
    return padded[1:-1, :-2] + padded[1:-1, 2:] + padded[:-2, 1:-1] + padded[2:, 1:-1] - 4.0 * field


def _initial_state(config: EcosystemConfig) -> EcosystemState:
    shape = (config.height, config.width)
    yy, xx = np.meshgrid(
        np.linspace(0.0, 1.0, config.height),
        np.linspace(0.0, 1.0, config.width),
        indexing="ij",
    )
    biomass = []
    for species in config.species:
        field = np.full(shape, species.initial_biomass, dtype=float)
        for region in species.seed_regions:
            mask = (xx - region.x) ** 2 + (yy - region.y) ** 2 <= region.radius**2
            field[mask] = region.biomass
        biomass.append(field)
    nutrients = []
    for nutrient in config.nutrients:
        field = np.full(shape, nutrient.initial, dtype=float)
        _apply_boundary(field, nutrient)
        nutrients.append(field)
    biomass = np.stack(biomass)
    nutrients = np.stack(nutrients)
    return EcosystemState(
        time_h=0.0,
        step=0,
        biomass=biomass,
        nutrients=nutrients,
        mutations=np.zeros((len(config.species), *shape), dtype=np.int64),
    )


def _apply_boundary(field: NDArray[np.float64], nutrient: NutrientConfig) -> None:
    """Apply fixed-concentration boundary edges in-place."""
    if nutrient.boundary_value is None:
        return
    value = nutrient.boundary_value
    if "top" in nutrient.boundary_edges:
        field[0, :] = value
    if "bottom" in nutrient.boundary_edges:
        field[-1, :] = value
    if "left" in nutrient.boundary_edges:
        field[:, 0] = value
    if "right" in nutrient.boundary_edges:
        field[:, -1] = value


def run(config: EcosystemConfig) -> EcosystemResult:
    """Run a deterministic ecosystem experiment and retain every frame."""
    max_diffusivity = max(n.diffusivity for n in config.nutrients)
    stability = config.timestep_h * max_diffusivity / config.cell_size_um**2
    if stability > 0.25:
        raise EcosystemError(
            f"timestep is unstable for explicit diffusion (CFL={stability:.3g}; maximum is 0.25)"
        )

    state = _initial_state(config)
    mutations = state.mutations.copy()
    seed_registry = SeedRegistry(config.seed)
    mutation_rng = seed_registry.stream("ecosystem.mutations")
    frames = [
        EcosystemFrame(0.0, state.biomass.copy(), state.nutrients.copy(), state.mutations.copy())
    ]
    mutation_interval = config.mutation_interval_h
    next_mutation = mutation_interval if mutation_interval is not None else np.inf

    for step in range(1, config.steps + 1):
        dt = min(config.timestep_h, config.duration_h - state.time_h)
        biomass = state.biomass.copy()
        nutrients = state.nutrients.copy()
        for nutrient_index, nutrient in enumerate(config.nutrients):
            nutrients[nutrient_index] += (
                dt
                * nutrient.diffusivity
                / config.cell_size_um**2
                * _laplacian(nutrients[nutrient_index])
            )
            _apply_boundary(nutrients[nutrient_index], nutrient)

        for species_index, species in enumerate(config.species):
            growth_factor = np.where(
                mutations[species_index] > 0, species.mutation_growth_multiplier, 1.0
            )
            if state.time_h + dt >= next_mutation:
                events = mutation_rng.random(growth_factor.shape) < species.mutation_probability
                new_events = events & (mutations[species_index] == 0)
                mutations[species_index] += new_events
                growth_factor = np.where(
                    new_events, species.mutation_growth_multiplier, growth_factor
                )

            limitation = np.ones_like(biomass[species_index])
            for nutrient_index, (half_saturation, yield_value) in enumerate(
                zip(species.half_saturation, species.yield_per_nutrient, strict=True)
            ):
                concentration = np.maximum(nutrients[nutrient_index], 0.0)
                limitation *= concentration / (half_saturation + concentration)
                nutrients[nutrient_index] -= (
                    dt
                    * species.maximum_growth_per_h
                    * limitation
                    * biomass[species_index]
                    / yield_value
                )
            biomass[species_index] *= np.exp(
                dt * species.maximum_growth_per_h * limitation * growth_factor
            )

        nutrients = np.maximum(nutrients, 0.0)
        biomass = np.maximum(biomass, 0.0)
        if not (np.all(np.isfinite(biomass)) and np.all(np.isfinite(nutrients))):
            raise EcosystemError(f"state became non-finite at step {step}")
        time_h = state.time_h + dt
        state = EcosystemState(time_h, step, biomass, nutrients, mutations.copy())
        frames.append(EcosystemFrame(time_h, biomass.copy(), nutrients.copy(), mutations.copy()))
        while next_mutation <= time_h + 1e-12:
            next_mutation += mutation_interval if mutation_interval is not None else np.inf
        if dt <= 0:
            break

    return EcosystemResult(config, tuple(frames), state)
