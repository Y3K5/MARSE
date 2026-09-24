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

from marse.additives import AdditiveEffect, apply_effect
from marse.core.seeds import SeedRegistry
from marse.ecosystem.providers import EcosystemProviders, NoBoundary
from marse.immune import ImmuneInteraction, MolecularNeutralizer, apply_immune_pressure
from marse.microbes.cardinal import cardinal_ph, cardinal_temperature
from marse.microbes.growth import monod
from marse.niche import Capability

__all__ = [
    "AdditiveConfig",
    "ConditionConfig",
    "EcosystemConfig",
    "EcosystemError",
    "EcosystemFrame",
    "EcosystemProviders",
    "EcosystemResult",
    "EcosystemState",
    "NutrientConfig",
    "PhenotypeConfig",
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


ConditionConfig = NutrientConfig


@dataclass(frozen=True, slots=True)
class AdditiveConfig:
    """A diffusing small-molecule field with first-order decay."""

    name: str
    initial: float
    diffusivity: float
    decay_per_h: float = 0.0
    boundary_value: float | None = None
    boundary_edges: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        NutrientConfig(
            self.name,
            self.initial,
            self.diffusivity,
            self.boundary_value,
            self.boundary_edges,
        )
        if self.decay_per_h < 0:
            raise EcosystemError("additive.decay_per_h must be non-negative")


@dataclass(frozen=True, slots=True)
class PhenotypeConfig:
    """A quorum-activated state with hysteresis and explicit multipliers."""

    name: str
    activation_threshold: float
    deactivation_threshold: float
    growth_multiplier: float = 1.0
    spreading_multiplier: float = 1.0
    minimum_dwell_h: float = 0.0

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise EcosystemError("phenotype.name: must not be empty")
        if self.activation_threshold < 0 or self.deactivation_threshold < 0:
            raise EcosystemError("phenotype thresholds must be non-negative")
        if self.deactivation_threshold > self.activation_threshold:
            raise EcosystemError("phenotype needs deactivation_threshold <= activation_threshold")
        if self.growth_multiplier <= 0 or self.spreading_multiplier < 0:
            raise EcosystemError("phenotype multipliers must be positive/non-negative")
        if self.minimum_dwell_h < 0:
            raise EcosystemError("phenotype.minimum_dwell_h must be non-negative")


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
    spreading_per_h: float = 0.0
    production_per_nutrient: tuple[float, ...] = ()
    capabilities: tuple[Capability, ...] = ()
    additive_effects: tuple[AdditiveEffect, ...] = ()
    chemotaxis_field: str | None = None
    chemotaxis_sensitivity: float = 0.0
    phenotypes: tuple[PhenotypeConfig, ...] = ()
    quorum_signal: str | None = None
    quorum_secretion_per_h: float = 0.0
    adhesion_per_h: float = 0.0
    detachment_per_h: float = 0.0
    adhesion_edges: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise EcosystemError("species.name: must not be empty")
        if self.initial_biomass < 0 or self.maximum_growth_per_h < 0:
            raise EcosystemError(f"species '{self.name}': biomass and growth must be non-negative")
        if len(self.half_saturation) == 0:
            raise EcosystemError(f"species '{self.name}': at least one nutrient is required")
        if len(self.half_saturation) != len(self.yield_per_nutrient):
            raise EcosystemError(f"species '{self.name}': nutrient parameter lengths differ")
        production = self.production_per_nutrient or (0.0,) * len(self.half_saturation)
        if len(production) != len(self.half_saturation):
            raise EcosystemError(f"species '{self.name}': production parameters differ")
        if any(value <= 0 for value in self.half_saturation + self.yield_per_nutrient):
            raise EcosystemError(f"species '{self.name}': nutrient parameters must be positive")
        if any(value < 0 for value in production):
            raise EcosystemError(f"species '{self.name}': production must be non-negative")
        if not 0 <= self.mutation_probability <= 1:
            raise EcosystemError(f"species '{self.name}': mutation probability must be in [0, 1]")
        if self.mutation_growth_multiplier <= 0:
            raise EcosystemError(f"species '{self.name}': mutation multiplier must be positive")
        if self.spreading_per_h < 0:
            raise EcosystemError(f"species '{self.name}': spreading rate must be non-negative")
        if not np.isfinite(self.chemotaxis_sensitivity):
            raise EcosystemError(f"species '{self.name}': chemotaxis sensitivity must be finite")
        if self.chemotaxis_field is not None and not self.chemotaxis_field.strip():
            raise EcosystemError(f"species '{self.name}': chemotaxis field must not be empty")
        if self.quorum_signal is not None and not self.quorum_signal.strip():
            raise EcosystemError(f"species '{self.name}': quorum signal must not be empty")
        if self.quorum_secretion_per_h < 0:
            raise EcosystemError(f"species '{self.name}': quorum secretion must be non-negative")
        if self.adhesion_per_h < 0 or self.detachment_per_h < 0:
            raise EcosystemError(
                f"species '{self.name}': adhesion and detachment must be non-negative"
            )
        if any(edge not in {"top", "bottom", "left", "right"} for edge in self.adhesion_edges):
            raise EcosystemError(f"species '{self.name}': adhesion edges are invalid")
        if len({phenotype.name for phenotype in self.phenotypes}) != len(self.phenotypes):
            raise EcosystemError(f"species '{self.name}': phenotype names must be unique")
        if len({effect.additive for effect in self.additive_effects}) != len(self.additive_effects):
            raise EcosystemError(f"species '{self.name}': additive effects must be unique")

    @property
    def production_coefficients(self) -> tuple[float, ...]:
        """Production per unit biomass growth for each nutrient field."""
        return self.production_per_nutrient or (0.0,) * len(self.half_saturation)


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
    carrying_capacity: float = 1.0
    competition_coefficients: tuple[tuple[float, ...], ...] | None = None
    conditions: tuple[ConditionConfig, ...] = ()
    additives: tuple[AdditiveConfig, ...] = ()
    immune_interactions: tuple[ImmuneInteraction, ...] = ()
    immune_neutralizers: tuple[MolecularNeutralizer, ...] = ()

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
        if self.carrying_capacity <= 0:
            raise EcosystemError("carrying_capacity must be positive")
        if self.competition_coefficients is not None:
            count = len(self.species)
            if len(self.competition_coefficients) != count or any(
                len(row) != count for row in self.competition_coefficients
            ):
                raise EcosystemError("competition_coefficients must be a species-by-species matrix")
            if any(value < 0 for row in self.competition_coefficients for value in row):
                raise EcosystemError("competition coefficients must be non-negative")
        condition_names = {condition.name for condition in self.conditions}
        additive_names = {additive.name for additive in self.additives}
        if len(condition_names) != len(self.conditions):
            raise EcosystemError("condition names must be unique")
        nutrient_names = {nutrient.name for nutrient in self.nutrients}
        if nutrient_names & condition_names:
            raise EcosystemError("nutrient and condition names must be distinct")
        if (nutrient_names | condition_names) & additive_names:
            raise EcosystemError("nutrient, condition, and additive names must be distinct")
        field_names = nutrient_names | condition_names
        if len(additive_names) != len(self.additives):
            raise EcosystemError("additive names must be unique")
        species_names = {species.name for species in self.species}
        for interaction in self.immune_interactions:
            if interaction.species not in species_names:
                raise EcosystemError(
                    f"immune interaction references unknown species '{interaction.species}'"
                )
            if interaction.effector not in additive_names:
                raise EcosystemError(
                    f"immune interaction references unknown effector '{interaction.effector}'"
                )
        for neutralizer in self.immune_neutralizers:
            if neutralizer.effector not in additive_names:
                raise EcosystemError(
                    f"immune neutralizer references unknown effector '{neutralizer.effector}'"
                )
            if neutralizer.molecule not in additive_names:
                raise EcosystemError(
                    f"immune neutralizer references unknown molecule '{neutralizer.molecule}'"
                )
        for species in self.species:
            if species.chemotaxis_field is not None and species.chemotaxis_field not in (
                field_names | additive_names
            ):
                raise EcosystemError(
                    f"species '{species.name}' chemotaxis references unknown field "
                    f"'{species.chemotaxis_field}'"
                )
            if species.quorum_signal is not None and species.quorum_signal not in additive_names:
                raise EcosystemError(
                    f"species '{species.name}' quorum signal references unknown additive "
                    f"'{species.quorum_signal}'"
                )
            for capability in species.capabilities:
                required = {capability.substrate}
                if capability.temperature_c is not None:
                    required.add("temperature_c")
                if capability.ph is not None:
                    required.add("ph")
                if capability.oxygen_half_saturation is not None:
                    required.add("oxygen")
                if not required.issubset(field_names):
                    missing = ", ".join(sorted(required - field_names))
                    raise EcosystemError(
                        f"species '{species.name}' capability '{capability.id}' "
                        f"requires missing condition field(s): {missing}"
                    )
            for effect in species.additive_effects:
                if effect.additive not in additive_names:
                    raise EcosystemError(
                        f"species '{species.name}' effect references unknown additive "
                        f"'{effect.additive}'"
                    )

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
    conditions: NDArray[np.float64]
    additives: NDArray[np.float64]
    mutations: NDArray[np.int64]
    phenotype_indices: NDArray[np.int64] | None = None
    phenotype_dwell_h: NDArray[np.float64] | None = None


@dataclass(frozen=True, slots=True)
class EcosystemFrame:
    """One visualization-ready snapshot with JSON-compatible metadata."""

    time_h: float
    biomass: NDArray[np.float64]
    nutrients: NDArray[np.float64]
    conditions: NDArray[np.float64]
    additives: NDArray[np.float64]
    mutations: NDArray[np.int64]
    niche_rates: NDArray[np.float64]
    niche_limiting_factors: tuple[NDArray[np.str_], ...]
    phenotype_indices: NDArray[np.int64] | None = None

    def to_dict(
        self,
        species_names: tuple[str, ...],
        nutrient_names: tuple[str, ...],
        carrying_capacity: float,
        condition_names: tuple[str, ...] = (),
        additive_names: tuple[str, ...] = (),
    ) -> dict[str, Any]:
        return {
            "time_h": self.time_h,
            "species": {name: self.biomass[i].tolist() for i, name in enumerate(species_names)},
            "nutrients": {
                name: self.nutrients[i].tolist() for i, name in enumerate(nutrient_names)
            },
            "conditions": {
                name: self.conditions[i].tolist() for i, name in enumerate(condition_names)
            },
            "additives": {
                name: self.additives[i].tolist() for i, name in enumerate(additive_names)
            },
            "niche": {
                "effective_growth_rate_per_h": {
                    name: self.niche_rates[i].tolist() for i, name in enumerate(species_names)
                },
                "limiting_factor": {
                    name: self.niche_limiting_factors[i].tolist()
                    for i, name in enumerate(species_names)
                },
            },
            "mutations": self.mutations.tolist(),
            "phenotypes": (
                self.phenotype_indices.tolist() if self.phenotype_indices is not None else []
            ),
            "statistics": {
                "species_total_biomass": {
                    name: float(self.biomass[i].sum()) for i, name in enumerate(species_names)
                },
                "species_occupied_cells": {
                    name: int(np.count_nonzero(self.biomass[i] > 0.01 * carrying_capacity))
                    for i, name in enumerate(species_names)
                },
                "nutrient_mean": {
                    name: float(self.nutrients[i].mean()) for i, name in enumerate(nutrient_names)
                },
                "mutation_count": int(np.count_nonzero(self.mutations)),
            },
        }


@dataclass(frozen=True, slots=True)
class EcosystemResult:
    config: EcosystemConfig
    frames: tuple[EcosystemFrame, ...]
    final_state: EcosystemState
    provider_versions: dict[str, str] | None = None

    def write_frames(self, path: str | Path) -> Path:
        destination = Path(path)
        payload = {
            "experiment_id": self.config.experiment_id,
            "width": self.config.width,
            "height": self.config.height,
            "cell_size_um": self.config.cell_size_um,
            "species": [s.name for s in self.config.species],
            "nutrients": [n.name for n in self.config.nutrients],
            "conditions": [c.name for c in self.config.conditions],
            "additives": [a.name for a in self.config.additives],
            "providers": self.provider_versions or {},
            "frames": [
                frame.to_dict(
                    tuple(s.name for s in self.config.species),
                    tuple(n.name for n in self.config.nutrients),
                    self.config.carrying_capacity,
                    tuple(c.name for c in self.config.conditions),
                    tuple(a.name for a in self.config.additives),
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


def _capabilities(raw: Any, where: str) -> tuple[Capability, ...]:
    if not isinstance(raw, list):
        raise EcosystemError(f"{where}: expected a list")
    parsed = []
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            raise EcosystemError(f"{where}[{index}]: expected an object")
        try:
            parsed.append(
                Capability(
                    id=str(item["id"]),
                    maximum_rate_per_h=_number(
                        item["maximum_rate_per_h"], f"{where}[{index}].maximum_rate_per_h"
                    ),
                    substrate=str(item["substrate"]),
                    half_saturation=_number(
                        item["half_saturation"], f"{where}[{index}].half_saturation"
                    ),
                    temperature_c=(
                        tuple(
                            _number(value, f"{where}[{index}].temperature_c[{j}]")
                            for j, value in enumerate(item["temperature_c"])
                        )
                        if item.get("temperature_c") is not None
                        else None
                    ),
                    ph=(
                        tuple(
                            _number(value, f"{where}[{index}].ph[{j}]")
                            for j, value in enumerate(item["ph"])
                        )
                        if item.get("ph") is not None
                        else None
                    ),
                    oxygen_half_saturation=(
                        _number(
                            item["oxygen_half_saturation"],
                            f"{where}[{index}].oxygen_half_saturation",
                        )
                        if item.get("oxygen_half_saturation") is not None
                        else None
                    ),
                    evidence_source=str(item.get("evidence_source", "")),
                    evidence_confidence=str(item.get("evidence_confidence", "")),
                )
            )
        except (KeyError, TypeError, ValueError) as error:
            raise EcosystemError(f"{where}[{index}]: invalid capability") from error
    return tuple(parsed)


def _additive_effects(raw: Any, where: str) -> tuple[AdditiveEffect, ...]:
    if not isinstance(raw, list):
        raise EcosystemError(f"{where}: expected a list")
    try:
        return tuple(
            AdditiveEffect(
                additive=str(item["additive"]),
                minimum_multiplier=_number(
                    item["minimum_multiplier"], f"{where}[{i}].minimum_multiplier"
                ),
                maximum_multiplier=_number(
                    item["maximum_multiplier"], f"{where}[{i}].maximum_multiplier"
                ),
                half_effect=_number(item["half_effect"], f"{where}[{i}].half_effect"),
                coefficient=_number(item.get("coefficient", 1.0), f"{where}[{i}].coefficient"),
                direction=str(item.get("direction", "increasing")),
            )
            if isinstance(item, dict)
            else (_ for _ in ()).throw(EcosystemError(f"{where}[{i}]: expected an object"))
            for i, item in enumerate(raw)
        )
    except (KeyError, TypeError, ValueError) as error:
        raise EcosystemError(f"{where}: invalid additive effect") from error


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
                production_per_nutrient=(
                    _tuple_numbers(
                        item["production_per_nutrient"],
                        f"species[{i}].production_per_nutrient",
                    )
                    if item.get("production_per_nutrient") is not None
                    else ()
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
                spreading_per_h=_number(
                    item.get("spreading_per_h", 0.0), f"species[{i}].spreading_per_h"
                ),
                capabilities=_capabilities(
                    item.get("capabilities", []), f"species[{i}].capabilities"
                ),
                additive_effects=_additive_effects(
                    item.get("additive_effects", []), f"species[{i}].additive_effects"
                ),
                chemotaxis_field=(
                    str(item["chemotaxis_field"])
                    if item.get("chemotaxis_field") is not None
                    else None
                ),
                chemotaxis_sensitivity=_number(
                    item.get("chemotaxis_sensitivity", 0.0),
                    f"species[{i}].chemotaxis_sensitivity",
                ),
                quorum_signal=(
                    str(item["quorum_signal"]) if item.get("quorum_signal") is not None else None
                ),
                quorum_secretion_per_h=_number(
                    item.get("quorum_secretion_per_h", 0.0),
                    f"species[{i}].quorum_secretion_per_h",
                ),
                adhesion_per_h=_number(
                    item.get("adhesion_per_h", 0.0), f"species[{i}].adhesion_per_h"
                ),
                detachment_per_h=_number(
                    item.get("detachment_per_h", 0.0), f"species[{i}].detachment_per_h"
                ),
                adhesion_edges=tuple(str(edge) for edge in item.get("adhesion_edges", ())),
                phenotypes=tuple(
                    PhenotypeConfig(
                        name=str(phenotype.get("name", "")),
                        activation_threshold=_number(
                            phenotype.get("activation_threshold"),
                            f"species[{i}].phenotypes[{j}].activation_threshold",
                        ),
                        deactivation_threshold=_number(
                            phenotype.get("deactivation_threshold"),
                            f"species[{i}].phenotypes[{j}].deactivation_threshold",
                        ),
                        growth_multiplier=_number(
                            phenotype.get("growth_multiplier", 1.0),
                            f"species[{i}].phenotypes[{j}].growth_multiplier",
                        ),
                        spreading_multiplier=_number(
                            phenotype.get("spreading_multiplier", 1.0),
                            f"species[{i}].phenotypes[{j}].spreading_multiplier",
                        ),
                        minimum_dwell_h=_number(
                            phenotype.get("minimum_dwell_h", 0.0),
                            f"species[{i}].phenotypes[{j}].minimum_dwell_h",
                        ),
                    )
                    if isinstance(phenotype, dict)
                    else (_ for _ in ()).throw(
                        EcosystemError(f"species[{i}].phenotypes[{j}]: expected an object")
                    )
                    for j, phenotype in enumerate(item.get("phenotypes", ()))
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
        carrying_capacity=_number(raw.get("carrying_capacity", 1.0), "carrying_capacity"),
        competition_coefficients=(
            tuple(
                tuple(
                    _number(value, f"competition_coefficients[{i}][{j}]")
                    for j, value in enumerate(row)
                )
                for i, row in enumerate(raw["competition_coefficients"])
            )
            if raw.get("competition_coefficients") is not None
            else None
        ),
        conditions=tuple(
            NutrientConfig(
                name=str(item.get("name", "")),
                initial=_number(item.get("initial"), f"conditions[{i}].initial"),
                diffusivity=_number(item.get("diffusivity"), f"conditions[{i}].diffusivity"),
                boundary_value=(
                    _number(item["boundary_value"], f"conditions[{i}].boundary_value")
                    if item.get("boundary_value") is not None
                    else None
                ),
                boundary_edges=tuple(str(edge) for edge in item.get("boundary_edges", ())),
            )
            if isinstance(item, dict)
            else (_ for _ in ()).throw(EcosystemError(f"conditions[{i}]: expected an object"))
            for i, item in enumerate(raw.get("conditions", ()))
        ),
        additives=tuple(
            AdditiveConfig(
                name=str(item.get("name", "")),
                initial=_number(item.get("initial"), f"additives[{i}].initial"),
                diffusivity=_number(item.get("diffusivity"), f"additives[{i}].diffusivity"),
                decay_per_h=_number(item.get("decay_per_h", 0.0), f"additives[{i}].decay_per_h"),
                boundary_value=(
                    _number(item["boundary_value"], f"additives[{i}].boundary_value")
                    if item.get("boundary_value") is not None
                    else None
                ),
                boundary_edges=tuple(str(edge) for edge in item.get("boundary_edges", ())),
            )
            if isinstance(item, dict)
            else (_ for _ in ()).throw(EcosystemError(f"additives[{i}]: expected an object"))
            for i, item in enumerate(raw.get("additives", ()))
        ),
        immune_interactions=tuple(
            ImmuneInteraction(
                species=str(item.get("species", "")),
                effector=str(item.get("effector", "")),
                maximum_kill_per_h=_number(
                    item.get("maximum_kill_per_h"),
                    f"immune_interactions[{i}].maximum_kill_per_h",
                ),
                half_effect=_number(
                    item.get("half_effect"),
                    f"immune_interactions[{i}].half_effect",
                ),
                hill_coefficient=_number(
                    item.get("hill_coefficient", 1.0),
                    f"immune_interactions[{i}].hill_coefficient",
                ),
                susceptibility=_number(
                    item.get("susceptibility", 1.0),
                    f"immune_interactions[{i}].susceptibility",
                ),
            )
            if isinstance(item, dict)
            else (_ for _ in ()).throw(
                EcosystemError(f"immune_interactions[{i}]: expected an object")
            )
            for i, item in enumerate(raw.get("immune_interactions", ()))
        ),
        immune_neutralizers=tuple(
            MolecularNeutralizer(
                effector=str(item.get("effector", "")),
                molecule=str(item.get("molecule", "")),
                neutralization_fraction=_number(
                    item.get("neutralization_fraction"),
                    f"immune_neutralizers[{i}].neutralization_fraction",
                ),
                half_effect=_number(
                    item.get("half_effect"),
                    f"immune_neutralizers[{i}].half_effect",
                ),
                hill_coefficient=_number(
                    item.get("hill_coefficient", 1.0),
                    f"immune_neutralizers[{i}].hill_coefficient",
                ),
            )
            if isinstance(item, dict)
            else (_ for _ in ()).throw(
                EcosystemError(f"immune_neutralizers[{i}]: expected an object")
            )
            for i, item in enumerate(raw.get("immune_neutralizers", ()))
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
    conditions = []
    for condition in config.conditions:
        field = np.full(shape, condition.initial, dtype=float)
        _apply_boundary(field, condition)
        conditions.append(field)
    condition_array = np.stack(conditions) if conditions else np.empty((0, *shape), dtype=float)
    additives = []
    for additive in config.additives:
        field = np.full(shape, additive.initial, dtype=float)
        _apply_boundary(field, additive)
        additives.append(field)
    additive_array = np.stack(additives) if additives else np.empty((0, *shape), dtype=float)
    return EcosystemState(
        time_h=0.0,
        step=0,
        biomass=biomass,
        nutrients=nutrients,
        conditions=condition_array,
        additives=additive_array,
        mutations=np.zeros((len(config.species), *shape), dtype=np.int64),
        phenotype_indices=np.zeros((len(config.species), *shape), dtype=np.int64),
        phenotype_dwell_h=np.zeros((len(config.species), *shape), dtype=float),
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


def _niche_maps(
    species: SpeciesConfig,
    nutrients: NDArray[np.float64],
    conditions: NDArray[np.float64],
    nutrient_names: tuple[str, ...],
    condition_names: tuple[str, ...],
) -> tuple[NDArray[np.float64], NDArray[np.str_]]:
    shape = nutrients.shape[1:]
    if not species.capabilities:
        return np.full(shape, species.maximum_growth_per_h), np.full(shape, "", dtype="<U1")
    fields = {name: nutrients[index] for index, name in enumerate(nutrient_names)} | {
        name: conditions[index] for index, name in enumerate(condition_names)
    }
    rates = []
    limiting = []
    for capability in species.capabilities:
        factors = {
            capability.substrate: np.asarray(
                monod(fields[capability.substrate], 1.0, capability.half_saturation)
            )
        }
        if capability.temperature_c is not None:
            factors["temperature"] = np.asarray(
                cardinal_temperature(fields["temperature_c"], *capability.temperature_c)
            )
        if capability.ph is not None:
            factors["ph"] = np.asarray(cardinal_ph(fields["ph"], *capability.ph))
        if capability.oxygen_half_saturation is not None:
            factors["oxygen"] = np.asarray(
                monod(fields["oxygen"], 1.0, capability.oxygen_half_saturation)
            )
        names = tuple(factors)
        factor_array = np.stack(tuple(np.clip(factors[name], 0.0, 1.0) for name in names))
        rates.append(capability.maximum_rate_per_h * np.prod(factor_array, axis=0))
        limiting.append(np.asarray(names)[np.argmin(factor_array, axis=0)])
    rate_array = np.stack(rates)
    best = np.argmax(rate_array, axis=0)
    return np.take_along_axis(rate_array, best[None, ...], axis=0)[0], np.choose(
        best, np.stack(limiting)
    )


def _additive_multiplier(
    species: SpeciesConfig,
    additives: NDArray[np.float64],
    additive_names: tuple[str, ...],
) -> NDArray[np.float64]:
    multiplier = np.ones(additives.shape[1:], dtype=float)
    fields = {name: additives[index] for index, name in enumerate(additive_names)}
    for effect in species.additive_effects:
        multiplier *= apply_effect(effect, fields[effect.additive])
    return multiplier


def _update_phenotypes(
    indices: NDArray[np.int64],
    dwell: NDArray[np.float64],
    biomass: NDArray[np.float64],
    species: SpeciesConfig,
    *,
    dt: float,
    carrying_capacity: float,
    quorum_signal: NDArray[np.float64] | None = None,
) -> None:
    """Apply quorum hysteresis in-place for one species."""
    if not species.phenotypes:
        return
    previous = indices.copy()
    dwell += dt
    signal = (
        biomass / carrying_capacity if quorum_signal is None else np.maximum(quorum_signal, 0.0)
    )
    for phenotype_index, phenotype in enumerate(species.phenotypes, start=1):
        activate = signal >= phenotype.activation_threshold
        deactivate = signal <= phenotype.deactivation_threshold
        eligible = dwell >= phenotype.minimum_dwell_h
        if phenotype_index == 1:
            indices[(previous == 0) & activate & eligible] = phenotype_index
        else:
            indices[(previous == phenotype_index - 1) & activate & eligible] = phenotype_index
        indices[(previous == phenotype_index) & deactivate & eligible] = 0
    dwell[indices != previous] = 0.0


def _chemotaxis_step(
    biomass: NDArray[np.float64],
    signal: NDArray[np.float64],
    *,
    sensitivity: float,
    dt: float,
    cell_size_um: float,
) -> NDArray[np.float64]:
    """Move biomass up a signal gradient with conservative edge no-fluxes."""
    gradient_x = np.diff(signal, axis=1) / cell_size_um
    gradient_y = np.diff(signal, axis=0) / cell_size_um
    flux_x = sensitivity * gradient_x * np.where(gradient_x >= 0.0, biomass[:, :-1], biomass[:, 1:])
    flux_y = sensitivity * gradient_y * np.where(gradient_y >= 0.0, biomass[:-1, :], biomass[1:, :])
    divergence_x = np.zeros_like(biomass)
    divergence_y = np.zeros_like(biomass)
    divergence_x[:, 1:-1] = (flux_x[:, 1:] - flux_x[:, :-1]) / cell_size_um
    divergence_x[:, 0] = flux_x[:, 0] / cell_size_um
    divergence_x[:, -1] = -flux_x[:, -1] / cell_size_um
    divergence_y[1:-1, :] = (flux_y[1:, :] - flux_y[:-1, :]) / cell_size_um
    divergence_y[0, :] = flux_y[0, :] / cell_size_um
    divergence_y[-1, :] = -flux_y[-1, :] / cell_size_um
    divergence = divergence_x + divergence_y
    return np.maximum(biomass - dt * divergence, 0.0)


def _adhesion_step(
    biomass: NDArray[np.float64],
    *,
    adhesion_per_h: float,
    detachment_per_h: float,
    edges: tuple[str, ...],
    dt: float,
) -> NDArray[np.float64]:
    """Transfer biomass toward selected surfaces and remove detached edge biomass."""
    if not edges or (adhesion_per_h == 0.0 and detachment_per_h == 0.0):
        return biomass
    updated = biomass.copy()
    for edge in edges:
        if edge == "top":
            surface, interior = (0, slice(None)), (1, slice(None))
        elif edge == "bottom":
            surface, interior = (-1, slice(None)), (-2, slice(None))
        elif edge == "left":
            surface, interior = (slice(None), 0), (slice(None), 1)
        else:
            surface, interior = (slice(None), -1), (slice(None), -2)
        transfer = np.minimum(
            biomass[interior] * dt * adhesion_per_h,
            biomass[interior],
        )
        updated[interior] -= transfer
        updated[surface] += transfer
        updated[surface] *= np.exp(-dt * detachment_per_h)
    return np.maximum(updated, 0.0)


def run(config: EcosystemConfig, providers: EcosystemProviders | None = None) -> EcosystemResult:
    """Run a deterministic ecosystem experiment and retain every frame."""
    providers = providers or EcosystemProviders()
    max_diffusivity = max(
        (field.diffusivity for field in (*config.nutrients, *config.conditions, *config.additives)),
        default=0.0,
    )
    stability = config.timestep_h * max_diffusivity / config.cell_size_um**2
    max_spreading = max(
        (
            s.spreading_per_h
            * max((phenotype.spreading_multiplier for phenotype in s.phenotypes), default=1.0)
            for s in config.species
        ),
        default=0.0,
    )
    spreading_stability = config.timestep_h * max_spreading / config.cell_size_um**2
    if max(stability, spreading_stability) > 0.25:
        raise EcosystemError(
            "timestep is unstable for explicit transport "
            f"(CFL={max(stability, spreading_stability):.3g}; maximum is 0.25)"
        )

    state = _initial_state(config)
    mutations = state.mutations.copy()
    phenotype_indices = state.phenotype_indices.copy()
    phenotype_dwell_h = state.phenotype_dwell_h.copy()
    nutrient_names = tuple(n.name for n in config.nutrients)
    condition_names = tuple(c.name for c in config.conditions)
    additive_names = tuple(a.name for a in config.additives)
    competition = (
        np.asarray(config.competition_coefficients, dtype=float)
        if config.competition_coefficients is not None
        else np.zeros((len(config.species), len(config.species)))
    )
    seed_registry = SeedRegistry(config.seed)
    mutation_rng = seed_registry.stream("ecosystem.mutations")
    frames = [
        EcosystemFrame(
            0.0,
            state.biomass.copy(),
            state.nutrients.copy(),
            state.conditions.copy(),
            state.additives.copy(),
            state.mutations.copy(),
            tuple(
                _niche_maps(
                    species,
                    state.nutrients,
                    state.conditions,
                    nutrient_names,
                    condition_names,
                )[0]
                for species in config.species
            ),
            tuple(
                _niche_maps(
                    species,
                    state.nutrients,
                    state.conditions,
                    nutrient_names,
                    condition_names,
                )[1]
                for species in config.species
            ),
            phenotype_indices.copy(),
        )
    ]
    mutation_interval = config.mutation_interval_h
    next_mutation = mutation_interval if mutation_interval is not None else np.inf

    for step in range(1, config.steps + 1):
        dt = min(config.timestep_h, config.duration_h - state.time_h)
        biomass = state.biomass.copy()
        nutrients = state.nutrients.copy()
        conditions = state.conditions.copy()
        additives = state.additives.copy()
        for condition_index, condition in enumerate(config.conditions):
            conditions[condition_index] = providers.transport.advance(
                conditions[condition_index],
                diffusivity=condition.diffusivity,
                dt=dt,
                cell_size_um=config.cell_size_um,
                boundary=condition,
            )
        for additive_index, additive in enumerate(config.additives):
            additives[additive_index] = providers.transport.advance(
                additives[additive_index],
                diffusivity=additive.diffusivity,
                dt=dt,
                cell_size_um=config.cell_size_um,
                decay_per_h=additive.decay_per_h,
                boundary=additive,
            )
        for nutrient_index, nutrient in enumerate(config.nutrients):
            nutrients[nutrient_index] = providers.transport.advance(
                nutrients[nutrient_index],
                diffusivity=nutrient.diffusivity,
                dt=dt,
                cell_size_um=config.cell_size_um,
                boundary=nutrient,
            )
        additive_indices = {name: index for index, name in enumerate(additive_names)}
        for species_index, species in enumerate(config.species):
            if species.quorum_signal is not None and species.quorum_secretion_per_h:
                additives[additive_indices[species.quorum_signal]] += (
                    dt * species.quorum_secretion_per_h * biomass[species_index]
                )
        field_map = {
            **{name: nutrients[index] for index, name in enumerate(nutrient_names)},
            **{name: conditions[index] for index, name in enumerate(condition_names)},
            **{name: additives[index] for index, name in enumerate(additive_names)},
        }

        for species_index, species in enumerate(config.species):
            _update_phenotypes(
                phenotype_indices[species_index],
                phenotype_dwell_h[species_index],
                biomass[species_index],
                species,
                dt=dt,
                carrying_capacity=config.carrying_capacity,
                quorum_signal=(
                    field_map[species.quorum_signal] if species.quorum_signal is not None else None
                ),
            )
            phenotype_index = phenotype_indices[species_index]
            active_phenotypes = [phenotype for phenotype in species.phenotypes]
            growth_multiplier = np.ones_like(biomass[species_index])
            spreading_multiplier = np.ones_like(biomass[species_index])
            for index, phenotype in enumerate(active_phenotypes, start=1):
                growth_multiplier = np.where(
                    phenotype_index == index, phenotype.growth_multiplier, growth_multiplier
                )
                spreading_multiplier = np.where(
                    phenotype_index == index, phenotype.spreading_multiplier, spreading_multiplier
                )
            if species.spreading_per_h:
                biomass[species_index] = providers.biomass_transport.advance(
                    biomass[species_index],
                    diffusivity=species.spreading_per_h * spreading_multiplier,
                    dt=dt,
                    cell_size_um=config.cell_size_um,
                    boundary=NoBoundary(),
                )
                biomass[species_index] = np.maximum(biomass[species_index], 0.0)
            if species.chemotaxis_field is not None and species.chemotaxis_sensitivity != 0.0:
                biomass[species_index] = _chemotaxis_step(
                    biomass[species_index],
                    field_map[species.chemotaxis_field],
                    sensitivity=species.chemotaxis_sensitivity,
                    dt=dt,
                    cell_size_um=config.cell_size_um,
                )
            biomass[species_index] = _adhesion_step(
                biomass[species_index],
                adhesion_per_h=species.adhesion_per_h,
                detachment_per_h=species.detachment_per_h,
                edges=species.adhesion_edges,
                dt=dt,
            )
            growth_factor = np.where(
                mutations[species_index] > 0, species.mutation_growth_multiplier, 1.0
            )
            growth_factor *= growth_multiplier
            if state.time_h + dt >= next_mutation:
                events = mutation_rng.random(growth_factor.shape) < species.mutation_probability
                new_events = events & (mutations[species_index] == 0)
                mutations[species_index] += new_events
                growth_factor = np.where(
                    new_events, species.mutation_growth_multiplier, growth_factor
                )

            nutrient_limitation = np.ones_like(biomass[species_index])
            for nutrient_index, (half_saturation, yield_value) in enumerate(
                zip(species.half_saturation, species.yield_per_nutrient, strict=True)
            ):
                concentration = np.maximum(nutrients[nutrient_index], 0.0)
                nutrient_limitation *= concentration / (half_saturation + concentration)
                nutrients[nutrient_index] -= (
                    dt
                    * species.maximum_growth_per_h
                    * nutrient_limitation
                    * biomass[species_index]
                    / yield_value
                )
            competition_pressure = np.zeros_like(nutrient_limitation)
            for competitor_index in range(len(config.species)):
                competition_pressure += (
                    competition[species_index, competitor_index]
                    * biomass[competitor_index]
                    / config.carrying_capacity
                )
            niche_rate, _ = _niche_maps(
                species, nutrients, conditions, nutrient_names, condition_names
            )
            niche_multiplier = (
                niche_rate / species.maximum_growth_per_h
                if species.maximum_growth_per_h > 0
                else np.zeros_like(niche_rate)
            )
            limitation = nutrient_limitation * niche_multiplier
            limitation *= _additive_multiplier(species, additives, additive_names)
            limitation /= 1.0 + competition_pressure
            growth_rate = (
                species.maximum_growth_per_h * limitation * growth_factor * biomass[species_index]
            )
            for nutrient_index, production in enumerate(species.production_coefficients):
                nutrients[nutrient_index] += dt * production * growth_rate
            biomass[species_index] *= np.exp(
                dt
                * species.maximum_growth_per_h
                * limitation
                * growth_factor
                * np.maximum(1.0 - biomass[species_index] / config.carrying_capacity, 0.0)
            )

        additive_fields = {name: additives[index] for index, name in enumerate(additive_names)}
        species_indices = {species.name: index for index, species in enumerate(config.species)}
        for interaction in config.immune_interactions:
            target_index = species_indices[interaction.species]
            pressure = apply_immune_pressure(
                biomass[target_index],
                effectors=additive_fields,
                molecules=additive_fields,
                interaction=interaction,
                neutralizers=config.immune_neutralizers,
                dt=dt,
            )
            biomass[target_index] = pressure.biomass

        nutrients = np.maximum(nutrients, 0.0)
        biomass = np.maximum(biomass, 0.0)
        if not (np.all(np.isfinite(biomass)) and np.all(np.isfinite(nutrients))):
            raise EcosystemError(f"state became non-finite at step {step}")
        time_h = state.time_h + dt
        state = EcosystemState(
            time_h,
            step,
            biomass,
            nutrients,
            conditions,
            additives,
            mutations.copy(),
            phenotype_indices.copy(),
            phenotype_dwell_h.copy(),
        )
        niche_maps = tuple(
            _niche_maps(species, nutrients, conditions, nutrient_names, condition_names)
            for species in config.species
        )
        frames.append(
            EcosystemFrame(
                time_h,
                biomass.copy(),
                nutrients.copy(),
                conditions.copy(),
                additives.copy(),
                mutations.copy(),
                tuple(result[0] for result in niche_maps),
                tuple(result[1] for result in niche_maps),
                phenotype_indices.copy(),
            )
        )
        while next_mutation <= time_h + 1e-12:
            next_mutation += mutation_interval if mutation_interval is not None else np.inf
        if dt <= 0:
            break

    return EcosystemResult(config, tuple(frames), state, providers.versions)
