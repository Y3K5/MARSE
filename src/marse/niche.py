"""Niche response functions, capability records, and condition scans.

This layer deliberately stops at an interpretable phenotype boundary:
environmental conditions modify evidence-backed capabilities, and capabilities
produce effective rates. It does not infer biology from arbitrary DNA strings.
"""

from __future__ import annotations

import csv
import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from marse.microbes.cardinal import cardinal_ph, cardinal_temperature
from marse.microbes.growth import monod

__all__ = [
    "Capability",
    "CapabilityEvaluation",
    "NicheError",
    "NicheScan",
    "NicheScanResult",
    "NicheSpecies",
    "evaluate_capability",
    "load_niche_scan",
    "run_niche_scan",
]


class NicheError(ValueError):
    """A niche or capability definition is invalid."""


@dataclass(frozen=True, slots=True)
class Capability:
    """A parameterized physiological capability with provenance."""

    id: str
    maximum_rate_per_h: float
    substrate: str
    half_saturation: float
    temperature_c: tuple[float, float, float] | None = None
    ph: tuple[float, float, float] | None = None
    oxygen_half_saturation: float | None = None
    evidence_source: str = ""
    evidence_confidence: str = ""

    def __post_init__(self) -> None:
        if not self.id.strip() or not self.substrate.strip():
            raise NicheError("capability id and substrate must not be empty")
        if self.maximum_rate_per_h < 0 or self.half_saturation <= 0:
            raise NicheError("capability rate must be non-negative and half_saturation positive")
        if self.oxygen_half_saturation is not None and self.oxygen_half_saturation <= 0:
            raise NicheError("oxygen_half_saturation must be positive")
        for name, values in (("temperature_c", self.temperature_c), ("ph", self.ph)):
            if values is not None and not values[0] < values[1] < values[2]:
                raise NicheError(f"{name} must satisfy min < optimum < max")


@dataclass(frozen=True, slots=True)
class NicheSpecies:
    """A species assembled from one or more capabilities."""

    name: str
    capabilities: tuple[Capability, ...]

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.capabilities:
            raise NicheError("species needs a name and at least one capability")


@dataclass(frozen=True, slots=True)
class CapabilityEvaluation:
    """Effective rate and the factors that explain it."""

    capability_id: str
    rate_per_h: float
    factors: dict[str, float]
    limiting_factor: str


@dataclass(frozen=True, slots=True)
class NicheScan:
    """A Cartesian scan over scalar environmental axes."""

    experiment_id: str
    axes: dict[str, tuple[float, ...]]
    species: tuple[NicheSpecies, ...]


@dataclass(frozen=True, slots=True)
class NicheScanResult:
    scan: NicheScan
    rows: tuple[dict[str, Any], ...]

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(
            json.dumps(
                {"experiment_id": self.scan.experiment_id, "rows": list(self.rows)},
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return destination

    def write_csv(self, path: str | Path) -> Path:
        destination = Path(path)
        fieldnames = list(self.rows[0]) if self.rows else []
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.rows)
        return destination


def _factor(value: float, name: str) -> float:
    if not np.isfinite(value):
        raise NicheError(f"{name} must be finite")
    return float(np.clip(value, 0.0, 1.0))


def evaluate_capability(
    capability: Capability, conditions: dict[str, float]
) -> CapabilityEvaluation:
    """Evaluate one capability and identify its strongest limiting factor."""
    def condition(name: str) -> float:
        if name not in conditions:
            raise NicheError(f"missing condition '{name}'")
        return float(conditions[name])

    substrate = condition(capability.substrate)
    factors: dict[str, float] = {
        capability.substrate: _factor(
            float(monod(max(substrate, 0.0), 1.0, capability.half_saturation)),
            capability.substrate,
        )
    }
    if capability.temperature_c is not None:
        factors["temperature"] = _factor(
            float(cardinal_temperature(condition("temperature_c"), *capability.temperature_c)),
            "temperature",
        )
    if capability.ph is not None:
        factors["ph"] = _factor(float(cardinal_ph(condition("ph"), *capability.ph)), "ph")
    if capability.oxygen_half_saturation is not None:
        factors["oxygen"] = _factor(
            float(monod(max(condition("oxygen"), 0.0), 1.0, capability.oxygen_half_saturation)),
            "oxygen",
        )
    limiting_factor = min(factors, key=factors.get)
    rate = capability.maximum_rate_per_h * float(np.prod(tuple(factors.values())))
    return CapabilityEvaluation(capability.id, rate, factors, limiting_factor)


def _capability_from_dict(raw: dict[str, Any], index: int) -> Capability:
    try:
        return Capability(
            id=str(raw["id"]),
            maximum_rate_per_h=float(raw["maximum_rate_per_h"]),
            substrate=str(raw["substrate"]),
            half_saturation=float(raw["half_saturation"]),
            temperature_c=(
                tuple(float(value) for value in raw["temperature_c"])
                if raw.get("temperature_c") is not None
                else None
            ),
            ph=(
                tuple(float(value) for value in raw["ph"])
                if raw.get("ph") is not None
                else None
            ),
            oxygen_half_saturation=(
                float(raw["oxygen_half_saturation"])
                if raw.get("oxygen_half_saturation") is not None
                else None
            ),
            evidence_source=str(raw.get("evidence_source", "")),
            evidence_confidence=str(raw.get("evidence_confidence", "")),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise NicheError(f"capabilities[{index}]: invalid definition") from error


def load_niche_scan(path: str | Path) -> NicheScan:
    """Load a niche scan definition from JSON."""
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        species = tuple(
            NicheSpecies(
                name=str(item["name"]),
                capabilities=tuple(
                    _capability_from_dict(capability, j)
                    for j, capability in enumerate(item["capabilities"])
                ),
            )
            for item in raw["species"]
        )
        axes = {
            str(name): tuple(float(value) for value in values)
            for name, values in raw["axes"].items()
        }
        if not axes or not species:
            raise NicheError("scan needs axes and species")
        return NicheScan(str(raw["experiment_id"]), axes, species)
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        raise NicheError(f"invalid niche scan: {error}") from error


def run_niche_scan(scan: NicheScan) -> NicheScanResult:
    """Evaluate every species capability at every Cartesian axis combination."""
    axis_names = tuple(scan.axes)
    required_axes = {
        capability.substrate
        for species in scan.species
        for capability in species.capabilities
    }
    required_axes.update(
        condition
        for species in scan.species
        for capability in species.capabilities
        for condition, configured in (
            ("temperature_c", capability.temperature_c),
            ("ph", capability.ph),
            ("oxygen", capability.oxygen_half_saturation),
        )
        if configured is not None
    )
    missing_axes = sorted(required_axes.difference(scan.axes))
    if missing_axes:
        joined = ", ".join(missing_axes)
        raise NicheError(f"scan axes missing required conditions: {joined}")
    rows: list[dict[str, Any]] = []
    for values in itertools.product(*(scan.axes[name] for name in axis_names)):
        conditions = dict(zip(axis_names, values, strict=True))
        for species in scan.species:
            evaluations = [evaluate_capability(capability, conditions) for capability in species.capabilities]
            best = max(evaluations, key=lambda evaluation: evaluation.rate_per_h)
            rows.append(
                {
                    **conditions,
                    "species": species.name,
                    "capability": best.capability_id,
                    "growth_rate_per_h": best.rate_per_h,
                    "limiting_factor": best.limiting_factor,
                    "viable": best.rate_per_h > 0.0,
                }
            )
    return NicheScanResult(scan, tuple(rows))
