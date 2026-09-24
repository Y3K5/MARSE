"""Reproducible ensembles of ecosystem simulations.

An ensemble separates the parameter space being explored from the state space
produced by each run. Scenarios are content-addressed, duplicate definitions
are removed before execution, and results retain compact signatures so large
scans do not require loading every frame into memory.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from marse.ecosystem.model import ecosystem_from_dict, run

__all__ = [
    "RunCatalog",
    "RunRecord",
    "Scenario",
    "ScenarioBatch",
    "run_batch",
]


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _checksum(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _set_path(document: dict[str, Any], path: str, value: Any) -> dict[str, Any]:
    result = json.loads(json.dumps(document))
    parts = path.split(".")
    target: Any = result
    for part in parts[:-1]:
        if isinstance(target, dict):
            if part not in target:
                raise ValueError(f"parameter path '{path}' does not exist")
            target = target[part]
        elif isinstance(target, list) and part.isdigit():
            index = int(part)
            if index >= len(target):
                raise ValueError(f"parameter path '{path}' does not exist")
            target = target[index]
        else:
            raise ValueError(f"parameter path '{path}' does not address a container")
    final = parts[-1]
    if isinstance(target, dict) and final in target:
        target[final] = value
    elif isinstance(target, list) and final.isdigit() and int(final) < len(target):
        target[int(final)] = value
    else:
        raise ValueError(f"parameter path '{path}' does not exist")
    return result


@dataclass(frozen=True, slots=True)
class Scenario:
    """One fully materialized ecosystem configuration."""

    config: dict[str, Any]
    scenario_id: str

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Scenario:
        if not isinstance(config, dict):
            raise ValueError("scenario configuration must be an object")
        return cls(config=json.loads(_canonical(config)), scenario_id=f"SCN-{_checksum(config)[:16]}")


@dataclass(frozen=True, slots=True)
class ScenarioBatch:
    """A deduplicated, reproducible collection of scenarios."""

    scenarios: tuple[Scenario, ...]

    @classmethod
    def grid(
        cls, base_config: dict[str, Any], parameters: dict[str, tuple[Any, ...] | list[Any]]
    ) -> ScenarioBatch:
        if not parameters:
            return cls((Scenario.from_config(base_config),))
        names = tuple(parameters)
        values = tuple(tuple(parameters[name]) for name in names)
        if any(not options for options in values):
            raise ValueError("parameter grids must not contain empty value lists")
        scenarios = {}
        for combination in itertools.product(*values):
            config = base_config
            for name, value in zip(names, combination, strict=True):
                config = _set_path(config, name, value)
            scenario = Scenario.from_config(config)
            scenarios[scenario.scenario_id] = scenario
        return cls(tuple(scenarios.values()))


@dataclass(frozen=True, slots=True)
class RunRecord:
    scenario_id: str
    experiment_id: str
    signature: dict[str, float]
    status: str = "completed"
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "experiment_id": self.experiment_id,
            "signature": self.signature,
            "status": self.status,
            "error": self.error,
        }


@dataclass(frozen=True, slots=True)
class RunCatalog:
    """Compact index of completed and failed scenario runs."""

    records: tuple[RunRecord, ...]

    def query(self, *, status: str | None = None) -> tuple[RunRecord, ...]:
        if status is None:
            return self.records
        return tuple(record for record in self.records if record.status == status)

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(
            json.dumps({"runs": [record.to_dict() for record in self.records]}, indent=2) + "\n",
            encoding="utf-8",
        )
        return destination


def _execute(scenario: Scenario) -> RunRecord:
    try:
        config = ecosystem_from_dict(scenario.config)
        result = run(config)
        final = result.final_state
        signature: dict[str, float] = {
            f"biomass:{species.name}": float(final.biomass[index].sum())
            for index, species in enumerate(config.species)
        }
        signature.update(
            {
                f"nutrient:{nutrient.name}": float(final.nutrients[index].mean())
                for index, nutrient in enumerate(config.nutrients)
            }
        )
        for index, species in enumerate(config.species):
            signature[f"occupied:{species.name}"] = float(
                (final.biomass[index] > 0.01 * config.carrying_capacity).sum()
            )
        return RunRecord(scenario.scenario_id, config.experiment_id, signature)
    except (TypeError, ValueError, OSError) as error:
        return RunRecord(
            scenario.scenario_id,
            str(scenario.config.get("experiment_id", "")),
            {},
            status="failed",
            error=f"{type(error).__name__}: {error}",
        )


def run_batch(batch: ScenarioBatch, *, workers: int = 1) -> RunCatalog:
    """Execute scenarios independently, preserving batch order."""
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1:
        records = tuple(_execute(scenario) for scenario in batch.scenarios)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            records = tuple(executor.map(_execute, batch.scenarios))
    return RunCatalog(records)
