"""Running a network: the settings that turn a reaction network into an experiment.

A version 2 document becomes runnable when every process has a rate and the
document states what to run: an ``experiment_id``, how long to run
(``duration_h``) in what steps (``timestep_h``), and the initial amounts
(``initial_mol_per_m3``). ``marse run`` then integrates it:

- in a closed, well-mixed box (:mod:`marse.core.well_mixed`), or
- in space, when the document also has a ``domain``
  (:mod:`marse.schemas.domain`, :mod:`marse.core.reactive_transport`).

The run's manifest records this configuration so that ``marse replay`` can
reproduce it. Every component starts at the amount given, or at zero if none is
given, and the zeros are written out, so a manifest shows the whole initial
state.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from marse.core.config import ConfigError
from marse.schemas._reading import load_json, plain
from marse.schemas.domain import Domain, read_domain
from marse.schemas.network import Network, _known, read_document

__all__ = [
    "ReactiveTransportConfig",
    "WellMixedConfig",
    "experiment_from_dict",
    "load_experiment",
]

DEFAULT_RELATIVE_TOLERANCE = 1e-6
DEFAULT_RELATIVE_TOLERANCE_IN_SPACE = 1e-4
"""Measured to give an error near 5e-6 of each component's peak (docs/theory.md, section 9.8).

The error estimate is conservative by a factor of 20 to 50, and in space each
step costs a multigrid solve, so the looser default buys the same accuracy the
well-mixed default gives, at a practical cost."""
DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3 = 1e-9
"""Picomolar: far below any concentration that matters, so traces cost no accuracy.

With a femtomolar default, the controller demanded relative accuracy of
products still near zero and took thousands of needless substeps."""

type Tolerance = float | dict[str, float]
"""One absolute tolerance for every component, or one for each, by name."""


def _per_component(tolerance: Tolerance, names: tuple[str, ...]) -> float | NDArray[np.float64]:
    if isinstance(tolerance, dict):
        return np.array([tolerance[n] for n in names], dtype=float)
    return tolerance


@dataclass(frozen=True, slots=True)
class WellMixedConfig:
    """A reaction network with its initial state and clock, ready to run.

    Satisfies :class:`marse.core.provenance.RunConfig`, so a run is recorded in,
    and rebuilt from, a manifest like any other.
    """

    experiment_id: str
    network: Network
    initial_mol_per_m3: dict[str, float]
    duration_h: float
    timestep_h: float
    record_interval_h: float
    seed: int = 0
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE
    absolute_tolerance_mol_per_m3: Tolerance = DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3

    @property
    def kind(self) -> str:
        return "well_mixed"

    def absolute_tolerances(self) -> float | NDArray[np.float64]:
        """The absolute tolerance the engine takes: one number, or one per component in order."""
        return _per_component(self.absolute_tolerance_mol_per_m3, self.network.component_names)

    @property
    def steps(self) -> int:
        """Number of timesteps; the last is shortened to land on ``duration_h``.

        The ratio is rounded before the ceiling is taken, as in
        :attr:`marse.core.config.ExperimentConfig.steps`, so that a duration
        that is a decimal multiple of the timestep gains no spurious step.
        """
        return math.ceil(round(self.duration_h / self.timestep_h, 9))

    @property
    def record_every(self) -> int:
        """Steps between recorded rows of the trajectory."""
        return max(1, round(self.record_interval_h / self.timestep_h))

    def to_dict(self) -> dict[str, Any]:
        """The whole configuration, defaults and zero initial amounts written out."""
        return self.network.to_dict() | {
            "experiment_id": self.experiment_id,
            "initial_mol_per_m3": dict(self.initial_mol_per_m3),
            "duration_h": self.duration_h,
            "timestep_h": self.timestep_h,
            "record_interval_h": self.record_interval_h,
            "relative_tolerance": self.relative_tolerance,
            "absolute_tolerance_mol_per_m3": _written(self.absolute_tolerance_mol_per_m3),
            "seed": self.seed,
        }


@dataclass(frozen=True, slots=True)
class ReactiveTransportConfig:
    """A reaction network in a box of voxels over a surface, with its initial state and clock."""

    experiment_id: str
    network: Network
    domain: Domain
    initial_mol_per_m3: dict[str, float]
    duration_h: float
    timestep_h: float
    record_interval_h: float
    seed: int = 0
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE_IN_SPACE
    absolute_tolerance_mol_per_m3: Tolerance = DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3

    @property
    def kind(self) -> str:
        return "reactive_transport"

    def absolute_tolerances(self) -> float | NDArray[np.float64]:
        """The absolute tolerance the engine takes: one number, or one per component in order."""
        return _per_component(self.absolute_tolerance_mol_per_m3, self.network.component_names)

    @property
    def steps(self) -> int:
        """Number of timesteps; the last is shortened to land on ``duration_h``."""
        return math.ceil(round(self.duration_h / self.timestep_h, 9))

    @property
    def record_every(self) -> int:
        """Steps between recorded frames and rows."""
        return max(1, round(self.record_interval_h / self.timestep_h))

    def to_dict(self) -> dict[str, Any]:
        """The whole configuration, the domain included, defaults and zeros written out."""
        return self.network.to_dict() | {
            "experiment_id": self.experiment_id,
            "initial_mol_per_m3": dict(self.initial_mol_per_m3),
            "duration_h": self.duration_h,
            "timestep_h": self.timestep_h,
            "record_interval_h": self.record_interval_h,
            "relative_tolerance": self.relative_tolerance,
            "absolute_tolerance_mol_per_m3": _written(self.absolute_tolerance_mol_per_m3),
            "seed": self.seed,
            "domain": self.domain.to_dict(),
        }


type RunConfig = WellMixedConfig | ReactiveTransportConfig


def _written(tolerance: Tolerance) -> float | dict[str, float]:
    return dict(tolerance) if isinstance(tolerance, dict) else tolerance


def experiment_from_dict(raw: Any) -> RunConfig:
    """Read a runnable version 2 document; see :mod:`marse.schemas.network` for the rest.

    A document with a ``domain`` runs in space; one without runs well mixed.
    """
    network, values = read_document(raw)
    missing = [k for k in ("experiment_id", "duration_h", "timestep_h") if k not in values]
    if missing:
        listed = ", ".join(f"'{m}'" for m in missing)
        raise ConfigError(f"experiment: missing {listed}, which a network needs to run")
    unrated = [p.name for p in network.processes if p.rate is None]
    if unrated:
        listed = ", ".join(f"'{n}'" for n in unrated)
        raise ConfigError(f"experiment: every process needs a rate to run; {listed} has none")
    experiment_id = values["experiment_id"].strip()
    if not experiment_id:
        raise ConfigError("experiment.experiment_id: must not be empty")
    duration = float(values["duration_h"])
    timestep = float(values["timestep_h"])
    if not duration > 0:
        raise ConfigError("experiment.duration_h: must be positive")
    if not 0 < timestep <= duration:
        raise ConfigError("experiment.timestep_h: must be positive and at most duration_h")
    record = float(values.get("record_interval_h", timestep))
    if not timestep <= record <= duration:
        raise ConfigError(
            "experiment.record_interval_h: must lie between timestep_h and duration_h"
        )
    given = values.get("initial_mol_per_m3", {})
    _known(list(given), {c.name: c for c in network.components}, "experiment.initial_mol_per_m3")
    negative = [name for name, amount in given.items() if amount < 0]
    if negative:
        listed = ", ".join(f"'{n}'" for n in negative)
        raise ConfigError(f"experiment.initial_mol_per_m3: {listed} must not be negative")
    seed = values.get("seed", 0)
    if seed < 0:
        raise ConfigError("experiment.seed: must not be negative")
    in_space = "domain" in values
    default_relative = (
        DEFAULT_RELATIVE_TOLERANCE_IN_SPACE if in_space else DEFAULT_RELATIVE_TOLERANCE
    )
    relative = float(values.get("relative_tolerance", default_relative))
    if not 0 < relative < 1:
        raise ConfigError("experiment.relative_tolerance: must lie between 0 and 1")
    absolute = _absolute_tolerance(values, network)
    initial = {
        c.name: float(plain(given[c.name])) if c.name in given else 0.0 for c in network.components
    }
    settings = {
        "experiment_id": experiment_id,
        "network": network,
        "initial_mol_per_m3": initial,
        "duration_h": duration,
        "timestep_h": timestep,
        "record_interval_h": record,
        "seed": seed,
        "relative_tolerance": relative,
        "absolute_tolerance_mol_per_m3": absolute,
    }
    if in_space:
        return ReactiveTransportConfig(domain=read_domain(values["domain"], network), **settings)
    return WellMixedConfig(**settings)


def _absolute_tolerance(values: dict[str, Any], network: Network) -> Tolerance:
    """One number, or, from an object by component, one per component, the rest at the default.

    The object is written back complete, every component named, so a manifest
    states the tolerance each was held to.
    """
    where = "experiment.absolute_tolerance_mol_per_m3"
    given = values.get("absolute_tolerance_mol_per_m3", DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3)
    if not isinstance(given, dict):
        if not given > 0:
            raise ConfigError(f"{where}: must be positive")
        return float(given)
    _known(list(given), {c.name: c for c in network.components}, where)
    for name, amount in given.items():
        if not amount > 0:
            raise ConfigError(f"{where}.{name}: must be positive")
    return {
        c.name: float(given[c.name]) if c.name in given else DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3
        for c in network.components
    }


def load_experiment(path: str | Path) -> RunConfig:
    """Read a runnable version 2 document from a JSON file."""
    return experiment_from_dict(load_json(path))
