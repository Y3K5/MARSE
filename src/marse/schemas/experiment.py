"""Running a network: the settings that turn a reaction network into an experiment.

A version 2 document becomes runnable when every process has a rate and the
document states what to run: an ``experiment_id``, how long to run
(``duration_h``) in what steps (``timestep_h``), and the initial amounts
(``initial_mol_per_m3``). ``marse run`` then integrates it in a closed,
well-mixed box (:mod:`marse.core.well_mixed`), and the run's manifest records
this configuration so that ``marse replay`` can reproduce it.

Every component starts at the amount given, or at zero if none is given, and
the zeros are written out, so a manifest shows the whole initial state.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from marse.core.config import ConfigError
from marse.schemas._reading import load_json, plain
from marse.schemas.network import Network, _known, read_document

__all__ = ["WellMixedConfig", "experiment_from_dict", "load_experiment"]

DEFAULT_RELATIVE_TOLERANCE = 1e-6
DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3 = 1e-9
"""Picomolar: far below any concentration that matters, so traces cost no accuracy.

With a femtomolar default, the controller demanded relative accuracy of
products still near zero and took thousands of needless substeps."""


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
    absolute_tolerance_mol_per_m3: float = DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3

    @property
    def kind(self) -> str:
        return "well_mixed"

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
            "absolute_tolerance_mol_per_m3": self.absolute_tolerance_mol_per_m3,
            "seed": self.seed,
        }


def experiment_from_dict(raw: Any) -> WellMixedConfig:
    """Read a runnable version 2 document; see :mod:`marse.schemas.network` for the rest."""
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
    relative = float(values.get("relative_tolerance", DEFAULT_RELATIVE_TOLERANCE))
    if not 0 < relative < 1:
        raise ConfigError("experiment.relative_tolerance: must lie between 0 and 1")
    absolute = float(
        values.get("absolute_tolerance_mol_per_m3", DEFAULT_ABSOLUTE_TOLERANCE_MOL_PER_M3)
    )
    if not absolute > 0:
        raise ConfigError("experiment.absolute_tolerance_mol_per_m3: must be positive")
    return WellMixedConfig(
        experiment_id=experiment_id,
        network=network,
        initial_mol_per_m3={
            c.name: float(plain(given[c.name])) if c.name in given else 0.0
            for c in network.components
        },
        duration_h=duration,
        timestep_h=timestep,
        record_interval_h=record,
        seed=seed,
        relative_tolerance=relative,
        absolute_tolerance_mol_per_m3=absolute,
    )


def load_experiment(path: str | Path) -> WellMixedConfig:
    """Read a runnable version 2 document from a JSON file."""
    return experiment_from_dict(load_json(path))
