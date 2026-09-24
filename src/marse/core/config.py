"""Validated experiment configuration.

Configuration is the boundary where a researcher's intent enters MARSE, so it
is also where units are fixed and mistakes are caught. Every numeric field
names its unit, every value is range-checked on construction, and an invalid
configuration raises before any computation starts rather than producing a
plausible-looking wrong answer.

Configurations are read from JSON. YAML support is a planned optional extra;
JSON keeps the core dependency-free and is unambiguous to checksum, which the
run manifest relies on (see :mod:`marse.core.provenance`).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "BiofilmConfig",
    "EnvironmentConfig",
    "ExperimentConfig",
    "OrganismConfig",
    "SubstrateConfig",
    "load_experiment",
]


class ConfigError(ValueError):
    """An experiment configuration is invalid. The message names the field."""


def _check(condition: bool, where: str, message: str) -> None:
    if not condition:
        raise ConfigError(f"{where}: {message}")


def _number(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{where}: expected a number, got {type(value).__name__}")
    number = float(value)
    if number != number or number in (float("inf"), float("-inf")):
        raise ConfigError(f"{where}: must be finite")
    return number


def _triple(value: Any, where: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ConfigError(f"{where}: expected three values [min, opt, max]")
    low, opt, high = (_number(v, f"{where}[{i}]") for i, v in enumerate(value))
    _check(low < opt < high, where, f"need min < opt < max, got {low} < {opt} < {high}")
    return low, opt, high


@dataclass(frozen=True, slots=True)
class EnvironmentConfig:
    """Physical and chemical context, held constant for the run."""

    temperature_c: float
    ph: float

    def __post_init__(self) -> None:
        _check(
            -20.0 <= self.temperature_c <= 150.0, "environment.temperature_c", "outside -20..150 C"
        )
        _check(0.0 <= self.ph <= 14.0, "environment.ph", "outside 0..14")


@dataclass(frozen=True, slots=True)
class SubstrateConfig:
    """The single limiting substrate shared by every organism in the run."""

    name: str
    initial_mm: float

    def __post_init__(self) -> None:
        _check(bool(self.name.strip()), "substrate.name", "must not be empty")
        _check(self.initial_mm >= 0.0, "substrate.initial_mm", "must be non-negative")


@dataclass(frozen=True, slots=True)
class OrganismConfig:
    """A species or strain, with every rate constant carrying its unit.

    ``cardinal_temperature_c`` and ``cardinal_ph`` are optional. When present,
    the maximum growth rate is scaled by the cardinal models
    (docs/theory.md, section 2); when absent, ``mu_opt_per_h`` is used directly
    and the manifest records that no environmental scaling was applied.
    """

    name: str
    initial_biomass_g_per_l: float
    mu_opt_per_h: float
    k_s_mm: float
    yield_g_per_mmol: float
    maintenance_mmol_per_g_per_h: float = 0.0
    decay_per_h: float = 0.0
    cardinal_temperature_c: tuple[float, float, float] | None = None
    cardinal_ph: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        where = f"organism '{self.name}'"
        _check(bool(self.name.strip()), "organism.name", "must not be empty")
        _check(
            self.initial_biomass_g_per_l > 0.0, where, "initial_biomass_g_per_l must be positive"
        )
        _check(self.mu_opt_per_h > 0.0, where, "mu_opt_per_h must be positive")
        _check(self.k_s_mm > 0.0, where, "k_s_mm must be positive")
        _check(self.yield_g_per_mmol > 0.0, where, "yield_g_per_mmol must be positive")
        _check(
            self.maintenance_mmol_per_g_per_h >= 0.0,
            where,
            "maintenance_mmol_per_g_per_h must be non-negative",
        )
        _check(self.decay_per_h >= 0.0, where, "decay_per_h must be non-negative")
        if self.cardinal_temperature_c is not None:
            low, opt, high = self.cardinal_temperature_c
            # The CTMI is only well posed when the optimum lies nearer the maximum.
            _check(
                opt >= (low + high) / 2,
                where,
                "cardinal_temperature_c needs opt >= (min + max) / 2 for the CTMI",
            )


@dataclass(frozen=True, slots=True)
class BiofilmConfig:
    """A flat biofilm of fixed biomass, solved to steady state.

    Adding this block changes what the experiment *is*. A batch run evolves a
    well-mixed culture over time; a biofilm run holds the biomass fixed and
    solves the depth profile the gradient supports, which has no time axis at
    all. The two are different computations, so ``duration_h`` and
    ``timestep_h`` are not required here and are refused if given, rather than
    being silently ignored.

    ``diffusivity_um2_per_h`` names its unit because the alternative is the
    mistake that matters most in this model: reference tables quote
    diffusivities per second while growth rates are per hour, and mixing them
    shortens the penetration depth sixtyfold while still looking plausible
    (docs/theory.md section 4.6).
    """

    thickness_um: float
    cells: int
    diffusivity_um2_per_h: float

    def __post_init__(self) -> None:
        _check(self.thickness_um > 0.0, "biofilm.thickness_um", "must be positive")
        _check(
            isinstance(self.cells, int) and self.cells >= 2,
            "biofilm.cells",
            "must be an integer of at least 2",
        )
        _check(
            self.diffusivity_um2_per_h > 0.0, "biofilm.diffusivity_um2_per_h", "must be positive"
        )


@dataclass(frozen=True, slots=True)
class ExperimentConfig:
    """A complete, runnable experiment: a batch culture or a biofilm profile."""

    experiment_id: str
    organisms: tuple[OrganismConfig, ...]
    environment: EnvironmentConfig
    substrate: SubstrateConfig
    seed: int
    duration_h: float | None = None
    timestep_h: float | None = None
    checkpoint_interval_h: float | None = None
    biofilm: BiofilmConfig | None = None
    description: str = ""
    assumptions: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        _check(bool(self.experiment_id.strip()), "experiment_id", "must not be empty")
        _check(len(self.organisms) > 0, "organisms", "at least one organism is required")
        names = [o.name for o in self.organisms]
        _check(len(names) == len(set(names)), "organisms", f"names must be unique, got {names}")
        _check(
            isinstance(self.seed, int) and self.seed >= 0, "seed", "must be a non-negative integer"
        )

        if self.biofilm is not None:
            for name in ("duration_h", "timestep_h", "checkpoint_interval_h"):
                _check(
                    getattr(self, name) is None,
                    name,
                    "does not apply to a biofilm profile, which is solved to steady state "
                    "and has no time axis; remove it",
                )
            return

        _check(self.duration_h is not None, "duration_h", "is required for a batch run")
        _check(self.timestep_h is not None, "timestep_h", "is required for a batch run")
        _check(self.duration_h > 0.0, "duration_h", "must be positive")
        _check(self.timestep_h > 0.0, "timestep_h", "must be positive")
        _check(
            self.timestep_h <= self.duration_h,
            "timestep_h",
            "must not exceed duration_h",
        )
        _check(
            self.checkpoint_interval_h is not None
            and self.checkpoint_interval_h >= self.timestep_h,
            "checkpoint_interval_h",
            "must be at least one timestep",
        )

    @property
    def kind(self) -> str:
        """``"biofilm_profile"`` or ``"batch"``. Recorded in the run manifest."""
        return "biofilm_profile" if self.biofilm is not None else "batch"

    @property
    def steps(self) -> int:
        """Number of timesteps; the last step is shortened to land on ``duration_h``.

        The ratio is rounded before taking the ceiling so that a duration which
        is an exact multiple of the timestep in decimal, such as 24 / 0.001,
        does not gain a spurious extra step from binary floating point.
        """
        if self.duration_h is None or self.timestep_h is None:
            raise ConfigError(f"{self.kind}: has no timestepping, so no step count")
        return math.ceil(round(self.duration_h / self.timestep_h, 9))

    def to_dict(self) -> dict[str, Any]:
        """Plain data, suitable for JSON and for checksumming."""
        return asdict(self)


def _organism_from_dict(raw: dict[str, Any], index: int) -> OrganismConfig:
    where = f"organisms[{index}]"
    if not isinstance(raw, dict):
        raise ConfigError(f"{where}: expected an object")
    known = {f for f in OrganismConfig.__dataclass_fields__}
    unknown = set(raw) - known
    if unknown:
        raise ConfigError(
            f"{where}: unknown field(s) {sorted(unknown)}; known fields are {sorted(known)}"
        )
    for required in (
        "name",
        "initial_biomass_g_per_l",
        "mu_opt_per_h",
        "k_s_mm",
        "yield_g_per_mmol",
    ):
        if required not in raw:
            raise ConfigError(f"{where}: missing required field '{required}'")
    return OrganismConfig(
        name=str(raw["name"]),
        initial_biomass_g_per_l=_number(
            raw["initial_biomass_g_per_l"], f"{where}.initial_biomass_g_per_l"
        ),
        mu_opt_per_h=_number(raw["mu_opt_per_h"], f"{where}.mu_opt_per_h"),
        k_s_mm=_number(raw["k_s_mm"], f"{where}.k_s_mm"),
        yield_g_per_mmol=_number(raw["yield_g_per_mmol"], f"{where}.yield_g_per_mmol"),
        maintenance_mmol_per_g_per_h=_number(
            raw.get("maintenance_mmol_per_g_per_h", 0.0), f"{where}.maintenance_mmol_per_g_per_h"
        ),
        decay_per_h=_number(raw.get("decay_per_h", 0.0), f"{where}.decay_per_h"),
        cardinal_temperature_c=(
            _triple(raw["cardinal_temperature_c"], f"{where}.cardinal_temperature_c")
            if raw.get("cardinal_temperature_c") is not None
            else None
        ),
        cardinal_ph=(
            _triple(raw["cardinal_ph"], f"{where}.cardinal_ph")
            if raw.get("cardinal_ph") is not None
            else None
        ),
    )


def _biofilm_from_dict(raw: Any) -> BiofilmConfig:
    if not isinstance(raw, dict):
        raise ConfigError("biofilm: expected an object")
    known = {f for f in BiofilmConfig.__dataclass_fields__}
    unknown = set(raw) - known
    if unknown:
        raise ConfigError(
            f"biofilm: unknown field(s) {sorted(unknown)}; known fields are {sorted(known)}"
        )
    for name in known:
        if name not in raw:
            raise ConfigError(f"biofilm: missing required field '{name}'")
    cells = raw["cells"]
    if isinstance(cells, bool) or not isinstance(cells, int):
        raise ConfigError("biofilm.cells: expected an integer")
    return BiofilmConfig(
        thickness_um=_number(raw["thickness_um"], "biofilm.thickness_um"),
        cells=cells,
        diffusivity_um2_per_h=_number(
            raw["diffusivity_um2_per_h"], "biofilm.diffusivity_um2_per_h"
        ),
    )


def experiment_from_dict(raw: dict[str, Any]) -> ExperimentConfig:
    """Build a validated :class:`ExperimentConfig` from plain data."""
    if not isinstance(raw, dict):
        raise ConfigError("experiment: expected a JSON object at the top level")
    known = {f for f in ExperimentConfig.__dataclass_fields__}
    unknown = set(raw) - known
    if unknown:
        raise ConfigError(
            f"experiment: unknown field(s) {sorted(unknown)}; known fields are {sorted(known)}"
        )

    # A round-tripped config carries every optional key with a null value, so
    # presence alone does not mean a field was set.
    def given(name: str) -> bool:
        return raw.get(name) is not None

    required = ["experiment_id", "organisms", "environment", "substrate", "seed"]
    if not given("biofilm"):  # a batch run needs a clock; a biofilm profile does not
        required += ["duration_h", "timestep_h"]
    for name in required:
        if not given(name):
            raise ConfigError(f"experiment: missing required field '{name}'")

    environment = raw["environment"]
    if not isinstance(environment, dict):
        raise ConfigError("environment: expected an object")
    substrate = raw["substrate"]
    if not isinstance(substrate, dict):
        raise ConfigError("substrate: expected an object")
    organisms = raw["organisms"]
    if not isinstance(organisms, list) or not organisms:
        raise ConfigError("organisms: expected a non-empty list")

    timestep = _number(raw["timestep_h"], "timestep_h") if given("timestep_h") else None
    checkpoint = raw["checkpoint_interval_h"] if given("checkpoint_interval_h") else timestep
    return ExperimentConfig(
        experiment_id=str(raw["experiment_id"]),
        organisms=tuple(_organism_from_dict(o, i) for i, o in enumerate(organisms)),
        environment=EnvironmentConfig(
            temperature_c=_number(environment.get("temperature_c"), "environment.temperature_c"),
            ph=_number(environment.get("ph"), "environment.ph"),
        ),
        substrate=SubstrateConfig(
            name=str(substrate.get("name", "")),
            initial_mm=_number(substrate.get("initial_mm"), "substrate.initial_mm"),
        ),
        duration_h=_number(raw["duration_h"], "duration_h") if given("duration_h") else None,
        timestep_h=timestep,
        checkpoint_interval_h=(
            _number(checkpoint, "checkpoint_interval_h") if checkpoint is not None else None
        ),
        biofilm=_biofilm_from_dict(raw["biofilm"]) if given("biofilm") else None,
        seed=int(raw["seed"]),
        description=str(raw.get("description", "")),
        assumptions=tuple(str(a) for a in raw.get("assumptions", ())),
    )


def load_experiment(path: str | Path) -> ExperimentConfig:
    """Read and validate an experiment configuration from a JSON file."""
    text = Path(path).read_text(encoding="utf-8")
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as error:
        raise ConfigError(
            f"{Path(path).name}: invalid JSON ({error.msg} at line {error.lineno})"
        ) from error
    return experiment_from_dict(raw)
