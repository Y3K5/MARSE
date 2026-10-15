"""Load one YAML benchmark case, and refuse it before a run if it is incomplete."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

_COMMANDS = frozenset({"run", "ecosystem"})
_KINDS = frozenset(
    {
        "manifest_present",
        "summary_present",
        "report_present",
        "unverified_disclosed",
        "replay_matches",
        "balance_residual",
        "array_reference",
    }
)
_ARRAY_METRICS = frozenset({"rmse", "max_absolute", "max_relative"})
_ID = "id must be lowercase letters, digits, and hyphens, and start with a letter or digit"


class CaseError(ValueError):
    """A case file is missing a field or names a path outside the repository."""


@dataclass(frozen=True, slots=True)
class MetricSpec:
    """One declared comparison. Arrays are scored later; this object only stores them."""

    id: str
    kind: str
    tolerance: float | None = None
    quantity: str | None = None
    metric: str | None = None
    unit: str | None = None
    expected: tuple[float, ...] | None = None
    simulated: tuple[float, ...] | None = None


@dataclass(frozen=True, slots=True)
class Case:
    """A case file after it has loaded. ``experiment`` is relative to the repository root."""

    id: str
    title: str
    adapter: str
    experiment: str
    command: str
    metrics: tuple[MetricSpec, ...]
    notes: str


def load_cases(directory: Path, root: Path) -> tuple[Case, ...]:
    """Every ``*.yml`` and ``*.yaml`` file in ``directory``, in id order."""
    directory = Path(directory)
    if not directory.is_dir():
        raise CaseError(f"case directory not found: {directory}")
    paths = sorted(directory.glob("*.yml")) + sorted(directory.glob("*.yaml"))
    cases = tuple(load_case(path, root) for path in paths)
    ids = [case.id for case in cases]
    if len(ids) != len(set(ids)):
        raise CaseError("case ids must be unique in one directory")
    return tuple(sorted(cases, key=lambda case: case.id))


def load_case(path: Path, root: Path) -> Case:
    """Read one case. The experiment file must already exist inside ``root``."""
    raw = _mapping(path)
    case_id = _token(raw, "id", path)
    if not _safe_id(case_id):
        raise CaseError(f"{path}: {_ID}")
    title = _text(raw, "title", path)
    adapter = _text(raw, "adapter", path)
    command = _text(raw, "command", path)
    if command not in _COMMANDS:
        raise CaseError(f"{path}: command must be 'run' or 'ecosystem'")
    notes = _text(raw, "notes", path)
    experiment = _experiment(raw.get("experiment"), path, root)
    metrics_raw = raw.get("metrics")
    if not isinstance(metrics_raw, list) or not metrics_raw:
        raise CaseError(f"{path}: metrics must be a non-empty list")
    metrics = tuple(_metric(item, path) for item in metrics_raw)
    metric_ids = [metric.id for metric in metrics]
    if len(metric_ids) != len(set(metric_ids)):
        raise CaseError(f"{path}: metric ids must be unique")
    return Case(case_id, title, adapter, experiment, command, metrics, notes)


def _mapping(path: Path) -> dict[str, Any]:
    try:
        loaded: Any = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise CaseError(f"{path}: {error}") from error
    if not isinstance(loaded, dict):
        raise CaseError(f"{path}: a case file must be a mapping")
    return loaded


def _text(raw: dict[str, Any], key: str, path: Path) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise CaseError(f"{path}: {key} is required")
    return value.strip()


def _token(raw: dict[str, Any], key: str, path: Path) -> str:
    value = _text(raw, key, path)
    if any(character.isspace() for character in value):
        raise CaseError(f"{path}: {key} must be a single token")
    return value


def _safe_id(value: str) -> bool:
    if not value or len(value) > 64:
        return False
    first = value[0]
    if not (first.isdigit() or (first.isascii() and first.islower())):
        return False
    allowed = set("abcdefghijklmnopqrstuvwxyz0123456789-")
    return all(character in allowed for character in value)


def _experiment(value: object, path: Path, root: Path) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CaseError(f"{path}: experiment is required")
    relative = Path(value)
    # A leading slash is absolute on POSIX and anchored on Windows. Either one
    # would leave the repository root, so both are refused.
    if relative.is_absolute() or relative.anchor or ".." in relative.parts:
        raise CaseError(f"{path}: experiment must stay inside the repository")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise CaseError(f"{path}: experiment file not found: {value}")
    return relative.as_posix()


def _metric(item: object, path: Path) -> MetricSpec:
    if not isinstance(item, dict):
        raise CaseError(f"{path}: each metric must be a mapping")
    metric_id = item.get("id")
    kind = item.get("kind")
    if not isinstance(metric_id, str) or not metric_id.strip():
        raise CaseError(f"{path}: a metric needs an id")
    if kind not in _KINDS:
        raise CaseError(f"{path}: unsupported metric kind {kind!r}")
    tolerance = _optional_float(item.get("tolerance"), path, "tolerance")
    if kind in {"balance_residual", "array_reference"} and tolerance is None:
        raise CaseError(f"{path}: metric '{metric_id}' needs a tolerance")
    quantity = _optional_text(item.get("quantity"))
    metric_name = _optional_text(item.get("metric"))
    unit = _optional_text(item.get("unit"))
    expected = _optional_floats(item.get("expected"), path, "expected")
    simulated = _optional_floats(item.get("simulated"), path, "simulated")
    if kind == "balance_residual" and not quantity:
        raise CaseError(f"{path}: metric '{metric_id}' needs a quantity")
    if kind == "array_reference":
        _check_arrays(path, str(metric_id), metric_name, unit, expected, simulated)
    return MetricSpec(
        str(metric_id).strip(),
        str(kind),
        tolerance,
        quantity,
        metric_name,
        unit,
        expected,
        simulated,
    )


def _check_arrays(
    path: Path,
    metric_id: str,
    metric_name: str | None,
    unit: str | None,
    expected: tuple[float, ...] | None,
    simulated: tuple[float, ...] | None,
) -> None:
    if metric_name not in _ARRAY_METRICS:
        raise CaseError(
            f"{path}: metric '{metric_id}' needs metric rmse, max_absolute, or max_relative"
        )
    if not unit:
        raise CaseError(f"{path}: metric '{metric_id}' needs a unit")
    if expected is None or simulated is None or len(expected) != len(simulated):
        raise CaseError(
            f"{path}: metric '{metric_id}' needs expected and simulated arrays of the same length"
        )


def _optional_text(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise CaseError("a metric field that is present must be text")
    return value.strip()


def _optional_float(value: object, path: Path, key: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise CaseError(f"{path}: {key} must be a number")
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise CaseError(f"{path}: {key} must be finite and non-negative")
    return number


def _optional_floats(value: object, path: Path, key: str) -> tuple[float, ...] | None:
    if value is None:
        return None
    if not isinstance(value, list) or not value:
        raise CaseError(f"{path}: {key} must be a non-empty list of numbers")
    numbers: list[float] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, int | float):
            raise CaseError(f"{path}: {key} must contain only finite numbers")
        number = float(item)
        if not math.isfinite(number):
            raise CaseError(f"{path}: {key} must contain only finite numbers")
        numbers.append(number)
    return tuple(numbers)
