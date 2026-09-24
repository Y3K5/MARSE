"""Evidence-backed culture conditions and measurement metadata.

This module records the context around biological observations; it does not
turn a literature value into a universal optimum.  Raw observations and fitted
kinetics remain separate so a simulation can report which values were measured
and which were inferred.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "AgarProtocol",
    "CultureDataset",
    "CultureError",
    "CultureRecord",
    "EvidenceGrade",
    "KineticObservation",
    "MeasurementMethod",
    "MediumRecipe",
    "SourceRecord",
    "dataset_checksum",
    "load_culture_dataset",
]


class CultureError(ValueError):
    """A culture evidence record is invalid or incomplete."""


EvidenceGrade = str
_EVIDENCE_GRADES = {"A", "B", "C"}
_ENDPOINTS = {"cfu_ml", "od", "colony_area", "dry_mass", "metabolite"}
_MODELS = {"monod", "logistic", "gompertz", "baranyi", "ratkowsky", "none"}


def _text(value: Any, where: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CultureError(f"{where}: must be a non-empty string")
    return value.strip()


def _number(value: Any, where: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CultureError(f"{where}: expected a number")
    result = float(value)
    if not math.isfinite(result):
        raise CultureError(f"{where}: must be finite")
    if minimum is not None and result < minimum:
        raise CultureError(f"{where}: must be at least {minimum}")
    return result


@dataclass(frozen=True, slots=True)
class SourceRecord:
    """A bibliographic or standards source supporting one or more records."""

    id: str
    title: str
    organization_or_authors: str
    year: int
    source_type: str
    evidence_grade: EvidenceGrade
    doi: str = ""
    url: str = ""

    def __post_init__(self) -> None:
        _text(self.id, "source.id")
        _text(self.title, f"source '{self.id}'.title")
        _text(self.organization_or_authors, f"source '{self.id}'.organization_or_authors")
        if not isinstance(self.year, int) or not 1800 <= self.year <= 2200:
            raise CultureError(f"source '{self.id}'.year: expected 1800..2200")
        _text(self.source_type, f"source '{self.id}'.source_type")
        if self.evidence_grade not in _EVIDENCE_GRADES:
            raise CultureError(f"source '{self.id}'.evidence_grade: expected A, B, or C")
        if not self.doi and not self.url:
            raise CultureError(f"source '{self.id}': doi or url is required")


@dataclass(frozen=True, slots=True)
class MediumRecipe:
    """A named medium whose composition can be referenced by a culture record."""

    id: str
    name: str
    components: dict[str, str]
    ph: float | None = None
    source_ids: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        _text(self.id, "medium.id")
        _text(self.name, f"medium '{self.id}'.name")
        if not self.components:
            raise CultureError(f"medium '{self.id}'.components: at least one component is required")
        if self.ph is not None and not 0.0 <= self.ph <= 14.0:
            raise CultureError(f"medium '{self.id}'.ph: expected 0..14")
        if not self.source_ids:
            raise CultureError(f"medium '{self.id}'.source_ids: at least one source is required")


@dataclass(frozen=True, slots=True)
class AgarProtocol:
    """Physical setup for a solid or semi-solid culture assay."""

    agar_concentration_g_per_l: float
    plate_depth_mm: float
    volume_ml: float
    inoculum_method: str
    atmosphere: str

    def __post_init__(self) -> None:
        _number(
            self.agar_concentration_g_per_l,
            "agar.agar_concentration_g_per_l",
            minimum=0.0,
        )
        _number(self.plate_depth_mm, "agar.plate_depth_mm", minimum=0.0)
        if self.volume_ml <= 0.0:
            raise CultureError("agar.volume_ml: must be positive")
        _text(self.inoculum_method, "agar.inoculum_method")
        _text(self.atmosphere, "agar.atmosphere")


@dataclass(frozen=True, slots=True)
class MeasurementMethod:
    """Definition of an endpoint and its calibration requirements."""

    id: str
    endpoint: str
    unit: str
    instrument_or_method: str
    calibration_id: str | None = None

    def __post_init__(self) -> None:
        _text(self.id, "measurement.id")
        if self.endpoint not in _ENDPOINTS:
            raise CultureError(f"measurement '{self.id}'.endpoint: unsupported endpoint")
        _text(self.unit, f"measurement '{self.id}'.unit")
        _text(self.instrument_or_method, f"measurement '{self.id}'.instrument_or_method")
        if self.endpoint == "od" and not self.calibration_id:
            raise CultureError(
                f"measurement '{self.id}': OD requires a condition-specific calibration_id"
            )


@dataclass(frozen=True, slots=True)
class KineticObservation:
    """A measured or fitted kinetic result with explicit uncertainty metadata."""

    measurement_id: str
    model: str
    growth_rate_per_h: float | None = None
    lag_h: float | None = None
    carrying_capacity: float | None = None
    uncertainty: float | None = None
    replicates: int = 1
    fitting_interval_h: tuple[float, float] | None = None
    response_transform: str = "none"
    error_model: str = "unspecified"

    def __post_init__(self) -> None:
        _text(self.measurement_id, "kinetics.measurement_id")
        if self.model not in _MODELS:
            raise CultureError(f"kinetics.model: unsupported model '{self.model}'")
        for name, value in (
            ("growth_rate_per_h", self.growth_rate_per_h),
            ("lag_h", self.lag_h),
            ("carrying_capacity", self.carrying_capacity),
            ("uncertainty", self.uncertainty),
        ):
            if value is not None:
                _number(value, f"kinetics.{name}", minimum=0.0)
        if not isinstance(self.replicates, int) or self.replicates < 1:
            raise CultureError("kinetics.replicates: must be a positive integer")
        if self.fitting_interval_h is not None:
            low, high = self.fitting_interval_h
            if low < 0.0 or high <= low:
                raise CultureError("kinetics.fitting_interval_h: need 0 <= start < end")


@dataclass(frozen=True, slots=True)
class CultureRecord:
    """One strain- and condition-specific culture observation."""

    id: str
    organism: str
    strain: str | None
    medium_id: str
    format: str
    temperature_c: float
    ph: float
    atmosphere: str
    inoculum: str
    incubation_h: float
    measurement_id: str
    source_ids: tuple[str, ...]
    agar: AgarProtocol | None = None
    kinetics: KineticObservation | None = None
    optimum_claim: bool = False
    assumptions: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        _text(self.id, "culture.id")
        organism = _text(self.organism, f"culture '{self.id}'.organism")
        if organism.lower() in {"all", "any", "universal", "*"}:
            raise CultureError(f"culture '{self.id}': universal organism claims are unsupported")
        if self.optimum_claim and (not self.strain or not self.medium_id):
            raise CultureError(
                f"culture '{self.id}': an optimum claim requires a strain and medium context"
            )
        _text(self.medium_id, f"culture '{self.id}'.medium_id")
        if self.format not in {"broth", "agar", "semi_solid", "biofilm"}:
            raise CultureError(f"culture '{self.id}'.format: unsupported culture format")
        if not -20.0 <= self.temperature_c <= 150.0:
            raise CultureError(f"culture '{self.id}'.temperature_c: outside -20..150 C")
        if not 0.0 <= self.ph <= 14.0:
            raise CultureError(f"culture '{self.id}'.ph: outside 0..14")
        _text(self.atmosphere, f"culture '{self.id}'.atmosphere")
        _text(self.inoculum, f"culture '{self.id}'.inoculum")
        if self.incubation_h <= 0.0:
            raise CultureError(f"culture '{self.id}'.incubation_h: must be positive")
        _text(self.measurement_id, f"culture '{self.id}'.measurement_id")
        if not self.source_ids:
            raise CultureError(f"culture '{self.id}'.source_ids: at least one source is required")
        if self.format in {"agar", "semi_solid"} and self.agar is None:
            raise CultureError(f"culture '{self.id}': agar protocol is required for {self.format}")


@dataclass(frozen=True, slots=True)
class CultureDataset:
    """Versioned, self-contained evidence set for culture simulations."""

    dataset_version: str
    sources: tuple[SourceRecord, ...]
    media: tuple[MediumRecipe, ...]
    measurements: tuple[MeasurementMethod, ...]
    cultures: tuple[CultureRecord, ...]

    def __post_init__(self) -> None:
        source_ids = {source.id for source in self.sources}
        medium_ids = {medium.id for medium in self.media}
        measurement_ids = {measurement.id for measurement in self.measurements}
        if len(source_ids) != len(self.sources):
            raise CultureError("sources: ids must be unique")
        if len(medium_ids) != len(self.media):
            raise CultureError("media: ids must be unique")
        if len(measurement_ids) != len(self.measurements):
            raise CultureError("measurements: ids must be unique")
        for medium in self.media:
            _references(medium.source_ids, source_ids, f"medium '{medium.id}'")
        for culture in self.cultures:
            _references(culture.source_ids, source_ids, f"culture '{culture.id}'")
            if culture.medium_id not in medium_ids:
                raise CultureError(f"culture '{culture.id}': unknown medium '{culture.medium_id}'")
            if culture.measurement_id not in measurement_ids:
                raise CultureError(
                    f"culture '{culture.id}': unknown measurement '{culture.measurement_id}'"
                )
            if (
                culture.kinetics is not None
                and culture.kinetics.measurement_id != culture.measurement_id
            ):
                raise CultureError(
                    f"culture '{culture.id}': kinetics measurement must match culture measurement"
                )
        for measurement in self.measurements:
            if measurement.calibration_id and measurement.calibration_id not in measurement_ids:
                raise CultureError(
                    f"measurement '{measurement.id}': unknown calibration "
                    f"'{measurement.calibration_id}'"
                )
        _text(self.dataset_version, "dataset_version")

    def to_dict(self) -> dict[str, Any]:
        return _as_json(self)

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return destination


def _references(references: tuple[str, ...], known: set[str], where: str) -> None:
    missing = sorted(set(references) - known)
    if missing:
        raise CultureError(f"{where}: unknown source ids {missing}")


def _as_json(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return {name: _as_json(getattr(value, name)) for name in value.__dataclass_fields__}
    if isinstance(value, tuple):
        return [_as_json(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _as_json(item) for key, item in value.items()}
    return value


def dataset_checksum(dataset: CultureDataset) -> str:
    """Return a stable SHA-256 checksum of the canonical evidence set."""
    canonical = json.dumps(dataset.to_dict(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _source(raw: dict[str, Any]) -> SourceRecord:
    return SourceRecord(**raw)


def _medium(raw: dict[str, Any]) -> MediumRecipe:
    return MediumRecipe(
        id=raw["id"],
        name=raw["name"],
        components=dict(raw["components"]),
        ph=raw.get("ph"),
        source_ids=tuple(raw.get("source_ids", ())),
    )


def _measurement(raw: dict[str, Any]) -> MeasurementMethod:
    return MeasurementMethod(**raw)


def _culture(raw: dict[str, Any]) -> CultureRecord:
    agar = raw.get("agar")
    kinetics = raw.get("kinetics")
    if kinetics is not None:
        kinetics = {
            **kinetics,
            "fitting_interval_h": (
                tuple(kinetics["fitting_interval_h"])
                if kinetics.get("fitting_interval_h") is not None
                else None
            ),
        }
    return CultureRecord(
        **{
            **raw,
            "source_ids": tuple(raw.get("source_ids", ())),
            "assumptions": tuple(raw.get("assumptions", ())),
            "agar": AgarProtocol(**agar) if agar is not None else None,
            "kinetics": KineticObservation(**kinetics) if kinetics is not None else None,
        }
    )


def load_culture_dataset(path: str | Path) -> CultureDataset:
    """Load and validate a versioned culture evidence dataset from JSON."""
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        return CultureDataset(
            dataset_version=raw["dataset_version"],
            sources=tuple(_source(item) for item in raw["sources"]),
            media=tuple(_medium(item) for item in raw["media"]),
            measurements=tuple(_measurement(item) for item in raw["measurements"]),
            cultures=tuple(_culture(item) for item in raw["cultures"]),
        )
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        if isinstance(error, CultureError):
            raise
        raise CultureError(f"invalid culture dataset: {error}") from error
