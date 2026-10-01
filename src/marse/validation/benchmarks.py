"""Source-linked benchmark cases and observed-vs-simulated comparisons.

Benchmark cases are validation contracts, not claims that a model is generally
true. Each case declares its source, units, comparison metric, and tolerance.
The registry starts with analytical limits and can later contain curated
organism experiments without changing the comparison machinery.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from marse.validation.analytical import exponential_growth, first_order_profile, logistic_growth

__all__ = [
    "BenchmarkCase",
    "BenchmarkError",
    "BenchmarkResult",
    "BenchmarkSuite",
    "analytical_benchmarks",
    "compare_benchmark",
]

Metric = str
_METRICS = {"rmse", "max_absolute", "max_relative"}


class BenchmarkError(ValueError):
    """A benchmark definition or comparison is invalid."""


def _finite_array(value: ArrayLike, where: str) -> NDArray[np.float64]:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0 or array.size == 0 or not np.all(np.isfinite(array)):
        raise BenchmarkError(f"{where}: expected a non-empty finite array")
    return array


@dataclass(frozen=True, slots=True)
class BenchmarkCase:
    """A reproducible comparison contract for one observable."""

    id: str
    title: str
    source_ids: tuple[str, ...]
    inputs: NDArray[np.float64]
    expected: NDArray[np.float64]
    unit: str
    metric: Metric = "rmse"
    tolerance: float = 1e-6
    relative_floor: float = 1e-12
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.id.strip() or not self.title.strip():
            raise BenchmarkError("benchmark id and title must not be empty")
        if not self.source_ids or any(not source.strip() for source in self.source_ids):
            raise BenchmarkError(f"benchmark '{self.id}': source_ids are required")
        inputs = _finite_array(self.inputs, f"benchmark '{self.id}'.inputs")
        expected = _finite_array(self.expected, f"benchmark '{self.id}'.expected")
        if inputs.shape != expected.shape:
            raise BenchmarkError(f"benchmark '{self.id}': inputs and expected shapes differ")
        if not self.unit.strip():
            raise BenchmarkError(f"benchmark '{self.id}': unit is required")
        if self.metric not in _METRICS:
            raise BenchmarkError(f"benchmark '{self.id}': unsupported metric '{self.metric}'")
        if self.tolerance < 0 or not np.isfinite(self.tolerance):
            raise BenchmarkError(
                f"benchmark '{self.id}': tolerance must be finite and non-negative"
            )
        if self.relative_floor <= 0 or not np.isfinite(self.relative_floor):
            raise BenchmarkError(f"benchmark '{self.id}': relative_floor must be positive")
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "expected", expected)


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """Comparison result with enough detail to audit a pass or failure."""

    benchmark_id: str
    metric: Metric
    score: float
    tolerance: float
    passed: bool
    residuals: NDArray[np.float64]

    def to_dict(self) -> dict[str, Any]:
        return {
            "benchmark_id": self.benchmark_id,
            "metric": self.metric,
            "score": self.score,
            "tolerance": self.tolerance,
            "passed": self.passed,
            "residuals": self.residuals.tolist(),
        }


@dataclass(frozen=True, slots=True)
class BenchmarkSuite:
    """A named collection of benchmark cases."""

    name: str
    cases: tuple[BenchmarkCase, ...]

    def __post_init__(self) -> None:
        ids = [case.id for case in self.cases]
        if not self.name.strip() or not self.cases:
            raise BenchmarkError("benchmark suite needs a name and at least one case")
        if len(ids) != len(set(ids)):
            raise BenchmarkError("benchmark ids must be unique within a suite")

    def evaluate(
        self, evaluator: Callable[[BenchmarkCase], ArrayLike]
    ) -> tuple[BenchmarkResult, ...]:
        return tuple(compare_benchmark(case, evaluator(case)) for case in self.cases)

    def write_json(self, path: str) -> None:
        payload = {
            "name": self.name,
            "cases": [
                {
                    "id": case.id,
                    "title": case.title,
                    "source_ids": list(case.source_ids),
                    "inputs": case.inputs.tolist(),
                    "expected": case.expected.tolist(),
                    "unit": case.unit,
                    "metric": case.metric,
                    "tolerance": case.tolerance,
                    "relative_floor": case.relative_floor,
                    "notes": case.notes,
                }
                for case in self.cases
            ],
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")


def compare_benchmark(case: BenchmarkCase, simulated: ArrayLike) -> BenchmarkResult:
    """Compare one simulation output against a benchmark's expected output."""
    actual = _finite_array(simulated, f"benchmark '{case.id}'.simulated")
    if actual.shape != case.expected.shape:
        raise BenchmarkError(
            f"benchmark '{case.id}': simulated shape {actual.shape} "
            f"does not match expected {case.expected.shape}"
        )
    residuals = actual - case.expected
    if case.metric == "rmse":
        score = float(np.sqrt(np.mean(residuals**2)))
    elif case.metric == "max_absolute":
        score = float(np.max(np.abs(residuals)))
    else:
        score = float(
            np.max(np.abs(residuals) / np.maximum(np.abs(case.expected), case.relative_floor))
        )
    return BenchmarkResult(
        case.id, case.metric, score, case.tolerance, score <= case.tolerance, residuals
    )


def analytical_benchmarks() -> BenchmarkSuite:
    """Return source-linked analytical limits used as the initial registry."""
    exponential_time = np.linspace(0.0, 10.0, 11)
    logistic_time = np.linspace(0.0, 12.0, 13)
    depth = np.linspace(0.0, 200.0, 21)
    return BenchmarkSuite(
        "analytical-v1",
        (
            BenchmarkCase(
                "exponential-growth-v1",
                "Exponential growth reference",
                ("analytical:exponential_growth",),
                exponential_time,
                np.asarray(exponential_growth(exponential_time, 0.01, 0.7)),
                "biomass_g_per_l",
                tolerance=1e-12,
                notes="Idealized unrestricted growth, not a biological universal.",
            ),
            BenchmarkCase(
                "logistic-growth-v1",
                "Logistic growth reference",
                ("analytical:logistic_growth",),
                logistic_time,
                np.asarray(logistic_growth(logistic_time, 0.01, 4.0, 1.3)),
                "biomass_g_per_l",
                tolerance=1e-12,
                notes="Idealized carrying-capacity limit.",
            ),
            BenchmarkCase(
                "first-order-oxygen-profile-v1",
                "First-order oxygen profile reference",
                ("analytical:first_order_profile",),
                depth,
                np.asarray(
                    first_order_profile(
                        depth,
                        thickness=200.0,
                        diffusivity=100.0,
                        surface=1.0,
                        max_uptake=0.5,
                        half_saturation=0.1,
                    )
                ),
                "concentration",
                tolerance=1e-12,
                notes=(
                    "Steady one-dimensional first-order limit. Depends on uptake only "
                    "through max_uptake / half_saturation (5 per h); a Monod model meets "
                    "it only when the solute stays far below half_saturation."
                ),
            ),
        ),
    )
