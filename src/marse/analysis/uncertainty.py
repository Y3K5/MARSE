"""Reproducible uncertainty propagation and screening sensitivity analysis.

This module operates on scalar parameter evaluators so it can wrap ecosystem
signatures, calibrated models, or other deterministic MARSE outputs without
coupling uncertainty logic to one simulator. It provides transparent sampling
and rank correlations, not a substitute for a posterior statistical model.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "ParameterRange",
    "SensitivityResult",
    "UncertaintyError",
    "UncertaintyResult",
    "propagate",
    "rank_sensitivity",
]


class UncertaintyError(ValueError):
    """An uncertainty or sensitivity request is invalid."""


@dataclass(frozen=True, slots=True)
class ParameterRange:
    """A finite uniform parameter interval."""

    name: str
    lower: float
    upper: float

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise UncertaintyError("parameter name must not be empty")
        if not np.isfinite((self.lower, self.upper)).all():
            raise UncertaintyError(f"parameter '{self.name}' bounds must be finite")
        if self.upper <= self.lower:
            raise UncertaintyError(f"parameter '{self.name}' needs lower < upper")


@dataclass(frozen=True, slots=True)
class UncertaintyResult:
    """Samples, evaluator outputs, and reproducible summary statistics."""

    parameters: tuple[str, ...]
    samples: NDArray[np.float64]
    outputs: NDArray[np.float64]
    seed: int
    summary: dict[str, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "parameters": list(self.parameters),
            "samples": self.samples.tolist(),
            "outputs": self.outputs.tolist(),
            "seed": self.seed,
            "summary": self.summary,
        }


@dataclass(frozen=True, slots=True)
class SensitivityResult:
    """Rank-based screening sensitivity for one scalar output."""

    output_name: str
    correlations: dict[str, float]
    sample_count: int
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "output_name": self.output_name,
            "correlations": self.correlations,
            "sample_count": self.sample_count,
            "warnings": list(self.warnings),
        }


def _validate_count(count: int) -> None:
    if isinstance(count, bool) or not isinstance(count, int) or count < 2:
        raise UncertaintyError("sample_count must be an integer of at least two")


def _rank(values: NDArray[np.float64]) -> NDArray[np.float64]:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    sorted_values = values[order]
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0
        start = end
    return ranks


def propagate(
    parameters: tuple[ParameterRange, ...],
    evaluator: Callable[[Mapping[str, float]], float],
    *,
    sample_count: int = 100,
    seed: int = 0,
    latin_hypercube: bool = True,
) -> UncertaintyResult:
    """Evaluate a deterministic model over reproducible bounded samples."""
    _validate_count(sample_count)
    if not parameters:
        raise UncertaintyError("at least one parameter range is required")
    names = tuple(parameter.name for parameter in parameters)
    if len(set(names)) != len(names):
        raise UncertaintyError("parameter names must be unique")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise UncertaintyError("seed must be a non-negative integer")
    rng = np.random.default_rng(seed)
    unit = rng.random((sample_count, len(parameters)))
    if latin_hypercube:
        for index in range(len(parameters)):
            unit[:, index] = (rng.permutation(sample_count) + unit[:, index]) / sample_count
    samples = np.empty_like(unit)
    for index, parameter in enumerate(parameters):
        samples[:, index] = parameter.lower + unit[:, index] * (parameter.upper - parameter.lower)
    outputs = np.array(
        [evaluator(dict(zip(names, row, strict=True))) for row in samples],
        dtype=float,
    )
    if outputs.shape != (sample_count,) or not np.all(np.isfinite(outputs)):
        raise UncertaintyError("evaluator must return one finite scalar per sample")
    summary = {
        "mean": float(outputs.mean()),
        "standard_deviation": float(outputs.std(ddof=1)),
        "minimum": float(outputs.min()),
        "maximum": float(outputs.max()),
        "q05": float(np.quantile(outputs, 0.05)),
        "q50": float(np.quantile(outputs, 0.50)),
        "q95": float(np.quantile(outputs, 0.95)),
    }
    return UncertaintyResult(names, samples, outputs, seed, summary)


def rank_sensitivity(
    result: UncertaintyResult, *, output_name: str = "output"
) -> SensitivityResult:
    """Compute Spearman-style rank correlations for screening sensitivity."""
    warnings: list[str] = []
    output_ranks = _rank(result.outputs)
    correlations: dict[str, float] = {}
    for index, name in enumerate(result.parameters):
        parameter = result.samples[:, index]
        if np.all(parameter == parameter[0]):
            correlations[name] = 0.0
            warnings.append(f"parameter '{name}' has no sampled variation")
            continue
        correlations[name] = float(np.corrcoef(_rank(parameter), output_ranks)[0, 1])
    if result.samples.shape[0] < 20:
        warnings.append("fewer than 20 samples: rank sensitivity is exploratory")
    return SensitivityResult(output_name, correlations, result.samples.shape[0], tuple(warnings))
