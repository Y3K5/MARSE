"""Score a finished run directory against the metrics a case declared.

The numbers come from ``manifest.json``, ``summary.json``, and ``report.html``.
Replay uses the same flattened outputs as ``marse replay``. An unverified
engine cannot produce a passing case; the metric rows are still recorded.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from marse.bench.cases import Case, MetricSpec
from marse.cli.replay import _comparable
from marse.core.provenance import Manifest
from marse.engines import registry
from marse.validation.benchmarks import BenchmarkCase, BenchmarkError, compare_benchmark

_SHA = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True, slots=True)
class MetricResult:
    """One metric after it has been scored."""

    id: str
    kind: str
    passed: bool
    score: float | None
    tolerance: float | None
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "detail": self.detail,
            "id": self.id,
            "kind": self.kind,
            "passed": self.passed,
            "score": self.score,
            "tolerance": self.tolerance,
        }


@dataclass(frozen=True, slots=True)
class Evaluation:
    """The case result. ``passed`` is false when the engine is unverified."""

    config_sha256: str
    metrics: tuple[MetricResult, ...]
    passed: bool
    reason: str


def evaluate(case: Case, output_dir: Path, *, unverified: bool) -> Evaluation:
    """Score ``output_dir``. Missing files fail the metric that needed them."""
    manifest = _read_json(output_dir / "manifest.json")
    summary = _read_json(output_dir / "summary.json")
    report_path = output_dir / "report.html"
    report = report_path.read_text(encoding="utf-8") if report_path.is_file() else None
    checksum = ""
    if manifest is not None and isinstance(manifest.get("config_sha256"), str):
        checksum = manifest["config_sha256"]
    results = tuple(
        _score(spec, output_dir, manifest, summary, report, checksum, unverified)
        for spec in case.metrics
    )
    if unverified:
        return Evaluation(checksum, results, False, "unverified engine")
    if all(result.passed for result in results):
        return Evaluation(checksum, results, True, "")
    failed = next(result.id for result in results if not result.passed)
    return Evaluation(checksum, results, False, f"metric failed: {failed}")


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        loaded: Any = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return loaded if isinstance(loaded, dict) else None


def _score(
    spec: MetricSpec,
    output_dir: Path,
    manifest: dict[str, Any] | None,
    summary: dict[str, Any] | None,
    report: str | None,
    checksum: str,
    unverified: bool,
) -> MetricResult:
    if spec.kind == "manifest_present":
        passed = manifest is not None and _SHA.fullmatch(checksum) is not None
        detail = "config_sha256 is present" if passed else "manifest is missing or unsigned"
        return _result(spec, passed, None, detail)
    if spec.kind == "summary_present":
        return _summary(spec, manifest, summary, checksum)
    if spec.kind == "report_present":
        passed = report is not None and bool(checksum) and checksum in report
        detail = "report names the checksum" if passed else "report.html is missing"
        return _result(spec, passed, None, detail)
    if spec.kind == "unverified_disclosed":
        return _disclosed(spec, summary, unverified)
    if spec.kind == "replay_matches":
        return _replay(spec, output_dir)
    if spec.kind == "balance_residual":
        return _balance(spec, manifest)
    return _arrays(spec)


def _summary(
    spec: MetricSpec,
    manifest: dict[str, Any] | None,
    summary: dict[str, Any] | None,
    checksum: str,
) -> MetricResult:
    if manifest is None or summary is None:
        return _result(spec, False, None, "summary.json or manifest.json is missing")
    same = (
        summary.get("layout_version") == 1
        and summary.get("experiment_id") == manifest.get("experiment_id")
        and summary.get("config_sha256") == checksum
        and bool(checksum)
    )
    detail = "summary matches the manifest" if same else "summary does not match the manifest"
    return _result(spec, same, None, detail)


def _disclosed(spec: MetricSpec, summary: dict[str, Any] | None, unverified: bool) -> MetricResult:
    flag = None if summary is None else summary.get("unverified")
    if not isinstance(flag, bool):
        return _result(
            spec, False, None, "summary.json does not say whether the engine is unverified"
        )
    passed = flag if unverified else flag is False
    if not passed:
        return _result(spec, False, None, "summary.json disagrees with the engine")
    detail = "unverified engine is disclosed" if unverified else "engine is not marked unverified"
    return _result(spec, True, None, detail)


def _replay(spec: MetricSpec, output_dir: Path) -> MetricResult:
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.is_file():
        return _result(spec, False, None, "manifest.json is missing")
    try:
        original = Manifest.read(manifest_path)
        engine = registry.for_kind(original.kind)
        completed = engine.replay(original, None)
        fresh = completed.result.manifest
        if fresh is None:
            return _result(spec, False, None, "replay produced no manifest")
        passed = _comparable(original.outputs) == _comparable(fresh.outputs)
    except (KeyError, TypeError, ValueError, OSError) as error:
        return _result(spec, False, None, str(error))
    detail = "replay reproduced the recorded outputs" if passed else "replay differed"
    return _result(spec, passed, None, detail)


def _balance(spec: MetricSpec, manifest: dict[str, Any] | None) -> MetricResult:
    if manifest is None:
        return _result(spec, False, None, "manifest.json is missing")
    outputs = manifest.get("outputs")
    balance = outputs.get("balance") if isinstance(outputs, dict) else None
    entry = balance.get(spec.quantity) if isinstance(balance, dict) else None
    if not isinstance(entry, dict) or "largest_relative_residual" not in entry:
        return _result(spec, False, None, f"no balance residual for {spec.quantity}")
    residual = float(entry["largest_relative_residual"])
    assert spec.tolerance is not None
    passed = residual <= spec.tolerance
    return _result(spec, passed, residual, f"{spec.quantity} residual {residual:.3e}")


def _arrays(spec: MetricSpec) -> MetricResult:
    if spec.expected is None or spec.simulated is None:
        return _result(spec, False, None, "expected and simulated arrays are required")
    if spec.metric is None or spec.unit is None or spec.tolerance is None:
        return _result(spec, False, None, "metric, unit, and tolerance are required")
    try:
        expected = np.asarray(spec.expected, dtype=float)
        benchmark = BenchmarkCase(
            spec.id,
            spec.id,
            (f"case:{spec.id}",),
            np.arange(expected.size, dtype=float),
            expected,
            spec.unit,
            metric=spec.metric,
            tolerance=spec.tolerance,
        )
        compared = compare_benchmark(benchmark, np.asarray(spec.simulated, dtype=float))
    except BenchmarkError as error:
        return _result(spec, False, None, str(error))
    return _result(spec, compared.passed, compared.score, f"{spec.metric} {compared.score:.3e}")


def _result(spec: MetricSpec, passed: bool, score: float | None, detail: str) -> MetricResult:
    return MetricResult(spec.id, spec.kind, passed, score, spec.tolerance, detail)
