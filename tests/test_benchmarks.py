"""Tests for source-linked benchmark definitions and comparisons."""

from pathlib import Path

import numpy as np
import pytest

from marse.validation.benchmarks import (
    BenchmarkCase,
    BenchmarkError,
    analytical_benchmarks,
    compare_benchmark,
)


def test_analytical_registry_is_source_linked_and_passes_exact_outputs(tmp_path: Path):
    suite = analytical_benchmarks()
    results = suite.evaluate(lambda case: case.expected)
    assert len(results) == 3
    assert all(result.passed for result in results)
    suite.write_json(tmp_path / "benchmarks.json")
    assert (tmp_path / "benchmarks.json").is_file()


def test_comparison_reports_residuals_and_tolerance_failure():
    case = BenchmarkCase(
        "test",
        "Test benchmark",
        ("source:test",),
        np.array([0.0, 1.0]),
        np.array([1.0, 2.0]),
        "unit",
        metric="max_absolute",
        tolerance=0.1,
    )
    result = compare_benchmark(case, [1.0, 2.3])
    assert not result.passed
    assert result.score == pytest.approx(0.3)
    np.testing.assert_allclose(result.residuals, [0.0, 0.3])


def test_invalid_benchmark_context_is_rejected():
    with pytest.raises(BenchmarkError, match="source_ids"):
        BenchmarkCase("id", "title", (), [0.0], [1.0], "unit")
    with pytest.raises(BenchmarkError, match="shapes"):
        BenchmarkCase("id", "title", ("source",), [0.0], [1.0, 2.0], "unit")
    case = BenchmarkCase("id", "title", ("source",), [0.0], [1.0], "unit")
    with pytest.raises(BenchmarkError, match="shape"):
        compare_benchmark(case, [1.0, 2.0])
