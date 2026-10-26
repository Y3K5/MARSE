"""Tests for source-linked benchmark definitions and comparisons."""

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from marse.spatial.diffusion import solve_steady_state
from marse.spatial.domain import Grid1D
from marse.validation.benchmarks import (
    BenchmarkCase,
    BenchmarkError,
    analytical_benchmarks,
    compare_benchmark,
)


def test_registry_machinery_accepts_an_exact_reference(tmp_path: Path):
    """Plumbing only: evaluation, comparison and serialisation work end to end.

    Each case is scored against its own expected values, so this test says
    nothing about whether any model agrees with any benchmark. That is what
    the simulator-scored tests below are for.
    """
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


def _first_order_case() -> BenchmarkCase:
    return next(c for c in analytical_benchmarks().cases if c.id == "first-order-oxygen-profile-v1")


def _solver_on_benchmark_depths(case: BenchmarkCase, cells: int) -> np.ndarray:
    """Run the validated Monod solver in the case's first-order limit.

    The reference depends on the uptake only through the first-order rate
    ``max_uptake / half_saturation`` (5 per hour here). Monod uptake reduces
    to it only when the solute stays far below the half-saturation constant,
    so the solver is run with ``K`` large enough that surface / K is 1e-8 and
    the rate held at the case's value. Run with the case's literal Monod
    values instead (surface / K = 10) and the profile is nearly zero order and
    rightly fails.
    """
    first_order_rate = 0.5 / 0.1
    half_saturation = 1e8
    grid = Grid1D(thickness=200.0, cells=cells)
    profile = solve_steady_state(
        grid,
        diffusivity=100.0,
        surface=1.0,
        max_uptake=first_order_rate * half_saturation,
        half_saturation=half_saturation,
    )
    # The surface is a boundary condition, not a solved node; supply it for z = 0.
    depths = np.concatenate(([0.0], grid.depths))
    concentration = np.concatenate(([1.0], profile.concentration))
    return np.interp(case.inputs, depths, concentration)


def test_first_order_benchmark_is_met_by_the_validated_solver():
    """The first registry case scored against simulator output, not itself.

    The registry's 1e-12 tolerance suits an exact reference only. A
    discretised solver is held instead to 1e-5 absolute at 0.1 um spacing,
    twice its measured error there. The order check below matters as much:
    a tolerance alone can be met by luck, a second-order error that falls
    fourfold per halving cannot.
    """
    case = dataclasses.replace(_first_order_case(), metric="max_absolute", tolerance=1e-5)
    result = compare_benchmark(case, _solver_on_benchmark_depths(case, cells=2000))
    assert result.passed, f"score {result.score:.3e} exceeds {case.tolerance:.0e}"

    coarse = compare_benchmark(case, _solver_on_benchmark_depths(case, cells=1000)).score
    fine = result.score
    order = np.log2(coarse / fine)
    assert 1.9 < order < 2.1, f"observed order {order:.2f}, expected second order"
