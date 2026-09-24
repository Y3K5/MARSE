"""Analytical references and source-linked benchmark contracts.

The validation suite is shipped with MARSE so anyone can re-run its analytical
and literature-backed cases. See docs/validation.md.
"""

from marse.validation.benchmarks import (
    BenchmarkCase,
    BenchmarkError,
    BenchmarkResult,
    BenchmarkSuite,
    analytical_benchmarks,
    compare_benchmark,
)

__all__ = [
    "BenchmarkCase",
    "BenchmarkError",
    "BenchmarkResult",
    "BenchmarkSuite",
    "analytical_benchmarks",
    "compare_benchmark",
]
