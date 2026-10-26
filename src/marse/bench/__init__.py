"""Benchmark cases, the in-process adapter, and scored reports.

A passing case has met the comparison it declares. It is not a claim that a
model is a calibrated description of a habitat. See
docs/benchmarks/phase2-plan.md.
"""

from marse.bench.cases import Case, CaseError, MetricSpec, load_cases

__all__ = ["Case", "CaseError", "MetricSpec", "load_cases"]
