"""The periodontal benchmark scores the column. It does not score version 1."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from marse.bench.cases import load_cases
from marse.cli import main

CASES = Path("examples/benchmarks/periodontal")
LEGACY = Path("examples/experiments/periodontal_pathogen_biofilm.json")
COLUMN = Path("examples/networks/periodontal_pathogen_v2.json")


def test_the_cases_declare_only_what_the_files_can_support():
    cases = {case.id: case for case in load_cases(CASES, Path(".").resolve())}
    column = cases["periodontal-v2"]
    legacy = cases["periodontal-v1"]
    assert column.command == "run"
    assert column.experiment == COLUMN.as_posix()
    assert {metric.kind for metric in column.metrics if metric.kind == "balance_residual"} == {
        "balance_residual"
    }
    assert {metric.quantity for metric in column.metrics if metric.kind == "balance_residual"} == {
        "carbon",
        "nitrogen",
        "electrons",
    }
    assert all(
        metric.tolerance == 1e-12 for metric in column.metrics if metric.kind == "balance_residual"
    )
    assert legacy.command == "ecosystem"
    assert legacy.experiment == LEGACY.as_posix()
    assert {metric.kind for metric in legacy.metrics} == {
        "manifest_present",
        "summary_present",
        "report_present",
        "unverified_disclosed",
    }
    assert "not the same experiment" in legacy.notes
    assert "not a pocket" in column.notes
    assert LEGACY.is_file()
    assert COLUMN.is_file()


def test_the_column_passes_and_version_one_cannot(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    output = tmp_path / "bench"
    code = main(["bench", "run", str(CASES), "-o", str(output), "--root", "."])
    captured = capsys.readouterr()
    assert code == 1
    assert "unverified" in captured.err
    assert "periodontal-v2  pass" in captured.out
    assert "periodontal-v1  fail" in captured.out

    column = json.loads((output / "periodontal-v2" / "bench.json").read_text(encoding="utf-8"))
    legacy = json.loads((output / "periodontal-v1" / "bench.json").read_text(encoding="utf-8"))
    assert column["passed"] is True
    assert column["engine"] == "reactive_transport"
    assert column["unverified"] is False
    assert column["container"] is None
    for metric in column["metrics"]:
        assert metric["passed"] is True
        if metric["kind"] == "balance_residual":
            assert metric["score"] <= metric["tolerance"] == 1e-12
    assert legacy["passed"] is False
    assert legacy["reason"] == "unverified engine"
    assert legacy["engine"] == "legacy_ecosystem_v1"
    assert legacy["unverified"] is True
    assert "balance_residual" not in {metric["kind"] for metric in legacy["metrics"]}
    assert "array_reference" not in {metric["kind"] for metric in legacy["metrics"]}
    assert legacy["metrics"][-1]["passed"] is True

    for case_id in ("periodontal-v2", "periodontal-v1"):
        record = (output / case_id / "bench.json").read_text(encoding="utf-8")
        page = (output / case_id / "bench.html").read_text(encoding="utf-8")
        assert str(tmp_path) not in record
        assert "not a result to cite" in page or "not a pocket" in page
    assert "not the same experiment" in (output / "periodontal-v1" / "bench.html").read_text(
        encoding="utf-8"
    )
    assert main(["bench", "report", str(output)]) == 1
