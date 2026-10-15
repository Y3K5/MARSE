"""Cases that record a missing score. They do not run a model."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from marse.bench.cases import CaseError, load_case
from marse.cli import main

CASES = Path("examples/benchmarks/limits")


def test_a_not_scored_metric_needs_a_statement(tmp_path: Path):
    marker = tmp_path / "not-run.json"
    marker.write_text("{}\n", encoding="utf-8")
    case = tmp_path / "bare.yml"
    case.write_text(
        "id: bare\ntitle: t\nadapter: marse\ncommand: none\nnotes: n\n"
        f"experiment: {marker.name}\nmetrics:\n  - id: limit\n    kind: not_scored\n",
        encoding="utf-8",
    )
    with pytest.raises(CaseError, match="statement"):
        load_case(case, tmp_path)


def test_peptidoglycan_and_antigen_selection_are_not_scored(tmp_path: Path):
    output = tmp_path / "bench"
    assert main(["bench", "run", str(CASES), "-o", str(output), "--root", "."]) == 1
    for case_id, phrase in (
        ("peptidoglycan", "No peptidoglycan reaction"),
        ("antigen-selection", "does not select an antigen"),
    ):
        payload = json.loads((output / case_id / "bench.json").read_text(encoding="utf-8"))
        assert payload["passed"] is False
        assert payload["engine"] == "not_run"
        assert payload["reason"] == "metric failed: limit"
        assert payload["metrics"][0]["passed"] is False
        assert phrase in payload["metrics"][0]["detail"]
        assert not (output / case_id / "manifest.json").exists()
    assert main(["bench", "report", str(output)]) == 1
