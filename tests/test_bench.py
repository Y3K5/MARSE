"""The benchmark harness scores a run directory. It does not add a biological claim."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from marse.bench.adapters import MarseAdapter, adapter_for
from marse.bench.cases import CaseError, load_case, load_cases
from marse.cli import main


def _batch(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "experiment_id": "bench-batch",
                "seed": 1,
                "duration_h": 0.01,
                "timestep_h": 0.01,
                "checkpoint_interval_h": 0.01,
                "environment": {"temperature_c": 37.0, "ph": 7.0},
                "substrate": {"name": "glucose", "initial_mm": 10.0},
                "organisms": [
                    {
                        "name": "solo",
                        "initial_biomass_g_per_l": 0.01,
                        "mu_opt_per_h": 0.8,
                        "k_s_mm": 0.4,
                        "yield_g_per_mmol": 0.08,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


def _case(directory: Path, name: str, metrics: str, experiment: str = "batch.json") -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / name).write_text(
        "\n".join(
            [
                "id: tiny-batch",
                "title: A short batch used to test the harness",
                "adapter: marse",
                f"experiment: {experiment}",
                "command: run",
                "notes: This checks the harness, not a habitat.",
                "metrics:",
                metrics,
            ]
        )
        + "\n",
        encoding="utf-8",
    )


_PASSING_METRICS = """
  - id: recorded
    kind: manifest_present
  - id: summary
    kind: summary_present
  - id: page
    kind: report_present
  - id: disclosed
    kind: unverified_disclosed
  - id: replay
    kind: replay_matches
  - id: reference
    kind: array_reference
    metric: max_absolute
    unit: biomass_g_per_l
    tolerance: 0.0
    expected: [1.0, 2.0]
    simulated: [1.0, 2.0]
"""


def test_a_case_missing_a_required_field_is_refused(tmp_path: Path):
    experiment = tmp_path / "batch.json"
    _batch(experiment)
    case = tmp_path / "broken.yml"
    case.write_text("id: tiny-batch\ntitle: Missing the rest\n", encoding="utf-8")
    with pytest.raises(CaseError, match="adapter"):
        load_case(case, tmp_path)
    metrics = "  - id: carbon\n    kind: balance_residual\n    quantity: C\n"
    _case(tmp_path / "cases", "no-tolerance.yml", metrics)
    with pytest.raises(CaseError, match="tolerance"):
        load_cases(tmp_path / "cases", tmp_path)
    outside = tmp_path / "outside.yml"
    outside.write_text(
        "id: tiny-batch\ntitle: t\nadapter: marse\ncommand: run\nnotes: n\n"
        "experiment: /tmp/secret.json\nmetrics:\n  - id: recorded\n    kind: manifest_present\n",
        encoding="utf-8",
    )
    with pytest.raises(CaseError, match="inside the repository"):
        load_case(outside, tmp_path)


def test_list_prints_cases_and_writes_nothing(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    _batch(tmp_path / "batch.json")
    _case(tmp_path / "cases", "tiny.yml", _PASSING_METRICS)
    assert main(["bench", "list", str(tmp_path / "cases"), "--root", str(tmp_path)]) == 0
    captured = capsys.readouterr()
    assert "tiny-batch  marse  batch.json" in captured.out
    assert not (tmp_path / "manifest.json").exists()
    assert list(tmp_path.rglob("bench.json")) == []


def test_run_writes_the_scored_report_for_a_tiny_batch(tmp_path: Path):
    _batch(tmp_path / "batch.json")
    _case(tmp_path / "cases", "tiny.yml", _PASSING_METRICS)
    output = tmp_path / "bench"
    code = main(
        ["bench", "run", str(tmp_path / "cases"), "-o", str(output), "--root", str(tmp_path)]
    )
    assert code == 0
    case_dir = output / "tiny-batch"
    for name in ("manifest.json", "summary.json", "report.html", "bench.json", "bench.html"):
        assert (case_dir / name).is_file()
    payload = json.loads((case_dir / "bench.json").read_text(encoding="utf-8"))
    assert payload["passed"] is True
    assert payload["container"] is None
    assert payload["engine"] == "legacy_batch"
    assert payload["unverified"] is False
    assert payload["manifest"] == "manifest.json"
    assert str(tmp_path) not in (case_dir / "bench.json").read_text(encoding="utf-8")
    page = (case_dir / "bench.html").read_text(encoding="utf-8")
    assert "does not show that the model is right" in page
    assert "report.html" in page
    assert adapter_for("marse").container == ""
    assert MarseAdapter.container == ""
    assert main(["bench", "report", str(output)]) == 0


def test_a_metric_outside_its_tolerance_fails_the_case(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    _batch(tmp_path / "batch.json")
    _case(
        tmp_path / "cases",
        "tiny.yml",
        "  - id: carbon\n    kind: balance_residual\n    quantity: C\n    tolerance: 1.0e-12\n",
    )
    output = tmp_path / "bench"
    code = main(
        ["bench", "run", str(tmp_path / "cases"), "-o", str(output), "--root", str(tmp_path)]
    )
    assert code == 1
    payload = json.loads((output / "tiny-batch" / "bench.json").read_text(encoding="utf-8"))
    assert payload["passed"] is False
    assert payload["metrics"][0]["passed"] is False
    assert main(["bench", "report", str(output)]) == 1
    assert "fail" in capsys.readouterr().out


def test_an_unknown_adapter_is_refused(tmp_path: Path):
    _batch(tmp_path / "batch.json")
    _case(tmp_path / "cases", "tiny.yml", "  - id: recorded\n    kind: manifest_present\n")
    path = tmp_path / "cases" / "tiny.yml"
    path.write_text(
        path.read_text(encoding="utf-8").replace("adapter: marse", "adapter: other"),
        encoding="utf-8",
    )
    with pytest.raises(CaseError, match="unknown adapter"):
        adapter_for("other")
    argv = [
        "bench",
        "run",
        str(tmp_path / "cases"),
        "-o",
        str(tmp_path / "out"),
        "--root",
        str(tmp_path),
    ]
    code = main(argv)
    assert code == 2


def test_yaml_that_does_not_parse_exits_2(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    case = tmp_path / "broken.yml"
    case.write_text("id: [unterminated\n", encoding="utf-8")
    with pytest.raises(CaseError, match=r"broken\.yml"):
        load_case(case, tmp_path)
    cases = tmp_path / "cases"
    cases.mkdir()
    (cases / "broken.yml").write_text("id: [unterminated\n", encoding="utf-8")
    code = main(["bench", "list", str(cases), "--root", str(tmp_path)])
    assert code == 2
    assert "invalid benchmark case" in capsys.readouterr().out


def test_a_directory_with_no_cases_is_not_a_pass(tmp_path: Path):
    empty = tmp_path / "empty"
    empty.mkdir()
    _batch(tmp_path / "batch.json")
    code = main(["bench", "run", str(empty), "-o", str(tmp_path / "out"), "--root", str(tmp_path)])
    assert code == 2
    assert not (tmp_path / "out").exists()


def test_an_unverified_engine_cannot_pass(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    experiment = "examples/experiments/two_species_ecosystem.json"
    cases = tmp_path / "cases"
    cases.mkdir()
    (cases / "eco.yml").write_text(
        "\n".join(
            [
                "id: eco",
                "title: Unverified engine",
                "adapter: marse",
                f"experiment: {experiment}",
                "command: ecosystem",
                "notes: This engine is unverified. The case must not pass.",
                "metrics:",
                "  - id: disclosed",
                "    kind: unverified_disclosed",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    root = Path(".").resolve()
    output = tmp_path / "out"
    code = main(["bench", "run", str(cases), "-o", str(output), "--root", str(root)])
    captured = capsys.readouterr()
    assert code == 1
    assert "unverified" in captured.err
    payload = json.loads((output / "eco" / "bench.json").read_text(encoding="utf-8"))
    assert payload["passed"] is False
    assert payload["reason"] == "unverified engine"
    assert payload["unverified"] is True
    assert payload["engine"] == "legacy_ecosystem_v1"
    assert payload["metrics"][0]["passed"] is True
    assert str(tmp_path) not in (output / "eco" / "bench.json").read_text(encoding="utf-8")


def test_bench_help_exits_zero():
    for argv in (
        ["bench", "--help"],
        ["bench", "list", "--help"],
        ["bench", "run", "--help"],
        ["bench", "report", "--help"],
    ):
        with pytest.raises(SystemExit) as caught:
            main(argv)
        assert caught.value.code == 0
