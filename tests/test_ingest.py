"""Table importers check synthetic rows. They do not claim a measurement."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from marse.cli import main
from marse.ingest import Feed, apply_fragment, importer_for, register, registered_kinds, unregister
from marse.ingest.errors import IngestError

EXAMPLES = Path("examples/ingest")


def test_the_builtin_kinds_are_registered():
    assert registered_kinds() == ("abundance", "environment", "metabolomics", "transcriptomics")


def test_a_caller_can_register_an_importer():
    class Probe:
        kind = "probe"

        def load(self, path: Path) -> Feed:
            raise IngestError(path.name)

    register(Probe())
    try:
        assert importer_for("probe").kind == "probe"
        with pytest.raises(IngestError, match="already registered"):
            register(Probe())
    finally:
        unregister("probe")
    assert "probe" not in registered_kinds()


def test_environment_feed_copies_only_temperature_and_ph(tmp_path: Path):
    code = main(
        [
            "ingest",
            "environment",
            str(EXAMPLES / "synthetic_environment.csv"),
            "-o",
            str(tmp_path / "feed.json"),
        ]
    )
    assert code == 0
    feed = json.loads((tmp_path / "feed.json").read_text(encoding="utf-8"))
    assert feed["format"] == "marse.ingest/1"
    assert feed["qc"]["passed"] is True
    assert feed["experiment_fragment"] == {"environment": {"temperature_c": 37.0, "ph": 7.0}}
    assert str(tmp_path) not in (tmp_path / "feed.json").read_text(encoding="utf-8")
    experiment = {"experiment_id": "synthetic", "environment": {"temperature_c": 37.0, "ph": 7.0}}
    assert apply_fragment(experiment, feed["experiment_fragment"])["experiment_id"] == "synthetic"
    with pytest.raises(IngestError, match="already has environment"):
        apply_fragment(
            {"environment": {"temperature_c": 20.0, "ph": 7.0}}, feed["experiment_fragment"]
        )


def test_abundance_and_transcripts_do_not_fill_an_experiment(tmp_path: Path):
    for kind, name in (
        ("abundance", "synthetic_abundance.csv"),
        ("transcriptomics", "synthetic_transcripts.csv"),
    ):
        output = tmp_path / f"{kind}.json"
        assert main(["ingest", kind, str(EXAMPLES / name), "-o", str(output)]) == 0
        feed = json.loads(output.read_text(encoding="utf-8"))
        assert feed["experiment_fragment"] is None
        assert feed["qc"]["passed"] is True


def test_one_metabolite_in_mol_per_m3_fills_one_substrate(tmp_path: Path):
    output = tmp_path / "feed.json"
    argv = ["ingest", "metabolomics", str(EXAMPLES / "synthetic_metabolite.csv"), "-o", str(output)]
    assert main(argv) == 0
    feed = json.loads(output.read_text(encoding="utf-8"))
    assert feed["experiment_fragment"] == {"substrate": {"name": "glucose", "initial_mm": 10.0}}
    merged = apply_fragment({"experiment_id": "synthetic"}, feed["experiment_fragment"])
    assert merged["substrate"]["initial_mm"] == 10.0


def test_two_metabolites_stay_in_the_table(tmp_path: Path):
    table = tmp_path / "two.csv"
    table.write_text(
        "sample,feature,value,unit\nsample-1,glucose,10,mmol_per_l\n"
        "sample-1,acetate,1,mmol_per_l\n",
        encoding="utf-8",
    )
    feed = importer_for("metabolomics").load(table)
    assert feed.experiment_fragment is None


def test_a_bad_abundance_table_is_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    table = tmp_path / "bad.csv"
    _write(
        table,
        ["sample", "taxon", "abundance", "unit"],
        [["sample-1", "taxon-a", "-1", "relative"]],
    )
    code = main(["ingest", "abundance", str(table), "-o", str(tmp_path / "out.json")])
    assert code == 2
    assert "invalid table" in capsys.readouterr().out
    assert not (tmp_path / "out.json").exists()
    _write(
        table,
        ["sample", "taxon", "abundance", "unit"],
        [["sample-1", "taxon-a", "0.2", "relative"], ["sample-1", "taxon-b", "0.2", "relative"]],
    )
    with pytest.raises(IngestError, match="sum to"):
        importer_for("abundance").load(table)


def test_unknown_environment_fields_and_units_are_refused(tmp_path: Path):
    table = tmp_path / "env.csv"
    _write(table, ["name", "value", "unit"], [["salinity", "35", "psu"]])
    with pytest.raises(IngestError, match="not an environment field"):
        importer_for("environment").load(table)
    _write(table, ["name", "value", "unit"], [["temperature_c", "37", "K"], ["ph", "7", "pH"]])
    with pytest.raises(IngestError, match="unit"):
        importer_for("environment").load(table)


def test_a_tab_separated_table_loads(tmp_path: Path):
    table = tmp_path / "wide.tsv"
    table.write_text(
        "sample\tfeature\tvalue\tunit\nsample-1\tfeature-a\t3\tcount\n",
        encoding="utf-8",
    )
    feed = importer_for("transcriptomics").load(table)
    assert feed.rows[0]["feature"] == "feature-a"
    with pytest.raises(IngestError, match="unknown importer"):
        importer_for("missing")


def _write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
