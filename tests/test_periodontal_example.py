import json

from marse.ecosystem import load_experiment, run


def test_periodontal_pathogen_example_runs_and_keeps_three_species():
    result = run(load_experiment("examples/experiments/periodontal_pathogen_biofilm.json"))

    assert [species.name for species in result.config.species] == [
        "T. denticola",
        "T. forsythia",
        "P. gingivalis",
    ]
    assert len(result.frames) == result.config.steps + 1
    assert all(frame.biomass.shape == (3, 40, 60) for frame in result.frames)
    assert all(frame.biomass.min() >= 0 for frame in result.frames)


def test_periodontal_variant_workflow_writes_compact_controls(tmp_path, monkeypatch):
    from examples.periodontal_variant_analysis import main

    monkeypatch.setattr(
        "sys.argv",
        ["periodontal_variant_analysis", "--output", str(tmp_path)],
    )
    main()
    assert (tmp_path / "summary.csv").is_file()
    assert (tmp_path / "manifest.json").is_file()
    assert len((tmp_path / "summary.csv").read_text(encoding="utf-8").splitlines()) == 73
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["provenance_status"] == "exploratory-placeholder"
    assert manifest["calibration_status"] == ["not-calibrated"]
    assert "community-interactions" in manifest["parameter_group_ids"]


def test_periodontal_controls_cover_pairwise_and_resource_controls(tmp_path, monkeypatch):
    from examples.periodontal_controls_analysis import main

    monkeypatch.setattr(
        "sys.argv",
        ["periodontal_controls_analysis", "--output", str(tmp_path)],
    )
    main()
    rows = (tmp_path / "summary.csv").read_text(encoding="utf-8").splitlines()
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert len(rows) == 169
    assert manifest["scenario_count"] == 168
    assert len(manifest["communities"]) == 7
