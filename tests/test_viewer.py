from marse.ecosystem import load_experiment, run, write_viewer


def test_viewer_contains_interactive_layers_and_controls(tmp_path):
    result = run(load_experiment("examples/experiments/two_species_ecosystem.json"))
    viewer = write_viewer(result, tmp_path / "viewer.html")
    html = viewer.read_text(encoding="utf-8")
    assert "phenotype state" in html
    assert "Keyboard: Space" in html
    assert "additive" in html
