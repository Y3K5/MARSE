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
