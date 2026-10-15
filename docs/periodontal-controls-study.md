# Periodontal community controls

**Unverified.** These controls use the same ecosystem v1 engine as
[periodontal-pathogen-study.md](periodontal-pathogen-study.md). The
`marse ecosystem` command warns on every run. Do not read a difference
between a control and the three-species case as a biological result until
the study runs on the version 2 engine.

`examples/periodontal_controls_analysis.py` is the next analysis step after
the exploratory three-species run. It is designed to answer whether a result
comes from community composition, the competition rule, or a limiting resource.

Run it locally:

```bash
python examples/periodontal_controls_analysis.py \
  --output /tmp/periodontal-controls
```

The workflow runs **168 compact scenarios**:

- seven communities: three single-species controls, three pairwise controls,
  and the mixed three-species community;
- reference competition and a no-competition counterfactual;
- reference, oxygen-limited, heme-limited, and peptide-limited resources;
- three deterministic replicate seeds;
- fixed reference temperature (37 C) and reference moisture.

It writes only `summary.csv` and `manifest.json`. The summary includes final
biomass, occupied cells, and cells occupied by at least two species. It does
not export browser trajectories, so the control sweep remains small enough to
repeat during development.

## Interpretation

The no-competition condition is a counterfactual, not a biological claim. The
resource perturbations halve one model field and are useful for sensitivity
ranking, but they are not calibrated representations of oxygen, heme, or
peptide depletion in a periodontal pocket. Pairwise and mixed results should
be compared against measured abundance or growth data before being described
as pathogen interactions.
