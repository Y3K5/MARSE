# Local Mac visualizer

MARSE's ecosystem viewer is dependency-free and runs well on an Apple Silicon
Mac. Use the repository environment to generate a simulation and open the
result in Safari or Chrome:

```bash
python -m marse.cli ecosystem examples/experiments/two_species_ecosystem.json \
  --output runs/two-species
open runs/two-species/viewer.html
```

The viewer renders species, nutrients, conditions, additives, phenotype states,
effective growth rates, and limiting factors. It supports playback speed,
single-frame stepping, keyboard controls, responsive canvas sizing, and
click-to-inspect local fields.

For a year-long experiment, do not export every integration step as a browser
frame. Use a larger numerical timestep where the stability and biological
timescales permit it, then export daily or weekly snapshots. Checkpointed
streaming execution is the next infrastructure step for very long runs.
