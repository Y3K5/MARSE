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
click-to-inspect local fields. The default Smooth fields mode uses browser
interpolation for continuous density gradients and can be switched to Cell grid
to inspect the exact numerical cells. Zoom enlarges the canvas for colony-scale
inspection without changing the simulation data, while opacity helps compare a
field with the dark spatial background.

Fields + particles adds deterministic population representatives sampled from
biomass density. Agents only displays explicit agent positions when a frame
contains them, and draws short movement trails between adjacent frames. The
current ecosystem engine exports continuous biomass fields by default, so these
particles are visual representatives rather than individually simulated cells.

A run's frames are streamed to disk as it goes, into `frames/` beside the
viewer: one NumPy `.npy` file per field, in single precision, plus an
`index.json` with each frame's time and step. Memory therefore stays flat
however long the run. By default about 200 evenly spaced frames are stored,
whatever the run's length; `--frame-interval-h` sets the spacing instead (for
example `--frame-interval-h 24` for daily frames). The viewer embeds up to 100
of the stored frames, always including the first and the last, with values
rounded to four significant digits of each field's largest value. The exact
result of a run is its final state, recorded in `manifest.json` and reproduced
by `marse replay`; the frames are for looking at how it got there.

Frames no longer limit how long a run can be. The step size still does: the
explicit transport scheme needs short steps, which is what the second stage of
the [roadmap](roadmap.md#order-of-work-correctness-first) replaces.
