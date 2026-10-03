"""A day of plaque, and how deep oxygen reaches into it.

Runs the two scenes of examples/environments/oral that Stage S2 added:

- the oxygen profile: 400 um of plaque under saliva at the air for half an
  hour, then a rinse of 10% sucrose. Prints how far into the plaque oxygen
  reaches, to 1% of saturation, every 3 minutes;
- the day: from 07:00, meals that are chewed, sweets, gum, brushing at 07:45
  and 22:00, and the tongue and cheeks wearing the surface. Prints the
  plaque's thickness every hour, and the lowest pH under it in the hour
  before.

Then checks, as the example's own test:

- P6, set before Stage S2 was built: under saliva, the plaque is anoxic below
  200 to 250 um, and after sucrose the oxygen reaches less far, 120 to
  180 um (von Ohle et al. 2010: about 220 and 150 um);
- each brushing takes 42% of the plaque off, and the plaque grows back
  between them, to within a sixth of where it started by the next morning;
- P1: the box, and the box and the mouth together, conserve every element
  over the day, to 1e-12, counting what the air, the meals, the brush, the
  swallows and the glands brought or took.

The respiration and growth on saliva are calibrated to the first check
(confidence C, docs/parameters.md). The day takes about two minutes.

Run with: python examples/plaque_day.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from marse.core.reactive_transport import run
from marse.schemas.experiment import experiment_from_dict

SCENES = Path(__file__).resolve().parent / "environments" / "oral"
EVERY_MIN = 3.0


def oxygen_reach() -> tuple[np.ndarray, np.ndarray, float]:
    """Minutes, and how far oxygen reaches into the plaque then; and when the rinse is taken."""
    raw = json.loads((SCENES / "oxygen_profile.json").read_text("utf-8"))
    config = experiment_from_dict(raw)
    names = config.network.component_names
    saturation = config.domain.air.saturation_mol_per_m3["oxygen"]
    surface = config.domain.plaque.maximum_um
    heights = config.domain.grid.heights_um()
    inside = heights < surface
    minutes, reach = [], []

    def frames(_: int, time_h: float, fields: np.ndarray) -> None:
        # From the plaque's surface down, the first voxel below 1% of saturation.
        oxygen = fields[names.index("oxygen")][inside][::-1]
        depth = (surface - heights[inside])[::-1]
        anoxic = np.flatnonzero(oxygen < 0.01 * saturation)
        minutes.append(time_h * 60.0)
        reach.append(float(depth[anoxic[0]]) if anoxic.size else surface)

    run(config, frames=frames)
    return np.array(minutes), np.array(reach), config.domain.diet[0].start_h * 60.0


def main() -> None:
    minutes, reach, rinse = oxygen_reach()
    print("How far oxygen reaches into 400 um of plaque, to 1% of saturation")
    for minute in np.arange(0.0, minutes[-1] + 1e-9, EVERY_MIN):
        i = int(np.argmin(np.abs(minutes - minute)))
        note = "  <- the sucrose rinse" if abs(minutes[i] - rinse) < 0.5 * EVERY_MIN else ""
        print(f"{minutes[i]:>6.0f} min {reach[i]:>7.1f} um{note}")
    under_saliva = reach[minutes <= rinse][-1]
    after_sucrose = reach[minutes > rinse].min()

    result = run(experiment_from_dict(json.loads((SCENES / "plaque_day.json").read_text("utf-8"))))
    times = result.times_h
    thickness = result.plaque["thickness_um"]
    ph = result.ph["substratum_mean"]
    print("\nThe day: the plaque's thickness, and the lowest pH under it in the hour before")
    for hour in range(25):
        i = int(np.argmin(np.abs(times - hour)))
        lowest = ph[(times > hour - 1 - 1e-9) & (times <= hour + 1e-9)].min()
        print(f"{(7 + hour) % 24:02d}:00 {thickness[i]:>7.1f} um   pH {lowest:.2f}")
    outputs = result.manifest.outputs
    box = max(entry["largest_relative_residual"] for entry in outputs["balance"].values())
    whole = max(
        entry["largest_relative_residual"] for entry in outputs["mouth"]["balance"].values()
    )
    brushings = [event.start_h for event in result.config.domain.hygiene]
    taken = []
    for at in brushings:
        before = thickness[np.flatnonzero(times <= at + 1e-9)[-1]]
        after = thickness[np.flatnonzero(times > at + 1e-9)[0]]
        taken.append(1.0 - after / before)
    print(
        f"\n{outputs['mouth']['swallows']} swallows; the air gave "
        f"{outputs['air']['exchanged_mol_per_m2']['oxygen']:.3f} mol of oxygen per m2; "
        f"conserved to {box:.1e} in the box and {whole:.1e} with the mouth"
    )

    checks = {
        "under saliva, anoxic below 200 to 250 um": 200.0 <= under_saliva <= 250.0,
        "after sucrose, oxygen reaches less far: 120 to 180 um": (
            120.0 <= after_sucrose <= 180.0 and after_sucrose < under_saliva
        ),
        # Within the five minutes after each brushing, the plaque also grows and wears a little.
        "each brushing takes about 42% of the plaque off": all(
            abs(share - 0.42) < 0.02 for share in taken
        ),
        "the plaque grows back between brushings": thickness[
            np.flatnonzero(times < brushings[1])[-1]
        ]
        > thickness[np.flatnonzero(times > brushings[0])[0]],
        "by the next morning, to within a sixth of where it started": (
            abs(thickness[-1] / thickness[0] - 1.0) < 1 / 6
        ),
        "the box conserves every element to 1e-12": box < 1e-12,
        "and so do the box and the mouth together": whole < 1e-12,
    }
    print()
    for check, passed in checks.items():
        print(f"{'pass' if passed else 'FAIL'}  {check}")
    if not all(checks.values()):
        raise SystemExit("the plaque's day failed its checks")


if __name__ == "__main__":
    main()
