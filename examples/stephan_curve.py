"""The Stephan curve, and what keeps plaque acid for longer.

Runs the three oral scenes of examples/environments/oral: a rinse of 10% sucrose
held for a minute (the Stephan curve), 100 mL of a drink of 10% sucrose sipped
over 20 minutes, and the same rinse followed by food left on the teeth. Prints
the pH at the substratum, under 150 um of plaque, every 3 minutes, then for
each scene the lowest pH and when, the minutes and the area below pH 5.5, and
when the pH is back above 6.

Then checks, as the example's own test:

- the rinse meets the criteria set for a Stephan curve before Stage S1 was
  built (docs/validation.md, "The Stephan curve"): the pH falls at least one
  unit, to a minimum of 4.5 to 5.5 within 5 to 20 minutes, is back above 6
  within an hour, and the plaque holds 10 to 60 mM more lactate at 7 minutes;
- the sipped drink and the food left on the teeth each keep the plaque below
  pH 5.5 for longer than the rinse alone, as a high-sugar eater's habits do.

The acid production, buffer, plaque thickness and film speed are calibrated
to the first criteria (confidence C, docs/parameters.md).

Run with: python examples/stephan_curve.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from marse.core.reactive_transport import ReactiveTransportResult, run
from marse.oral import CRITICAL_PH, StephanCurve
from marse.schemas.experiment import experiment_from_dict

SCENES = Path(__file__).resolve().parent / "environments" / "oral"
NAMES = {"stephan_rinse.json": "rinse", "sipping.json": "sipping", "pocket.json": "pocket"}
EVERY_MIN = 3.0


def lactate_rise_at(minutes: float, name: str) -> tuple[ReactiveTransportResult, float]:
    """Run a scene, and how much the plaque's mean lactate has risen at ``minutes``."""
    config = experiment_from_dict(json.loads((SCENES / name).read_text("utf-8")))
    names = config.network.component_names
    start = config.domain.initial_state(names, config.initial_mol_per_m3, config.seed)
    plaque = start[names.index("bacteria")] > 0
    times, lactate = [], []

    def frames(_: int, time_h: float, fields: np.ndarray) -> None:
        times.append(time_h * 60.0)
        lactate.append(float(fields[names.index("lactate")][plaque].mean()))

    result = run(config, frames=frames)
    return result, float(np.interp(minutes, times, lactate)) - lactate[0]


def main() -> None:
    results, rises = {}, {}
    for file, label in NAMES.items():
        results[label], rises[label] = lactate_rise_at(7.0, file)
    times = results["rinse"].times_h
    print("pH at the substratum")
    print(f"{'min':>6}" + "".join(f"{label:>10}" for label in results))
    for minute in np.arange(0.0, times[-1] * 60.0 + 1e-9, EVERY_MIN):
        row = [np.interp(minute, times * 60.0, r.ph["substratum_mean"]) for r in results.values()]
        print(f"{minute:>6.0f}" + "".join(f"{v:>10.2f}" for v in row))

    curves = {
        label: StephanCurve.of(result.times_h, result.ph["substratum_mean"])
        for label, result in results.items()
    }
    print(f"\n{'':<10}{'lowest':>8}{'at min':>8}{'min <5.5':>10}{'area <5.5':>11}{'>6 at min':>11}")
    for label, curve in curves.items():
        back = "never" if curve.back_above_6_min is None else f"{curve.back_above_6_min:.1f}"
        print(
            f"{label:<10}{curve.minimum:>8.2f}{curve.minimum_at_min:>8.1f}"
            f"{curve.minutes_below_critical:>10.1f}{curve.area_below_critical:>11.1f}{back:>11}"
        )
    print(
        f"\n(area below pH {CRITICAL_PH} in pH x minutes; lactate in the plaque at 7 min: "
        f"+{rises['rinse']:.1f} mM after the rinse)"
    )

    rinse = curves["rinse"]
    checks = {
        "the rinse's pH falls at least one unit": rinse.start - rinse.minimum >= 1.0,
        "to a minimum of 4.5 to 5.5": 4.5 <= rinse.minimum <= 5.5,
        "within 5 to 20 minutes": 5.0 <= rinse.minimum_at_min <= 20.0,
        "back above 6 within an hour": (
            rinse.back_above_6_min is not None and rinse.back_above_6_min <= 60.0
        ),
        "with 10 to 60 mM more lactate at 7 minutes": 10.0 <= rises["rinse"] <= 60.0,
        "sipping keeps the plaque below 5.5 for longer": (
            curves["sipping"].minutes_below_critical > rinse.minutes_below_critical
        ),
        "and so does food left on the teeth": (
            curves["pocket"].minutes_below_critical > rinse.minutes_below_critical
        ),
    }
    print()
    for check, passed in checks.items():
        print(f"{'pass' if passed else 'FAIL'}  {check}")
    if not all(checks.values()):
        raise SystemExit("the Stephan curve's checks failed")


if __name__ == "__main__":
    main()
