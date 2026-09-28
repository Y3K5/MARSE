"""The oral scenes: a rinse, a sipped drink and food left on the teeth, and the Stephan curve.

Two of the criteria set before Stage S1 was built are checked here, as slow
tests, on the rinse of examples/environments/oral/stephan_rinse.json:

- G5: the pH at the substratum falls at least one unit, to a minimum of 4.5 to
  5.5 within 5 to 20 minutes, is back above 6 within an hour, and the plaque
  holds 10 to 60 mM more lactate at 7 minutes;
- G7: a two-hour curve in 100 voxels runs in under 30 seconds.

So are the directions the literature reports, each against the same rinse:
sipping and food left on the teeth keep the plaque acid for longer; low
salivary flow (Lingstrom and Birkhed 1993) and a slower salivary film
(Macpherson and Dawes 1991) deepen the fall and slow the return; thicker plaque
keeps it acid for longer (Dawes and Dibdin 1986); and more fixed buffer makes
the fall shallower and the return slower (Dibdin 1990).
"""

import copy
import json
import math
import time
from pathlib import Path

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.chemistry import ChargeBalance
from marse.cli import main
from marse.oral import CRITICAL_PH, StephanCurve, area_below, back_above, minutes_below
from marse.schemas import experiment_from_dict

SCENES = Path(__file__).resolve().parents[1] / "examples" / "environments" / "oral"
FILES = ("stephan_rinse.json", "sipping.json", "pocket.json")


def load(name):
    return json.loads((SCENES / name).read_text("utf-8"))


def run(raw):
    """A scene's run, and the plaque's mean lactate at every record."""
    config = experiment_from_dict(raw)
    names = config.network.component_names
    start = config.domain.initial_state(names, config.initial_mol_per_m3, config.seed)
    plaque = start[names.index("bacteria")] > 0
    lactate = []

    def frames(_, __, fields):
        lactate.append(float(fields[names.index("lactate")][plaque].mean()))

    result = reactive_transport.run(config, frames=frames)
    return result, np.array(lactate)


def curve(result):
    return StephanCurve.of(result.times_h, result.ph["substratum_mean"])


def later(a, b):
    """Whether return time a is later than b; a curve that never returns is latest."""
    return (math.inf if a is None else a) > (math.inf if b is None else b)


# --- the scenes ------------------------------------------------------------------------


@pytest.mark.parametrize("name", FILES)
def test_each_scene_starts_neutral_at_the_ph_saliva_was_measured_at(name, capsys):
    raw = load(name)
    config = experiment_from_dict(raw)
    balance = ChargeBalance.of(config.network)
    names = config.network.component_names
    mouth = config.domain.mouth

    def ph(composition):
        return float(balance.ph(np.array([composition.get(n, 0.0) for n in names])))

    assert ph(mouth.saliva_mol_per_m3) == pytest.approx(6.8, abs=1e-3)  # Bardow et al. 2000
    assert ph(mouth.stimulated_saliva_mol_per_m3) == pytest.approx(7.2, abs=1e-3)
    box = balance.ph(config.domain.initial_state(names, config.initial_mol_per_m3, config.seed))
    np.testing.assert_allclose(box, 6.8, atol=1e-3)
    assert main(["check", str(SCENES / name)]) == 0
    assert "Bardow et al. (2000)" in " ".join(capsys.readouterr().out.split())


def test_a_short_rinse_conserves_every_element():
    raw = load("stephan_rinse.json")
    raw["duration_h"] = 0.05
    result, _ = run(raw)
    outputs = result.manifest.outputs
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-12, quantity


# --- the measures of a curve ------------------------------------------------------------------


def test_the_measures_of_a_curve_are_exact_for_the_lines_between_records():
    times_h = np.array([0.0, 10.0, 20.0, 40.0]) / 60.0
    ph = np.array([7.0, 5.0, 5.0, 7.0])
    # Below 5.5 from 7.5 to 25 minutes; the area below it is two triangles and a rectangle.
    assert minutes_below(times_h, ph, 5.5) == pytest.approx(17.5)
    assert area_below(times_h, ph, 5.5) == pytest.approx(0.5 * 0.5 * 2.5 + 0.5 * 10 + 0.5 * 0.5 * 5)
    assert back_above(times_h, ph, 6.0) == pytest.approx(30.0)
    assert back_above(times_h, ph, 7.5) is None
    summary = StephanCurve.of(times_h, ph)
    assert (summary.start, summary.minimum, summary.minimum_at_min) == (7.0, 5.0, 10.0)
    assert summary.minutes_below_critical == minutes_below(times_h, ph, CRITICAL_PH)
    with pytest.raises(ValueError, match="must increase"):
        minutes_below(times_h[::-1], ph, 5.5)


# --- G5, G7 and the directions (slow) --------------------------------------------------------


@pytest.fixture(scope="module")
def rinse():
    return run(load("stephan_rinse.json"))


@pytest.mark.slow
def test_g5_the_rinse_gives_a_stephan_curve(rinse):
    result, lactate = rinse
    shape = curve(result)
    rise = float(np.interp(7.0, result.times_h * 60.0, lactate)) - lactate[0]
    assert shape.start - shape.minimum >= 1.0
    assert 4.5 <= shape.minimum <= 5.5
    assert 5.0 <= shape.minimum_at_min <= 20.0
    assert shape.back_above_6_min is not None
    assert shape.back_above_6_min <= 60.0
    assert 10.0 <= rise <= 60.0
    # The film above the plaque is never more acid than the plaque itself.
    assert result.ph["box_min"].min() == pytest.approx(result.ph["substratum_min"].min(), abs=0.01)


@pytest.mark.slow
def test_g7_a_two_hour_curve_in_100_voxels_runs_in_under_30_seconds():
    raw = load("stephan_rinse.json")
    raw["duration_h"] = 2.0
    assert raw["domain"]["voxels"] == [100]
    started = time.perf_counter()
    reactive_transport.run(experiment_from_dict(raw))
    assert time.perf_counter() - started < 30.0


def _vary(change):
    raw = load("stephan_rinse.json")
    change(raw)
    return curve(run(raw)[0])


@pytest.mark.slow
@pytest.mark.parametrize("name", ["sipping.json", "pocket.json"])
def test_a_continuous_presence_of_sugar_keeps_the_plaque_acid_for_longer(name, rinse):
    base = curve(rinse[0])
    shape = curve(run(load(name))[0])
    assert shape.minutes_below_critical > 1.4 * base.minutes_below_critical
    assert later(shape.back_above_6_min, base.back_above_6_min)


@pytest.mark.slow
def test_low_salivary_flow_deepens_the_fall_and_slows_the_return(rinse):
    def low(raw):
        raw["domain"]["mouth"].update(
            unstimulated_flow_ml_per_min=0.1, stimulated_flow_ml_per_min=1.0
        )

    base, shape = curve(rinse[0]), _vary(low)
    assert shape.minimum < base.minimum
    assert later(shape.back_above_6_min, base.back_above_6_min)


@pytest.mark.slow
def test_a_slower_film_deepens_the_fall_and_slows_the_return(rinse):
    base = curve(rinse[0])
    shape = _vary(lambda raw: raw["domain"]["film"].update(velocity_mm_per_min=2))
    assert shape.minimum < base.minimum - 0.2
    assert later(shape.back_above_6_min, base.back_above_6_min)


@pytest.mark.slow
def test_thicker_plaque_keeps_it_acid_for_longer(rinse):
    def thicker(raw):
        domain = raw["domain"]
        domain["voxels"] = [160]  # 300 um of plaque under the same 100 um of film
        for colony in domain["colonies"]:
            colony["radius_um"] = 300

    base, shape = curve(rinse[0]), _vary(thicker)
    assert shape.minutes_below_critical > 1.5 * base.minutes_below_critical
    assert later(shape.back_above_6_min, base.back_above_6_min)


@pytest.mark.slow
def test_more_fixed_buffer_makes_the_fall_shallower_and_the_return_slower(rinse):
    def buffered(raw):
        colonies = {c["component"]: c for c in raw["domain"]["colonies"]}
        scale = 1.5
        for name in ("carboxyl_groups", "bound_potassium"):
            colonies[name]["concentration_mol_per_m3"] *= scale

    base, shape = curve(rinse[0]), _vary(buffered)
    assert shape.minimum > base.minimum
    assert later(shape.back_above_6_min, base.back_above_6_min)


def test_the_scenes_differ_only_in_what_is_eaten():
    rinse, sipping, pocket = (copy.deepcopy(load(name)) for name in FILES)
    for raw in (rinse, sipping, pocket):
        raw.pop("description")
        raw.pop("experiment_id")
        raw["domain"].pop("diet")
    assert sipping == rinse
    pocket["components"] = [c for c in pocket["components"] if c["name"] != "food_sugar"]
    pocket["processes"] = [p for p in pocket["processes"] if p["name"] != "dissolving"]
    assert pocket == rinse
