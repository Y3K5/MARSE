"""The scenes of Stage S2: oxygen in thick plaque, and a day of plaque.

The criteria set before Stage S2 was built that these check, on the scenes of
examples/environments/oral (docs/validation.md, "A day of plaque"):

- P6: under saliva, 400 um of plaque is anoxic below about 220 um, and after
  a sucrose rinse oxygen reaches less far, about 150 um (von Ohle et al. 2010);
- P1: over a day of meals, sweets, gum, two brushings, wear and the air, the
  box, and the box and the mouth together, conserve every element to 1e-12;
- P8: that day, in a column of 100 voxels, runs in minutes.

The day's brushings each take 42% of the plaque, and the plaque grows back
between them, on sugar and on saliva.
"""

import json
import time
from pathlib import Path

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.chemistry import ChargeBalance
from marse.cli import main
from marse.schemas import experiment_from_dict

SCENES = Path(__file__).resolve().parents[1] / "examples" / "environments" / "oral"
FILES = ("oxygen_profile.json", "plaque_day.json")


def load(name):
    return json.loads((SCENES / name).read_text("utf-8"))


@pytest.mark.parametrize("name", FILES)
def test_each_scene_breathes_spreads_and_starts_neutral(name, capsys):
    config = experiment_from_dict(load(name))
    domain = config.domain
    names = config.network.component_names
    balance = ChargeBalance.of(config.network)
    saliva = np.array([domain.mouth.saliva_mol_per_m3.get(n, 0.0) for n in names])
    assert float(balance.ph(saliva)) == pytest.approx(6.8, abs=1e-3)  # Bardow et al. 2000
    assert domain.air.saturation_mol_per_m3 == {"oxygen": 0.21}
    assert config.network.component("bacteria").density_mol_per_m3 == 800
    assert domain.plaque is domain.spreading
    assert domain.spreading.mechanism == "packed"
    assert main(["check", str(SCENES / name)]) == 0
    assert "the film rides on its surface" in capsys.readouterr().out


def test_the_scenes_share_one_network_built_on_the_stephan_curves():
    profile, day = load("oxygen_profile.json"), load("plaque_day.json")
    assert profile["components"] == day["components"]
    assert profile["processes"] == day["processes"]
    rinse = load("stephan_rinse.json")
    # The day's bacteria grow and pack, so they also state an oxygen role and a density.
    added = {"oxygen_role", "density_mol_per_m3"}
    components = [{k: v for k, v in c.items() if k not in added} for c in day["components"]]
    assert all(c in components for c in rinse["components"])
    assert all(p in day["processes"] for p in rinse["processes"])


def _reach(raw):
    """Minutes, and how far into the plaque oxygen stays above 1% of saturation."""
    config = experiment_from_dict(raw)
    names = config.network.component_names
    surface = config.domain.plaque.maximum_um
    heights = config.domain.grid.heights_um()
    inside = heights < surface
    minutes, reach = [], []

    def frames(_, time_h, fields):
        oxygen = fields[names.index("oxygen")][inside][::-1]
        depth = (surface - heights[inside])[::-1]
        anoxic = np.flatnonzero(oxygen < 0.01 * 0.21)
        minutes.append(time_h * 60.0)
        reach.append(float(depth[anoxic[0]]) if anoxic.size else surface)

    result = reactive_transport.run(config, frames=frames)
    return np.array(minutes), np.array(reach), result


@pytest.mark.slow
def test_p6_oxygen_reaches_as_deep_into_plaque_as_von_ohle_measured():
    minutes, reach, result = _reach(load("oxygen_profile.json"))
    under_saliva = reach[minutes <= 30.0 + 1e-9][-1]
    after_sucrose = reach[minutes > 30.0 + 1e-9].min()
    assert under_saliva == pytest.approx(216.2, abs=1.0)  # about 220 um
    assert after_sucrose == pytest.approx(151.2, abs=1.0)  # about 150 um
    assert reach[-1] > 200.0  # back as the sugar is cleared
    outputs = result.manifest.outputs
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-12, quantity


@pytest.mark.slow
def test_p1_and_p8_a_day_of_plaque_balances_and_runs_in_minutes():
    started = time.perf_counter()
    result = reactive_transport.run(experiment_from_dict(load("plaque_day.json")))
    elapsed = time.perf_counter() - started
    assert elapsed < 600.0  # P8: about two minutes
    outputs = result.manifest.outputs
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-12, quantity  # P1
    mouth = outputs["mouth"]
    assert mouth["intakes"] == 6
    assert mouth["cleanings"] == 2
    assert mouth["swallows"] > 2000
    times, thickness = result.times_h, result.plaque["thickness_um"]
    for at in (0.75, 15.0):  # each brushing takes 42%, less what grows in the five minutes after
        before = thickness[np.flatnonzero(times <= at + 1e-9)[-1]]
        after = thickness[np.flatnonzero(times > at + 1e-9)[0]]
        assert 1.0 - after / before == pytest.approx(0.42, abs=0.02)
    morning = thickness[np.flatnonzero(times > 1.0)[0]]
    evening = thickness[np.flatnonzero(times < 15.0)[-1]]
    assert evening > 2 * morning  # grown back, on sugar and on saliva
    assert thickness[-1] == pytest.approx(thickness[0], rel=1 / 6)
    assert outputs["air"]["exchanged_mol_per_m2"]["oxygen"] > 0.5
