"""Oxygen from the air: through the top face of a column, and into the film under the mouth.

The criterion set before Stage S2 was built that these check (docs/validation.md,
"Oxygen from the air"):

- P5: oxygen entering through the surface, held at saturation by the air,
  fills a slab as the series solution says, and under zero-order uptake
  reaches sqrt(2 D C_s / k0) into it, both to 2% after the grid is refined.

And conservation: both ledgers count what the air gave and took, so they
close to rounding, under the mouth as in a column.
"""

import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

import marse.core.reactive_transport as reactive_transport
from marse.cli import main
from marse.core.config import ConfigError
from marse.schemas import experiment_from_dict

ROOT = Path(__file__).resolve().parents[1]
RINSE = ROOT / "examples" / "environments" / "oral" / "stephan_rinse.json"
D = 2.0e-9  # m2/s
D_UM2_PER_H = D * 1e12 * 3600.0
SATURATION = 0.2  # mol/m3


def _column(height_um, voxel_um, *, uptake_per_h=0.0, hours=0.02, records_h=None):
    """A column at the air, with oxygen taken up at a constant rate, zero order, where it is."""
    voxels = round(height_um / voxel_um)
    colonies = []
    if uptake_per_h:
        colonies = [
            {
                "component": "biomass",
                "center_um": [],
                "radius_um": height_um,
                "concentration_mol_per_m3": 1.0,
            }
        ]
    return {
        "schema_version": 2,
        "experiment_id": "column-at-the-air",
        "components": [
            {"name": "oxygen", "phase": "dissolved", "formula": "O2"},
            {"name": "donor", "phase": "dissolved", "formula": "CH2O"},
            {"name": "carbon_dioxide", "phase": "dissolved", "formula": "CO2"},
            {"name": "biomass", "phase": "particulate", "formula": "CH1.8O0.5N0.2"},
        ],
        "processes": [
            {
                "name": "respiration",
                "kind": "reaction",
                "stoichiometry_mol_per_mol": {"donor": -1, "oxygen": -1, "carbon_dioxide": 1},
                "rate": {
                    "maximum_per_h": uptake_per_h,
                    "proportional_to": "biomass",
                    "factors": [
                        {"component": "oxygen", "form": "monod", "half_saturation_mol_per_m3": 1e-5}
                    ],
                    "assumed_in_excess": ["donor"],
                },
            }
        ],
        "initial_mol_per_m3": {"donor": 1000.0},
        "duration_h": hours,
        "timestep_h": records_h or hours,
        "relative_tolerance": 1e-5,
        "absolute_tolerance_mol_per_m3": 1e-9,
        "domain": {
            "voxels": [voxels],
            "voxel_um": voxel_um,
            "diffusivity_m2_per_s": {"oxygen": D, "donor": 1e-9, "carbon_dioxide": 1e-9},
            "colonies": colonies,
            "air": {"saturation_mol_per_m3": {"oxygen": SATURATION}},
        },
    }


def _filling(height_um, depth_from_bottom_um, t_h, terms=200):
    """A slab closed below and held at saturation above, filling from nothing (Crank 1975)."""
    total = 0.0
    for n in range(terms):
        m = 2 * n + 1
        total += (
            4.0
            / (m * math.pi)
            * (-1) ** n
            * np.cos(m * math.pi * depth_from_bottom_um / (2 * height_um))
            * math.exp(-(m**2) * math.pi**2 * D_UM2_PER_H * t_h / (4 * height_um**2))
        )
    return SATURATION * (1.0 - total)


def _worst_filling_error(voxel_um, height_um=200.0):
    """The largest error, as a share of saturation, at four times as the slab fills."""
    seen = {}

    def frames(index, t_h, fields):
        seen[round(t_h, 6)] = fields[0].copy()

    raw = _column(height_um, voxel_um, hours=0.02, records_h=0.0005)
    result = reactive_transport.run(experiment_from_dict(raw), frames=frames)
    z = result.config.domain.grid.heights_um()
    return max(
        float(np.max(np.abs(seen[t] - _filling(height_um, z, t)))) / SATURATION
        for t in (0.0005, 0.002, 0.005, 0.02)
    )


def test_p5_oxygen_fills_a_slab_from_the_air_as_the_series_says():
    coarse, fine = _worst_filling_error(5.0), _worst_filling_error(2.5)
    assert fine < 1e-4  # 5.5e-5 of saturation: far inside the criterion's 2%
    assert coarse / fine == pytest.approx(4.0, rel=0.1)  # second order in the voxel


def _penetration(voxel_um, reach_um=150.0, height_um=300.0):
    uptake = 2 * D_UM2_PER_H * SATURATION / reach_um**2  # mol/m3 per h, at 1 mol/m3 of biomass
    # Steady within seconds: diffusion across 150 um takes reach^2 / D, 11 s.
    raw = _column(height_um, voxel_um, uptake_per_h=uptake, hours=0.05)
    raw["absolute_tolerance_mol_per_m3"] = 1e-8
    result = reactive_transport.run(experiment_from_dict(raw))
    depth = height_um - result.config.domain.grid.heights_um()
    exact = np.where(depth < reach_um, SATURATION * (1 - depth / reach_um) ** 2, 0.0)
    oxygen = result.final_state[0]
    # Where it falls below 1% of saturation, between voxel centres, from the surface down.
    below = int(np.argmax(oxygen[::-1] < 0.01 * SATURATION))
    c, d = oxygen[::-1], depth[::-1]
    found = d[below - 1] + (d[below] - d[below - 1]) * (c[below - 1] - 0.01 * SATURATION) / (
        c[below - 1] - c[below]
    )
    return float(np.max(np.abs(oxygen - exact))) / SATURATION, found, result


def test_p5_zero_order_uptake_lets_oxygen_reach_as_far_as_theory_says():
    error, found, result = _penetration(5.0)
    reach = 150.0 * (1 - math.sqrt(0.01))  # where (1 - x / delta)^2 is 1%: 135 um
    assert found == pytest.approx(reach, rel=0.005)  # 135.3 um
    # Uptake is Monod with K = 1e-5, not quite zero order where oxygen runs out.
    assert error < 1e-3
    outputs = result.manifest.outputs
    for quantity, entry in outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-12, quantity
    # At steady state the air gives what the column takes up: all of it, through the surface.
    taken = outputs["totals_mol_per_m2"]["carbon_dioxide"]
    given = outputs["air"]["exchanged_mol_per_m2"]["oxygen"]
    assert given == pytest.approx(taken + outputs["totals_mol_per_m2"]["oxygen"], rel=1e-12)


def test_what_the_air_gives_is_booked_and_the_column_conserves_it():
    result = reactive_transport.run(experiment_from_dict(_column(200.0, 5.0, hours=0.01)))
    outputs = result.manifest.outputs
    held = outputs["totals_mol_per_m2"]["oxygen"]
    assert held > 0
    assert outputs["air"]["exchanged_mol_per_m2"]["oxygen"] == pytest.approx(held, rel=1e-13)
    assert outputs["imported_mol_per_m2"]["oxygen"] == pytest.approx(held, rel=1e-13)
    assert outputs["air"]["saturation_mol_per_m3"] == {"oxygen": SATURATION}
    # Nothing else crosses the top face: the donor stays where it was.
    assert outputs["imported_mol_per_m2"]["donor"] == 0.0
    for quantity, entry in outputs["balance"].items():
        assert entry["largest_relative_residual"] < 1e-13, quantity


# --- under the mouth -------------------------------------------------------------------------


def _breathing_rinse(duration_h=0.25, uptake=0.004):
    """S1's rinse with oxygen: the plaque respires sugar and the film's surface is at the air."""
    raw = json.loads(RINSE.read_text("utf-8"))
    raw["duration_h"] = duration_h
    raw["components"].append({"name": "oxygen", "phase": "dissolved", "formula": "O2"})
    raw["processes"].append(
        {
            "name": "respiration",
            "kind": "reaction",
            "stoichiometry_mol_per_mol": {"sugar": -1, "oxygen": -6, "carbonate": 6},
            "rate": {
                "maximum_per_h": uptake,
                "proportional_to": "bacteria",
                "factors": [
                    {"component": "sugar", "form": "monod", "half_saturation_mol_per_m3": 1},
                    {"component": "oxygen", "form": "monod", "half_saturation_mol_per_m3": 0.005},
                ],
            },
        }
    )
    raw["domain"]["diffusivity_m2_per_s"]["oxygen"] = 1.13e-9
    raw["domain"]["air"] = {"saturation_mol_per_m3": {"oxygen": 0.21}}
    return raw


def test_under_the_mouth_both_ledgers_count_what_the_air_exchanged():
    result = reactive_transport.run(experiment_from_dict(_breathing_rinse()))
    outputs = result.manifest.outputs
    for balance in (outputs["balance"], outputs["mouth"]["balance"]):
        for quantity, entry in balance.items():
            assert entry["largest_relative_residual"] < 1e-13, quantity
    # The air holds the mouth's oxygen at saturation, through swallows and the rinse.
    np.testing.assert_allclose(result.mouth["oxygen_mol_per_m3"], 0.21, rtol=1e-12)
    # It gave what the plaque burned, less what the box and the mouth gained.
    assert outputs["air"]["exchanged_mol_per_m2"]["oxygen"] > 0
    # Oxygen comes from the film's surface, at saturation, and falls into the plaque.
    oxygen = result.final_state[result.component_names.index("oxygen")]
    assert oxygen[-1] == pytest.approx(0.21, rel=0.01)
    assert np.all(np.diff(oxygen) > 0)
    assert oxygen[0] < oxygen[-1]


def test_a_rinse_that_holds_oxygen_of_its_own_is_brought_back_to_saturation():
    raw = _breathing_rinse(duration_h=0.05, uptake=0.0)
    raw["domain"]["diet"][0]["composition_mol_per_m3"]["oxygen"] = 0.05
    result = reactive_transport.run(experiment_from_dict(raw))
    np.testing.assert_allclose(result.mouth["oxygen_mol_per_m3"], 0.21, rtol=1e-12)
    for quantity, entry in result.manifest.outputs["mouth"]["balance"].items():
        assert entry["largest_relative_residual"] < 1e-13, quantity


# --- the schema --------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (lambda a: a.update(saturation_mol_per_m3={}), "at least one gas"),
        (lambda a: a.update(saturation_mol_per_m3={"nitrogen": 0.5}), "nitrogen"),
        (lambda a: a.update(saturation_mol_per_m3={"biomass": 0.5}), "particulate"),
        (lambda a: a.update(saturation_mol_per_m3={"oxygen": -0.1}), "must not be negative"),
        (lambda a: a.update(pressure_atm=1), "pressure_atm"),
    ],
)
def test_impossible_air_is_refused(change, match):
    raw = _column(100.0, 5.0)
    change(raw["domain"]["air"])
    with pytest.raises(ConfigError, match=match):
        experiment_from_dict(raw)


def test_carbon_dioxide_with_its_total_and_ions_stay_out_of_the_air():
    raw = _breathing_rinse()
    raw["domain"]["air"]["saturation_mol_per_m3"] = {"carbonate": 0.01}
    with pytest.raises(ConfigError, match="acid-base total"):
        experiment_from_dict(raw)
    raw["domain"]["air"]["saturation_mol_per_m3"] = {"sodium": 0.01}
    with pytest.raises(ConfigError, match="an ion"):
        experiment_from_dict(raw)


def test_a_gas_the_air_holds_is_left_out_of_the_saliva_and_the_bulk():
    raw = _breathing_rinse()
    raw["domain"]["mouth"]["saliva_mol_per_m3"]["oxygen"] = 0.1
    with pytest.raises(ConfigError, match=r"leave it out of mouth\.saliva_mol_per_m3"):
        experiment_from_dict(raw)
    raw = _column(100.0, 5.0)
    raw["domain"]["bulk_mol_per_m3"] = {"oxygen": 0.2}
    with pytest.raises(ConfigError, match="borders no bulk liquid"):
        experiment_from_dict(raw)


def test_the_air_reads_back_as_written_and_a_run_at_it_replays(tmp_path, capsys):
    raw = _column(100.0, 5.0, hours=0.01)
    config = experiment_from_dict(raw)
    written = config.domain.to_dict()
    assert written["air"] == {"saturation_mol_per_m3": {"oxygen": SATURATION}}
    again = copy.deepcopy(raw)
    again["domain"] = written
    assert experiment_from_dict(again).domain == config.domain
    path = tmp_path / "column.json"
    path.write_text(json.dumps(raw), "utf-8")
    assert main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "air         the top face is at the air: oxygen held at 0.2 mol per m3" in out
    run = tmp_path / "run"
    assert main(["run", str(path), "-o", str(run)]) == 0
    assert "air         gave oxygen" in capsys.readouterr().out
    assert main(["replay", str(run / "manifest.json")]) == 0
    assert "reproduced the recorded results exactly" in capsys.readouterr().out
