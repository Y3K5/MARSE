"""Spreading: biomass that outgrows its voxel pushes the excess on, in a column.

Stage 2d, increment 2d.2 (docs/stage-2d-plan.md). The criteria D4 to D13 were
set before it was built; docs/validation.md ("Spreading in a column") records
what was measured. The scenes:

- **a labelled film**: one species in two neutral labels, a below b, growing
  without limit at a known rate. Every material point moves from z to
  z e^(mu t), so the labels' boundary has a closed form;
- **a fed film**: one species on glucose that diffuses in from the bulk
  liquid. Once the film is deeper than the glucose reaches, it thickens at
  Y J / rho, set by the flux J that the validated steady solver (case V3)
  gives for the same column.
"""

import json
import math
from contextlib import contextmanager
from itertools import pairwise
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from marse.cli import main
from marse.core.config import ConfigError
from marse.core.provenance import Manifest
from marse.core.reactive_transport import run
from marse.core.simulation import ConservationError
from marse.schemas import experiment_from_dict
from marse.schemas.domain import DEFAULT_SPREADING_INTERVAL_H, Domain
from marse.spatial.diffusion import solve_steady_state
from marse.spatial.domain import Grid1D
from marse.spatial.spreading import (
    ContinuumSpreading,
    SpreadingError,
    _column_pressure,
    volume_fraction,
)

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "networks" / "spreading_column.json"
LN2 = math.log(2.0)
RHO = 1000.0  # packing density of every species here, C-mol per m3


def species(name: str, density: float | None = RHO) -> dict:
    raw = {
        "name": name,
        "phase": "particulate",
        "formula": "CH1.8O0.5",
        "oxygen_role": "aerotolerant",
    }
    if density is not None:
        raw["density_mol_per_m3"] = density
    return raw


SOLUTES = [
    {"name": "glucose", "phase": "dissolved", "formula": "C6H12O6"},
    {"name": "ethanol", "phase": "dissolved", "formula": "C2H6O"},
    {"name": "carbon_dioxide", "phase": "dissolved", "formula": "CO2"},
]
DIFFUSIVITY = {"glucose": 6.7e-10, "ethanol": 1.2e-9, "carbon_dioxide": 1.9e-9}


def fermentation(name: str, rate: dict) -> dict:
    return {
        "name": f"{name}_growth",
        "kind": "growth",
        "biomass": name,
        "substrate": "glucose",
        "yield_mol_per_mol": 1.5,
        "balanced_by": ["ethanol", "carbon_dioxide"],
        "rate": rate,
    }


def labelled_film(
    voxel_um: float = 2.0,
    *,
    mu: float = LN2,
    duration_h: float = 3.0,
    interval_h: float | None = 0.05,
    height_um: float = 200.0,
    film_um: float = 20.0,
    timestep_h: float = 0.25,
) -> dict:
    """Labels a and b of one species, a film of film_um, growing without limit at mu.

    Glucose is held so high that it never limits, which the processes state as
    an assumption, so every voxel of biomass grows at exactly mu. The film starts
    as b; :func:`lower_layer` relabels its lower part a.
    """
    unlimited = {"maximum_per_h": mu, "assumed_in_excess": ["glucose"]}
    raw = {
        "schema_version": 2,
        "experiment_id": "labelled-film",
        "components": [*map(dict, SOLUTES), species("a"), species("b")],
        "processes": [fermentation("a", unlimited), fermentation("b", unlimited)],
        "initial_mol_per_m3": {"glucose": 1000.0},
        "duration_h": duration_h,
        "timestep_h": timestep_h,
        "domain": {
            "voxels": [round(height_um / voxel_um)],
            "voxel_um": voxel_um,
            "bulk_mol_per_m3": {"glucose": 1000.0},
            "diffusivity_m2_per_s": dict(DIFFUSIVITY),
            "colonies": [
                {
                    "component": "b",
                    "center_um": [],
                    "radius_um": film_um,
                    "concentration_mol_per_m3": RHO,
                }
            ],
            "spreading": {"mechanism": "continuum"},
        },
    }
    if interval_h is not None:
        raw["domain"]["spreading"]["interval_h"] = interval_h
    return raw


@contextmanager
def lower_layer(name: str, below_um: float, other: str):
    """Start with component ``name`` in place of ``other`` below ``below_um``.

    Colonies are hemispheres from the substratum, so two of them cannot make a
    layer of b on a layer of a; the test sets the starting state directly.
    """
    original = Domain.initial_state

    def initial_state(self, names, uniform, seed):
        state = original(self, names, uniform, seed)
        below = (np.arange(state.shape[-1]) + 0.5) * self.voxel_um < below_um
        state[names.index(name)][below] = state[names.index(other)][below]
        state[names.index(other)][below] = 0.0
        return state

    with mock.patch.object(Domain, "initial_state", initial_state):
        yield


def labels(result) -> tuple[np.ndarray, np.ndarray, float]:
    names = result.component_names
    a, b = result.final_state[names.index("a")], result.final_state[names.index("b")]
    return a, b, result.config.domain.voxel_um


MU, K, Y, SB = 1.0, 0.05, 1.5, 0.02
"""The fed film: glucose at 20 uM above, consumed with K = 50 uM, mu 1/h, yield 1.5."""


def fed_film(
    film_um: float,
    *,
    voxel_um: float = 1.0,
    liquid_um: float = 50.0,
    duration_h: float = 1.0,
    interval_h: float = 0.05,
    timestep_h: float = 0.1,
    relative_tolerance: float | None = None,
) -> dict:
    monod = {"component": "glucose", "form": "monod", "half_saturation_mol_per_m3": K}
    raw = {
        "schema_version": 2,
        "experiment_id": "fed-film",
        "components": [*map(dict, SOLUTES), species("bug")],
        "processes": [fermentation("bug", {"maximum_per_h": MU, "factors": [monod]})],
        "initial_mol_per_m3": {"glucose": SB},
        "duration_h": duration_h,
        "timestep_h": timestep_h,
        "domain": {
            "voxels": [round((film_um + liquid_um) / voxel_um)],
            "voxel_um": voxel_um,
            "bulk_mol_per_m3": {"glucose": SB},
            "diffusivity_m2_per_s": dict(DIFFUSIVITY),
            "colonies": [
                {
                    "component": "bug",
                    "center_um": [],
                    "radius_um": film_um,
                    "concentration_mol_per_m3": RHO,
                }
            ],
            "spreading": {"mechanism": "continuum", "interval_h": interval_h},
        },
    }
    if relative_tolerance is not None:
        raw["relative_tolerance"] = relative_tolerance
    return raw


def thickening_by_v3(phi: np.ndarray, voxel_um: float, refine: int = 4) -> float:
    """Y J / rho in um per hour, J from the validated steady solver (case V3).

    V3 solves the same column, from the bulk held at its top face down to the
    substratum, with each voxel's biomass as uptake capacity, on nodes four
    times finer than the voxels.
    """
    height = phi.size * voxel_um
    grid = Grid1D(height, phi.size * refine)
    voxel = np.clip(((height - grid.depths) // voxel_um).astype(int), 0, phi.size - 1)
    capacity = MU * RHO * phi[voxel] / Y / 3600.0  # glucose, mol per m3 per s
    profile = solve_steady_state(
        grid,
        diffusivity=DIFFUSIVITY["glucose"] * 1e12,
        surface=SB,
        max_uptake=capacity,
        half_saturation=K,
    )
    return Y * profile.total_uptake(capacity, K) * 3600.0 / RHO


# --- the schema: densities, the spreading block and what it refuses (D12) -------------


def test_a_density_on_a_dissolved_component_is_refused():
    raw = labelled_film()
    raw["components"][0]["density_mol_per_m3"] = 100
    with pytest.raises(ConfigError, match="'glucose' is dissolved; only a particulate"):
        experiment_from_dict(raw)


def test_a_density_must_be_positive():
    raw = labelled_film()
    raw["components"][3]["density_mol_per_m3"] = 0
    with pytest.raises(ConfigError, match="density_mol_per_m3: must be positive"):
        experiment_from_dict(raw)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda r: r["components"].__setitem__(3, species("a", None)), "give it a density"),
        (lambda r: r["domain"]["spreading"].update(interval_h=0), "interval_h: must be positive"),
        (lambda r: r["domain"]["spreading"].update(mechanism="shoving"), "expected one of"),
        (lambda r: r["domain"].update(voxels=[4, 100]), "arrives with Stage 2d, increment 2d.3"),
        (lambda r: r["domain"].update(voxels=[1]), "at least two voxels"),
    ],
)
def test_a_domain_that_cannot_spread_is_refused_with_the_reason(change, message):
    raw = labelled_film()
    change(raw)
    if raw["domain"]["voxels"] == [4, 100]:
        raw["domain"]["colonies"][0]["center_um"] = [4.0]
    with pytest.raises(ConfigError, match=message):
        experiment_from_dict(raw)


def test_spreading_under_a_salivary_film_is_refused_until_it_is_supported():
    raw = json.loads((ROOT / "examples/environments/oral/stephan_rinse.json").read_text("utf-8"))
    bacteria = next(c for c in raw["components"] if c["name"] == "bacteria")
    bacteria["density_mol_per_m3"] = 800.0
    raw["domain"]["spreading"] = {"mechanism": "continuum"}
    with pytest.raises(
        ConfigError, match="continuum spreading under a salivary film is not supported"
    ):
        experiment_from_dict(raw)


def test_a_start_that_overfills_a_voxel_is_refused_naming_where():
    raw = labelled_film()
    raw["domain"]["colonies"].append(
        {"component": "a", "center_um": [], "radius_um": 10.0, "concentration_mol_per_m3": 1.0}
    )
    with pytest.raises(ConfigError, match=r"fills 1\.001 of the voxel at 1 um"):
        experiment_from_dict(raw)


def test_two_species_share_the_room_at_the_start():
    raw = labelled_film()
    raw["domain"]["colonies"][0]["concentration_mol_per_m3"] = RHO / 2
    raw["domain"]["colonies"].append(
        {"component": "a", "center_um": [], "radius_um": 20.0, "concentration_mol_per_m3": RHO / 2}
    )
    experiment_from_dict(raw)  # half and half: exactly full
    raw["domain"]["colonies"][1]["concentration_mol_per_m3"] = RHO / 2 * 1.01
    with pytest.raises(ConfigError, match="more than it can hold"):
        experiment_from_dict(raw)


def test_biomass_that_starts_in_the_top_layer_is_refused():
    raw = labelled_film(height_um=20.0)
    with pytest.raises(ConfigError, match="biomass starts in the top layer"):
        experiment_from_dict(raw)


def test_the_default_interval_is_the_measured_one():
    config = experiment_from_dict(labelled_film(interval_h=None))
    assert config.domain.spreading.interval_h == DEFAULT_SPREADING_INTERVAL_H == 0.25


def test_a_domain_that_spreads_survives_the_trip_through_json():
    config = experiment_from_dict(labelled_film())
    written = json.loads(json.dumps(config.to_dict()))
    assert written["domain"]["spreading"] == {
        "mechanism": "continuum",
        "interval_h": 0.05,
        "carried": [],
    }
    assert written["components"][3]["density_mol_per_m3"] == RHO
    assert experiment_from_dict(written) == config


def test_a_domain_that_does_not_spread_writes_nothing_new():
    raw = labelled_film()
    del raw["domain"]["spreading"]
    written = experiment_from_dict(raw).to_dict()
    assert "spreading" not in written["domain"]


# --- the column's pressure and sweep, alone ---------------------------------------------


def column(*layers: tuple[float, float], voxels: int = 12) -> tuple[np.ndarray, np.ndarray]:
    """A state of two moving components and one that takes no room, from the bottom up."""
    state = np.zeros((3, voxels))
    for v, (a, b) in enumerate(layers):
        state[0, v], state[1, v] = a * RHO, b * RHO
    state[2] = np.linspace(1.0, 2.0, voxels)  # a solute, which spreading must not touch
    return state, np.array([RHO, RHO, 0.0])


@pytest.mark.numerical
@pytest.mark.parametrize("seed", range(20))
def test_the_column_pressure_drives_the_wanner_gujer_displacement(seed):
    # D6: in a film full from the substratum, the flux through each face is the
    # excess summed below it, the discrete u(z) = integral of the growth.
    rng = np.random.default_rng(seed)
    filled = int(rng.integers(3, 40))
    full = np.zeros(filled + 10, dtype=bool)
    full[:filled] = True
    excess = np.zeros(full.size)
    excess[:filled] = rng.uniform(0.0, 0.3, filled)
    pressure = _column_pressure(full, excess)
    flux = pressure[:-1] - pressure[1:]
    np.testing.assert_allclose(flux[:filled], np.cumsum(excess[:filled]), rtol=1e-12, atol=0)
    assert np.all(flux[filled:] == 0)


def test_a_floating_full_region_spreads_both_ways_and_conserves():
    state, densities = column((0.5, 0), (0, 0), (1.2, 0), (1.0, 0), (0.9, 0))
    spread, stats = ContinuumSpreading((12,), densities).spread(state)
    phi = volume_fraction(spread, densities)
    assert phi.max() <= 1.0 + 1e-12
    assert phi[1] > 0  # some went down into the gap
    np.testing.assert_allclose(spread[:2].sum(axis=1), state[:2].sum(axis=1), rtol=1e-15)
    assert stats.rounds >= 1


@pytest.mark.invariance
@pytest.mark.parametrize("seed", range(30))
def test_spreading_conserves_stays_positive_and_fills_no_voxel_beyond_its_room(seed):
    # D4, D5 and D10 for the mechanism alone: random overfilled columns, mixed
    # species, gaps, growth up to fourfold in a step.
    rng = np.random.default_rng(seed)
    voxels = 40
    layers = []
    for _ in range(int(rng.integers(1, 20))):
        total = rng.uniform(0.0, 4.0) if rng.random() < 0.8 else 0.0
        share = rng.random()
        layers.append((total * share, total * (1 - share)))
    state, densities = column(*layers, voxels=voxels)
    spread, _ = ContinuumSpreading((voxels,), densities).spread(state)
    assert np.all(spread >= 0)
    assert volume_fraction(spread, densities).max() <= 1.0 + 1e-12
    np.testing.assert_allclose(spread[:2].sum(axis=1), state[:2].sum(axis=1), rtol=1e-14)
    assert np.array_equal(spread[2], state[2])  # what takes no room stays where it was


def test_a_layer_stays_a_layer_as_it_is_pushed():
    # Ordered transport: a below b, and the bottom voxel triples. The material
    # that leaves a voxel through a face is the material nearest that face, so
    # a stays below b, mixed only in the one voxel that straddles them.
    state, densities = column((3.0, 0), (1.0, 0), (0, 1.0), (0, 1.0), (0, 0.5))
    spread, _ = ContinuumSpreading((12,), densities).spread(state)
    a, b = spread[0] / RHO, spread[1] / RHO
    np.testing.assert_allclose(a[:4], [1, 1, 1, 1], rtol=1e-14)
    np.testing.assert_allclose(b[4:6], [1, 1], rtol=1e-14)
    assert a[4:].max() == 0
    assert b[:4].max() == 0


def test_a_full_box_cannot_spread():
    state, densities = column(*[(1.0, 0)] * 11 + [(1.1, 0)], voxels=12)
    with pytest.raises(SpreadingError, match="the box is full"):
        ContinuumSpreading((12,), densities).spread(state)


# --- the engine: conservation, capacity and positivity (D4, D5, D10) ---------------------


def test_a_run_that_spreads_records_its_structure_and_balance():
    with lower_layer("a", 10.0, "b"):
        result = run(experiment_from_dict(labelled_film(duration_h=1.0)))
    outputs = result.manifest.outputs
    assert result.manifest.models["spreading"] == "continuum_pressure_v1"
    assert outputs["spreading"]["spreads"] == 20
    assert np.all(result.structure["largest_volume_fraction"] <= 1.0 + 1e-12)
    # Unlimited growth for an hour doubles the film, to the integrator's tolerance.
    assert result.structure["biovolume_um3_per_um2"][-1] == pytest.approx(40.0, rel=1e-4)
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    assert worst <= 1e-12


@pytest.mark.parametrize("interval_h", [0.01, 1.0, 10.0])
def test_any_interval_conserves_and_stays_positive(interval_h):
    # D10: a thousand spreads, ten, or one after ten hours of growth in place.
    with lower_layer("a", 10.0, "b"):
        config = experiment_from_dict(
            labelled_film(
                mu=0.05,
                duration_h=10.0,
                timestep_h=10.0 if interval_h == 10.0 else 1.0,
                interval_h=interval_h,
                height_um=60.0,
            )
        )
        result = run(config)
    outputs = result.manifest.outputs
    assert outputs["spreading"]["spreads"] == round(10.0 / interval_h)
    assert result.final_state.min() >= 0
    assert max(q["largest_relative_residual"] for q in outputs["balance"].values()) <= 1e-12
    assert result.structure["biovolume_um3_per_um2"][-1] == pytest.approx(
        20.0 * math.exp(0.5), rel=1e-4
    )


@pytest.mark.slow
def test_ten_thousand_spreads_conserve_to_the_ledger_tolerance():
    # D4, as set: 10^4 spreading steps, every one checked against the ledger.
    with lower_layer("a", 10.0, "b"):
        config = experiment_from_dict(
            labelled_film(
                mu=0.05, duration_h=10.0, timestep_h=1.0, interval_h=0.001, height_um=60.0
            )
        )
        result = run(config)
    outputs = result.manifest.outputs
    assert outputs["spreading"]["spreads"] == 10_000
    assert max(q["largest_relative_residual"] for q in outputs["balance"].values()) <= 1e-12
    assert result.final_state.min() >= 0


def test_a_leak_in_spreading_is_caught_at_the_first_spread(monkeypatch):
    honest = ContinuumSpreading.spread

    def leaky(self, state):
        spread, stats = honest(self, state)
        spread[self.moving] *= 1.0 - 1e-6
        return spread, stats

    monkeypatch.setattr(ContinuumSpreading, "spread", leaky)
    with pytest.raises(ConservationError, match=r"changed the total of a component .* t = 0\.05 h"):
        run(experiment_from_dict(labelled_film(duration_h=1.0)))


def test_a_mechanism_that_overfills_a_voxel_is_caught(monkeypatch):
    def piling(self, state):
        piled = np.array(state)
        piled[self.moving, 1] += piled[self.moving, 2]
        piled[self.moving, 2] = 0.0
        return piled, None

    monkeypatch.setattr(ContinuumSpreading, "spread", piling)
    with pytest.raises(ConservationError, match="no voxel may hold more than all of it"):
        run(experiment_from_dict(labelled_film(duration_h=1.0)))


def test_a_mechanism_that_moves_a_solute_is_caught(monkeypatch):
    honest = ContinuumSpreading.spread

    def stirring(self, state):
        spread, stats = honest(self, state)
        spread[~self.moving] = spread[~self.moving][..., ::-1]
        return spread, stats

    monkeypatch.setattr(ContinuumSpreading, "spread", stirring)
    with pytest.raises(
        ConservationError, match="changed a component that neither takes room nor is carried"
    ):
        run(experiment_from_dict(labelled_film(duration_h=1.0)))


def test_a_film_that_reaches_the_top_layer_stops_the_run_naming_the_time():
    # D12: 20 um doubling every hour in a 60 um box reaches the top layer
    # between 1 and 2 hours.
    with pytest.raises(SpreadingError, match=r"at t = 1\.\d+ h the biofilm reached the top layer"):
        run(experiment_from_dict(labelled_film(height_um=60.0, duration_h=3.0)))


# --- against closed forms and the validated solver (D6 to D9) --------------------------


@pytest.mark.numerical
def test_a_labelled_band_moves_as_the_film_stretches():
    # D8: a fills 0 to 10 um and b 10 to 20 um; three doublings move the
    # boundary to 80 um. Its volume is exact, so the band must sit where a's
    # volume per area ends, within a voxel.
    with lower_layer("a", 10.0, "b"):
        result = run(experiment_from_dict(labelled_film(1.25)))
    a, b, voxel = labels(result)
    share = np.divide(a, a + b, out=np.zeros_like(a), where=(a + b) > 0)
    k = int(np.flatnonzero(share < 0.5)[0])
    centre = (k - 0.5) * voxel  # the centre of voxel k - 1
    band = centre + (share[k - 1] - 0.5) / (share[k - 1] - share[k]) * voxel
    exact = a.sum() * voxel / RHO
    assert exact == pytest.approx(80.0, rel=1e-4)
    assert abs(band - exact) < voxel


@pytest.mark.slow
@pytest.mark.numerical
def test_the_labelled_band_converges_as_the_voxels_shrink():
    # D8: the centroid of a, whose exact value is half its height, converges
    # at better than first order.
    errors = []
    for voxel in (2.5, 1.25, 0.625):
        with lower_layer("a", 10.0, "b"):
            result = run(experiment_from_dict(labelled_film(voxel)))
        a, _, _ = labels(result)
        heights = (np.arange(a.size) + 0.5) * voxel
        centroid = float((a * heights).sum() / a.sum())
        errors.append(abs(centroid - a.sum() * voxel / RHO / 2))
    orders = [math.log2(e0 / e1) for e0, e1 in pairwise(errors)]
    assert errors[-1] < 0.625
    assert min(orders) >= 0.9, (errors, orders)


@pytest.mark.numerical
def test_a_deep_film_thickens_at_the_rate_its_flux_sets_whatever_its_depth():
    # D9: films of 100 and 200 um, under the same 50 um of liquid, are both far
    # deeper than glucose reaches (about 25 um). Only their top layers grow, so
    # they thicken at the same rate, Y J / rho, with J the flux V3 gives.
    rates = []
    for film in (100.0, 200.0):
        config = experiment_from_dict(fed_film(film))
        result = run(config)
        biovolume = result.structure["biovolume_um3_per_um2"]
        rate = (biovolume[-1] - biovolume[-2]) / (result.times_h[-1] - result.times_h[-2])
        phi = volume_fraction(result.final_state, config.network.densities())
        reference = thickening_by_v3(phi, 1.0)
        assert rate == pytest.approx(reference, rel=0.01)
        rates.append(rate)
    assert rates[1] == pytest.approx(rates[0], rel=0.02)


@pytest.mark.slow
@pytest.mark.numerical
def test_the_splitting_converges_at_first_order_in_the_interval():
    # D7: the biovolume after an hour, against a run spread every 0.003125 h.
    def biovolume(interval_h: float) -> float:
        config = experiment_from_dict(
            fed_film(100.0, interval_h=interval_h, relative_tolerance=1e-7)
        )
        return float(run(config).structure["biovolume_um3_per_um2"][-1])

    reference = biovolume(0.003125)
    errors = [abs(biovolume(dt) - reference) for dt in (0.1, 0.05, 0.025)]
    orders = [math.log2(e0 / e1) for e0, e1 in pairwise(errors)]
    assert min(orders) >= 0.9, (errors, orders)


# --- reproducibility and the command line (D11, D13) -------------------------------------


def test_marse_check_reports_the_packing_and_the_spreading(capsys):
    assert main(["check", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert (
        "spreading: continuum, every 0.25 h; packing densities, mol per m3: heterotroph 1000" in out
    )
    assert "the fullest voxel starts 1 full" in out


def test_marse_run_writes_the_structure_and_replays_exactly(tmp_path, capsys):
    raw = json.loads(EXAMPLE.read_text("utf-8"))
    raw["duration_h"] = 1.0
    path = tmp_path / "column.json"
    path.write_text(json.dumps(raw), "utf-8")
    assert main(["run", str(path), "-o", str(tmp_path / "out")]) == 0
    out = capsys.readouterr().out
    assert "spreading   continuum, every 0.25 h: 4 spreads" in out
    structure = (tmp_path / "out" / "structure.csv").read_text("utf-8").splitlines()
    assert structure[0] == (
        "time_h,biovolume_um3_per_um2,maximum_thickness_um,largest_volume_fraction"
    )
    assert structure[1].startswith("0.000000,40,40,1")
    manifest = tmp_path / "out" / "manifest.json"
    assert Manifest.read(manifest).models["spreading"] == "continuum_pressure_v1"
    assert main(["replay", str(manifest)]) == 0
    assert "replay reproduced the recorded results exactly" in capsys.readouterr().out
