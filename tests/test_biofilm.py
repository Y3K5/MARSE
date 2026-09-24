"""Growth coupled to the solute gradient: the Phase 3 spatial result.

The behaviour being checked is stratification. A biofilm thinner than the
penetration depth should behave like a well-mixed culture, and a thicker one
should confine its growth to a surface layer whose thickness stops increasing,
so that the biofilm-averaged growth rate falls as more biomass is added below.
"""

from __future__ import annotations

import numpy as np
import pytest

from marse.biofilm.biomass import BiofilmState, Population, solve_growth_profile
from marse.spatial.domain import Grid1D
from marse.spatial.solutes import (
    OXYGEN_MOLAR_MASS,
    oxygen_diffusivity_um2_per_h,
    oxygen_saturation_mg_per_l,
)

# Oxygen at 37 C, in the hours-and-micrometres units MARSE works in.
DIFFUSIVITY = 0.43 * float(oxygen_diffusivity_um2_per_h(37.0))
SURFACE = float(oxygen_saturation_mg_per_l(37.0)) / OXYGEN_MOLAR_MASS
DENSITY = 25.0

AEROBE = Population("aerobe", mu_max=0.3, half_saturation=1e-3, yield_per_substrate=0.08)


def biofilm(thickness: float, populations=(AEROBE,), densities=(DENSITY,)) -> BiofilmState:
    grid = Grid1D(thickness, max(100, int(thickness * 3)))
    return BiofilmState.uniform(grid, tuple(populations), list(densities))


def profile(thickness: float, **kwargs):
    return solve_growth_profile(
        biofilm(thickness, **kwargs), diffusivity=DIFFUSIVITY, surface=SURFACE
    )


# --- the state ---------------------------------------------------------------


def test_uniform_state_has_the_expected_shape_and_totals():
    state = biofilm(100.0)
    assert state.density.shape == (1, state.grid.cells)
    np.testing.assert_allclose(state.total_density, DENSITY)
    assert state.areal_density() == pytest.approx(DENSITY * 100.0)
    assert state.names == ("aerobe",)


def test_invalid_states_are_rejected():
    grid = Grid1D(100.0, 10)
    with pytest.raises(ValueError, match="shape"):
        BiofilmState(grid, (AEROBE,), np.ones((1, 5)))
    with pytest.raises(ValueError, match="non-negative"):
        BiofilmState(grid, (AEROBE,), -np.ones((1, 10)))
    with pytest.raises(ValueError, match="unique"):
        BiofilmState.uniform(grid, (AEROBE, AEROBE), [1.0, 1.0])
    with pytest.raises(ValueError, match="at least one population"):
        BiofilmState(grid, (), np.zeros((0, 10)))


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        ({"mu_max": 0.0}, "mu_max"),
        ({"half_saturation": -1.0}, "half_saturation"),
        ({"yield_per_substrate": 0.0}, "yield"),
        ({"name": "  "}, "name"),
    ],
)
def test_invalid_populations_are_rejected(kwargs, expected):
    args = {
        "name": "x",
        "mu_max": 0.3,
        "half_saturation": 1e-3,
        "yield_per_substrate": 0.08,
    } | kwargs
    with pytest.raises(ValueError, match=expected):
        Population(**args)


def test_populations_with_different_affinities_are_refused_for_now():
    a = Population("a", mu_max=0.3, half_saturation=1e-3, yield_per_substrate=0.08)
    b = Population("b", mu_max=0.3, half_saturation=1e-2, yield_per_substrate=0.08)
    with pytest.raises(ValueError, match="share a half_saturation"):
        solve_growth_profile(
            biofilm(100.0, populations=(a, b), densities=(10.0, 10.0)),
            diffusivity=DIFFUSIVITY,
            surface=SURFACE,
        )


# --- the well-mixed limit ----------------------------------------------------


def test_a_thin_biofilm_reduces_to_the_well_mixed_answer():
    """With no gradient to speak of, every depth grows at the surface rate."""
    result = profile(10.0)
    surface_rate = float(result.growth_rate[0][0])
    assert result.mean_growth_rate()[0] == pytest.approx(surface_rate, rel=1e-3)
    assert np.ptp(result.growth_rate[0]) / surface_rate < 1e-2
    assert result.solute.penetration_depth() == result.state.grid.thickness


def test_the_surface_grows_at_close_to_the_unlimited_rate():
    result = profile(200.0)
    # Oxygen at the surface is far above the half-saturation constant.
    assert float(result.growth_rate[0][0]) == pytest.approx(AEROBE.mu_max, rel=0.02)


# --- stratification ----------------------------------------------------------


def test_a_thick_biofilm_confines_growth_to_a_surface_layer():
    result = profile(400.0)
    assert result.active_zone() < 0.5 * 400.0
    assert result.mean_growth_rate()[0] < 0.5 * float(result.growth_rate[0][0])
    # The base is starved however much biomass sits there.
    assert float(result.growth_rate[0][-1]) < 1e-3 * AEROBE.mu_max


def test_the_active_zone_stops_growing_once_the_film_exceeds_it():
    """The signature of stratification: added depth adds no active biomass."""
    zones = [profile(t).active_zone() for t in (200.0, 400.0, 800.0)]
    assert max(zones) - min(zones) < 0.05 * min(zones), zones


def test_mean_growth_falls_inversely_with_thickness_once_stratified():
    """A fixed active layer inside a growing total means mu_mean ~ 1/thickness."""
    thicknesses = [200.0, 400.0, 800.0]
    means = [float(profile(t).mean_growth_rate()[0]) for t in thicknesses]
    for i in range(len(means) - 1):
        assert means[i + 1] == pytest.approx(means[i] / 2.0, rel=0.15), means


def test_growth_decreases_with_depth():
    result = profile(300.0)
    rates = result.growth_rate[0]
    assert np.all(np.diff(rates) <= 1e-12)
    assert rates[0] > rates[-1]


def test_more_biomass_means_a_shallower_active_zone():
    """Denser biofilms consume their supply faster, so it reaches less far."""
    zones = []
    for density in (5.0, 25.0, 100.0):
        state = biofilm(400.0, densities=(density,))
        zones.append(
            solve_growth_profile(state, diffusivity=DIFFUSIVITY, surface=SURFACE).active_zone()
        )
    assert zones == sorted(zones, reverse=True), zones


def test_an_empty_biofilm_neither_consumes_nor_grows():
    state = biofilm(200.0, densities=(0.0,))
    result = solve_growth_profile(state, diffusivity=DIFFUSIVITY, surface=SURFACE)
    np.testing.assert_allclose(result.solute.concentration, SURFACE, rtol=1e-9)
    np.testing.assert_allclose(result.production, 0.0)
    assert result.active_zone() == 0.0
    assert result.active_fraction == 0.0


# --- several populations -----------------------------------------------------


def test_populations_sharing_a_gradient_keep_their_relative_rates():
    """Same affinity and same local oxygen: the ratio of rates is the ratio of mu_max."""
    fast = Population("fast", mu_max=0.5, half_saturation=1e-3, yield_per_substrate=0.08)
    slow = Population("slow", mu_max=0.2, half_saturation=1e-3, yield_per_substrate=0.08)
    result = profile(150.0, populations=(fast, slow), densities=(12.0, 12.0))
    means = result.mean_growth_rate()
    assert means[0] / means[1] == pytest.approx(fast.mu_max / slow.mu_max, rel=1e-6)


def test_a_second_population_deepens_the_oxygen_demand():
    alone = profile(300.0, densities=(12.0,))
    other = Population("other", mu_max=0.3, half_saturation=1e-3, yield_per_substrate=0.08)
    together = profile(300.0, populations=(AEROBE, other), densities=(12.0, 12.0))
    assert together.solute.penetration_depth() < alone.solute.penetration_depth()


def test_production_is_growth_times_biomass():
    result = profile(200.0, populations=(AEROBE,), densities=(DENSITY,))
    np.testing.assert_allclose(result.production, result.growth_rate * DENSITY)


def test_active_zone_threshold_is_validated():
    result = profile(200.0)
    for bad in (0.0, 1.0, -1.0):
        with pytest.raises(ValueError, match="fraction_of_max"):
            result.active_zone(fraction_of_max=bad)


# --- the unit trap -----------------------------------------------------------


def test_the_hours_conversion_matches_the_per_second_value():
    per_second = float(oxygen_diffusivity_um2_per_h(37.0)) / 3600.0
    from marse.spatial.solutes import oxygen_diffusivity_m2_per_s

    assert per_second == pytest.approx(float(oxygen_diffusivity_m2_per_s(37.0)) * 1e12)


def test_mixing_per_second_diffusivity_with_per_hour_rates_is_visibly_wrong():
    """Documents the failure mode the docstring warns about.

    Using a per-second diffusivity alongside per-hour growth rates overstates
    demand 3600-fold and shortens the penetration depth sixtyfold, to a value
    that still looks plausible. The guard against it is the unit conversion in
    `spatial.solutes`, not a runtime check, so this test pins the size of the
    error rather than any error being raised.
    """
    correct = profile(400.0).solute.penetration_depth()
    wrong = solve_growth_profile(
        biofilm(400.0), diffusivity=DIFFUSIVITY / 3600.0, surface=SURFACE
    ).solute.penetration_depth()
    assert correct / wrong == pytest.approx(60.0, rel=0.1)
