"""Unit tests for growth kinetics, cardinal models and oxygen properties."""

import numpy as np
import pytest

from marse.microbes.cardinal import cardinal_ph, cardinal_temperature, ratkowsky
from marse.microbes.growth import (
    baranyi_roberts,
    doubling_time,
    haldane,
    monod,
    product_formation_rate,
    specific_growth_rate,
    substrate_uptake_rate,
)
from marse.spatial.solutes import (
    OXYGEN_MOLAR_MASS,
    oxygen_diffusivity_m2_per_s,
    oxygen_saturation_mg_per_l,
)

# --- primary growth models -------------------------------------------------------


def test_doubling_time_and_growth_rate_are_inverse():
    assert doubling_time(specific_growth_rate(20.0)) == pytest.approx(20.0)
    assert specific_growth_rate(1 / 3) == pytest.approx(2.0794, abs=1e-4)  # 20 min, per hour


def test_monod_shape():
    assert monod(0.0, 2.0, 0.5) == 0.0
    assert monod(0.5, 2.0, 0.5) == pytest.approx(1.0)  # half-saturation
    assert monod(1e9, 2.0, 0.5) == pytest.approx(2.0)
    assert monod(-1e-12, 2.0, 0.5) == 0.0  # numerical undershoot counts as zero


def test_monod_is_elementwise_and_broadcasts():
    rates = monod(np.array([[0.0, 0.5], [1.5, 1e9]]), mu_max=2.0, k_s=np.array([0.5, 0.5]))
    np.testing.assert_allclose(rates, [[0.0, 1.0], [1.5, 2.0]], rtol=1e-8)


@pytest.mark.parametrize(("mu_max", "k_s"), [(-1.0, 0.5), (1.0, 0.0)])
def test_monod_rejects_invalid_parameters(mu_max, k_s):
    with pytest.raises(ValueError, match="must be"):
        monod(1.0, mu_max, k_s)


def test_haldane_peaks_at_geometric_mean_and_tends_to_monod():
    s = np.linspace(0.01, 50, 50_000)
    rates = haldane(s, 1.0, k_s=0.5, k_i=8.0)
    assert s[np.argmax(rates)] == pytest.approx(np.sqrt(0.5 * 8.0), rel=1e-3)
    assert haldane(2.0, 1.0, 0.5, 1e12) == pytest.approx(monod(2.0, 1.0, 0.5))


def test_pirt_uptake_and_luedeking_piret_production():
    assert substrate_uptake_rate(0.5, 2.0, yield_coefficient=0.5) == pytest.approx(2.0)
    assert substrate_uptake_rate(0.0, 2.0, 0.5, maintenance=0.1) == pytest.approx(0.2)
    assert product_formation_rate(0.5, 2.0, alpha=0.3, beta=0.1) == pytest.approx(0.5)


def test_baranyi_roberts_limits():
    t = np.array([0.0, 60.0])
    y = baranyi_roberts(t, y0=-2.0, y_max=5.0, mu_max=1.2, lag=3.0)
    assert y[0] == pytest.approx(-2.0, abs=1e-12)
    assert y[1] == pytest.approx(5.0, abs=1e-9)


def test_baranyi_roberts_exponential_phase_is_shifted_by_the_lag():
    mu, lag = 0.8, 2.5
    # The shift is exact asymptotically: well after the lag, far below the maximum.
    t = np.linspace(12.0, 16.0, 5)
    y = baranyi_roberts(t, y0=0.0, y_max=40.0, mu_max=mu, lag=lag)
    np.testing.assert_allclose(y, mu * (t - lag), atol=1e-3)


def test_baranyi_roberts_is_monotonic_and_without_lag_is_exponential():
    t = np.linspace(0, 30, 301)
    y = baranyi_roberts(t, y0=0.0, y_max=12.0, mu_max=1.0, lag=4.0)
    assert np.all(np.diff(y) >= 0)
    early = np.linspace(0, 2, 5)
    no_lag = baranyi_roberts(early, y0=0.0, y_max=50.0, mu_max=1.0, lag=0.0)
    np.testing.assert_allclose(no_lag, early, atol=1e-9)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"t": -1.0},
        {"mu_max": 0.0},
        {"lag": -1.0},
        {"y_max": -5.0},
        {"curvature": 0.0},
    ],
)
def test_baranyi_roberts_rejects_invalid_input(kwargs):
    args = {"t": 1.0, "y0": 0.0, "y_max": 5.0, "mu_max": 1.0, "lag": 1.0} | kwargs
    with pytest.raises(ValueError, match="must"):
        baranyi_roberts(**args)


# --- secondary (cardinal) models -------------------------------------------------

T_MIN, T_OPT, T_MAX = 6.0, 40.0, 47.0


def test_ctmi_equals_one_at_optimum_and_zero_at_and_beyond_the_limits():
    assert cardinal_temperature(T_OPT, T_MIN, T_OPT, T_MAX) == pytest.approx(1.0)
    gammas = cardinal_temperature([0.0, T_MIN, T_MAX, 60.0], T_MIN, T_OPT, T_MAX)
    np.testing.assert_array_equal(gammas, 0.0)


def test_ctmi_has_its_maximum_at_the_optimum():
    temps = np.linspace(T_MIN, T_MAX, 41_001)
    gammas = cardinal_temperature(temps, T_MIN, T_OPT, T_MAX)
    assert temps[np.argmax(gammas)] == pytest.approx(T_OPT, abs=1e-3)
    assert gammas.max() <= 1.0 + 1e-12
    assert np.all(gammas >= 0)


def test_ctmi_is_the_same_in_celsius_and_kelvin():
    kelvin = cardinal_temperature(310.15, T_MIN + 273.15, T_OPT + 273.15, T_MAX + 273.15)
    assert kelvin == pytest.approx(cardinal_temperature(37.0, T_MIN, T_OPT, T_MAX))


@pytest.mark.parametrize(("t_min", "t_opt", "t_max"), [(10, 5, 40), (5, 20, 40)])
def test_ctmi_rejects_invalid_or_ill_posed_cardinal_values(t_min, t_opt, t_max):
    with pytest.raises(ValueError, match=r"CTMI|cardinal"):
        cardinal_temperature(25.0, t_min, t_opt, t_max)


def test_cardinal_ph_model():
    assert cardinal_ph(7.0, 4.0, 7.0, 9.0) == pytest.approx(1.0)
    np.testing.assert_array_equal(cardinal_ph([3.0, 4.0, 9.0, 10.0], 4.0, 7.0, 9.0), 0.0)
    values = np.linspace(4.0, 9.0, 5001)
    gammas = cardinal_ph(values, 4.0, 7.0, 9.0)
    assert values[np.argmax(gammas)] == pytest.approx(7.0, abs=1e-3)
    assert gammas.max() <= 1.0 + 1e-12


def test_ratkowsky_square_root_model():
    b, t_min, t_max, c = 0.03, 5.0, 48.0, 0.3
    np.testing.assert_array_equal(ratkowsky([t_min, t_max, 50.0], b, t_min, t_max, c), 0.0)
    # Far below the maximum, sqrt(mu) grows linearly with temperature (Ratkowsky 1982).
    low = np.array([10.0, 15.0, 20.0])
    np.testing.assert_allclose(
        np.sqrt(ratkowsky(low, b, t_min, t_max, c)), b * (low - t_min), rtol=1e-3
    )


# --- oxygen ------------------------------------------------------------------------


def test_oxygen_saturation_matches_standard_table_values():
    # Air-saturated fresh water at 1 atm: 9.09 mg/L at 20 C and 8.26 mg/L at 25 C.
    assert oxygen_saturation_mg_per_l(20.0) == pytest.approx(9.09, abs=0.01)
    assert oxygen_saturation_mg_per_l(25.0) == pytest.approx(8.26, abs=0.01)
    at_37 = oxygen_saturation_mg_per_l(37.0) / OXYGEN_MOLAR_MASS  # mM
    assert at_37 == pytest.approx(0.210, abs=0.001)
    assert np.all(np.diff(oxygen_saturation_mg_per_l(np.arange(0, 41))) < 0)


def test_oxygen_diffusivity_rises_with_temperature():
    assert oxygen_diffusivity_m2_per_s(25.0) == pytest.approx(2.0e-9, rel=0.02)
    assert np.all(np.diff(oxygen_diffusivity_m2_per_s(np.arange(0, 96))) > 0)


@pytest.mark.parametrize("function", [oxygen_saturation_mg_per_l, oxygen_diffusivity_m2_per_s])
def test_oxygen_correlations_refuse_to_extrapolate(function):
    with pytest.raises(ValueError, match="valid from"):
        function(-5.0)
    with pytest.raises(ValueError, match="valid from"):
        function(120.0)
