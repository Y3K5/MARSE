"""Validation case V3: the reaction-diffusion solver, independently of any biology.

The solver is checked against three analytical limits that bracket real Monod
uptake — none, first order and zero order — plus a flux balance that ties the
boundary condition to the interior, and a convergence study confirming the
discretisation is second order. A solver that passes these is transporting
correctly, whatever the biology later asks of it.
"""

from __future__ import annotations

import numpy as np
import pytest

from marse.spatial.diffusion import ConvergenceError, solve_steady_state
from marse.spatial.domain import Grid1D
from marse.validation.analytical import (
    first_order_profile,
    first_order_surface_flux,
    zero_order_penetration_depth,
)

pytestmark = pytest.mark.numerical

# Oxygen in a dense biofilm at 37 C (docs/parameters.md): um, s, mM.
DIFFUSIVITY = 1128.0
SURFACE = 0.210
THICKNESS = 200.0


# --- the grid ----------------------------------------------------------------


def test_grid_geometry():
    grid = Grid1D(thickness=100.0, cells=4)
    assert grid.dx == 25.0
    np.testing.assert_allclose(grid.depths, [25.0, 50.0, 75.0, 100.0])
    assert grid.depths[-1] == grid.thickness


@pytest.mark.parametrize(("thickness", "cells"), [(0.0, 10), (-1.0, 10), (10.0, 1), (10.0, 2.5)])
def test_invalid_grids_are_rejected(thickness, cells):
    with pytest.raises(ValueError, match=r"thickness|cells"):
        Grid1D(thickness=thickness, cells=cells)


def test_integration_includes_the_surface():
    """Omitting the boundary value would bias the integral by half a spacing."""
    grid = Grid1D(thickness=1.0, cells=4)
    # A constant field of 1 integrates to the thickness.
    assert grid.integrate(1.0, np.ones(4)) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="node values"):
        grid.integrate(1.0, np.ones(3))


# --- limit 1: no uptake ------------------------------------------------------


def test_without_uptake_the_profile_is_uniform():
    profile = solve_steady_state(
        Grid1D(THICKNESS, 200),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=0.0,
        half_saturation=1.0,
    )
    np.testing.assert_allclose(profile.concentration, SURFACE, rtol=1e-12)
    assert profile.penetration_depth() == THICKNESS
    assert profile.anoxic_fraction == 0.0


# --- limit 2: first order (K much larger than C) -----------------------------


def test_first_order_limit_matches_the_analytical_profile():
    half_saturation, max_uptake = 1e4, 5.0  # K >> C, so uptake is ~linear in C
    profile = solve_steady_state(
        Grid1D(THICKNESS, 400),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=max_uptake,
        half_saturation=half_saturation,
    )
    exact = first_order_profile(
        profile.depths,
        thickness=THICKNESS,
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=max_uptake,
        half_saturation=half_saturation,
    )
    assert np.max(np.abs(profile.concentration - exact)) / SURFACE < 1e-5

    expected_flux = first_order_surface_flux(
        thickness=THICKNESS,
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=max_uptake,
        half_saturation=half_saturation,
    )
    assert profile.surface_flux(DIFFUSIVITY) == pytest.approx(expected_flux, rel=1e-3)


def test_first_order_profile_is_stable_for_a_very_thick_slab():
    """cosh(L/lambda) overflows if written naively; the reference must not."""
    values = first_order_profile(
        [0.0, 5.0, 5000.0],
        thickness=5000.0,
        diffusivity=1.0,
        surface=1.0,
        max_uptake=1.0,
        half_saturation=1e-4,  # lambda = 0.01, so L/lambda = 500_000
    )
    assert np.all(np.isfinite(values))
    assert values[0] == pytest.approx(1.0)
    assert values[-1] >= 0.0


# --- limit 3: zero order (K much smaller than C) -----------------------------


def test_zero_order_limit_matches_the_parabolic_profile_and_its_depth():
    half_saturation, max_uptake = 1e-6, 0.19
    profile = solve_steady_state(
        Grid1D(THICKNESS, 2000),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=max_uptake,
        half_saturation=half_saturation,
    )
    delta = zero_order_penetration_depth(DIFFUSIVITY, SURFACE, max_uptake)

    z = profile.depths
    expected = np.where(z < delta, SURFACE * (1.0 - z / delta) ** 2, 0.0)
    assert np.max(np.abs(profile.concentration - expected)) / SURFACE < 0.01

    # The measured depth uses a 1% threshold, and C = 0.01 C0 where
    # (1 - z/delta)^2 = 0.01, i.e. at 0.9 delta. The two definitions differ by
    # that factor by construction, not by error.
    assert profile.penetration_depth() == pytest.approx(0.9 * delta, rel=0.02)
    assert profile.penetration_depth(threshold_fraction=0.25) == pytest.approx(
        0.5 * delta, rel=0.02
    )


def test_a_thick_biofilm_is_stratified():
    """The result of theory.md section 5.1: an active surface layer over an anoxic interior."""
    profile = solve_steady_state(
        Grid1D(THICKNESS, 800),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=0.19,
        half_saturation=1e-4,
    )
    assert 0.6 < profile.anoxic_fraction < 0.95
    assert profile.penetration_depth() < THICKNESS
    assert profile.concentration[-1] < 1e-3 * SURFACE


# --- conservation and convergence -------------------------------------------


@pytest.mark.parametrize("half_saturation", [1e-4, 0.01, 0.5, 100.0])
def test_surface_flux_equals_total_uptake(half_saturation):
    """Every molecule entering the slab must be consumed inside it."""
    max_uptake = 0.05
    profile = solve_steady_state(
        Grid1D(THICKNESS, 1600),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=max_uptake,
        half_saturation=half_saturation,
    )
    flux = profile.surface_flux(DIFFUSIVITY)
    uptake = profile.total_uptake(max_uptake, half_saturation)
    assert flux == pytest.approx(uptake, rel=1e-4)


def test_the_discretisation_is_second_order():
    """Halving the spacing should cut the error roughly fourfold.

    The analytical solution is the *first-order* limit of Monod uptake, so the
    comparison carries two errors: the discretisation error, which shrinks with
    refinement, and the linearisation error from ``C / (K + C)`` not being
    exactly ``C / K``, which does not. The second must be made negligible or it
    becomes a floor that hides the convergence entirely — at ``K = 1e4`` here,
    ``C / K`` is 2e-5 and the measured order collapses to zero. ``K = 1e8``
    puts it far below the grid error. The upper grid size is capped for the
    same reason at the other end: past ~200 cells the error reaches rounding
    noise.
    """
    half_saturation, max_uptake = 1e8, 5e4
    errors = []
    for cells in (25, 50, 100, 200):
        profile = solve_steady_state(
            Grid1D(THICKNESS, cells),
            diffusivity=DIFFUSIVITY,
            surface=SURFACE,
            max_uptake=max_uptake,
            half_saturation=half_saturation,
        )
        exact = first_order_profile(
            profile.depths,
            thickness=THICKNESS,
            diffusivity=DIFFUSIVITY,
            surface=SURFACE,
            max_uptake=max_uptake,
            half_saturation=half_saturation,
        )
        errors.append(float(np.max(np.abs(profile.concentration - exact))))

    orders = [np.log2(errors[i] / errors[i + 1]) for i in range(len(errors) - 1)]
    assert all(o > 1.8 for o in orders), f"convergence orders {orders} are not second order"


def test_the_profile_decreases_with_depth():
    profile = solve_steady_state(
        Grid1D(THICKNESS, 400),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=0.05,
        half_saturation=0.01,
    )
    assert profile.concentration[0] <= SURFACE
    assert np.all(np.diff(profile.concentration) <= 1e-12)
    assert np.all(profile.concentration >= 0.0)


def test_more_demand_means_less_penetration():
    depths = [
        solve_steady_state(
            Grid1D(THICKNESS, 800),
            diffusivity=DIFFUSIVITY,
            surface=SURFACE,
            max_uptake=rate,
            half_saturation=1e-4,
        ).penetration_depth()
        for rate in (0.02, 0.05, 0.19, 0.5)
    ]
    assert depths == sorted(depths, reverse=True)


# --- input handling ----------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        ({"diffusivity": 0.0}, "diffusivity"),
        ({"surface": -1.0}, "surface"),
        ({"max_uptake": -1.0}, "max_uptake"),
        ({"half_saturation": 0.0}, "half_saturation"),
        ({"tolerance": 0.0}, "tolerance"),
    ],
)
def test_invalid_parameters_are_rejected(kwargs, expected):
    args = {
        "diffusivity": DIFFUSIVITY,
        "surface": SURFACE,
        "max_uptake": 0.05,
        "half_saturation": 0.01,
    } | kwargs
    with pytest.raises(ValueError, match=expected):
        solve_steady_state(Grid1D(THICKNESS, 50), **args)


def test_failure_to_converge_is_reported():
    with pytest.raises(ConvergenceError, match="did not converge"):
        solve_steady_state(
            Grid1D(THICKNESS, 200),
            diffusivity=DIFFUSIVITY,
            surface=SURFACE,
            max_uptake=0.19,
            half_saturation=1e-8,
            max_iterations=1,
        )


def test_penetration_depth_threshold_is_validated():
    profile = solve_steady_state(
        Grid1D(THICKNESS, 50),
        diffusivity=DIFFUSIVITY,
        surface=SURFACE,
        max_uptake=0.05,
        half_saturation=0.01,
    )
    for bad in (0.0, 1.0, -0.5):
        with pytest.raises(ValueError, match="threshold_fraction"):
            profile.penetration_depth(threshold_fraction=bad)
