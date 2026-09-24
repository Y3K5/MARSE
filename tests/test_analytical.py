"""Checks of the analytical reference solutions against independent numerical ones.

A reference that is wrong would silently validate a wrong simulation, so each
closed form is compared with a direct numerical solution of the same problem.
"""

import numpy as np
import pytest

from marse.microbes.growth import monod
from marse.validation.analytical import (
    batch_final_biomass,
    chemostat_break_even,
    exponential_growth,
    logistic_growth,
    monod_batch_time,
    point_source_diffusion_2d,
    zero_order_penetration_depth,
)

pytestmark = pytest.mark.numerical


def rk4(rhs, y0, t_end, steps):
    """Classic fourth-order Runge-Kutta for an autonomous system; returns (t, y)."""
    h = t_end / steps
    y = np.asarray(y0, dtype=float)
    ys = [y]
    for _ in range(steps):
        k1 = rhs(y)
        k2 = rhs(y + h / 2 * k1)
        k3 = rhs(y + h / 2 * k2)
        k4 = rhs(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ys.append(y)
    return np.linspace(0.0, t_end, steps + 1), np.array(ys)


def solve_tridiagonal(lower, diag, upper, rhs):
    """Thomas algorithm for a tridiagonal system (lower[0] and upper[-1] unused)."""
    n = len(diag)
    c, d = np.zeros(n), np.zeros(n)
    c[0], d[0] = upper[0] / diag[0], rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / m
    x = np.zeros(n)
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def test_exponential_and_logistic_growth_match_numerical_integration():
    t, y = rk4(lambda x: 0.7 * x, [0.01], 10.0, 2000)
    np.testing.assert_allclose(exponential_growth(t, 0.01, 0.7), y[:, 0], rtol=1e-9)
    t, y = rk4(lambda x: 1.3 * x * (1 - x / 4.0), [0.01], 15.0, 3000)
    np.testing.assert_allclose(logistic_growth(t, 0.01, 4.0, 1.3), y[:, 0], rtol=1e-9)


@pytest.mark.parametrize("k_s", [1e-3, 0.5, 5.0])
def test_integrated_monod_batch_solution(k_s):
    mu_max, yield_coefficient, s0, x0 = 1.0, 0.5, 5.0, 0.01

    def rhs(state):
        s, x = state
        growth = monod(s, mu_max, k_s) * x
        return np.array([-growth / yield_coefficient, growth])

    t, states = rk4(rhs, [s0, x0], 25.0, 12_500)
    s, x = states[:, 0], states[:, 1]
    # Compare over 99% of the substrate; in the last 1% the fixed-step RK4 reference
    # itself loses accuracy once S nears K_s and the dynamics speed up.
    usable = s > 1e-2 * s0
    predicted = monod_batch_time(s[usable], s0, x0, mu_max, k_s, yield_coefficient)
    np.testing.assert_allclose(predicted, t[usable], rtol=1e-7, atol=1e-9)
    # Mass balance holds throughout: biomass made = yield x substrate used.
    np.testing.assert_allclose(x + yield_coefficient * s, x0 + yield_coefficient * s0, rtol=1e-12)
    assert abs(s[-1]) < 1e-3 * s0
    assert x[-1] == pytest.approx(batch_final_biomass(s0, x0, yield_coefficient), rel=1e-3)


def test_zero_order_penetration_depth_matches_reaction_diffusion_steady_state():
    # D c'' = k c / (K + c) with K << c0: essentially zero-order uptake until c runs out.
    diffusivity, c0, k, half_saturation = 1.0, 1.0, 2.0, 1e-4
    depth = zero_order_penetration_depth(diffusivity, c0, k)  # = 1.0
    length, n = 2.0, 4000
    dx = length / n
    x = np.linspace(dx, length, n)
    c = np.maximum(c0 * (1 - x / depth), 0.0) ** 2 + 1e-9  # start near the answer
    off = diffusivity / dx**2
    for _ in range(100):  # Newton iterations on the discretized equations
        left = np.concatenate([[c0], c[:-1]])
        right = np.concatenate([c[1:], [c[-2]]])  # zero-flux end
        residual = off * (left - 2 * c + right) - k * c / (half_saturation + c)
        slope = k * half_saturation / (half_saturation + c) ** 2
        lower = np.full(n, off)
        lower[-1] = 2 * off
        step = solve_tridiagonal(lower, -2 * off - slope, np.full(n, off), -residual)
        c = np.maximum(c + step, 1e-14)
        if np.max(np.abs(step)) < 1e-12:
            break
    else:
        pytest.fail("Newton iteration did not converge")
    # All uptake is fed by the surface flux, so flux / k is the depth that is active.
    surface_flux = diffusivity * (3 * c0 - 4 * c[0] + c[1]) / (2 * dx)
    assert surface_flux / k == pytest.approx(depth, rel=5e-3)
    expected = np.where(x < depth, c0 * (1 - x / depth) ** 2, 0.0)
    assert np.max(np.abs(c - expected)) < 0.01


def test_point_source_conserves_mass_and_solves_the_diffusion_equation():
    diffusivity, amount, t, dt = 600.0, 1.0, 2.0, 1e-3  # e.g. glucose in um^2/s
    h = 1.0
    grid = np.arange(-400.0, 400.0 + h, h)
    xx, yy = np.meshgrid(grid, grid)
    r = np.hypot(xx, yy)
    c = point_source_diffusion_2d(r, t, amount, diffusivity)
    assert c.sum() * h * h == pytest.approx(amount, rel=1e-6)
    dcdt = (
        point_source_diffusion_2d(r, t + dt, amount, diffusivity)
        - point_source_diffusion_2d(r, t - dt, amount, diffusivity)
    ) / (2 * dt)
    laplacian = (c[:-2, 1:-1] + c[2:, 1:-1] + c[1:-1, :-2] + c[1:-1, 2:] - 4 * c[1:-1, 1:-1]) / h**2
    scale = np.max(np.abs(dcdt))
    assert np.max(np.abs(dcdt[1:-1, 1:-1] - diffusivity * laplacian)) < 1e-3 * scale


@pytest.mark.parametrize(("dilution", "winner"), [(0.3, "gleaner"), (0.58, "fast_grower")])
def test_chemostat_competition_is_won_by_the_lowest_break_even_level(dilution, winner):
    species = {"fast_grower": (1.0, 0.5), "gleaner": (0.6, 0.05)}  # (mu_max, K_s)
    s_in, yield_coefficient = 10.0, 0.5
    (mu_a, k_a), (mu_b, k_b) = species["fast_grower"], species["gleaner"]

    def rhs(state):
        s, a, b = state
        growth_a, growth_b = monod(s, mu_a, k_a) * a, monod(s, mu_b, k_b) * b
        return np.array(
            [
                dilution * (s_in - s) - (growth_a + growth_b) / yield_coefficient,
                growth_a - dilution * a,
                growth_b - dilution * b,
            ]
        )

    _, states = rk4(rhs, [s_in, 0.01, 0.01], 400.0, 20_000)
    s, fast, gleaner = states[-1]
    levels = {name: chemostat_break_even(dilution, *p) for name, p in species.items()}
    assert min(levels, key=levels.get) == winner
    survivor, loser = (gleaner, fast) if winner == "gleaner" else (fast, gleaner)
    assert survivor > 100 * loser
    assert s == pytest.approx(levels[winner], rel=1e-3)
    assert survivor == pytest.approx(yield_coefficient * (s_in - levels[winner]), rel=1e-3)


def test_break_even_is_infinite_when_the_species_washes_out():
    assert chemostat_break_even(0.5, mu_max=0.4, k_s=1.0) == float("inf")
