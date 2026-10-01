"""Surfaces, part one: how cells reach a surface, what it is made of, and how they bind.

The transport correlations are checked against a numerical solution of the
boundary-layer problem they come from, the substratum's material map against
its geometry, and the closed-form adhesion kinetics against an independent
numerical integration. docs/validation.md, "Adhesion to surfaces".
"""

import math

import numpy as np
import pytest

from marse.spatial.colloids import (
    LEVEQUE,
    film_wall_shear_rate_per_s,
    leveque_transfer_um_per_s,
    stokes_einstein_um2_per_s,
)
from marse.spatial.grid import Grid
from marse.spatial.surface import Patch, Substratum
from marse.validation.analytical import adhesion_kinetics

# --- transport to the surface -----------------------------------------------------------


def _wall_flux_by_marching(points=800, height=14.0, ratio=1.02):
    """March y dc/dx = d2c/dy2 from a uniform inlet to x = 1, with D = gamma = 1.

    The concentration is held at zero on the wall and at one far from it. The
    marching is implicit in x with geometrically growing steps, on a grid in y
    crowded towards the wall. It returns the wall gradient at x = 1, which the
    similarity solution puts at the Leveque constant.
    """
    y = height * (np.arange(1, points + 1) / (points + 1)) ** 1.5
    full = np.concatenate(([0.0], y, [height]))
    below, above = np.diff(full)[:-1], np.diff(full)[1:]
    lower = 2.0 / (below * (below + above))
    upper = 2.0 / (above * (below + above))
    c = np.ones(points)
    x, dx = 0.0, 1e-9
    while x < 1.0:
        dx = min(dx, 1.0 - x)
        diagonal = -(lower + upper) - y / dx
        rhs = -y / dx * c
        rhs[-1] -= upper[-1]  # c = 1 beyond the top point
        b, d = diagonal.copy(), rhs.copy()
        for i in range(1, points):  # the Thomas algorithm
            w = lower[i] / b[i - 1]
            b[i] -= w * upper[i - 1]
            d[i] -= w * d[i - 1]
        c = np.empty(points)
        c[-1] = d[-1] / b[-1]
        for i in range(points - 2, -1, -1):
            c[i] = (d[i] - upper[i] * c[i + 1]) / b[i]
        x += dx
        dx *= ratio
    y1, y2, c1, c2 = y[0], y[1], c[0], c[1]
    return (c1 * y2**2 - c2 * y1**2) / (y1 * y2 * (y2 - y1))


def test_the_leveque_constant_is_the_wall_flux_of_the_boundary_layer_problem():
    assert math.isclose(LEVEQUE, 1.0 / (math.gamma(4 / 3) * 9 ** (1 / 3)), rel_tol=1e-15)
    assert abs(LEVEQUE - 0.538366) < 1e-6
    # An independent route: march the convection-diffusion equation itself.
    assert _wall_flux_by_marching() == pytest.approx(LEVEQUE, rel=1e-2)


def test_the_transfer_velocity_scales_as_leveque_says():
    base = leveque_transfer_um_per_s(0.5, 15.0, 1e4)
    # A streptococcus-sized cell at 15 per s, 1 cm downstream: about 0.039 um per s.
    assert base == pytest.approx(LEVEQUE * (0.5**2 * 15.0 / 1e4) ** (1 / 3), rel=1e-15)
    assert base == pytest.approx(0.0388, rel=1e-2)
    assert leveque_transfer_um_per_s(4.0, 15.0, 1e4) == pytest.approx(base * 8 ** (2 / 3))
    assert leveque_transfer_um_per_s(0.5, 120.0, 1e4) == pytest.approx(base * 2)
    assert leveque_transfer_um_per_s(0.5, 15.0, 8e4) == pytest.approx(base / 2)
    with pytest.raises(ValueError, match="shear rate must be positive"):
        leveque_transfer_um_per_s(0.5, 0.0, 1e4)


def test_stokes_einstein_gives_the_known_diffusivity_of_a_micron_sphere():
    # A 1 um sphere in water at 25 C (0.890 mPa s) diffuses at about 0.49 um^2 per s.
    d = stokes_einstein_um2_per_s(1.0, 25.0, 0.890)
    assert d == pytest.approx(0.4907, rel=1e-3)
    assert stokes_einstein_um2_per_s(2.0, 25.0, 0.890) == pytest.approx(d / 2)
    assert stokes_einstein_um2_per_s(1.0, 25.0, 1.780) == pytest.approx(d / 2)
    hotter = stokes_einstein_um2_per_s(1.0, 37.0, 0.890)
    assert hotter == pytest.approx(d * (37 + 273.15) / (25 + 273.15))
    with pytest.raises(ValueError, match="above absolute zero"):
        stokes_einstein_um2_per_s(1.0, -300.0, 0.890)


def test_the_salivary_film_shears_its_surface_at_under_six_per_second():
    # 70 to 100 um thick (Collins and Dawes 1987), 0.8 to 7.6 mm per min (Dawes et al. 1989).
    slowest = film_wall_shear_rate_per_s(800 / 60, 100.0)
    fastest = film_wall_shear_rate_per_s(7600 / 60, 70.0)
    assert slowest == pytest.approx(0.4)
    assert fastest == pytest.approx(5.429, rel=1e-3)


# --- the substratum ----------------------------------------------------------------------


def test_patches_give_every_face_exactly_one_material():
    grid = Grid((8, 4, 6), 5.0)  # 40 x 20 um of surface
    substratum = Substratum(
        grid,
        [
            Patch("enamel", (0.0, 10.0, 0.0, 20.0)),
            Patch("titanium", (10.0, 20.0, 0.0, 20.0)),
            Patch("zirconia", (20.0, 30.0, 0.0, 20.0)),
            Patch("pmma", (30.0, 40.0, 0.0, 20.0)),
        ],
    )
    assert substratum.materials == ("enamel", "titanium", "zirconia", "pmma")
    assert substratum.index.shape == (8, 4)
    assert np.all(substratum.index[:2] == 0)
    assert np.all(substratum.index[6:] == 3)
    areas = [substratum.area_um2(m) for m in substratum.materials]
    assert areas == [200.0] * 4
    assert sum(areas) == pytest.approx(grid.footprint_um2)
    values = substratum.per_face({"enamel": 1.0, "titanium": 2.0, "zirconia": 3.0, "pmma": 4.0})
    assert np.all(values[4:6] == 3.0)
    with pytest.raises(ValueError, match="no value for material"):
        substratum.per_face({"enamel": 1.0})


def test_a_material_may_repeat_across_patches_and_a_column_has_one_face():
    slice_ = Substratum(
        Grid((6, 4), 2.0),
        [Patch("glass", (0.0, 4.0)), Patch("coated", (4.0, 8.0)), Patch("glass", (8.0, 12.0))],
    )
    assert slice_.materials == ("glass", "coated")
    assert slice_.index.tolist() == [0, 0, 1, 1, 0, 0]
    column = Substratum(Grid((10,), 5.0), [Patch("enamel", ())])
    assert column.index.shape == ()
    assert column.area_um2("enamel") == pytest.approx(Grid((10,), 5.0).footprint_um2)


@pytest.mark.parametrize(
    ("patches", "message"),
    [
        ([Patch("a", (0.0, 25.0, 0.0, 20.0)), Patch("b", (20.0, 40.0, 0.0, 20.0))], "overlaps"),
        ([Patch("a", (0.0, 30.0, 0.0, 20.0))], "belong to no patch"),
        ([Patch("a", (0.0, 40.0))], "needs 4 bounds"),
        ([Patch("a", (0.0, 50.0, 0.0, 20.0))], "must satisfy"),
        ([Patch("a", (0.0, 40.0, 0.0, 20.0)), Patch("b", (1.0, 2.0, 1.0, 2.0))], "covers no face"),
        ([], "at least one patch"),
    ],
)
def test_impossible_patterns_are_refused_with_the_reason(patches, message):
    with pytest.raises(ValueError, match=message):
        Substratum(Grid((8, 4, 6), 5.0), patches)


# --- binding kinetics ----------------------------------------------------------------------


def _integrate(arrival, jamming, detachment, locking, times, start=(0.0, 0.0), h=1e-4):
    """Classical Runge-Kutta, sampled at each of ``times`` (increasing)."""

    def f(x):
        binding = arrival * (1.0 - (x[0] + x[1]) / jamming)
        return np.array([binding - (detachment + locking) * x[0], locking * x[0]])

    x, now, samples = np.array(start, dtype=float), 0.0, []
    for target in times:
        for _ in range(round((target - now) / h)):
            k1 = f(x)
            k2 = f(x + h / 2 * k1)
            k3 = f(x + h / 2 * k2)
            k4 = f(x + h * k3)
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        now = target
        samples.append(x.copy())
    return samples


@pytest.mark.parametrize(
    ("detachment", "locking", "start"),
    [
        (3.0, 20.0, (0.0, 0.0)),  # the general case
        (3.0, 0.0, (5.0, 40.0)),  # no locking: Langmuir binding and detachment
        (0.0, 0.12, (0.0, 0.0)),  # a repeated eigenvalue: j0 / n_J equals k_lock
        (50.0, 0.5, (100.0, 200.0)),  # a started surface
    ],
)
def test_the_kinetics_match_an_independent_integration(detachment, locking, start):
    arrival, jamming, times = 120.0, 1000.0, (0.5, 3.0, 12.0)
    samples = _integrate(arrival, jamming, detachment, locking, times, start)
    for span, expected in zip(times, samples, strict=True):
        got = adhesion_kinetics(span, arrival, jamming, detachment, locking, *start)
        np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-9 * jamming)


def test_bound_cells_approach_the_jamming_limit_and_never_pass_it():
    t = np.linspace(0.0, 200.0, 2001)
    reversible, locked = adhesion_kinetics(t, 120.0, 1000.0, 3.0, 20.0)
    total = reversible + locked
    assert np.all(total <= 1000.0 * (1 + 1e-12))
    assert np.all(np.diff(total) >= -1e-9)
    assert locked[-1] == pytest.approx(1000.0, rel=1e-9)
    settled, _ = adhesion_kinetics(1e6, 120.0, 1000.0, 3.0, 0.0)
    assert settled == pytest.approx(120.0 / (120.0 / 1000.0 + 3.0))  # n = j0 / (j0 / n_J + k_off)
