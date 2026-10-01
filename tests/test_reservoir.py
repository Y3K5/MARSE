"""A well-mixed pool bordering the box: the engine that carries the mouth's saliva.

The exchange between a film and its pool is checked against its closed form,
the bordered linear system against the full operator, and the flux-form
update against the pool's and the box's ledgers.
"""

import math
from itertools import pairwise

import numpy as np
import pytest

from marse.core.reservoir import ReservoirPath, ReservoirTransport
from marse.spatial.grid import Grid
from marse.spatial.transport import Diffusion


def _no_reactions(components):
    return (
        np.zeros((0, components)),
        lambda c: np.zeros((0, *c.shape[1:])),
        lambda c: np.zeros((0, components, *c.shape[1:])),
    )


def _steady(thickness, span_h, secreted, start=0.0):
    return ReservoirPath(
        start_h=start,
        times_h=np.array([0.0, span_h]),
        thickness_um=np.array([thickness, thickness]),
        growth_um_per_h=np.zeros(2),
        secreted_mol_per_m3=np.asarray(secreted, dtype=float),
    )


# --- the path ----------------------------------------------------------------------------


def test_the_path_follows_a_cubic_exactly_and_books_what_enters_as_it_grows():
    times = np.linspace(0.0, 0.1, 5)

    def cubic(t):
        return 700.0 + 3e3 * t + 2e4 * t**2 - 1e5 * t**3

    def slope(t):
        return 3e3 + 4e4 * t - 3e5 * t**2

    path = ReservoirPath(
        start_h=2.0,
        times_h=times,
        thickness_um=cubic(times),
        growth_um_per_h=slope(times),
        secreted_mol_per_m3=np.array([5.0, 0.0]),
        drink_um_per_h=100.0,
        drink_mol_per_m3=np.array([0.0, 300.0]),
        supply_per_h=np.array([0.0, 2.0]),
    )
    for t in np.linspace(0.0, 0.1, 17):
        value, rate = path.at(2.0 + t)
        assert value == pytest.approx(cubic(t), rel=1e-13)
        assert rate == pytest.approx(slope(t), rel=1e-12)
    grown = cubic(0.1) - cubic(0.0)
    np.testing.assert_allclose(
        path.added(2.0, 2.1), [5.0 * (grown - 10.0), 300.0 * 10.0 + 0.2], rtol=1e-13
    )
    pieces = sum(path.added(2.0 + a, 2.0 + b) for a, b in pairwise(times))
    np.testing.assert_allclose(pieces, path.added(2.0, 2.1), rtol=1e-13)


# --- the exchange ------------------------------------------------------------------------


def test_a_film_and_its_pool_relax_as_the_closed_form_says():
    """A film of 100 um under a pool of 700 um: c - m decays at k (1 + delta / H).

    Two voxels of 50 um, mixed by a diffusivity fast enough to keep them equal.
    """
    stoichiometry, rates, jacobian = _no_reactions(1)
    grid = Grid((2,), 50.0)
    engine = ReservoirTransport(
        Diffusion(grid, np.array([1e12]), np.zeros(1), closed_top=True),
        stoichiometry,
        rates,
        jacobian,
        exchanged=[0],
        reference_um=700.0,
    )
    k, delta, thickness, c0, m0, span = 30.0, 100.0, 700.0, 10.0, 2.0, 0.02
    engine.exchange_per_h = np.array([k, k])
    engine.path = _steady(thickness, span, [0.0])
    y = engine.pack(np.array([[c0, c0]]), np.array([m0]))
    y, imports, _ = engine.integrate(
        y, span, relative_tolerance=1e-7, absolute_tolerance=1e-12, first_step=1e-5
    )
    box, pool = engine.unpack(y)
    mean = (c0 * delta + m0 * thickness) / (delta + thickness)
    decay = math.exp(-k * (1 + delta / thickness) * span)
    assert box[0, 0] == pytest.approx(mean + (c0 - mean) * decay, rel=1e-6)
    assert pool[0] == pytest.approx(mean + (m0 - mean) * decay, rel=1e-6)
    assert box[0, 0] * delta + pool[0] * thickness == pytest.approx(
        c0 * delta + m0 * thickness, rel=1e-14
    )
    assert box[0, 1] == pytest.approx(box[0, 0], rel=1e-9)
    assert imports[0] * 50.0 == pytest.approx((box[0].mean() - c0) * delta, rel=1e-12)


def _reacting_column(film_voxels=4, dims=1):
    """A column of 10 voxels: a catalyst turns A into B below a film renewed from the pool."""
    shape = (3,) * (dims - 1) + (10,)
    grid = Grid(shape, 10.0)
    k1 = 40.0

    def rates(c):
        return (k1 * np.maximum(c[2], 0.0) * np.maximum(c[0], 0.0))[None]

    def jacobian(c):
        jac = np.zeros((1, 3, *c.shape[1:]))
        jac[0, 0] = k1 * np.maximum(c[2], 0.0) * (c[0] >= 0)
        jac[0, 2] = k1 * np.maximum(c[0], 0.0) * (c[2] >= 0)
        return jac

    engine = ReservoirTransport(
        Diffusion(grid, np.array([1e6, 5e5, 0.0]), np.zeros(3), closed_top=True),
        np.array([[-1.0, 1.0, 0.0]]),
        rates,
        jacobian,
        exchanged=[0, 1],
        reference_um=700.0,
    )
    profile = np.zeros(shape)
    profile[..., -film_voxels:] = 50.0
    engine.exchange_per_h = profile
    box = np.zeros((3, *shape))
    box[2, ..., :6] = 1.0
    return engine, box


@pytest.mark.parametrize("dims", [1, 2])
def test_box_and_pool_conserve_what_enters_to_rounding(dims):
    engine, box = _reacting_column(dims=dims)
    y = engine.pack(box, np.array([3.0, 0.0, 0.0]))
    start = engine._per_area(box) + engine.unpack(y)[1] * engine.reference_um
    entered = np.zeros(3)
    t, step = 0.0, None
    for span_s, (h0, h1) in [
        (20.0, (700.0, 850.0)),
        (30.0, (850.0, 1000.0)),
        (40.0, (650.0, 900.0)),
    ]:
        span = span_s / 3600.0
        growth = (h1 - h0) / span
        engine.path = ReservoirPath(
            start_h=t,
            times_h=np.array([0.0, span]),
            thickness_um=np.array([h0, h1]),
            growth_um_per_h=np.array([growth, growth]),
            secreted_mol_per_m3=np.array([5.0, 0.0, 0.0]),
        )
        entered += engine.path.added(t, t + span)
        y, _, stats = engine.integrate(  # the balance holds to rounding at any tolerance
            y, span, relative_tolerance=1e-3, absolute_tolerance=1e-8, first_step=step, start_h=t
        )
        step, t = stats.next_step, t + span
        # a swallow between spans takes the pool's liquid, so it is booked here as leaving
        box, pool = engine.unpack(y)
        if h1 > 950.0:
            kept = 650.0 / h1
            entered -= pool * engine.reference_um * (1.0 - kept)
            y = engine.pack(box, pool * kept)
    box, pool = engine.unpack(y)
    now = engine._per_area(box) + pool * engine.reference_um
    carbon = now[0] + now[1]  # A becomes B
    assert carbon == pytest.approx(start[0] + start[1] + entered[0] + entered[1], rel=1e-13)
    assert now[2] == pytest.approx(start[2], rel=1e-14)  # the catalyst stays in the box
    assert np.all(box >= 0.0)
    assert np.all(pool >= 0.0)
    assert box[1].max() > 0.0  # B was made, and some reached the pool
    assert pool[1] > 0.0


@pytest.mark.parametrize("dims", [1, 2])
def test_the_bordered_system_is_the_jacobian_and_is_solved_exactly(dims):
    engine, box = _reacting_column(dims=dims)
    rng = np.random.default_rng(dims)
    box = box + rng.uniform(0.5, 2.0, size=box.shape) * (box == 0) * np.array([1, 1, 0]).reshape(
        (3,) + (1,) * dims
    )
    y = engine.pack(box, np.array([3.0, 0.5, 0.0]))
    engine.path = _steady(800.0, 0.01, [5.0, 0.0, 0.0])
    a = 0.002
    system = engine._system(y, a)
    x = rng.normal(size=y.shape)
    x[2] = 0.0  # the catalyst does not move; its pool entry is inert
    step = 1e-6
    linear = (engine.evaluate(y + step * x).rate - engine.evaluate(y - step * x).rate) / (2 * step)
    np.testing.assert_allclose(system.apply(x), x - a * linear, rtol=1e-6, atol=1e-8)
    b = rng.normal(size=y.shape)
    solved, _ = system.solve(b, scale=np.ones_like(b), tolerance=1e-12)
    np.testing.assert_allclose(system.apply(solved), b, rtol=0, atol=1e-8)


def test_a_pool_that_would_give_more_than_it_holds_is_scaled_and_still_conserves():
    engine, box = _reacting_column()
    engine.path = _steady(700.0, 0.01, [0.0, 0.0, 0.0])
    pool = np.array([1e-6, 0.0, 0.0])
    count, first = 2, engine._reacting
    extents = np.zeros((engine.stoichiometry.shape[0], 10))
    extents[first, -4:] = 1.0  # asks the pool for far more A than it holds
    transfers = (np.zeros((3, 10)),)
    new_box, new_pool, extents, _, _ = engine._conserve(
        box, pool, np.zeros(3), transfers, extents, None
    )
    assert new_pool[0] >= 0.0
    given = engine._per_area(extents[first : first + count])[0]
    assert given <= pool[0] * engine.reference_um
    assert engine._per_area(new_box)[0] + new_pool[0] * engine.reference_um == pytest.approx(
        pool[0] * engine.reference_um, rel=1e-12
    )


def test_a_box_open_at_the_top_is_refused():
    stoichiometry, rates, jacobian = _no_reactions(1)
    with pytest.raises(ValueError, match="close its top face"):
        ReservoirTransport(
            Diffusion(Grid((4,), 10.0), np.array([1e6]), np.zeros(1)),
            stoichiometry,
            rates,
            jacobian,
            exchanged=[0],
            reference_um=700.0,
        )
