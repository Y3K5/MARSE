"""The positive, conservative integrator: order, positivity, balance and selectivity.

Each property in docs/theory.md, section 9.6, is checked on problems small
enough to reason about by hand, and then on random networks.
"""

import numpy as np
import pytest

from marse.core.integrators import integrate, positive_conservative_step

# Linear decay A -> B at rate k [A]: one process, two species.
DECAY = np.array([[-1.0, 1.0]])


def decay_rates(k):
    return lambda c: np.array([k * c[0]])


def test_an_unlimited_euler_step_is_plain_explicit_euler():
    y, limited = positive_conservative_step(np.array([2.0, 0.0]), 0.1, decay_rates(3.0), DECAY, 1)
    assert not limited
    np.testing.assert_allclose(y, [2.0 * (1 - 0.3), 2.0 * 0.3], rtol=1e-15)


def test_a_step_far_too_large_is_limited_and_stays_positive():
    # Explicit Euler would take A to 2 (1 - 100) = -198. The limiter keeps it
    # positive, conserves A + B, and moves (almost) all of A.
    y, limited = positive_conservative_step(np.array([2.0, 0.0]), 100.0, decay_rates(1.0), DECAY, 1)
    assert limited
    assert y[0] > 0
    assert y[0] < 1e-10
    assert y.sum() == pytest.approx(2.0, rel=1e-15)


@pytest.mark.numerical
def test_heun_is_second_order_on_linear_decay():
    k, span = 1.3, 2.0
    exact = 2.0 * np.exp(-k * span)
    errors = []
    for n in (20, 40, 80, 160):
        y = np.array([2.0, 0.0])
        for _ in range(n):
            y, _ = positive_conservative_step(y, span / n, decay_rates(k), DECAY)
        errors.append(abs(y[0] - exact))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.9), orders


# An intermediate that starts at zero and is used up fast:
# S -> L at r1 = X S/(0.5+S); L -> P at r2 = 50 X L/(0.01+L). Species S, L, X, P.
CHAIN = np.array([[-1.0, 1.0, 0.0, 0.0], [0.0, -1.0, 0.0, 1.0]])


def chain_rates(c):
    return np.array([c[2] * c[0] / (0.5 + c[0]), 50.0 * c[2] * c[1] / (0.01 + c[1])])


def test_an_intermediate_starting_at_zero_does_not_stall_the_step():
    start = np.array([5.0, 0.0, 1.0, 0.0])
    y, _ = positive_conservative_step(start, 0.5, chain_rates, CHAIN)
    assert np.all(y >= 0)
    assert y[3] > 0.2  # product formed: the intermediate kept flowing
    assert y[[0, 1, 3]].sum() == pytest.approx(5.0, rel=1e-14)


def test_limiting_slows_only_the_processes_that_consume_the_scarce_species():
    # Two independent pairs: A -> B, fast and overdrawn; C -> D, slow. The C -> D
    # result must be exactly what it is without the fast pair.
    both = np.array([[-1.0, 1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 1.0]])

    def rates(c):
        return np.array([1000.0 * c[0], 0.2 * c[2]])

    y, limited = positive_conservative_step(np.array([1.0, 0.0, 3.0, 0.0]), 0.5, rates, both)
    alone, _ = positive_conservative_step(
        np.array([3.0, 0.0]), 0.5, lambda c: np.array([0.2 * c[0]]), DECAY
    )
    assert limited
    np.testing.assert_array_equal(y[2:], alone)


def test_cells_are_independent_columns():
    columns = np.array([[5.0, 1.0], [0.0, 0.2], [1.0, 3.0], [0.0, 0.5]])
    together, _ = positive_conservative_step(columns, 0.05, chain_rates, CHAIN)
    for j in range(2):
        alone, _ = positive_conservative_step(columns[:, j], 0.05, chain_rates, CHAIN)
        np.testing.assert_array_equal(together[:, j], alone)


def test_negative_rates_are_refused():
    with pytest.raises(ValueError, match="must not be negative"):
        positive_conservative_step(np.array([1.0, 1.0]), 0.1, lambda c: np.array([-1.0]), DECAY)


def reference(start, span, rates, stoichiometry, n=8000):
    """Classical fourth-order Runge-Kutta, far more accurate than the tolerances tested."""

    def slope(y):
        return stoichiometry.T @ rates(y)

    y, h = start.copy(), span / n
    for _ in range(n):
        k1 = slope(y)
        k2 = slope(y + h / 2 * k1)
        k3 = slope(y + h / 2 * k2)
        k4 = slope(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


@pytest.mark.numerical
def test_adaptive_substeps_meet_the_tolerance():
    start = np.array([5.0, 0.0, 1.0, 0.0])
    exact = reference(start, 2.0, chain_rates, CHAIN)
    for tolerance in (1e-3, 1e-5):
        y, stats = integrate(
            start,
            2.0,
            chain_rates,
            CHAIN,
            first_step=0.1,
            relative_tolerance=tolerance,
            absolute_tolerance=1e-9,
        )
        assert np.max(np.abs(y - exact)) < 20 * tolerance * np.max(start)
        assert stats.accepted > 0


def test_adaptive_integration_is_deterministic():
    start = np.array([5.0, 0.0, 1.0, 0.0])
    kwargs = {"first_step": 0.1, "relative_tolerance": 1e-6, "absolute_tolerance": 1e-9}
    first, _ = integrate(start, 0.5, chain_rates, CHAIN, **kwargs)
    second, _ = integrate(start, 0.5, chain_rates, CHAIN, **kwargs)
    np.testing.assert_array_equal(first, second)


def random_system(seed):
    """Random processes over six species, each row conserving the column sums."""
    rng = np.random.default_rng(seed)
    processes = []
    for _ in range(8):
        row = np.zeros(6)
        reactants = rng.choice(6, size=2, replace=False)
        product = rng.choice([i for i in range(6) if i not in reactants])
        amounts = rng.uniform(0.2, 2.0, size=2)
        row[reactants] = -amounts
        row[product] = amounts.sum()  # conserves the plain sum of the species
        processes.append(row)
    stoichiometry = np.array(processes)
    speed = rng.uniform(0.1, 50.0, size=8)

    def rates(c):
        consumed = stoichiometry < 0
        limit = np.where(consumed, c[None, :] / (0.1 + c[None, :]), 1.0).prod(axis=1)
        return speed * limit

    return stoichiometry, rates, rng.uniform(0.0, 3.0, size=6)


@pytest.mark.parametrize("seed", range(6))
def test_huge_steps_on_random_systems_stay_positive_and_conserve(seed):
    stoichiometry, rates, y = random_system(seed)
    total = y.sum()
    for _ in range(20):
        y, _ = positive_conservative_step(y, 1e3, rates, stoichiometry)
        assert np.all(y >= 0)
    assert y.sum() == pytest.approx(total, rel=1e-13)


def random_cyclic_system(seed):
    """Random processes over 3 to 8 species, including cycles and exact zeros."""
    rng = np.random.default_rng(seed)
    species, processes = int(rng.integers(3, 9)), int(rng.integers(2, 12))
    rows = []
    for _ in range(processes):
        row = np.zeros(species)
        reactants = rng.choice(species, size=int(rng.integers(1, 3)), replace=False)
        others = [i for i in range(species) if i not in reactants]
        products = rng.choice(
            others, size=int(rng.integers(1, min(2, len(others)) + 1)), replace=False
        )
        amounts = rng.uniform(0.1, 3.0, size=len(reactants))
        row[reactants] = -amounts
        row[products] += amounts.sum() * rng.dirichlet(np.ones(len(products)))
        rows.append(row)
    stoichiometry = np.array(rows)
    speed = rng.uniform(0.01, 1e3, size=processes)
    cells = int(rng.integers(1, 4))
    start = rng.uniform(0, 2, size=(species, cells)) * (rng.random((species, cells)) > 0.3)

    def rates(c):
        limit = np.where((stoichiometry < 0)[:, :, None], c[None] / (0.05 + c[None]), 1.0)
        return speed[:, None] * limit.prod(axis=1)

    return stoichiometry, rates, start


def test_three_thousand_random_cyclic_systems_stay_positive_and_conserve():
    worst = 0.0
    for seed in range(3000):
        stoichiometry, rates, start = random_cyclic_system(seed)
        total = start.sum(axis=0)
        for dt in (1e-3, 1.0, 1e3, 1e6):
            y, _ = positive_conservative_step(start, dt, rates, stoichiometry)
            assert np.all(y >= 0), (seed, dt)
            drift = np.abs(y.sum(axis=0) - total) / np.maximum(total, 1e-300)
            worst = max(worst, float(drift.max()))
    assert worst < 1e-14
