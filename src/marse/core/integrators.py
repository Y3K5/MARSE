"""Positive, conservative time integration for reaction networks.

A reaction network changes as dc/dt = N^T r(c): each process p runs at a rate
r_p >= 0, and its row of N conserves carbon, nitrogen and electrons
(docs/theory.md, section 3.6). The integrator here keeps every one of those
balances to rounding, however large the step. It never produces a negative
concentration, and it does so without clipping, which would create matter.

**The scheme** (docs/theory.md, section 9.6). Heun's method is written, as Shu
and Osher (1988) did, as a convex combination of two forward-Euler steps:

    y1 = E_h(y),    y2 = (y + E_h(y1)) / 2.

Each Euler step E_h is limited process by process. If the step would take more
of a species than it holds, the extents of the processes consuming it are
scaled down, with everything else those processes do, just enough to keep the
species positive. Three properties follow:

- **conservative**: scaling a whole process row cannot change a balance the
  row already conserves;
- **positive**: each limited Euler step is, and y2 is an average of positive
  states;
- **selective**: only the processes that consume a depleted species slow
  down, and every other process runs at its full rate.

A single factor for the whole increment, as in the BBKS schemes (Bruggeman et
al. 2007), is also positive and conservative. But then one depleted species
slows every process in the box, and a species that starts at zero can halt the
run. MARSE limits per process for that reason.

**Accuracy.** Without limiting the method is second order. The Euler result y1
and y2 form an embedded pair whose difference estimates the local error, and
:func:`integrate` adapts the substep to meet a tolerance. The timestep a user
sets is therefore an upper bound, not a promise of accuracy. The procedure is
deterministic: the same inputs give the same substeps and the same result,
bit for bit.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["IntegrationStats", "integrate", "positive_conservative_step"]

type Rates = Callable[[NDArray[np.float64]], NDArray[np.float64]]

_MARGIN = 1.0 - 1e-12
"""A species drawn down to its limit keeps a 1e-12 share of what was available.

The margin sits far above the rounding in the update, so a limited species stays
strictly positive, and far below any tolerance, so it changes no result."""


def _limited_euler(
    y: NDArray[np.float64],
    h: float,
    rates: Rates,
    stoichiometry: NDArray[np.float64],
    shape: tuple[int, ...],
) -> tuple[NDArray[np.float64], bool]:
    """One forward-Euler step, each process's extent limited so no species goes negative.

    ``y`` has shape (species, cells). A species is overdrawn when the step would
    take more of it than it holds. Every process that consumes an overdrawn
    species is scaled by that species' share, and a process is limited by the
    scarcest species it consumes. Afterwards no species loses more than it
    held. A single pass therefore suffices, and every species keeps at least
    1e-12 of its stock.

    What the step produces of a species is deliberately not counted toward what
    may be consumed. Counting it would let a cycle of processes pass
    arbitrarily large amounts through a species within one step. The large
    cancelling flows would then cost conservation its precision. An
    intermediate made and used within a substep still flows through the second
    Heun stage, and the error control sizes the substeps.
    """
    processes = stoichiometry.shape[0]
    rate = np.asarray(rates(y.reshape(shape)), dtype=float).reshape(processes, -1)
    if not np.all(np.isfinite(rate)):
        raise FloatingPointError("a process rate became non-finite")
    if np.any(rate < 0):
        raise ValueError("process rates must not be negative")
    extent = h * rate
    consumed = np.maximum(-stoichiometry, 0.0)  # processes x species
    consumption = consumed.T @ extent  # species x cells
    overdrawn = consumption > y * _MARGIN
    limited = bool(overdrawn.any())
    if limited:
        share = np.where(overdrawn, y * _MARGIN / np.where(overdrawn, consumption, 1.0), 1.0)
        extent = extent * np.where((consumed > 0)[:, :, None], share[None, :, :], 1.0).min(axis=1)
    updated = y + stoichiometry.T @ extent
    if np.any(updated < 0):  # pragma: no cover - the margin keeps rounding above zero
        raise ArithmeticError("a limited step still produced a negative concentration")
    return updated, limited


def _pair(
    y: NDArray[np.float64], h: float, rates: Rates, stoichiometry: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64], bool]:
    """The first-order (Euler) and second-order (Heun) results of one substep."""
    shape = y.shape
    flat = y.reshape(shape[0], -1)
    # Traces far below any tolerance may underflow to zero; that loses nothing.
    with np.errstate(under="ignore"):
        euler, first = _limited_euler(flat, h, rates, stoichiometry, shape)
        ahead, second = _limited_euler(euler, h, rates, stoichiometry, shape)
        heun = 0.5 * flat + 0.5 * ahead
    return euler.reshape(shape), heun.reshape(shape), first or second


def positive_conservative_step(
    y: NDArray[np.float64],
    dt: float,
    rates: Rates,
    stoichiometry: NDArray[np.float64],
    order: int = 2,
) -> tuple[NDArray[np.float64], bool]:
    """One fixed step of size ``dt``: the new state, and whether limiting engaged.

    ``y`` has shape (species,) or (species, *cells); ``rates(y)`` returns
    (processes, *cells) and ``stoichiometry`` is (processes, species).
    ``order=1`` returns the limited Euler step alone.
    """
    if order not in (1, 2):
        raise ValueError("order must be 1 or 2")
    euler, heun, limited = _pair(y, dt, rates, stoichiometry)
    return (euler if order == 1 else heun), limited


@dataclass(frozen=True, slots=True)
class IntegrationStats:
    """How :func:`integrate` got there: substeps taken, retried and limited."""

    accepted: int
    rejected: int
    limited: int
    next_step: float


def integrate(
    y: NDArray[np.float64],
    span: float,
    rates: Rates,
    stoichiometry: NDArray[np.float64],
    *,
    first_step: float,
    relative_tolerance: float,
    absolute_tolerance: float,
    peak: NDArray[np.float64] | None = None,
    guard: Callable[[NDArray[np.float64]], None] | None = None,
) -> tuple[NDArray[np.float64], IntegrationStats]:
    """Advance ``y`` by ``span`` in adaptive substeps, trying ``first_step`` first.

    A substep is accepted when every species' error estimate is within
    ``absolute_tolerance + relative_tolerance * magnitude``. The magnitude is the
    larger of the species' current value and ``peak``, the largest value it has
    reached, which the caller tracks. A species that has run down to a trace is
    then held to the accuracy its own scale deserves. It is not resolved to ever
    smaller absolute errors while it no longer matters.

    ``guard``, if given, sees both results of every substep attempted. It may
    raise to stop the integration as soon as an assumption behind the rates
    fails.
    """
    if not span > 0 or not first_step > 0:
        raise ValueError("span and first_step must be positive")
    reference = np.abs(y) if peak is None else peak
    elapsed, h = 0.0, first_step
    accepted = rejected = limited = 0
    while span - elapsed > 1e-12 * span:
        step = min(h, span - elapsed)
        clipped = step < h  # shortened only to land on the end of the span
        low, high, was_limited = _pair(y, step, rates, stoichiometry)
        if guard is not None:
            guard(low)
            guard(high)
        magnitude = np.maximum(np.maximum(np.abs(y), np.abs(high)), reference)
        scale = absolute_tolerance + relative_tolerance * magnitude
        with np.errstate(under="ignore"):
            error = float(np.max(np.abs(high - low) / scale))
        grow = 0.9 / np.sqrt(error) if error > 0 else 2.0
        if error <= 1.0:
            y = high
            elapsed += step
            accepted += 1
            limited += was_limited
            if not clipped:  # a clipped step says nothing about the step size to try next
                h = step * min(2.0, max(0.2, grow))
        else:
            rejected += 1
            h = step * min(0.5, max(0.1, grow))
            if h <= 1e-14 * span:
                raise ArithmeticError(
                    "the adaptive substep shrank below 1e-14 of the step; the chemistry "
                    "is too stiff for these tolerances"
                )
    return y, IntegrationStats(accepted, rejected, limited, h)
