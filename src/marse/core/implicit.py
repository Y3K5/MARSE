"""Implicit, conservative, positive integration of reactions and transport together.

On a grid, concentrations change as

    dc/dt = divergence of the face fluxes + N^T r(c),

diffusion (:mod:`marse.spatial.transport`) and reactions
(:mod:`marse.microbes.kinetics`) at once. Diffusion at physical diffusivities is
stiff: an explicit step on 2 um voxels would have to be a quarter of a
millisecond (docs/theory.md, section 9.1). Splitting reactions from transport
fails differently: at practical steps it starves a biofilm, because each step
supplies one refill of the box instead of a continuous flux. A quasi-steady
treatment fails in exactly the transients MARSE is built to study. So
everything is integrated together, implicitly (docs/theory.md, section 9.8):

- **The scheme** is the two-stage SDIRK method of Alexander (1977), with
  gamma = 1 - 1/sqrt(2):

      Y1 = y + gamma h f(Y1),
      Y2 = y + h ((1 - gamma) f(Y1) + gamma f(Y2)).

  It is L-stable and stiffly accurate, and second order: long steps land on
  the quasi-steady state, and short steps follow a transient.
- **Conservation for any solver tolerance.** The new state is not taken from
  the stage solution. It is rebuilt in flux form, from the face transfers and
  process extents the two stages imply,

      y+ = y + transfers + N^T extents,

  and each transfer leaves one voxel and enters another, while each process
  row conserves carbon, nitrogen and electrons. So the balance holds to
  rounding whatever the Newton or multigrid tolerance.
- **Positivity by limiting.** No scheme above first order is positive by
  itself at every step size. Where the update would leave a value negative,
  the transfers out of that voxel and the processes consuming there are
  scaled down, each as a whole, net of what arrives. A single-pass limiter that
  counts only the stock is the guaranteed fallback. Nothing is clipped. The
  update sums what arrives in each voxel and what leaves it separately, so a
  voxel that loses nothing can only gain, even in rounding. Traces below the
  smallest normal number, where no relative margin survives rounding, stop
  giving instead of being scaled.
- **Newton's method, projected.** A value well above the absolute tolerance
  may fall at most to a tenth of itself per iteration, which keeps Newton out
  of the kink of the rate laws at zero while a large transient is resolved. A
  value already below the tolerance may go negative, because the stage
  solution itself can undershoot there. A Newton iteration that does not
  converge rejects the step, which is retried shorter.
- **Error control.** The difference between the second-order result and a
  first-order one, h gamma (f(Y2) - f(Y1)), filtered through the stage matrix
  so that stiff components do not dominate it, sets the substep. It is
  deterministic, so a run replays bit for bit. A run's first step comes from
  the rates and their change over a trial Euler step (Hairer, Norsett and
  Wanner 1993, section II.4): colonies placed in fresh liquid start far from
  their quasi-steady state, and a step that tried to cross that in one go
  would only fail.

The linear systems are solved by :mod:`marse.spatial.multigrid`. One matrix,
with the Jacobian at the start of the step, serves both stages, and it is
rebuilt only when Newton converges slowly.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.spatial.multigrid import ImplicitSystem
from marse.spatial.transport import Diffusion, divergence

__all__ = ["GAMMA", "NewtonFailure", "ReactionTransport", "StepStats"]

GAMMA = 1.0 - 1.0 / math.sqrt(2.0)
_MARGIN = 1.0 - 1e-12
_TINY = float(np.finfo(float).tiny)  # below it, rounding swallows any relative margin
_NEWTON_TOLERANCE = 1e-3  # corrections, in units of the error tolerance
_NEWTON_ITERATIONS = 20
_LINEAR_TOLERANCE = 1e-4

type Field = NDArray[np.float64]
type Transfers = tuple[Field, ...]


def _magnitude(*fields: Field) -> Field:
    """Each component's largest absolute value anywhere in the box, shaped to broadcast.

    Errors are measured against a component's own scale across the box, not
    against the value in each voxel: a trace of lactate seeping into the liquid
    far from the colonies is held to the accuracy lactate needs where it
    matters, not resolved ever more finely where it does not.
    """
    largest = np.max([np.abs(f).reshape(f.shape[0], -1).max(axis=1) for f in fields], axis=0)
    return largest.reshape((-1,) + (1,) * (fields[0].ndim - 1))


class NewtonFailure(ArithmeticError):
    """Newton's method did not converge within its iterations; the step is retried shorter."""


@dataclass(frozen=True, slots=True)
class StepStats:
    """How :meth:`ReactionTransport.integrate` got there."""

    accepted: int
    rejected: int
    limited: int
    newton_failures: int
    newton_iterations: int
    next_step: float


@dataclass(frozen=True, slots=True)
class _Stage:
    rate: Field
    fluxes: Transfers
    reactions: Field


class ReactionTransport:
    """A reaction network in a box of voxels, with the bulk liquid held above it."""

    def __init__(
        self,
        diffusion: Diffusion,
        stoichiometry: NDArray[np.float64],
        rates: Callable[[Field], Field],
        jacobian: Callable[[Field], Field],
    ) -> None:
        self.diffusion = diffusion
        self.stoichiometry = np.asarray(stoichiometry, dtype=float)  # processes x components
        self.consumed = np.maximum(-self.stoichiometry, 0.0)
        self.produced = np.maximum(self.stoichiometry, 0.0)
        self.rates = rates
        self.jacobian = jacobian
        self.shape = diffusion.grid.shape
        self.spacing = diffusion.spacing

    # -- the right-hand side --------------------------------------------------------------

    def _react(self, extents: Field) -> Field:
        return np.einsum("pj,p...->j...", self.stoichiometry, extents)

    def evaluate(self, c: Field) -> _Stage:
        fluxes = self.diffusion.fluxes(c)
        reactions = self.rates(c)
        return _Stage(divergence(fluxes, self.spacing) + self._react(reactions), fluxes, reactions)

    def _system(self, c: Field, a: float) -> ImplicitSystem:
        blocks = np.einsum("pj,pk...->jk...", self.stoichiometry, self.jacobian(c))
        return ImplicitSystem(
            self.shape, self.spacing, self.diffusion.diffusivity_um2_per_h, a, blocks
        )

    # -- one stage: Y = base + a f(Y) ----------------------------------------------------

    def _stage(
        self,
        base: Field,
        a: float,
        guess: Field,
        system: ImplicitSystem,
        reference: Field,
        atol: float,
        rtol: float,
    ) -> tuple[Field, ImplicitSystem, int]:
        y = guess.copy()
        previous = math.inf
        for iteration in range(1, _NEWTON_ITERATIONS + 1):
            residual = y - base - a * self.evaluate(y).rate
            scale = atol + rtol * np.maximum(reference, _magnitude(y))
            correction, _ = system.solve(-residual, scale=scale, tolerance=_LINEAR_TOLERANCE)
            size = float(np.max(np.abs(correction) / scale))
            if size < _NEWTON_TOLERANCE:
                return y + correction, system, iteration
            if size > 0.5 * previous:  # slow: the Jacobian has gone stale
                system = self._system(y, a)
            previous = size
            trial = y + correction
            y = np.where((y > atol) & (trial < 0.1 * y), 0.1 * y, trial)
        raise NewtonFailure(
            f"Newton's method did not converge in {_NEWTON_ITERATIONS} iterations "
            f"(last correction {size:.3g} of the tolerance)"
        )

    # -- the conservative, positive update ------------------------------------------------

    def _flows(self, transfers: Transfers, extents: Field) -> tuple[Field, Field]:
        """What arrives in each voxel and what leaves it, per component: both non-negative.

        Arrivals are inflow through the faces and production; departures are
        outflow and consumption. Each transfer departs from one voxel and
        arrives in its neighbour, so the two conserve exactly what the
        divergence does. Kept apart, they make the update positive by
        construction wherever nothing leaves: a voxel that loses nothing can
        only gain, even in rounding.
        """
        into = np.einsum("pj,p...->j...", self.produced, extents)
        out = np.einsum("pj,p...->j...", self.consumed, extents)
        dims = len(self.shape)
        for axis, upper in enumerate(transfers):
            h = self.spacing[axis]
            forward = np.maximum(upper, 0.0) / h  # leaves through the upper face
            backward = np.maximum(-upper, 0.0) / h  # enters through the upper face
            out += forward
            into += backward
            if axis < dims - 1:  # periodic: the upper face of i is the lower face of i + 1
                into += np.roll(forward, 1, axis=axis + 1)
                out += np.roll(backward, 1, axis=axis + 1)
            else:  # the substratum below the first voxel carries nothing
                into[..., 1:] += forward[..., :-1]
                out[..., 1:] += backward[..., :-1]
        return into, out

    def _apply(self, y: Field, transfers: Transfers, extents: Field) -> tuple[Field, Field, Field]:
        """The new state, with what arrived and what left."""
        into, out = self._flows(transfers, extents)
        return (y - out) + into, into, out

    def _scale(self, transfers: Transfers, extents: Field, share: Field) -> tuple[Transfers, Field]:
        """Scale what leaves voxel i by share[i]: transfers by their source, processes whole."""
        dims = len(self.shape)
        scaled = []
        for axis, upper in enumerate(transfers):
            if axis < dims - 1:
                source_below = np.roll(share, -1, axis=axis + 1)
            else:  # the bulk above the top face is never short
                source_below = np.concatenate(
                    (share[..., 1:], np.ones_like(share[..., :1])), axis=-1
                )
            scaled.append(np.where(upper > 0, upper * share, upper * source_below))
        limit = np.where(
            (self.consumed > 0).reshape(self.consumed.shape + (1,) * dims),
            share[None],
            1.0,
        ).min(axis=1)
        return tuple(scaled), extents * limit

    def _limited_update(
        self, y: Field, transfers: Transfers, extents: Field
    ) -> tuple[Field, Transfers, Field, int]:
        new, into, out = self._apply(y, transfers, extents)
        if np.all(new >= 0):
            return new, transfers, extents, 0
        for rounds in range(1, 51):
            # A voxel is short only if more leaves than it holds and receives, so out > 0.
            short = new < 0
            available = y + into
            share = np.where(short, available * _MARGIN / np.where(short, out, 1.0), 1.0)
            # Below the smallest normal number no relative margin survives rounding:
            # scaled by 4/6, a transfer of one unit in the last place rounds back up
            # to one. A voxel holding that little stops giving instead.
            share = np.where(short & (available < _TINY), 0.0, share)
            transfers, extents = self._scale(transfers, extents, np.minimum(share, 1.0))
            new, into, out = self._apply(y, transfers, extents)
            if np.all(new >= 0):
                return new, transfers, extents, rounds
        over = out > y * _MARGIN  # the stock-only limiter
        share = np.where(over, y * _MARGIN / np.where(over, out, 1.0), 1.0)
        while True:
            transfers, extents = self._scale(transfers, extents, share)
            new, _, _ = self._apply(y, transfers, extents)
            if np.all(new >= 0):
                return new, transfers, extents, 51
            # Whatever rounding leaves below zero stops giving: a voxel that gives
            # nothing can only gain, so each pass settles more voxels for good.
            share = np.where(new < 0, 0.0, 1.0)

    # -- one step and an adaptive span ----------------------------------------------------

    def step(
        self, y: Field, h: float, reference: Field, atol: float, rtol: float
    ) -> tuple[Field, NDArray[np.float64], Field, int, int]:
        """One SDIRK2 step: new state, imports, error estimate, limiter rounds, Newton steps."""
        # Traces far below any tolerance may underflow to zero; that loses nothing.
        with np.errstate(under="ignore"):
            return self._step(y, h, reference, atol, rtol)

    def _step(
        self, y: Field, h: float, reference: Field, atol: float, rtol: float
    ) -> tuple[Field, NDArray[np.float64], Field, int, int]:
        g = GAMMA
        reference = _magnitude(reference)
        system = self._system(y, g * h)
        stage1, system, n1 = self._stage(y, g * h, y, system, reference, atol, rtol)
        s1 = self.evaluate(stage1)
        base = y + (1 - g) * h * s1.rate
        stage2, system, n2 = self._stage(base, g * h, stage1, system, reference, atol, rtol)
        s2 = self.evaluate(stage2)
        transfers = tuple(
            h * ((1 - g) * f1 + g * f2) for f1, f2 in zip(s1.fluxes, s2.fluxes, strict=True)
        )
        extents = h * ((1 - g) * s1.reactions + g * s2.reactions)
        new, transfers, extents, rounds = self._limited_update(y, transfers, extents)
        scale = atol + rtol * np.maximum(reference, _magnitude(y, new))
        estimate, _ = system.solve(g * h * (s2.rate - s1.rate), scale=scale, tolerance=1e-2)
        top = -transfers[-1][..., -1]  # into the box through each top face
        imports = top.reshape(top.shape[0], -1).sum(axis=1) / self.spacing[-1]
        return new, imports, estimate, rounds, n1 + n2

    def starting_step(self, y: Field, reference: Field, atol: float, rtol: float) -> float:
        """A first step from the rates and their change over a trial Euler step.

        The algorithm of Hairer, Norsett and Wanner (1993, section II.4) for a
        second-order method, in the error norm of :meth:`integrate`: the Euler
        increment stays a hundredth of the state, and the local error, from the
        change in the rates, a hundredth of the tolerance.
        """
        with np.errstate(under="ignore"):
            return self._starting_step(y, reference, atol, rtol)

    def _starting_step(self, y: Field, reference: Field, atol: float, rtol: float) -> float:
        scale = atol + rtol * np.maximum(_magnitude(reference), _magnitude(y))
        rate = self.evaluate(y).rate
        size = float(np.max(np.abs(y) / scale))
        speed = float(np.max(np.abs(rate) / scale))
        trial = 1e-6 if size < 1e-5 or speed < 1e-5 else 0.01 * size / speed
        change = self.evaluate(y + trial * rate).rate - rate
        curvature = float(np.max(np.abs(change) / scale)) / trial
        if max(speed, curvature) <= 1e-15:
            return max(1e-6, trial * 1e-3)
        return min(100 * trial, (0.01 / max(speed, curvature)) ** (1 / 3))

    def integrate(
        self,
        y: Field,
        span: float,
        *,
        relative_tolerance: float,
        absolute_tolerance: float,
        first_step: float | None = None,
        peak: Field | None = None,
    ) -> tuple[Field, NDArray[np.float64], StepStats]:
        """Advance ``y`` by ``span`` in adaptive steps; imports are summed over the span.

        Without ``first_step``, the first step comes from :meth:`starting_step`.
        Imports are in the units the ledger sums: concentration summed over
        voxels, per component.
        """
        reference = _magnitude(y if peak is None else peak)
        if first_step is None:
            first_step = self.starting_step(y, reference, absolute_tolerance, relative_tolerance)
        if not span > 0 or not first_step > 0:
            raise ValueError("span and first_step must be positive")
        elapsed, h = 0.0, first_step
        imports = np.zeros(y.shape[0])
        accepted = rejected = limited = failures = newton = 0
        while span - elapsed > 1e-12 * span:
            step = min(h, span - elapsed)
            clipped = step < h
            try:
                new, entered, estimate, rounds, iterations = self.step(
                    y, step, reference, absolute_tolerance, relative_tolerance
                )
            except NewtonFailure:
                failures += 1
                rejected += 1
                h = step * 0.25
                if h <= 1e-14 * span:
                    raise
                continue
            newton += iterations
            scale = absolute_tolerance + relative_tolerance * np.maximum(
                reference, _magnitude(y, new)
            )
            with np.errstate(under="ignore"):
                error = float(np.max(np.abs(estimate) / scale))
            grow = 0.9 / math.sqrt(error) if error > 0 else 2.0
            if error <= 1.0:
                y = new
                imports += entered
                elapsed += step
                accepted += 1
                limited += rounds > 0
                if not clipped:
                    h = step * min(2.0, max(0.2, grow))
            else:
                rejected += 1
                h = step * min(0.5, max(0.1, grow))
                if h <= 1e-14 * span:
                    raise ArithmeticError(
                        "the adaptive step shrank below 1e-14 of the span; the system is too "
                        "stiff for these tolerances"
                    )
        return y, imports, StepStats(accepted, rejected, limited, failures, newton, h)
