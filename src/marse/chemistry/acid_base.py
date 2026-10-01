"""pH from electroneutrality, over acid-base totals and ions of fixed charge.

A network sets a pH when one of its components is an acid-base total, such as
lactic acid and lactate together (docs/networks.md, "Acids, bases and pH").
Protons are never components. In every voxel the concentration of hydrogen
ions, h, is the one that makes the liquid electrically neutral:

    F(h) = h - Kw / h + sum_k T_k zbar_k(h) + sum_s z_s c_s = 0,

where T_k is a total, zbar_k(h) the mean charge of its protonation states at
h, and c_s an ion of fixed charge z_s (docs/theory.md, section 3.8). A total
with pKa values pK_1 < ... < pK_n, written in its most protonated form of
charge z_0, is a fraction

    alpha_j(h) = (K_1 ... K_j / h^j) / sum_i (K_1 ... K_i / h^i)

in the form that has lost j protons, so zbar = z_0 - sum_j j alpha_j. Every
zbar rises with h, and so does F, from minus infinity to plus infinity: the
root exists and is unique. It is found by Newton's method in log h, kept
inside a bracket that shrinks with every evaluation, so it always converges.
A voxel is done when its Newton step falls below 1e-13 in log h, or when the
balance is as small as rounding can make it, a few units in the last place of
its largest charge terms. Near the root the second test is the one that
counts: there a step of 1e-13 can be below what rounding lets the balance
resolve, and Newton's method would step between two neighbouring values of h
for ever.

Rates that depend on pH need how h responds to each component. The implicit
function theorem gives it exactly: dh/dc_k = -(dF/dc_k) / (dF/dh), with
dF/dc_k the charge a component carries at h, and dF/dh = 1 + Kw/h^2 + sum_k
T_k var_k / h, var_k being the variance of the number of protons lost.
Concentrations below zero count as zero, as they do in the rates, and have no
slope.

Concentrations are in mol/m3, which is mmol/L. pKa and pKw are the
conventional logarithms of constants in mol/L.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.schemas.network import Network

__all__ = ["DEFAULT_PKW", "ChargeBalance", "ph_of_hydrogen"]

DEFAULT_PKW = 14.0
"""Water's ion product at 25 C; at 37 C it is about 13.6."""

_LN10 = math.log(10.0)
_LN_MOL_PER_M3 = math.log(1000.0)  # a constant in mol/L, expressed in mol/m3
_LOWEST = math.log(1e-14)  # pH 17, in mol/m3
_HIGHEST = math.log(1e4)  # pH -1
_NEUTRAL = math.log(1e-4)  # pH 7, where every solve starts
_ITERATIONS = 200
_CONVERGED = 1e-13  # in log h: h to about thirteen digits
_FLOOR = 16 * float(np.finfo(float).eps)  # the balance's rounding, per unit of its charge terms

type Field = NDArray[np.float64]


def ph_of_hydrogen(h: Field) -> Field:
    """pH from the hydrogen ion concentration in mol/m3."""
    return -np.log10(np.asarray(h) / 1000.0)


@dataclass(frozen=True, slots=True)
class _Total:
    component: int
    log_weight: NDArray[np.float64]  # ln(K_1 ... K_j) in mol/m3 units, j = 0..n
    top_charge: float  # the charge of the most protonated form


@dataclass(frozen=True, slots=True)
class ChargeBalance:
    """Every charge a network's components carry, and water's ion product."""

    totals: tuple[_Total, ...]
    ions: tuple[tuple[int, float], ...]  # component index and its fixed charge
    kw: float  # (mol/m3)^2
    components: int

    @classmethod
    def of(cls, network: Network) -> ChargeBalance:
        """The charge balance of a network's components."""
        totals, ions = [], []
        for index, component in enumerate(network.components):
            charge = float(component.formula.charge)
            if component.pka:
                ka = [-float(k) * _LN10 + _LN_MOL_PER_M3 for k in component.pka]
                weight = np.concatenate(([0.0], np.cumsum(ka)))
                totals.append(_Total(index, weight, charge))
            elif charge:
                ions.append((index, charge))
        pkw = DEFAULT_PKW if network.pkw is None else float(network.pkw)
        return cls(tuple(totals), tuple(ions), 10.0 ** (-pkw) * 1e6, len(network.components))

    # -- the balance and its slopes -------------------------------------------------------

    def _speciation(self, total: _Total, log_h: Field) -> tuple[Field, Field]:
        """The mean charge of a total at h, and the variance of the protons it has lost."""
        n = total.log_weight.size
        j = np.arange(n, dtype=float).reshape((n,) + (1,) * log_h.ndim)
        exponent = total.log_weight.reshape(j.shape) - j * log_h
        exponent = exponent - exponent.max(axis=0)
        weight = np.exp(exponent)
        alpha = weight / weight.sum(axis=0)
        lost = (j * alpha).sum(axis=0)
        variance = np.maximum((j * j * alpha).sum(axis=0) - lost * lost, 0.0)
        return total.top_charge - lost, variance

    def _balance(self, c: Field, h: Field) -> tuple[Field, Field, Field]:
        """F(h), dF/dh, and the sum of the charge terms' sizes, which bounds F's rounding."""
        log_h = np.log(h)
        water = self.kw / h
        net = h - water
        slope = 1.0 + water / h
        size = h + water
        for total in self.totals:
            charge, variance = self._speciation(total, log_h)
            amount = c[total.component]
            net = net + amount * charge
            slope = slope + amount * variance / h
            size = size + amount * np.abs(charge)
        for index, charge in self.ions:
            net = net + charge * c[index]
            size = size + abs(charge) * c[index]
        return net, slope, size

    def residual(self, c: Field, h: Field) -> tuple[Field, Field]:
        """F(h), the net charge per volume, and dF/dh, for concentrations c (clamped at zero)."""
        c = np.maximum(np.asarray(c, dtype=float), 0.0)
        net, slope, _ = self._balance(c, np.asarray(h, dtype=float))
        return net, slope

    def hydrogen(self, c: Field) -> Field:
        """h in mol/m3 wherever c is given: the root of the charge balance.

        Newton's method in log h from pH 7, inside a bracket [pH -1, pH 17]
        that every evaluation narrows. A step that would leave the bracket, or
        that is not at most half the step before it, bisects the bracket
        instead, as in the safeguarded Newton method of Press et al. (2007,
        section 9.4): Newton's method can otherwise cycle, each step landing on
        the other end of an unchanging bracket. A voxel stops where the balance
        is at its rounding floor, and the solve ends when every voxel has
        stopped or moves by no more than 1e-13 in log h.
        """
        c = np.maximum(np.asarray(c, dtype=float), 0.0)
        shape = c.shape[1:]
        low = np.full(shape, _LOWEST)
        high = np.full(shape, _HIGHEST)
        x = np.full(shape, _NEUTRAL)
        last = np.full(shape, _HIGHEST - _LOWEST)  # the step before, which a step must halve
        for _ in range(_ITERATIONS):
            h = np.exp(x)
            net, slope, size = self._balance(c, h)
            floor = np.abs(net) <= _FLOOR * size
            low = np.where(net < 0, x, low)
            high = np.where(net > 0, x, high)
            new = x - net / (slope * h)  # dF/d(log h) = h dF/dh
            step = np.abs(new - x)
            # Strictly outside: at convergence a step can round onto a bound, and a step of
            # rounding size need not halve the one before.
            outside = (new < low) | (new > high) | ~np.isfinite(new)
            slow = (step > 0.5 * last) & (step > _CONVERGED)
            new = np.where(outside | slow, 0.5 * (low + high), new)
            new = np.where(floor, x, new)
            done = bool(np.all(floor | (np.abs(new - x) <= _CONVERGED)))
            last = np.abs(new - x)
            x = new
            if done:
                return np.exp(x)
        raise ArithmeticError(  # pragma: no cover - the bracket makes this unreachable
            f"the charge balance did not converge in {_ITERATIONS} iterations"
        )

    def ph(self, c: Field) -> Field:
        return ph_of_hydrogen(self.hydrogen(c))

    def dissociated(self, component: int, h: Field) -> tuple[Field, Field]:
        """The protons each unit of an acid-base total has lost at h, and its slope in h.

        For a monoprotic acid, Ka / (Ka + h). The slope is minus the variance of
        the protons lost over h, as for the balance's own slope.
        """
        total = next((t for t in self.totals if t.component == component), None)
        if total is None:
            raise ValueError(f"component {component} is not an acid-base total")
        h = np.asarray(h, dtype=float)
        charge, variance = self._speciation(total, np.log(h))
        return total.top_charge - charge, -variance / h

    def hydrogen_slopes(self, c: Field, h: Field) -> Field:
        """dh/dc for every component, shape (components, *cells): zero below zero, as in rates."""
        c = np.asarray(c, dtype=float)
        h = np.asarray(h, dtype=float)
        _, slope, _ = self._balance(np.maximum(c, 0.0), h)
        log_h = np.log(h)
        slopes = np.zeros(c.shape)
        live = c >= 0
        for total in self.totals:
            charge, _ = self._speciation(total, log_h)
            slopes[total.component] = np.where(live[total.component], -charge / slope, 0.0)
        for index, charge in self.ions:
            slopes[index] = np.where(live[index], -charge / slope, 0.0)
        return slopes
