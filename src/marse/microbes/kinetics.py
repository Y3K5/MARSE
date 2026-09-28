"""Process rates of a reaction network, for one box or many cells at once.

Each process runs at

    r = maximum_per_h * c[proportional_to] * product of its factors,

in mol of the process's reference per m3 per hour (docs/theory.md, section
3.7). The factors are the dimensionless switching functions of
:mod:`marse.microbes.growth`: Monod limitation, non-competitive inhibition and
Haldane substrate inhibition, and the cardinal pH model of
:mod:`marse.microbes.cardinal`, of the pH that the network's charges set
(:mod:`marse.chemistry.acid_base`, docs/theory.md, section 3.8).

Rates are evaluated for concentrations of shape (components, *cells) and
returned with shape (processes, *cells), so the same code serves a well-mixed
box (no cell axes) and a grid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.chemistry.acid_base import ChargeBalance, ph_of_hydrogen
from marse.microbes.cardinal import cardinal_ph, cardinal_ph_slope
from marse.microbes.growth import haldane, monod, noncompetitive_inhibition
from marse.schemas.network import Network

__all__ = ["RateTerms", "compile_rates", "process_rates", "rate_jacobian"]


@dataclass(frozen=True, slots=True)
class _CompiledFactor:
    process: int
    component: int  # -1 for a pH factor
    form: str
    half_saturation: float
    inhibition: float
    cardinal: tuple[float, float, float] = (0.0, 0.0, 0.0)  # ph_min, ph_optimum, ph_max


@dataclass(frozen=True, slots=True)
class RateTerms:
    """A network's rate laws as index arrays, ready for repeated evaluation.

    ``charge_balance`` is set when a rate has a pH factor; the pH is then
    solved once per evaluation, for all of them.
    """

    maximum_per_h: NDArray[np.float64]
    proportional_to: NDArray[np.intp]
    factors: tuple[_CompiledFactor, ...]
    charge_balance: ChargeBalance | None = None


def compile_rates(network: Network) -> RateTerms:
    """Index every process's rate law against the network's component order."""
    index = {name: i for i, name in enumerate(network.component_names)}
    unrated = [p.name for p in network.processes if p.rate is None]
    if unrated:
        raise ValueError(f"processes without a rate cannot run: {', '.join(unrated)}")
    factors = []
    for p, process in enumerate(network.processes):
        assert process.rate is not None  # checked above
        for factor in process.rate.factors:
            if factor.form == "ph":
                assert factor.ph_min is not None
                assert factor.ph_optimum is not None
                assert factor.ph_max is not None
                cardinal = (float(factor.ph_min), float(factor.ph_optimum), float(factor.ph_max))
                factors.append(_CompiledFactor(p, -1, "ph", 0.0, 0.0, cardinal))
                continue
            assert factor.component is not None
            factors.append(
                _CompiledFactor(
                    process=p,
                    component=index[factor.component],
                    form=factor.form,
                    half_saturation=float(factor.half_saturation_mol_per_m3 or 0),
                    inhibition=float(factor.inhibition_mol_per_m3 or 0),
                )
            )
    return RateTerms(
        maximum_per_h=np.array(
            [float(p.rate.maximum_per_h) for p in network.processes if p.rate], dtype=float
        ),
        proportional_to=np.array(
            [index[p.rate.proportional_to] for p in network.processes if p.rate], dtype=np.intp
        ),
        factors=tuple(factors),
        charge_balance=(
            ChargeBalance.of(network) if any(f.form == "ph" for f in factors) else None
        ),
    )


def process_rates(terms: RateTerms, concentrations: NDArray[np.float64]) -> NDArray[np.float64]:
    """Every process's rate, shape (processes, *cells), from concentrations (components, *cells)."""
    cells = concentrations.shape[1:]
    maximum = terms.maximum_per_h.reshape((-1,) + (1,) * len(cells))
    rates = maximum * concentrations[terms.proportional_to]
    ph = None
    for f in terms.factors:
        if f.form == "ph":
            if ph is None:
                assert terms.charge_balance is not None
                ph = terms.charge_balance.ph(concentrations)
            rates[f.process] *= cardinal_ph(ph, *f.cardinal)
            continue
        c = concentrations[f.component]
        if f.form == "monod":
            rates[f.process] *= monod(c, 1.0, f.half_saturation)
        elif f.form == "inhibition":
            rates[f.process] *= noncompetitive_inhibition(c, f.inhibition)
        else:
            rates[f.process] *= haldane(c, 1.0, f.half_saturation, f.inhibition)
    return rates


def _factor(form: str, c: NDArray[np.float64], k: float, ki: float) -> tuple[NDArray, NDArray]:
    """A switching function and its derivative, at concentrations already clamped at zero."""
    if form == "monod":
        return c / (k + c), k / (k + c) ** 2
    if form == "inhibition":
        return ki / (ki + c), -ki / (ki + c) ** 2
    denominator = k + c + c * c / ki
    return c / denominator, (k - c * c / ki) / denominator**2


def rate_jacobian(terms: RateTerms, concentrations: NDArray[np.float64]) -> NDArray[np.float64]:
    """How every rate responds to every component: shape (processes, components, *cells).

    The derivative of k c_a prod f_i(c_j) by the product rule, of the rates as
    the engines evaluate them: at concentrations clamped at zero. A negative
    value, which the rates see as zero, therefore has no slope. Taking the slope
    at zero instead would tell Newton's method that consumption still responds
    there, and it would creep towards a negative stage value by a few percent
    per iteration. An implicit step (docs/theory.md, section 9.8) linearises the
    reactions with it.

    A pH factor responds to every charged component at once, through the
    hydrogen ion concentration h: d gamma / d c_k = gamma'(pH) (dpH/dh)
    (dh/dc_k), with dh/dc_k from the implicit function theorem
    (:meth:`~marse.chemistry.acid_base.ChargeBalance.hydrogen_slopes`).
    """
    c = np.maximum(concentrations, 0.0)
    live = concentrations >= 0  # slopes of the clamp: one above zero, none below
    processes = terms.maximum_per_h.size
    cells = c.shape[1:]
    jacobian = np.zeros((processes, c.shape[0], *cells))
    by_process: list[list[tuple[int, NDArray, NDArray]]] = [[] for _ in range(processes)]
    through_h = None  # d pH / d c_k for every component, once per evaluation
    for f in terms.factors:
        if f.form == "ph":
            if through_h is None:
                assert terms.charge_balance is not None
                h = terms.charge_balance.hydrogen(c)
                ph = ph_of_hydrogen(h)
                slopes = terms.charge_balance.hydrogen_slopes(concentrations, h)
                through_h = -slopes / (h * math.log(10.0))
            value = np.asarray(cardinal_ph(ph, *f.cardinal))
            slope = np.asarray(cardinal_ph_slope(ph, *f.cardinal))
            by_process[f.process].append((-1, value, slope))
            continue
        value, slope = _factor(f.form, c[f.component], f.half_saturation, f.inhibition)
        by_process[f.process].append((f.component, value, slope))
    for p in range(processes):
        k, a = terms.maximum_per_h[p], terms.proportional_to[p]
        factors = by_process[p]
        product = np.ones(cells)
        for _, value, _ in factors:
            product = product * value
        jacobian[p, a] += k * product * live[a]
        for i, (component, _, slope) in enumerate(factors):
            others = np.ones(cells)
            for m, (_, value, _) in enumerate(factors):
                if m != i:
                    others = others * value
            if component < 0:  # the pH factor, through every charged component
                assert through_h is not None
                jacobian[p] += (k * c[a] * slope * others)[None] * through_h
            else:
                jacobian[p, component] += k * c[a] * slope * others * live[component]
    return jacobian
