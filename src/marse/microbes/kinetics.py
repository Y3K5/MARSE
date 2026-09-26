"""Process rates of a reaction network, for one box or many cells at once.

Each process runs at

    r = maximum_per_h * c[proportional_to] * product of its factors,

in mol of the process's reference per m3 per hour (docs/theory.md, section
3.7). The factors are the dimensionless switching functions of
:mod:`marse.microbes.growth`: Monod limitation, non-competitive inhibition and
Haldane substrate inhibition.

Rates are evaluated for concentrations of shape (components, *cells) and
returned with shape (processes, *cells), so the same code serves a well-mixed
box (no cell axes) and a grid.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.microbes.growth import haldane, monod, noncompetitive_inhibition
from marse.schemas.network import Network

__all__ = ["RateTerms", "compile_rates", "process_rates"]


@dataclass(frozen=True, slots=True)
class _CompiledFactor:
    process: int
    component: int
    form: str
    half_saturation: float
    inhibition: float


@dataclass(frozen=True, slots=True)
class RateTerms:
    """A network's rate laws as index arrays, ready for repeated evaluation."""

    maximum_per_h: NDArray[np.float64]
    proportional_to: NDArray[np.intp]
    factors: tuple[_CompiledFactor, ...]


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
    )


def process_rates(terms: RateTerms, concentrations: NDArray[np.float64]) -> NDArray[np.float64]:
    """Every process's rate, shape (processes, *cells), from concentrations (components, *cells)."""
    cells = concentrations.shape[1:]
    maximum = terms.maximum_per_h.reshape((-1,) + (1,) * len(cells))
    rates = maximum * concentrations[terms.proportional_to]
    for f in terms.factors:
        c = concentrations[f.component]
        if f.form == "monod":
            rates[f.process] *= monod(c, 1.0, f.half_saturation)
        elif f.form == "inhibition":
            rates[f.process] *= noncompetitive_inhibition(c, f.inhibition)
        else:
            rates[f.process] *= haldane(c, 1.0, f.half_saturation, f.inhibition)
    return rates
