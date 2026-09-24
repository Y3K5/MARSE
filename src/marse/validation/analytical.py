"""Closed-form solutions that serve as ground truth for validation cases.

Each function is exact for an idealized problem named in docs/validation.md.
The test suite checks every one against an independent numerical solution, so
the references can be trusted before any simulation is compared with them.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from marse._numeric import FloatOrArray, as_array, require, unwrap

__all__ = [
    "batch_final_biomass",
    "chemostat_break_even",
    "exponential_growth",
    "logistic_growth",
    "monod_batch_time",
    "point_source_diffusion_2d",
    "zero_order_penetration_depth",
]


def exponential_growth(t: ArrayLike, x0: float, mu: float) -> FloatOrArray:
    """Unrestricted growth, x(t) = x0 exp(mu t). Case V1."""
    return unwrap(x0 * np.exp(mu * as_array(t)))


def logistic_growth(t: ArrayLike, x0: float, x_max: float, mu: float) -> FloatOrArray:
    """Growth towards a carrying capacity, dx/dt = mu x (1 - x / x_max). Case V1."""
    require(0 < x0 <= x_max, "need 0 < x0 <= x_max")
    return unwrap(x_max / (1.0 + (x_max / x0 - 1.0) * np.exp(-mu * as_array(t))))


def batch_final_biomass(s0: float, x0: float, yield_coefficient: float) -> float:
    """Biomass once a batch culture has used all its substrate, x0 + Y s0. Case V2."""
    return x0 + yield_coefficient * s0


def monod_batch_time(
    substrate: ArrayLike,
    s0: float,
    x0: float,
    mu_max: float,
    k_s: float,
    yield_coefficient: float,
) -> FloatOrArray:
    """Time at which a Monod batch culture has drawn its substrate down to ``substrate``.

    Integrates dX/dt = mu_max S X / (K_s + S) with dS/dt = -(1 / Y) dX/dt, whose mass
    balance X + Y S = x0 + Y s0 gives the implicit solution

        mu_max t = (K_s / C) ln(s0 / S) + (1 + K_s / C) ln(X / x0),

    with C = s0 + x0 / Y and X = x0 + Y (s0 - S). No maintenance, no death. Case V2.
    """
    s = as_array(substrate)
    require(np.all((s > 0) & (s <= s0)), "substrate must lie in (0, s0]")
    require(min(x0, mu_max, k_s, yield_coefficient) > 0, "parameters must be positive")
    total = s0 + x0 / yield_coefficient
    biomass = x0 + yield_coefficient * (s0 - s)
    t = (k_s / total) * np.log(s0 / s) + (1.0 + k_s / total) * np.log(biomass / x0)
    return unwrap(t / mu_max)


def zero_order_penetration_depth(
    diffusivity: ArrayLike, boundary_concentration: ArrayLike, uptake_rate: ArrayLike
) -> FloatOrArray:
    """Depth at which a solute consumed at a constant volumetric rate runs out.

    For D c'' = k in a slab with c = c0 at the surface, the concentration falls
    parabolically to zero at depth sqrt(2 D c0 / k). This is the classic estimate
    of how far oxygen penetrates a biofilm. Case V3.
    """
    d, c0, k = as_array(diffusivity), as_array(boundary_concentration), as_array(uptake_rate)
    require(np.all(d > 0) and np.all(c0 >= 0) and np.all(k > 0), "invalid transport values")
    return unwrap(np.sqrt(2.0 * d * c0 / k))


def point_source_diffusion_2d(
    r: ArrayLike, t: float, amount: float, diffusivity: float
) -> FloatOrArray:
    """Concentration a distance r from a point release on a plane after time t > 0.

    c(r, t) = M / (4 pi D t) exp(-r^2 / (4 D t)); the reference solution for the
    diffusion-only case. Case V3.
    """
    require(t > 0 and diffusivity > 0, "t and diffusivity must be positive")
    spread = 4.0 * diffusivity * t
    return unwrap(amount / (np.pi * spread) * np.exp(-(as_array(r) ** 2) / spread))


def chemostat_break_even(dilution_rate: float, mu_max: float, k_s: float) -> float:
    """Substrate level at which a Monod species' growth exactly balances washout.

    lambda = K_s D / (mu_max - D), or infinity when mu_max <= D (the species washes
    out). With one limiting substrate, the species with the lowest lambda excludes
    the others (Hsu, Hubbell and Waltman 1977). Case V5.
    """
    require(dilution_rate > 0 and k_s > 0 and mu_max >= 0, "invalid chemostat parameters")
    if mu_max <= dilution_rate:
        return float("inf")
    return k_s * dilution_rate / (mu_max - dilution_rate)
