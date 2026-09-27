"""Closed-form solutions that serve as ground truth for validation cases.

Each function is exact for an idealized problem named in docs/validation.md.
The test suite checks every one against an independent numerical solution, so
the references can be trusted before any simulation is compared with them.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from marse._numeric import FloatOrArray, as_array, require, unwrap

__all__ = [
    "adhesion_kinetics",
    "batch_final_biomass",
    "chemostat_break_even",
    "cosine_mode_rate",
    "exponential_growth",
    "first_order_profile",
    "first_order_surface_flux",
    "logistic_growth",
    "monod_batch_time",
    "point_source_diffusion_2d",
    "point_source_diffusion_3d",
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


def _first_order_length(diffusivity: float, max_uptake: float, half_saturation: float) -> float:
    """Reaction-diffusion length lambda = sqrt(D K / k_max) for first-order uptake."""
    require(diffusivity > 0, "diffusivity must be positive")
    require(max_uptake > 0, "max_uptake must be positive")
    require(half_saturation > 0, "half_saturation must be positive")
    return float(np.sqrt(diffusivity * half_saturation / max_uptake))


def first_order_profile(
    depth: ArrayLike,
    *,
    thickness: float,
    diffusivity: float,
    surface: float,
    max_uptake: float,
    half_saturation: float,
) -> FloatOrArray:
    """Steady profile of a solute consumed in proportion to its concentration.

    The opposite limit to :func:`zero_order_penetration_depth`. When the solute
    is everywhere far below the half-saturation constant, Monod uptake becomes
    first order, ``R = (k_max / K) C``, and the steady state solves
    ``D C'' = (k_max / K) C``. With ``C = surface`` at the top and no flux at
    the base of a slab of the given ``thickness``,

        C(z) = C0 cosh((L - z) / lambda) / cosh(L / lambda),
        lambda = sqrt(D K / k_max).

    Unlike the zero-order case the solute never runs out completely: it decays
    towards the base rather than reaching a penetration depth. Together the two
    limits bracket the behaviour of a real Monod uptake profile, which is why
    MARSE checks its solver against both. Case V3.
    """
    require(thickness > 0, "thickness must be positive")
    require(surface >= 0, "surface concentration must be non-negative")
    z = as_array(depth)
    require(np.all((z >= 0) & (z <= thickness)), "depth must lie within the slab")
    scale = _first_order_length(diffusivity, max_uptake, half_saturation)
    # cosh(a)/cosh(b) written so that a thick slab (b >> 1) cannot overflow.
    a, b = (thickness - z) / scale, thickness / scale
    ratio = np.exp(a - b) * (1.0 + np.exp(-2.0 * a)) / (1.0 + np.exp(-2.0 * b))
    return unwrap(surface * ratio)


def first_order_surface_flux(
    *,
    thickness: float,
    diffusivity: float,
    surface: float,
    max_uptake: float,
    half_saturation: float,
) -> float:
    """Flux into the slab for :func:`first_order_profile`, D C0 tanh(L/lambda) / lambda.

    At steady state this must equal the total uptake integrated over the slab,
    which is the strongest single check on a numerical solver: it ties the
    boundary condition to the interior chemistry.
    """
    require(thickness > 0, "thickness must be positive")
    scale = _first_order_length(diffusivity, max_uptake, half_saturation)
    return float(diffusivity * surface * np.tanh(thickness / scale) / scale)


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


def point_source_diffusion_3d(
    r: ArrayLike, t: float, amount: float, diffusivity: float
) -> FloatOrArray:
    """Concentration a distance r from a point release in unbounded space after time t > 0.

    c(r, t) = M / (4 pi D t)^(3/2) exp(-r^2 / (4 D t)). Around a point, three
    dimensions dilute a release as 1/r^3 of its spread where two dimensions
    dilute it as 1/r^2, which is why a two-dimensional model overstates how far
    a metabolite reaches. The reference for the three-dimensional transport
    solver (docs/validation.md).
    """
    require(t > 0 and diffusivity > 0, "t and diffusivity must be positive")
    spread = 4.0 * diffusivity * t
    return unwrap(amount / (np.pi * spread) ** 1.5 * np.exp(-(as_array(r) ** 2) / spread))


def cosine_mode_rate(
    wavenumbers: ArrayLike, diffusivity: float, spacing: float | None = None
) -> float:
    """Decay rate of a product of cosines under diffusion, per hour if D is per hour.

    A field cos(k_x x) cos(k_y y) cos(k_z z) decays as exp(rate t). In the
    continuum the rate is -D (k_x^2 + k_y^2 + k_z^2). On a grid of spacing h,
    with the field sampled at voxel centres, the finite-volume operator of
    :mod:`marse.spatial.transport` has the same field as an exact eigenvector,
    with rate -D sum (4 / h^2) sin^2(k h / 2). That holds for periodic lateral
    axes (k = 2 pi m / L), and in height, with no flux through the substratum
    and zero at the top face, for k = (2 n + 1) pi / (2 L). The two rates differ
    at second order in h.
    """
    k = as_array(wavenumbers)
    require(diffusivity >= 0, "diffusivity must not be negative")
    if spacing is None:
        return float(-diffusivity * np.sum(k**2))
    require(spacing > 0, "spacing must be positive")
    return float(-diffusivity * np.sum(4.0 / spacing**2 * np.sin(k * spacing / 2.0) ** 2))


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


def adhesion_kinetics(
    t: ArrayLike,
    arrival: float,
    jamming: float,
    detachment: float,
    locking: float,
    reversible0: float = 0.0,
    locked0: float = 0.0,
) -> tuple[FloatOrArray, FloatOrArray]:
    """Cells reversibly bound, R, and locked, L, per unit area of a surface over time.

    Cells arrive and bind at ``arrival`` per area and time on a bare surface,
    in proportion to the free area: bound cells block binding linearly, and
    ``jamming`` cells per area block it completely. Bound cells detach at
    ``detachment`` or lock at ``locking``, per unit time:

        dR/dt = j0 (1 - (R + L) / n_J) - (k_off + k_lock) R,
        dL/dt = k_lock R.

    The system is linear as long as R + L < n_J, and the dynamics never leave
    that region. With locking, every cell ends up locked and the surface
    jammed, (R, L) -> (0, n_J), along the two decaying modes of the matrix.
    Without it, L stays where it started and R settles to the Langmuir balance
    between binding and detachment. The reference for the adhesion exchange
    of the spatial engine (docs/validation.md).
    """
    time = as_array(t)
    require(arrival > 0 and jamming > 0, "arrival and jamming must be positive")
    require(detachment >= 0 and locking >= 0, "rates must not be negative")
    require(
        reversible0 >= 0 and locked0 >= 0 and reversible0 + locked0 <= jamming,
        "the starting cells must fit on the surface",
    )
    a = arrival / jamming
    if locking == 0:
        rate = a + detachment
        settled = (arrival - a * locked0) / rate
        reversible = settled + (reversible0 - settled) * np.exp(-rate * time)
        return unwrap(reversible), unwrap(np.full_like(time, locked0))
    matrix = np.array([[-(a + detachment + locking), -a], [locking, 0.0]])
    fixed = np.array([0.0, jamming])
    offset = np.array([reversible0, locked0]) - fixed
    half_trace = 0.5 * np.trace(matrix)
    # The discriminant (tr/2)^2 - det, written as a sum of non-negative terms so
    # that rounding cannot push a repeated eigenvalue's zero below zero.
    spread = math.sqrt(
        (0.5 * (a - locking)) ** 2 + 0.25 * detachment * (detachment + 2.0 * (a + locking))
    )
    identity = np.eye(2)
    if spread <= 1e-12 * abs(half_trace):  # repeated eigenvalue
        lam = half_trace
        decay = np.exp(lam * time)[..., None, None]
        propagator = decay * (identity + time[..., None, None] * (matrix - lam * identity))
    else:
        fast, slow = half_trace - spread, half_trace + spread
        e_fast = np.exp(fast * time)[..., None, None]
        e_slow = np.exp(slow * time)[..., None, None]
        propagator = (e_slow * (matrix - fast * identity) - e_fast * (matrix - slow * identity)) / (
            slow - fast
        )
    state = fixed + propagator @ offset
    return unwrap(state[..., 0]), unwrap(state[..., 1])
