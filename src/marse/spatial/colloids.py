"""How cells in a suspension reach a surface: the transport half of adhesion.

A cell suspended in the liquid reaches the substratum by Brownian diffusion,
carried towards the wall by the flow past it. Close to a wall any flow is a
simple shear, u = gamma y, and for particles much smaller than the depleted
layer the steady flux onto a surface that captures what arrives, a distance x
downstream of where the capture begins, is the Leveque (1928) solution. In
colloid deposition it is called the Smoluchowski-Levich approximation:

    j = 0.538 c (D^2 gamma / x)^(1/3),    0.538 = 1 / (Gamma(4/3) 9^(1/3)),

so a cell arrives at the transfer velocity k = j / c (docs/theory.md, section
6.4). Cells are Brownian spheres, with the Stokes-Einstein diffusivity

    D = k_B T / (3 pi eta d).

A flow chamber sets gamma directly. A thin film of liquid flowing over a
surface with its top free, such as the salivary film on a tooth, has a
half-parabolic profile, and its wall shear rate is three times its mean
velocity over its thickness.

Sedimentation is left out: the scenes this serves have surfaces that are
vertical or face down (docs/environments.md).
"""

from __future__ import annotations

import math

from marse._numeric import require

__all__ = [
    "BOLTZMANN_J_PER_K",
    "LEVEQUE",
    "film_wall_shear_rate_per_s",
    "leveque_transfer_um_per_s",
    "stokes_einstein_um2_per_s",
]

BOLTZMANN_J_PER_K = 1.380649e-23  # exact in the SI since 2019
LEVEQUE = 1.0 / (math.gamma(4.0 / 3.0) * 9.0 ** (1.0 / 3.0))
"""The Leveque constant, 0.5384: the wall flux of the similarity solution."""

_ZERO_CELSIUS_K = 273.15


def stokes_einstein_um2_per_s(
    diameter_um: float, temperature_c: float, viscosity_mpa_s: float
) -> float:
    """The Brownian diffusivity of a sphere, in um^2 per s: k_B T / (3 pi eta d)."""
    require(diameter_um > 0, "diameter_um must be positive")
    require(viscosity_mpa_s > 0, "viscosity_mpa_s must be positive")
    kelvin = temperature_c + _ZERO_CELSIUS_K
    require(kelvin > 0, "temperature_c must be above absolute zero")
    viscosity = viscosity_mpa_s * 1e-3  # Pa s
    diameter = diameter_um * 1e-6  # m
    return BOLTZMANN_J_PER_K * kelvin / (3.0 * math.pi * viscosity * diameter) * 1e12


def leveque_transfer_um_per_s(
    diffusivity_um2_per_s: float, wall_shear_rate_per_s: float, distance_um: float
) -> float:
    """The velocity at which a suspension is delivered to a capturing wall, in um per s.

    Multiplied by the concentration in the bulk it gives the flux onto the
    wall, a distance ``distance_um`` downstream of where capture begins.
    """
    require(diffusivity_um2_per_s > 0, "the diffusivity must be positive")
    require(wall_shear_rate_per_s > 0, "the wall shear rate must be positive")
    require(distance_um > 0, "the distance downstream must be positive")
    return LEVEQUE * (diffusivity_um2_per_s**2 * wall_shear_rate_per_s / distance_um) ** (1 / 3)


def film_wall_shear_rate_per_s(mean_velocity_um_per_s: float, thickness_um: float) -> float:
    """The wall shear rate of a thin film with a free top: 3 u / delta."""
    require(mean_velocity_um_per_s > 0, "the film's velocity must be positive")
    require(thickness_um > 0, "the film's thickness must be positive")
    return 3.0 * mean_velocity_um_per_s / thickness_um
