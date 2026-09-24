"""Physical properties of dissolved oxygen in water (docs/theory.md, section 2.5).

They set the oxygen concentration at an air interface and how fast oxygen
diffuses, both of which change markedly with temperature.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from marse._numeric import FloatOrArray, as_array, require, unwrap

__all__ = [
    "OXYGEN_MOLAR_MASS",
    "SECONDS_PER_HOUR",
    "oxygen_diffusivity_m2_per_s",
    "oxygen_diffusivity_um2_per_h",
    "oxygen_saturation_mg_per_l",
]

SECONDS_PER_HOUR = 3600.0
"""Diffusivities are published per second; MARSE works in hours internally."""

OXYGEN_MOLAR_MASS = 31.998
"""Molar mass of O2 in g/mol; divide mg/L by it to get mM."""


def oxygen_saturation_mg_per_l(temperature_c: ArrayLike) -> FloatOrArray:
    """Oxygen in fresh water at equilibrium with air at 1 atm, in mg/L.

    Benson & Krause (1984), as used for dissolved-oxygen tables; valid from 0 to
    40 degrees Celsius. Dissolved salts lower the value (culture media hold less
    oxygen than pure water).
    """
    temp = as_array(temperature_c)
    require(np.all((temp >= 0) & (temp <= 40)), "valid from 0 to 40 degrees Celsius")
    kelvin = temp + 273.15
    ln_c = (
        -139.34411
        + 1.575701e5 / kelvin
        - 6.642308e7 / kelvin**2
        + 1.243800e10 / kelvin**3
        - 8.621949e11 / kelvin**4
    )
    return unwrap(np.exp(ln_c))


def oxygen_diffusivity_m2_per_s(temperature_c: ArrayLike) -> FloatOrArray:
    """Diffusion coefficient of oxygen in water, in m^2/s.

    Interpolation formula of Han & Bartels (1996),
    log10(D / cm^2 s^-1) = -4.410 + 773.8 / T - (506.4 / T)^2 with T in kelvin,
    measured from 0 to 95 degrees Celsius.
    """
    temp = as_array(temperature_c)
    require(np.all((temp >= 0) & (temp <= 95)), "valid from 0 to 95 degrees Celsius")
    kelvin = temp + 273.15
    log10_cm2_per_s = -4.410 + 773.8 / kelvin - (506.4 / kelvin) ** 2
    return unwrap(1e-4 * 10.0**log10_cm2_per_s)


def oxygen_diffusivity_um2_per_h(temperature_c: ArrayLike) -> FloatOrArray:
    """Oxygen diffusivity in um^2/h, the unit the rest of MARSE works in.

    The same quantity as :func:`oxygen_diffusivity_m2_per_s`, converted. It
    exists because growth rates are per hour while diffusivities are published
    per second, and combining the two unconverted understates how far a solute
    penetrates by a factor of sixty while still producing a plausible-looking
    number.
    """
    return unwrap(as_array(oxygen_diffusivity_m2_per_s(temperature_c)) * 1e12 * SECONDS_PER_HOUR)
