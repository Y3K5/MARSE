"""Secondary growth models: how temperature and pH scale the maximum growth rate.

The cardinal models return a factor gamma between 0 and 1, which multiplies the
optimal rate following the gamma concept of Zwietering et al. (1992):
mu_max(T, pH) = mu_opt * gamma_T(T) * gamma_pH(pH). The factors assume that
temperature and pH act independently, an approximation that breaks down near the
growth limits (docs/theory.md, section 2).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from marse._numeric import FloatOrArray, as_array, require, unwrap

__all__ = ["cardinal_ph", "cardinal_temperature", "ratkowsky"]


def cardinal_temperature(
    temperature: ArrayLike, t_min: float, t_opt: float, t_max: float
) -> FloatOrArray:
    """Cardinal temperature model with inflection (CTMI) of Rosso et al. (1993).

    Returns gamma_T: 1 at ``t_opt``, falling to 0 at ``t_min`` and ``t_max`` and
    staying 0 outside them. Only temperature differences enter, so degrees Celsius
    and kelvin give the same result. The formula is well defined only when the
    optimum lies nearer the maximum, t_opt >= (t_min + t_max) / 2, as it does for
    typical microbial data; other inputs are rejected.
    """
    require(t_min < t_opt < t_max, "cardinal temperatures must satisfy t_min < t_opt < t_max")
    require(
        t_opt >= (t_min + t_max) / 2,
        "the CTMI needs t_opt >= (t_min + t_max) / 2; fit a different model for this organism",
    )
    temp = as_array(temperature)
    inside = (temp > t_min) & (temp < t_max)
    x = np.where(inside, temp, t_opt)  # evaluate only where the formula applies
    numerator = (x - t_max) * (x - t_min) ** 2
    denominator = (t_opt - t_min) * (
        (t_opt - t_min) * (x - t_opt) - (t_opt - t_max) * (t_opt + t_min - 2 * x)
    )
    return unwrap(np.where(inside, numerator / denominator, 0.0))


def cardinal_ph(ph: ArrayLike, ph_min: float, ph_opt: float, ph_max: float) -> FloatOrArray:
    """Cardinal pH model (CPM) of Rosso et al. (1995).

    Returns gamma_pH: 1 at ``ph_opt``, falling to 0 at ``ph_min`` and ``ph_max`` and
    staying 0 outside them.
    """
    require(ph_min < ph_opt < ph_max, "cardinal pH values must satisfy ph_min < ph_opt < ph_max")
    value = as_array(ph)
    inside = (value > ph_min) & (value < ph_max)
    x = np.where(inside, value, ph_opt)
    numerator = (x - ph_min) * (x - ph_max)
    return unwrap(np.where(inside, numerator / (numerator - (x - ph_opt) ** 2), 0.0))


def ratkowsky(
    temperature: ArrayLike, b: float, t_min: float, t_max: float, c: float
) -> FloatOrArray:
    """Square-root model over the whole growth range (Ratkowsky et al. 1983).

    Returns the growth rate itself, mu = [b (T - T_min) (1 - exp(c (T - T_max)))]^2,
    between ``t_min`` and ``t_max`` and 0 outside. Below the optimum it reduces to
    the linear square-root relation of Ratkowsky et al. (1982).
    """
    require(t_min < t_max, "t_min must be below t_max")
    require(b > 0 and c > 0, "b and c must be positive")
    temp = as_array(temperature)
    inside = (temp > t_min) & (temp < t_max)
    root = b * (temp - t_min) * (1.0 - np.exp(c * (temp - t_max)))
    return unwrap(np.where(inside, root**2, 0.0))
