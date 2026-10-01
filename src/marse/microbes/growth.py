"""Primary growth models: how a population grows under fixed conditions.

Every function is pure and element-wise: pass floats or NumPy arrays (parameters
broadcast too) and get the same shape back. Units are the caller's choice but must
be consistent, for example hours for time and mM for concentrations. The equations,
their assumptions and their sources are in docs/theory.md, section 1.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from marse._numeric import FloatOrArray, as_array, require, unwrap

__all__ = [
    "baranyi_roberts",
    "doubling_time",
    "haldane",
    "minimum_substrate_concentration",
    "monod",
    "net_growth_rate",
    "noncompetitive_inhibition",
    "product_formation_rate",
    "specific_growth_rate",
    "substrate_uptake_rate",
]


def doubling_time(mu: ArrayLike) -> FloatOrArray:
    """Doubling time ln(2)/mu of exponential growth at specific growth rate mu > 0."""
    rate = as_array(mu)
    require(np.all(rate > 0), "mu must be positive")
    return unwrap(np.log(2.0) / rate)


def specific_growth_rate(doubling: ArrayLike) -> FloatOrArray:
    """Specific growth rate ln(2)/t_d for a doubling time t_d > 0."""
    td = as_array(doubling)
    require(np.all(td > 0), "doubling time must be positive")
    return unwrap(np.log(2.0) / td)


def monod(substrate: ArrayLike, mu_max: ArrayLike, k_s: ArrayLike) -> FloatOrArray:
    """Monod (1949) specific growth rate, mu = mu_max * S / (K_s + S).

    Negative concentrations, which can appear as numerical undershoot, count as zero.
    For growth on two essential substrates (for example carbon and oxygen), multiply
    the terms: ``mu_max * monod(s, 1, k_s) * monod(o2, 1, k_o2)``.
    """
    s = np.maximum(as_array(substrate), 0.0)
    mu_max, k_s = as_array(mu_max), as_array(k_s)
    require(np.all(mu_max >= 0), "mu_max must be non-negative")
    require(np.all(k_s > 0), "k_s must be positive")
    return unwrap(mu_max * s / (k_s + s))


def haldane(
    substrate: ArrayLike, mu_max: ArrayLike, k_s: ArrayLike, k_i: ArrayLike
) -> FloatOrArray:
    """Substrate-inhibited growth (Andrews 1968), mu = mu_max S / (K_s + S + S^2 / K_i).

    The rate peaks at S = sqrt(K_s K_i) and approaches Monod kinetics as K_i grows.
    """
    s = np.maximum(as_array(substrate), 0.0)
    mu_max, k_s, k_i = as_array(mu_max), as_array(k_s), as_array(k_i)
    require(np.all(mu_max >= 0), "mu_max must be non-negative")
    require(np.all(k_s > 0) and np.all(k_i > 0), "k_s and k_i must be positive")
    return unwrap(mu_max * s / (k_s + s + s**2 / k_i))


def noncompetitive_inhibition(concentration: ArrayLike, k_i: ArrayLike) -> FloatOrArray:
    """Non-competitive inhibition factor, K_i / (K_i + C).

    1 without the inhibitor and 1/2 at C = K_i. Multiplying a rate by it is how
    an anaerobe's growth is slowed by oxygen. Negative concentrations count as
    zero.
    """
    c = np.maximum(as_array(concentration), 0.0)
    k_i = as_array(k_i)
    require(np.all(k_i > 0), "k_i must be positive")
    return unwrap(k_i / (k_i + c))


def net_growth_rate(
    substrate: ArrayLike, mu_max: ArrayLike, k_s: ArrayLike, decay: ArrayLike
) -> FloatOrArray:
    """Monod growth less maintenance and decay, mu_net = mu_max S / (K_s + S) - b.

    Unlike bare Monod kinetics this is *negative* at low substrate and crosses
    zero at a finite concentration, so a population stops growing where a real
    one does instead of creeping upward forever. That crossing point is
    :func:`minimum_substrate_concentration`. See docs/theory.md, section 1.3, and
    docs/modeling-landscape.md, section 3.4, for why the threshold matters
    whenever conditions approach those found in natural environments.
    """
    b = as_array(decay)
    require(np.all(b >= 0), "decay must be non-negative")
    return unwrap(as_array(monod(substrate, mu_max, k_s)) - b)


def minimum_substrate_concentration(
    mu_max: ArrayLike, k_s: ArrayLike, decay: ArrayLike
) -> FloatOrArray:
    """Substrate concentration below which there is no net growth (``S_min``).

    Solving ``net_growth_rate = 0`` gives S_min = K_s b / (mu_max - b); it is
    infinite when maintenance and decay match or exceed the maximum growth rate,
    meaning the population cannot sustain itself at any concentration. This is
    the same algebra as
    :func:`marse.validation.analytical.chemostat_break_even`, with maintenance
    in place of the dilution rate.
    """
    mu_max, k_s, b = as_array(mu_max), as_array(k_s), as_array(decay)
    require(np.all(mu_max >= 0), "mu_max must be non-negative")
    require(np.all(k_s > 0), "k_s must be positive")
    require(np.all(b >= 0), "decay must be non-negative")
    with np.errstate(divide="ignore", invalid="ignore"):
        threshold = k_s * b / (mu_max - b)
    return unwrap(np.where(b >= mu_max, np.inf, threshold))


def substrate_uptake_rate(
    mu: ArrayLike, biomass: ArrayLike, yield_coefficient: ArrayLike, maintenance: ArrayLike = 0.0
) -> FloatOrArray:
    """Volumetric substrate consumption (mu / Y + m) X, following Pirt (1965).

    Y is the true growth yield (biomass per substrate) and m the maintenance
    coefficient (substrate per biomass per time), consumed even without growth.
    """
    y, m = as_array(yield_coefficient), as_array(maintenance)
    require(np.all(y > 0), "yield_coefficient must be positive")
    require(np.all(m >= 0), "maintenance must be non-negative")
    return unwrap((as_array(mu) / y + m) * as_array(biomass))


def product_formation_rate(
    mu: ArrayLike, biomass: ArrayLike, alpha: ArrayLike, beta: ArrayLike = 0.0
) -> FloatOrArray:
    """Luedeking & Piret (1959) production rate (alpha mu + beta) X.

    alpha is the growth-associated and beta the non-growth-associated coefficient;
    MARSE uses this form for extracellular matrix and for metabolites that other
    species feed on.
    """
    return unwrap((as_array(alpha) * as_array(mu) + as_array(beta)) * as_array(biomass))


def baranyi_roberts(
    t: ArrayLike,
    y0: float,
    y_max: float,
    mu_max: float,
    lag: float,
    curvature: float = 1.0,
) -> FloatOrArray:
    """Growth curve of Baranyi & Roberts (1994) under constant conditions.

    Returns y(t) = ln N(t), the natural logarithm of population density, rising from
    ``y0`` after a lag, growing exponentially at ``mu_max`` and levelling off at
    ``y_max``. The lag enters through h0 = mu_max * lag, the physiological state of
    the inoculum; ``curvature`` (m in the original paper) shapes the approach to
    stationary phase. Times must be non-negative.
    """
    time = as_array(t)
    require(np.all(time >= 0), "t must be non-negative")
    require(mu_max > 0, "mu_max must be positive")
    require(lag >= 0, "lag must be non-negative")
    require(y_max > y0, "y_max must exceed y0")
    require(curvature > 0, "curvature must be positive")
    h0 = mu_max * lag
    decay = np.exp(-mu_max * time)
    with np.errstate(divide="ignore"):  # log1p(-1) = -inf at t = 0 is intended
        # A(t) = t + ln(e^(-mu t) + e^(-h0) - e^(-mu t - h0)) / mu, in a stable form.
        adjusted = time + np.logaddexp(-mu_max * time, -h0 + np.log1p(-decay)) / mu_max
    span = curvature * (y_max - y0)
    # y = y0 + mu A - ln(1 + (e^(m mu A) - 1) / e^(m (y_max - y0))) / m, in a stable form.
    braking = np.logaddexp(curvature * mu_max * adjusted - span, np.log1p(-np.exp(-span)))
    return unwrap(y0 + mu_max * adjusted - braking / curvature)
