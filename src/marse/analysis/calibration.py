"""Replicate-aware calibration of simple microbial growth curves.

The fitter is intentionally conservative and dependency-light. It supports
exponential, logistic, Gompertz, and Baranyi curves, reports residuals, and
refuses datasets that cannot identify the requested model. It is a calibration
aid, not a replacement for experimental design or statistical review.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from marse.microbes.growth import baranyi_roberts

__all__ = ["CalibrationError", "GrowthCurve", "GrowthFit", "fit_growth_curve"]

ModelName = Literal["exponential", "logistic", "gompertz", "baranyi"]


class CalibrationError(ValueError):
    """Growth-curve data cannot support the requested calibration."""


@dataclass(frozen=True, slots=True)
class GrowthCurve:
    """A time series with rows representing independent replicates."""

    time_h: NDArray[np.float64]
    observations: NDArray[np.float64]
    endpoint: str = "population"
    transform: str = "log"

    def __post_init__(self) -> None:
        time = np.asarray(self.time_h, dtype=float)
        values = np.asarray(self.observations, dtype=float)
        if time.ndim != 1 or time.size < 3:
            raise CalibrationError("time_h: at least three time points are required")
        if values.ndim == 1:
            values = values[None, :]
        if values.ndim != 2 or values.shape[1] != time.size:
            raise CalibrationError("observations must have shape (replicates, time points)")
        if np.any(~np.isfinite(time)) or np.any(~np.isfinite(values)):
            raise CalibrationError("time and observations must be finite")
        if np.any(time < 0) or np.any(np.diff(time) <= 0):
            raise CalibrationError("time_h must be strictly increasing and non-negative")
        if np.any(values <= 0):
            raise CalibrationError("growth observations must be positive")
        if self.transform not in {"log", "identity"}:
            raise CalibrationError("transform must be 'log' or 'identity'")
        object.__setattr__(self, "time_h", time)
        object.__setattr__(self, "observations", values)


@dataclass(frozen=True, slots=True)
class GrowthFit:
    """Fitted parameters and diagnostics, with units implied by the curve."""

    model: ModelName
    parameters: dict[str, float]
    predicted: NDArray[np.float64]
    residuals: NDArray[np.float64]
    rmse: float
    replicate_count: int
    identifiable: bool
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "model": self.model,
            "parameters": self.parameters,
            "rmse": self.rmse,
            "replicate_count": self.replicate_count,
            "identifiable": self.identifiable,
            "warnings": list(self.warnings),
        }


def _target(curve: GrowthCurve) -> NDArray[np.float64]:
    return np.log(curve.observations) if curve.transform == "log" else curve.observations


def _fit_exponential(curve: GrowthCurve) -> tuple[dict[str, float], NDArray[np.float64]]:
    target = _target(curve)
    mean = target.mean(axis=0)
    slope, intercept = np.polyfit(curve.time_h, mean, 1)
    if slope <= 0:
        raise CalibrationError("exponential fit found a non-positive growth rate")
    prediction = intercept + slope * curve.time_h
    return {"mu_per_h": float(slope), "intercept": float(intercept)}, prediction


def _curve_prediction(
    model: ModelName, time: NDArray[np.float64], p: NDArray[np.float64]
) -> NDArray:
    if model == "logistic":
        y0, ymax, mu = p
        return ymax - np.log1p((np.exp(ymax - y0) - 1.0) * np.exp(-mu * time))
    if model == "gompertz":
        y0, ymax, mu, lag = p
        return y0 + (ymax - y0) * np.exp(-np.exp(1.0 + mu * (lag - time) / (ymax - y0)))
    if model == "baranyi":
        y0, ymax, mu, lag = p
        return np.asarray(baranyi_roberts(time, y0, ymax, mu, lag))
    raise CalibrationError(f"unsupported nonlinear model '{model}'")


def _fit_nonlinear(curve: GrowthCurve, model: ModelName) -> tuple[dict[str, float], NDArray]:
    target = _target(curve)
    mean = target.mean(axis=0)
    y0 = float(mean[0])
    ymax = float(mean.max())
    if ymax <= y0:
        raise CalibrationError("curve has no measurable growth span")
    span = ymax - y0
    if model == "logistic":
        initial = np.array([y0, ymax + 0.1 * span, 1.0], dtype=float)
        bounds = np.array([[y0 - span, y0 + 1e-9], [ymax, ymax + 20 * span], [1e-6, 20.0]])
    else:
        initial = np.array([y0, ymax + 0.1 * span, 1.0, max(curve.time_h[1], 0.1)], dtype=float)
        bounds = np.array(
            [
                [y0 - span, y0 + 1e-9],
                [ymax, ymax + 20 * span],
                [1e-6, 20.0],
                [0.0, max(curve.time_h[-1] * 2, 1.0)],
            ]
        )
    params = initial

    def loss(candidate: NDArray) -> float:
        try:
            prediction = _curve_prediction(model, curve.time_h, candidate)
        except (FloatingPointError, ValueError):
            return float("inf")
        if not np.all(np.isfinite(prediction)):
            return float("inf")
        return float(np.mean((prediction - mean) ** 2))

    best = loss(params)
    steps = (bounds[:, 1] - bounds[:, 0]) / 4.0
    for _ in range(48):
        improved = False
        for index in range(params.size):
            for direction in (-1.0, 1.0):
                candidate = params.copy()
                candidate[index] = np.clip(
                    candidate[index] + direction * steps[index], *bounds[index]
                )
                score = loss(candidate)
                if score < best:
                    params, best, improved = candidate, score, True
        steps *= 0.75
        if not improved and np.max(steps) < 1e-5:
            break
    prediction = _curve_prediction(model, curve.time_h, params)
    names = (
        ("y0", "ymax", "mu_per_h")
        if model == "logistic"
        else (
            "y0",
            "ymax",
            "mu_per_h",
            "lag_h",
        )
    )
    return dict(zip(names, (float(value) for value in params), strict=True)), prediction


def fit_growth_curve(curve: GrowthCurve, model: ModelName = "logistic") -> GrowthFit:
    """Fit a growth curve to replicate means and return pointwise diagnostics."""
    if curve.observations.shape[0] < 2:
        warnings = ("only one replicate supplied; uncertainty cannot be estimated",)
    else:
        warnings = ()
    if model == "exponential":
        parameters, predicted = _fit_exponential(curve)
    else:
        parameters, predicted = _fit_nonlinear(curve, model)
    target = _target(curve)
    residuals = target - predicted
    rmse = float(np.sqrt(np.mean(residuals**2)))
    identifiable = curve.time_h.size >= (2 if model == "exponential" else 4)
    if not identifiable:
        warnings += ("too few time points for the requested parameterization",)
    return GrowthFit(
        model,
        parameters,
        np.broadcast_to(predicted, target.shape).copy(),
        residuals,
        rmse,
        curve.observations.shape[0],
        identifiable,
        warnings,
    )
