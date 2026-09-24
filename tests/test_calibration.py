"""Tests for conservative growth-curve calibration."""

import numpy as np
import pytest

from marse.calibration import CalibrationError, GrowthCurve, fit_growth_curve


def test_exponential_fit_recovers_rate_from_replicates():
    time = np.linspace(0.0, 4.0, 9)
    expected = np.exp(-1.0 + 0.7 * time)
    curve = GrowthCurve(time, np.vstack((expected, expected * 1.05)))
    fit = fit_growth_curve(curve, "exponential")
    assert fit.parameters["mu_per_h"] == pytest.approx(0.7, rel=1e-3)
    assert fit.replicate_count == 2
    assert fit.identifiable
    assert fit.rmse < 0.04


def test_logistic_fit_reports_saturation_and_residuals():
    time = np.linspace(0.0, 8.0, 17)
    y0, ymax, mu = -4.0, 2.0, 0.9
    values = ymax - np.log1p((np.exp(ymax - y0) - 1.0) * np.exp(-mu * time))
    fit = fit_growth_curve(GrowthCurve(time, np.exp(np.vstack((values, values)))), "logistic")
    assert fit.parameters["mu_per_h"] == pytest.approx(mu, rel=0.1)
    assert fit.parameters["ymax"] > fit.parameters["y0"]
    assert fit.rmse < 0.1


def test_baranyi_and_gompertz_are_supported():
    time = np.linspace(0.0, 10.0, 21)
    observations = np.exp(
        np.array(
            [
                -4.0 + 6.0 * np.exp(-np.exp(1.0 + 0.8 * (2.0 - time) / 6.0)),
                -4.0 + 6.0 * np.exp(-np.exp(1.0 + 0.8 * (2.0 - time) / 6.0)),
            ]
        )
    )
    assert fit_growth_curve(GrowthCurve(time, observations), "gompertz").identifiable
    assert fit_growth_curve(GrowthCurve(time, observations), "baranyi").identifiable


def test_invalid_and_underdetermined_data_are_explicit():
    with pytest.raises(CalibrationError, match="strictly increasing"):
        GrowthCurve(np.array([0.0, 1.0, 1.0]), np.ones((2, 3)))
    time = np.array([0.0, 1.0, 2.0])
    fit = fit_growth_curve(GrowthCurve(time, np.array([[1.0, 2.0, 4.0]])), "logistic")
    assert not fit.identifiable
    assert any("one replicate" in warning for warning in fit.warnings)
