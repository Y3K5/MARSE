"""Tests for reproducible uncertainty propagation and sensitivity."""

import numpy as np
import pytest

from marse.analysis.uncertainty import (
    ParameterRange,
    UncertaintyError,
    propagate,
    rank_sensitivity,
)


def test_latin_hypercube_propagation_is_reproducible_and_summarized():
    parameters = (ParameterRange("growth", 0.5, 1.5), ParameterRange("dose", 0.0, 2.0))

    def evaluate(values):
        return 2.0 * values["growth"] - values["dose"]

    first = propagate(parameters, evaluate, sample_count=64, seed=12)
    second = propagate(parameters, evaluate, sample_count=64, seed=12)
    np.testing.assert_array_equal(first.samples, second.samples)
    np.testing.assert_array_equal(first.outputs, second.outputs)
    assert first.summary["q05"] <= first.summary["q50"] <= first.summary["q95"]


def test_rank_sensitivity_identifies_monotonic_driver():
    parameters = (ParameterRange("strong", 0.0, 1.0), ParameterRange("weak", 0.0, 1.0))
    result = propagate(
        parameters,
        lambda values: 10.0 * values["strong"] + 0.01 * values["weak"],
        sample_count=100,
        seed=4,
    )
    sensitivity = rank_sensitivity(result)
    assert sensitivity.correlations["strong"] > sensitivity.correlations["weak"]
    assert sensitivity.correlations["strong"] > 0.95


def test_invalid_ranges_counts_and_nonfinite_outputs_are_rejected():
    with pytest.raises(UncertaintyError, match="lower < upper"):
        ParameterRange("bad", 1.0, 1.0)
    with pytest.raises(UncertaintyError, match="at least two"):
        propagate((ParameterRange("x", 0.0, 1.0),), lambda _: 1.0, sample_count=1)
    with pytest.raises(UncertaintyError, match="finite"):
        propagate(
            (ParameterRange("x", 0.0, 1.0),),
            lambda _: float("nan"),
            sample_count=4,
        )


def test_small_sensitivity_samples_are_flagged():
    result = propagate(
        (ParameterRange("x", 0.0, 1.0),),
        lambda values: values["x"],
        sample_count=4,
        seed=2,
    )
    sensitivity = rank_sensitivity(result)
    assert any("exploratory" in warning for warning in sensitivity.warnings)
