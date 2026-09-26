"""Rate laws and run settings in the version 2 schema.

A rate may not let a process consume what its rate does not depend on (unless
the assumption is stated), may not apply one component's factor twice (known
defect 4's pattern), and every constant carries its unit. The run settings are
refused unless they describe a run that can happen.
"""

import copy
import json
import re
from pathlib import Path

import pytest

from marse.core.config import ConfigError
from marse.schemas import experiment_from_dict, network_from_dict

EXAMPLE = (
    Path(__file__).resolve().parents[1] / "examples" / "networks" / "glucose_cross_feeding.json"
)


def example() -> dict:
    return json.loads(EXAMPLE.read_text("utf-8"))


def monod(component: str, k: float = 0.1) -> dict:
    return {"component": component, "form": "monod", "half_saturation_mol_per_m3": k}


def with_rate(process: int, rate: dict) -> dict:
    raw = example()
    raw["processes"][process]["rate"] = rate
    return raw


def refused(raw: dict, message: str, reader=network_from_dict) -> None:
    with pytest.raises(ConfigError, match=re.escape(message)):
        reader(raw)


def test_the_example_rates_are_read_exactly_and_round_trip():
    raw = example()
    config = experiment_from_dict(raw)
    rate = config.network.process("heterotroph_growth_on_glucose").rate
    assert rate.proportional_to == "heterotroph"  # the default for growth
    assert [f.component for f in rate.factors] == ["glucose", "oxygen", "ammonium"]
    again = experiment_from_dict(json.loads(json.dumps(config.to_dict())))
    assert again == config
    assert again.to_dict() == config.to_dict()


def test_a_substrate_limited_twice_is_refused():
    # Known defect 4: the v1 engine applied a capability's Monod term on top of
    # the substrate's own. Version 2 refuses a second factor for one component.
    rate = {"maximum_per_h": 0.4, "factors": [monod("glucose"), monod("glucose", 0.3)]}
    refused(with_rate(0, rate), "'glucose' appears more than once")


def test_consuming_what_the_rate_ignores_is_refused():
    # Aerobic growth consumes oxygen and ammonium; a rate on glucose alone would
    # keep consuming them after they run out.
    rate = {"maximum_per_h": 0.4, "factors": [monod("glucose")]}
    refused(with_rate(0, rate), "consumes 'oxygen', but its rate does not depend on it")


def test_an_explicit_assumption_of_excess_is_accepted_and_recorded():
    rate = {
        "maximum_per_h": 0.4,
        "factors": [monod("glucose"), monod("oxygen", 0.005)],
        "assumed_in_excess": ["ammonium"],
    }
    network = network_from_dict(with_rate(0, rate))
    written = network.process("heterotroph_growth_on_glucose").to_dict()["rate"]
    assert written["assumed_in_excess"] == ["ammonium"]


@pytest.mark.parametrize(
    ("rate", "message"),
    [
        (
            {
                "maximum_per_h": 0.4,
                "factors": [monod("glucose"), monod("oxygen")],
                "assumed_in_excess": ["ammonium", "lactate"],
            },
            "does not consume 'lactate'",
        ),
        (
            {
                "maximum_per_h": 0.4,
                "factors": [monod("glucose"), monod("oxygen"), monod("ammonium")],
                "assumed_in_excess": ["ammonium"],
            },
            "'ammonium' already limits this rate",
        ),
        ({"maximum_per_h": -0.1}, "must not be negative"),
        ({"maximum_per_h": 0.4, "proportional_to": "heterotrof"}, "did you mean 'heterotroph'"),
        (
            {"maximum_per_h": 0.4, "factors": [{"component": "glucose", "form": "monod"}]},
            "a monod factor needs half_saturation_mol_per_m3",
        ),
        (
            {
                "maximum_per_h": 0.4,
                "factors": [
                    {
                        "component": "glucose",
                        "form": "inhibition",
                        "half_saturation_mol_per_m3": 1.0,
                    }
                ],
            },
            "half_saturation_mol_per_m3 does not apply to an inhibition factor",
        ),
        (
            {
                "maximum_per_h": 0.4,
                "factors": [
                    {"component": "glucose", "form": "haldane", "half_saturation_mol_per_m3": 1.0}
                ],
            },
            "a haldane factor needs inhibition_mol_per_m3",
        ),
        ({"maximum_per_h": 0.4, "factors": [monod("glucose", 0.0)]}, "must be positive"),
        (
            {
                "maximum_per_h": 0.4,
                "factors": [
                    {"component": "glucose", "form": "switch", "half_saturation_mol_per_m3": 1.0}
                ],
            },
            "expected one of 'monod', 'inhibition'",
        ),
        (
            {"maximum_per_h": 0.4, "factors": [{**monod("glucose"), "half_saturation_mm": 1}]},
            "this field is 'half_saturation_mol_per_m3'",
        ),
        ({"maximum_per_s": 0.4}, "this field is 'maximum_per_h', in per hour"),
    ],
)
def test_malformed_rates_are_refused(rate, message):
    refused(with_rate(0, rate), message)


def test_a_reaction_rate_names_what_it_is_proportional_to():
    rate = {"maximum_per_h": 0.01, "factors": [monod("oxygen", 0.005)]}
    refused(with_rate(3, rate), "a reaction's rate needs proportional_to")


def test_an_inhibition_factor_does_not_count_as_depending_on_a_substrate():
    # Inhibition by oxygen does not slow consumption of oxygen as it runs out.
    rate = {
        "maximum_per_h": 0.4,
        "factors": [
            monod("glucose"),
            {"component": "oxygen", "form": "inhibition", "inhibition_mol_per_m3": 0.01},
            monod("ammonium"),
        ],
    }
    refused(with_rate(0, rate), "consumes 'oxygen', but its rate does not depend on it")


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"duration_h": None}, "missing 'duration_h'"),
        ({"experiment_id": "  "}, "must not be empty"),
        ({"timestep_h": 60.0}, "at most duration_h"),
        ({"record_interval_h": 0.1}, "between timestep_h and duration_h"),
        ({"initial_mol_per_m3": {"glucos": 1.0}}, "did you mean 'glucose'"),
        ({"initial_mol_per_m3": {"glucose": -1.0}}, "'glucose' must not be negative"),
        ({"relative_tolerance": 1.5}, "must lie between 0 and 1"),
        ({"absolute_tolerance_mol_per_m3": 0.0}, "must be positive"),
        ({"seed": -1}, "must not be negative"),
        ({"duration_s": 3.0}, "this field is 'duration_h', in hours"),
    ],
)
def test_run_settings_that_describe_no_possible_run_are_refused(change, message):
    raw = example()
    for key, value in change.items():
        if value is None:
            raw.pop(key)
        else:
            raw[key] = value
    refused(raw, message, experiment_from_dict)


def test_a_network_without_rates_is_checked_but_cannot_run():
    raw = example()
    stripped = copy.deepcopy(raw)
    del stripped["processes"][1]["rate"]
    network_from_dict(stripped)  # a valid network
    refused(
        stripped,
        "every process needs a rate to run; 'fermenter_growth_on_glucose'",
        experiment_from_dict,
    )


def test_unlisted_components_start_at_zero_and_the_zeros_are_written_out():
    config = experiment_from_dict(example())
    assert config.initial_mol_per_m3["lactate"] == 0.0
    assert config.to_dict()["initial_mol_per_m3"]["carbon_dioxide"] == 0.0
    assert config.to_dict()["relative_tolerance"] == 1e-6
