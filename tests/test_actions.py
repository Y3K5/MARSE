import pytest

from marse.experimental.host.actions import ActionError, ActionRule, ResourceBudget, try_action


def test_action_consumes_only_declared_resource():
    result = try_action(
        ActionRule("attack", "energy", 2.0),
        ResourceBudget((("energy", 3.0), ("oxygen", 1.0))),
    )
    assert result.executed
    assert result.budget.get("energy") == pytest.approx(1.0)
    assert result.budget.get("oxygen") == pytest.approx(1.0)


def test_action_requires_capability_and_sufficient_resource():
    rule = ActionRule("efflux", "energy", 1.0, prerequisite="drug_resistance")
    budget = ResourceBudget((("energy", 0.5),))
    missing = try_action(rule, budget)
    assert not missing.executed
    assert "prerequisite" in missing.reason
    insufficient = try_action(
        rule,
        ResourceBudget((("energy", 0.5),)),
        available_capabilities=frozenset({"drug_resistance"}),
    )
    assert not insufficient.executed
    assert "insufficient" in insufficient.reason


def test_invalid_action_cost_is_rejected():
    with pytest.raises(ActionError, match="cost"):
        ActionRule("invalid", "energy", -1.0)
