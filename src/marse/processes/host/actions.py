"""Bounded action costs and prerequisites for game-like simulations."""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["ActionError", "ActionResult", "ActionRule", "ResourceBudget", "try_action"]


class ActionError(ValueError):
    """An action rule or resource budget is invalid."""


@dataclass(frozen=True, slots=True)
class ActionRule:
    """One typed action with a resource cost and optional prerequisite."""

    name: str
    resource: str
    cost: float
    prerequisite: str | None = None

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.resource.strip():
            raise ActionError("action name and resource must not be empty")
        if self.cost < 0:
            raise ActionError("action cost must be non-negative")
        if self.prerequisite is not None and not self.prerequisite.strip():
            raise ActionError("action prerequisite must not be empty")


@dataclass(frozen=True, slots=True)
class ResourceBudget:
    """Immutable resource balances consumed by actions."""

    values: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        names = [name for name, _ in self.values]
        if len(set(names)) != len(names) or any(not name.strip() for name in names):
            raise ActionError("resource names must be unique and non-empty")
        if any(value < 0 for _, value in self.values):
            raise ActionError("resource balances must be non-negative")

    def get(self, resource: str) -> float:
        return dict(self.values).get(resource, 0.0)


@dataclass(frozen=True, slots=True)
class ActionResult:
    """Whether an action executed, with the resulting budget and reason."""

    executed: bool
    budget: ResourceBudget
    reason: str


def try_action(
    rule: ActionRule,
    budget: ResourceBudget,
    *,
    available_capabilities: frozenset[str] = frozenset(),
) -> ActionResult:
    """Execute an action only when its capability and resources are available."""
    if rule.prerequisite is not None and rule.prerequisite not in available_capabilities:
        return ActionResult(False, budget, f"missing prerequisite '{rule.prerequisite}'")
    current = budget.get(rule.resource)
    if current < rule.cost:
        return ActionResult(False, budget, f"insufficient resource '{rule.resource}'")
    updated = dict(budget.values)
    updated[rule.resource] = current - rule.cost
    return ActionResult(
        True,
        ResourceBudget(tuple(sorted(updated.items()))),
        "executed",
    )
