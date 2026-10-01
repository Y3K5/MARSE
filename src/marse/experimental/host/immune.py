"""Bounded immune and molecular interaction primitives.

These functions operate on explicit concentration fields and biomass arrays.
They are intentionally phenomenological: they do not infer receptor biology,
cell trajectories, or clinical outcomes from sparse parameters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from marse.experimental.host.actions import ActionRule, ResourceBudget, try_action
from marse.microbes.additives import hill_response

__all__ = [
    "ImmuneAgent",
    "ImmuneAgentStepResult",
    "ImmuneCellType",
    "ImmuneError",
    "ImmuneInteraction",
    "ImmunePressureResult",
    "MolecularNeutralizer",
    "apply_immune_pressure",
    "step_immune_agents",
]


class ImmuneError(ValueError):
    """An immune or molecular interaction is invalid."""


@dataclass(frozen=True, slots=True)
class ImmuneCellType:
    """Configurable action and resource rules for one immune-cell archetype."""

    name: str
    action: str
    effector: str
    attack_per_h: float = 0.0
    movement_per_h: float = 0.0
    secretion_per_h: float = 0.0

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.effector.strip():
            raise ImmuneError("immune cell name and effector must not be empty")
        if self.action not in {"move_toward", "attack", "secrete"}:
            raise ImmuneError("immune action must be move_toward, attack, or secrete")
        if min(self.attack_per_h, self.movement_per_h, self.secretion_per_h) < 0:
            raise ImmuneError("immune action rates must be non-negative")


@dataclass(frozen=True, slots=True)
class ImmuneAgent:
    """A grid-local immune cell with an explicit action budget."""

    cell_type: ImmuneCellType
    y: int
    x: int
    energy: float = 1.0

    def __post_init__(self) -> None:
        if self.y < 0 or self.x < 0 or self.energy < 0:
            raise ImmuneError("immune agent coordinates and energy must be non-negative")


@dataclass(frozen=True, slots=True)
class ImmuneAgentStepResult:
    """Updated immune agents, target biomass, and secreted effector fields."""

    agents: tuple[ImmuneAgent, ...]
    biomass: NDArray[np.float64]
    effectors: dict[str, NDArray[np.float64]]
    budgets: tuple[ResourceBudget, ...] | None = None
    blocked_actions: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ImmuneInteraction:
    """Effector-mediated loss for one named species."""

    species: str
    effector: str
    maximum_kill_per_h: float
    half_effect: float
    hill_coefficient: float = 1.0
    susceptibility: float = 1.0

    def __post_init__(self) -> None:
        if not self.species.strip() or not self.effector.strip():
            raise ImmuneError("species and effector names must not be empty")
        if self.maximum_kill_per_h < 0 or self.half_effect <= 0:
            raise ImmuneError("kill rate must be non-negative and half_effect positive")
        if self.hill_coefficient <= 0 or self.susceptibility < 0:
            raise ImmuneError("hill coefficient must be positive and susceptibility non-negative")


@dataclass(frozen=True, slots=True)
class MolecularNeutralizer:
    """Reduce an effector field by a competing molecular concentration."""

    effector: str
    molecule: str
    neutralization_fraction: float
    half_effect: float
    hill_coefficient: float = 1.0

    def __post_init__(self) -> None:
        if not self.effector.strip() or not self.molecule.strip():
            raise ImmuneError("effector and molecule names must not be empty")
        if not 0 <= self.neutralization_fraction <= 1:
            raise ImmuneError("neutralization_fraction must be in [0, 1]")
        if self.half_effect <= 0 or self.hill_coefficient <= 0:
            raise ImmuneError("neutralizer half_effect and hill_coefficient must be positive")


@dataclass(frozen=True, slots=True)
class ImmunePressureResult:
    """Biomass after pressure and the removed biomass field."""

    biomass: NDArray[np.float64]
    removed: NDArray[np.float64]


def apply_immune_pressure(
    biomass: ArrayLike,
    *,
    effectors: dict[str, ArrayLike],
    molecules: dict[str, ArrayLike] | None = None,
    interaction: ImmuneInteraction,
    neutralizers: tuple[MolecularNeutralizer, ...] = (),
    dt: float,
) -> ImmunePressureResult:
    """Apply deterministic, bounded effector killing to one biomass field."""
    if dt < 0:
        raise ImmuneError("dt must be non-negative")
    target = np.maximum(np.asarray(biomass, dtype=float), 0.0)
    if interaction.effector not in effectors:
        raise ImmuneError(f"missing effector field '{interaction.effector}'")
    effector = np.maximum(np.asarray(effectors[interaction.effector], dtype=float), 0.0)
    if effector.shape != target.shape:
        raise ImmuneError("effector and biomass fields must have matching shapes")
    effective_effector = effector.copy()
    for neutralizer in neutralizers:
        if neutralizer.effector != interaction.effector:
            continue
        if molecules is None or neutralizer.molecule not in molecules:
            raise ImmuneError(f"missing neutralizer molecule '{neutralizer.molecule}'")
        molecule = np.maximum(np.asarray(molecules[neutralizer.molecule], dtype=float), 0.0)
        if molecule.shape != target.shape:
            raise ImmuneError("molecule and biomass fields must have matching shapes")
        occupancy = hill_response(molecule, neutralizer.half_effect, neutralizer.hill_coefficient)
        effective_effector *= 1.0 - neutralizer.neutralization_fraction * occupancy
    kill_rate = (
        interaction.maximum_kill_per_h
        * interaction.susceptibility
        * hill_response(effective_effector, interaction.half_effect, interaction.hill_coefficient)
    )
    updated = target * np.exp(-dt * kill_rate)
    return ImmunePressureResult(updated, target - updated)


def step_immune_agents(
    agents: tuple[ImmuneAgent, ...],
    biomass: ArrayLike,
    *,
    dt: float,
    effectors: dict[str, ArrayLike] | None = None,
    budgets: tuple[ResourceBudget, ...] | None = None,
    capabilities: tuple[frozenset[str], ...] | None = None,
) -> ImmuneAgentStepResult:
    """Execute one deterministic action step for grid-local immune agents.

    ``move_toward`` follows the steepest adjacent biomass gradient, ``attack``
    removes biomass at the agent cell, and ``secrete`` adds an effector field.
    These are action primitives, not claims about a particular cell lineage.
    """
    if dt < 0:
        raise ImmuneError("dt must be non-negative")
    updated_biomass = np.maximum(np.asarray(biomass, dtype=float), 0.0).copy()
    if updated_biomass.ndim != 2:
        raise ImmuneError("biomass must be a 2D field")
    fields = {
        name: np.maximum(np.asarray(field, dtype=float), 0.0).copy()
        for name, field in (effectors or {}).items()
    }
    if budgets is not None and len(budgets) != len(agents):
        raise ImmuneError("one resource budget is required per immune agent")
    if capabilities is not None and len(capabilities) != len(agents):
        raise ImmuneError("one capability set is required per immune agent")
    next_budgets = list(budgets) if budgets is not None else None
    blocked: list[str] = []
    if any(field.shape != updated_biomass.shape for field in fields.values()):
        raise ImmuneError("effector and biomass fields must have matching shapes")
    moved: list[ImmuneAgent] = []
    height, width = updated_biomass.shape
    for index, agent in enumerate(agents):
        if agent.y >= height or agent.x >= width:
            raise ImmuneError("immune agent coordinates exceed biomass field")
        cell_type = agent.cell_type
        y, x = agent.y, agent.x
        budget = next_budgets[index] if next_budgets is not None else None
        available = capabilities[index] if capabilities is not None else frozenset()
        if cell_type.action == "move_toward":
            if budget is not None:
                result = try_action(
                    ActionRule("move", "energy", dt * cell_type.movement_per_h),
                    budget,
                    available_capabilities=available,
                )
                if not result.executed:
                    blocked.append("move")
                    moved.append(agent)
                    continue
                next_budgets[index] = result.budget
            candidates = [
                (y, x),
                *((y + dy, x + dx) for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1))),
            ]
            valid = [(cy, cx) for cy, cx in candidates if 0 <= cy < height and 0 <= cx < width]
            y, x = max(valid, key=lambda position: updated_biomass[position])
        elif cell_type.action == "attack":
            if budget is not None:
                result = try_action(
                    ActionRule("attack", "energy", dt * cell_type.attack_per_h),
                    budget,
                    available_capabilities=available,
                )
                if not result.executed:
                    blocked.append("attack")
                    moved.append(agent)
                    continue
                next_budgets[index] = result.budget
            updated_biomass[y, x] *= np.exp(-dt * cell_type.attack_per_h * agent.energy)
        elif cell_type.action == "secrete":
            if budget is not None:
                result = try_action(
                    ActionRule("secrete", "energy", dt * cell_type.secretion_per_h),
                    budget,
                    available_capabilities=available,
                )
                if not result.executed:
                    blocked.append("secrete")
                    moved.append(agent)
                    continue
                next_budgets[index] = result.budget
            field = fields.setdefault(cell_type.effector, np.zeros_like(updated_biomass))
            field[y, x] += dt * cell_type.secretion_per_h * agent.energy
        moved.append(ImmuneAgent(cell_type, y, x, agent.energy))
    return ImmuneAgentStepResult(
        tuple(moved),
        updated_biomass,
        fields,
        tuple(next_budgets) if next_budgets is not None else None,
        tuple(blocked),
    )
