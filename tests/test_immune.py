import numpy as np
import pytest

from marse.actions import ResourceBudget
from marse.immune import (
    ImmuneAgent,
    ImmuneCellType,
    ImmuneError,
    ImmuneInteraction,
    MolecularNeutralizer,
    apply_immune_pressure,
    step_immune_agents,
)


def test_effector_pressure_reduces_biomass_and_tracks_removed_mass():
    result = apply_immune_pressure(
        np.ones((2, 2)),
        effectors={"macrophage_signal": np.full((2, 2), 2.0)},
        interaction=ImmuneInteraction("target", "macrophage_signal", 1.0, 1.0),
        dt=1.0,
    )
    assert np.all(result.biomass < 1.0)
    np.testing.assert_allclose(result.removed, 1.0 - result.biomass)


def test_molecular_neutralization_reduces_pressure():
    interaction = ImmuneInteraction("target", "effector", 2.0, 1.0)
    without = apply_immune_pressure(
        [1.0], effectors={"effector": [2.0]}, interaction=interaction, dt=1.0
    )
    with_neutralizer = apply_immune_pressure(
        [1.0],
        effectors={"effector": [2.0]},
        molecules={"drug": [100.0]},
        interaction=interaction,
        neutralizers=(MolecularNeutralizer("effector", "drug", 1.0, 1.0),),
        dt=1.0,
    )
    assert with_neutralizer.biomass[0] > without.biomass[0]


def test_missing_effector_is_explicitly_rejected():
    with pytest.raises(ImmuneError, match="missing effector"):
        apply_immune_pressure(
            [1.0],
            effectors={},
            interaction=ImmuneInteraction("target", "effector", 1.0, 1.0),
            dt=1.0,
        )


def test_immune_actions_move_attack_and_secrete_explicitly():
    mover = ImmuneCellType("macrophage", "move_toward", "il1", movement_per_h=1.0)
    moved = step_immune_agents(
        (ImmuneAgent(mover, 1, 1),),
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 0.0]]),
        dt=0.1,
    )
    assert moved.agents[0].x == 2

    attacker = ImmuneCellType("neutrophil", "attack", "ros", attack_per_h=2.0)
    attacked = step_immune_agents((ImmuneAgent(attacker, 0, 0),), [[1.0]], dt=1.0)
    assert attacked.biomass[0, 0] < 1.0

    secretor = ImmuneCellType("tcell", "secrete", "ifng", secretion_per_h=3.0)
    secreted = step_immune_agents((ImmuneAgent(secretor, 0, 0),), [[0.0]], dt=0.5)
    assert secreted.effectors["ifng"][0, 0] == pytest.approx(1.5)


def test_immune_actions_consume_budgets_and_report_blocked_actions():
    attacker = ImmuneCellType("neutrophil", "attack", "ros", attack_per_h=2.0)
    result = step_immune_agents(
        (ImmuneAgent(attacker, 0, 0),),
        [[1.0]],
        dt=1.0,
        budgets=(ResourceBudget((("energy", 1.0),)),),
    )
    assert result.blocked_actions == ("attack",)
    assert result.biomass[0, 0] == pytest.approx(1.0)
    assert result.budgets is not None
