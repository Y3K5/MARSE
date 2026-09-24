import numpy as np
import pytest

from marse.immune import (
    ImmuneError,
    ImmuneInteraction,
    MolecularNeutralizer,
    apply_immune_pressure,
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
