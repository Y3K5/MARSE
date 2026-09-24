"""End-to-end composition test for genotype mapping and immune pressure."""

import numpy as np

from marse.ecosystem import AdditiveConfig, EcosystemConfig, NutrientConfig, SpeciesConfig, run
from marse.genotype import GenotypeRule, map_genotype_to_capabilities
from marse.immune import ImmuneInteraction
from marse.niche import Capability


def test_resistance_allele_changes_ecosystem_survival_under_effector_pressure():
    base_capability = Capability("efflux", 1.0, "food", 0.2)
    mapped = map_genotype_to_capabilities(
        (base_capability,),
        {"resistance_locus": "R"},
        (GenotypeRule("resistance_locus", "R", "efflux", maximum_rate_multiplier=4.0),),
    )
    resistance_factor = (
        mapped.capabilities[0].maximum_rate_per_h / base_capability.maximum_rate_per_h
    )
    species = (
        SpeciesConfig("susceptible", 0.1, 0.0, (0.2,), (1.0,)),
        SpeciesConfig("resistant", 0.1, 0.0, (0.2,), (1.0,)),
    )
    result = run(
        EcosystemConfig(
            experiment_id="genotype-immune-integration",
            width=4,
            height=4,
            cell_size_um=10.0,
            duration_h=0.2,
            timestep_h=0.01,
            seed=11,
            nutrients=(NutrientConfig("food", 1.0, 0.0),),
            additives=(AdditiveConfig("effector", 2.0, 0.0),),
            species=species,
            immune_interactions=(
                ImmuneInteraction("susceptible", "effector", 2.0, 1.0, susceptibility=1.0),
                ImmuneInteraction(
                    "resistant",
                    "effector",
                    2.0,
                    1.0,
                    susceptibility=1.0 / resistance_factor,
                ),
            ),
        )
    )
    susceptible, resistant = result.final_state.biomass.sum(axis=(1, 2))
    assert resistant > susceptible
    assert np.all(np.isfinite(result.final_state.biomass))
