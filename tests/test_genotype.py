import pytest

from marse.genotype import GenotypeError, GenotypeRule, map_genotype_to_capabilities
from marse.niche import Capability


def test_declared_allele_modifies_only_matching_capability():
    capabilities = (
        Capability("respiration", 2.0, "carbon", 0.5),
        Capability("fermentation", 1.0, "carbon", 0.2),
    )
    result = map_genotype_to_capabilities(
        capabilities,
        {"oxyR": "loss_of_function"},
        (
            GenotypeRule(
                "oxyR",
                "loss_of_function",
                "respiration",
                maximum_rate_multiplier=0.25,
                evidence_source="paper:example",
            ),
        ),
    )
    assert result.capabilities[0].maximum_rate_per_h == pytest.approx(0.5)
    assert result.capabilities[1] == capabilities[1]
    assert result.applied_rules[0].evidence_source == "paper:example"


def test_unknown_capability_rule_is_rejected_only_when_allele_matches():
    with pytest.raises(GenotypeError, match="unknown capability"):
        map_genotype_to_capabilities(
            (Capability("respiration", 1.0, "carbon", 0.5),),
            {"gene": "A"},
            (GenotypeRule("gene", "A", "missing"),),
        )
