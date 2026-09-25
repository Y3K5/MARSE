"""Explicit, evidence-bounded genotype-to-capability mappings.

This module maps declared alleles to parameter modifiers. It does not infer
phenotypes from arbitrary DNA sequences or claim a universal genotype model.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from marse.microbes.niche import Capability

__all__ = [
    "CapabilityMapResult",
    "GenotypeError",
    "GenotypeRule",
    "map_genotype_to_capabilities",
]


class GenotypeError(ValueError):
    """A declared genotype mapping is invalid."""


@dataclass(frozen=True, slots=True)
class GenotypeRule:
    """Modify one capability when a locus has a declared allele."""

    locus: str
    allele: str
    capability_id: str
    maximum_rate_multiplier: float = 1.0
    half_saturation_multiplier: float = 1.0
    evidence_source: str = ""

    def __post_init__(self) -> None:
        if not self.locus.strip() or not self.allele.strip() or not self.capability_id.strip():
            raise GenotypeError("genotype rule identifiers must not be empty")
        if self.maximum_rate_multiplier < 0 or self.half_saturation_multiplier <= 0:
            raise GenotypeError(
                "genotype multipliers must be non-negative and half-saturation positive"
            )


@dataclass(frozen=True, slots=True)
class CapabilityMapResult:
    """Mapped capabilities plus the rules that were actually applied."""

    capabilities: tuple[Capability, ...]
    applied_rules: tuple[GenotypeRule, ...]


def map_genotype_to_capabilities(
    capabilities: tuple[Capability, ...],
    genotype: dict[str, str],
    rules: tuple[GenotypeRule, ...],
) -> CapabilityMapResult:
    """Apply matching declared rules without inventing unobserved traits."""
    by_id = {capability.id: capability for capability in capabilities}
    if len(by_id) != len(capabilities):
        raise GenotypeError("capability IDs must be unique")
    mapped = dict(by_id)
    applied: list[GenotypeRule] = []
    for rule in rules:
        if genotype.get(rule.locus) != rule.allele:
            continue
        if rule.capability_id not in mapped:
            raise GenotypeError(f"rule references unknown capability '{rule.capability_id}'")
        capability = mapped[rule.capability_id]
        mapped[rule.capability_id] = replace(
            capability,
            maximum_rate_per_h=capability.maximum_rate_per_h * rule.maximum_rate_multiplier,
            half_saturation=capability.half_saturation * rule.half_saturation_multiplier,
        )
        applied.append(rule)
    return CapabilityMapResult(
        tuple(mapped[capability.id] for capability in capabilities), tuple(applied)
    )
