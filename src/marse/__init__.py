"""MARSE: Microbial Adaptability Resource Simulation Engine.

An open, modular framework for reproducible spatial simulation of microbial
populations and biofilms. MARSE is under development and has no stable
release yet; docs/roadmap.md says what exists today, docs/validation.md what
has been verified, and the roadmap what comes next.
"""

__version__ = "0.1.0.dev0"

from marse.genotype import (
    CapabilityMapResult,
    GenotypeError,
    GenotypeRule,
    map_genotype_to_capabilities,
)

__all__ = [
    "CapabilityMapResult",
    "GenotypeError",
    "GenotypeRule",
    "__version__",
    "map_genotype_to_capabilities",
]
