"""MARSE: Microbial Adaptability Resource Engine.

An open, modular framework for reproducible spatial simulation of microbial
populations and biofilms. MARSE is in Phase 0 (specification); see
docs/roadmap.md for what exists today and what comes next.
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
