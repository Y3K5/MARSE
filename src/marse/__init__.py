"""MARSE: Microbial Adaptability Resource Simulation Engine.

An open, modular framework for reproducible spatial simulation of microbial
populations and biofilms. It is under development and has no stable release yet;
docs/roadmap.md says what exists today, docs/validation.md what has been
verified, and the roadmap what comes next.
"""

__version__ = "0.1.0.dev0"

from marse._compat import install as _install_compat
from marse.processes.genotype import (
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

_install_compat()
