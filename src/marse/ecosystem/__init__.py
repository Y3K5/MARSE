"""Deterministic 2D ecosystem simulation and visualization frames.

Unverified. See the module docstring of :mod:`marse.ecosystem.model` and
docs/validation.md. The command line warns; importing this package does not,
so library tests can treat warnings as errors.
"""

from marse.ecosystem.model import (
    AdditiveConfig,
    ConditionConfig,
    EcosystemConfig,
    EcosystemError,
    EcosystemFrame,
    EcosystemResult,
    EcosystemState,
    NutrientConfig,
    PhenotypeConfig,
    SeedRegion,
    SpeciesConfig,
    load_experiment,
    run,
)
from marse.ecosystem.providers import EcosystemProviders, ExplicitTransportProvider, NoBoundary
from marse.ecosystem.viewer import write_viewer

__all__ = [
    "AdditiveConfig",
    "ConditionConfig",
    "EcosystemConfig",
    "EcosystemError",
    "EcosystemFrame",
    "EcosystemProviders",
    "EcosystemResult",
    "EcosystemState",
    "ExplicitTransportProvider",
    "NoBoundary",
    "NutrientConfig",
    "PhenotypeConfig",
    "SeedRegion",
    "SpeciesConfig",
    "load_experiment",
    "run",
    "write_viewer",
]
