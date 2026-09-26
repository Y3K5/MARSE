"""Validated definitions of what a MARSE simulation is made of.

Configuration schema version 2 is built here, one part at a time
(docs/roadmap.md, Stage 2). The first part is the reaction network
(docs/networks.md): components with a chemical composition, and processes that
are proven, as they are read, to conserve carbon, nitrogen and electrons
exactly. Rates, initial amounts and a clock make a network runnable
(:mod:`marse.schemas.experiment`); transport follows with the spatial engine.
"""

from marse.schemas.experiment import WellMixedConfig, experiment_from_dict
from marse.schemas.formula import QUANTITIES, Formula, parse_formula
from marse.schemas.network import (
    Component,
    ContinuityError,
    Factor,
    Growth,
    Network,
    Process,
    RateLaw,
    load_network,
    network_from_dict,
)

__all__ = [
    "QUANTITIES",
    "Component",
    "ContinuityError",
    "Factor",
    "Formula",
    "Growth",
    "Network",
    "Process",
    "RateLaw",
    "WellMixedConfig",
    "experiment_from_dict",
    "load_network",
    "network_from_dict",
    "parse_formula",
]
