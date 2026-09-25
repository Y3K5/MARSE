"""Validated definitions of what a MARSE simulation is made of.

Configuration schema version 2 is built here, one part at a time
(docs/roadmap.md, Stage 2). The first part is the reaction network
(docs/networks.md): components with a chemical composition, and processes that
are proven, as they are read, to conserve carbon, nitrogen and electrons
exactly. The rest of the experiment (initial state, rates, transport and time)
follows as the engine that runs these networks is built.
"""

from marse.schemas.formula import QUANTITIES, Formula, parse_formula
from marse.schemas.network import (
    Component,
    ContinuityError,
    Growth,
    Network,
    Process,
    load_network,
    network_from_dict,
)

__all__ = [
    "QUANTITIES",
    "Component",
    "ContinuityError",
    "Formula",
    "Growth",
    "Network",
    "Process",
    "load_network",
    "network_from_dict",
    "parse_formula",
]
