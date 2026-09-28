"""The mouth: its saliva, the film over the teeth, and what is eaten (docs/environments.md)."""

from marse.oral.diet import NOTHING, Diet, Event, Inflow
from marse.oral.film import film_layers, renewal_per_h
from marse.oral.mouth import OralFluid, Stretch

__all__ = [
    "NOTHING",
    "Diet",
    "Event",
    "Inflow",
    "OralFluid",
    "Stretch",
    "film_layers",
    "renewal_per_h",
]
