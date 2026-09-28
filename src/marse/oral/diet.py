"""What is eaten and drunk: the intakes of a run, and the food they leave on the teeth.

A high-sugar eater differs from a low-sugar eater in more than the amount of
sugar. Sugar stays in the mouth for longer, sipped in drinks or sucked from
sweets, and food retained on the teeth keeps it concentrated where the plaque
is (Kashket, Zhang and Van Houte 1996). A diet lists intakes, each from a start
for a duration (docs/theory.md, section 4.9):

- **a rinse** adds its volume at once. The mouth holds it without swallowing,
  and at the end expels everything above its resting volume. The Stephan
  curve is the response to a sugar rinse;
- **a drink** flows in steadily over the duration, and is swallowed as the
  mouth fills;
- **a food** releases its amounts into the saliva steadily over the
  duration, as a sweet sucked slowly does, and adds no liquid.

While an intake is in the mouth, the film over the plaque is mixed with the
mouth's liquid (Dibdin 1990), at a rate the intake may state. Food it leaves
on the teeth is placed in the film when it ends, and the network's own
processes release what dissolves from it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.oral.film import film_layers
from marse.schemas.domain import Film, Intake
from marse.spatial.grid import Grid
from marse.spatial.surface import faces_within

__all__ = ["NOTHING", "Diet", "Event", "Inflow"]

_ML = 1e-6  # m3
_UM = 1e-6  # m
_MMOL = 1e-3  # mol


@dataclass(frozen=True, slots=True)
class Inflow:
    """What enters the mouth while an intake lasts, besides saliva.

    ``liquid_m3_per_s`` of a drink of ``liquid_mol_per_m3``;
    ``released_mol_per_s`` from a food, without liquid; ``held`` while a rinse
    is held, when the mouth does not swallow; and ``mixing_per_h``, how fast
    the film mixes with the mouth's liquid in each of its voxels.
    """

    liquid_m3_per_s: float = 0.0
    liquid_mol_per_m3: NDArray[np.float64] | None = None
    released_mol_per_s: NDArray[np.float64] | None = None
    held: bool = False
    mixing_per_h: float = 0.0

    def stimulus_mol_per_s(self, stimulus: int | None) -> float:
        """How fast the stimulus enters the mouth, which tastes it and secretes faster."""
        rate = 0.0
        if stimulus is None:
            return rate
        if self.liquid_mol_per_m3 is not None:
            rate += self.liquid_m3_per_s * float(self.liquid_mol_per_m3[stimulus])
        if self.released_mol_per_s is not None:
            rate += float(self.released_mol_per_s[stimulus])
        return rate


NOTHING = Inflow()
"""Between intakes, only saliva enters the mouth."""


@dataclass(frozen=True, slots=True)
class Event:
    """An intake's start or end."""

    time_h: float
    intake: Intake
    starts: bool


class Diet:
    """The intakes of a run as events in order, each start then its end.

    The reader has checked that the intakes come in order, one at a time, so
    the events do too. :meth:`due` moves through them as the run reaches
    them; :attr:`current` is the intake in the mouth between events.
    """

    def __init__(self, intakes: Sequence[Intake], names: tuple[str, ...]) -> None:
        self.names = names
        self._events = [
            Event(time, intake, starts)
            for intake in intakes
            for time, starts in ((intake.start_h, True), (intake.end_h, False))
        ]
        self._next = 0
        self.current: Intake | None = None
        self.taken = 0  # intakes started so far

    def next_h(self) -> float:
        """When the next event falls, or infinity after the last."""
        return self._events[self._next].time_h if self._next < len(self._events) else math.inf

    def due(self, now_h: float) -> list[Event]:
        """The events that fall by ``now_h``, to rounding, in order; the run applies them now."""
        tolerance = 1e-12 * max(1.0, now_h)
        due = []
        while self._next < len(self._events) and (
            self._events[self._next].time_h <= now_h + tolerance
        ):
            event = self._events[self._next]
            self.current = event.intake if event.starts else None
            self.taken += event.starts
            due.append(event)
            self._next += 1
        return due

    def vector(self, amounts: dict[str, float] | None) -> NDArray[np.float64]:
        """Amounts by name as a vector over every component, zero where unnamed."""
        given = amounts or {}
        return np.array([given.get(n, 0.0) for n in self.names])

    def inflow(self) -> Inflow:
        """What the intake in the mouth brings in, per second, or nothing between intakes."""
        intake = self.current
        if intake is None:
            return NOTHING
        seconds = intake.duration_min * 60.0
        mixing = intake.mixing_per_s * 3600.0
        if intake.kind == "rinse":
            return Inflow(held=True, mixing_per_h=mixing)
        if intake.kind == "drink":
            assert intake.volume_ml is not None
            return Inflow(
                liquid_m3_per_s=intake.volume_ml * _ML / seconds,
                liquid_mol_per_m3=self.vector(intake.composition_mol_per_m3),
                mixing_per_h=mixing,
            )
        return Inflow(
            released_mol_per_s=self.vector(intake.released_mmol) * _MMOL / seconds,
            mixing_per_h=mixing,
        )

    def rinse(self, intake: Intake) -> tuple[float, NDArray[np.float64]]:
        """A rinse's volume, in m3, and the amounts it brings, in mol."""
        assert intake.volume_ml is not None
        volume = intake.volume_ml * _ML
        return volume, self.vector(intake.composition_mol_per_m3) * volume

    def pocket(self, intake: Intake, grid: Grid, film: Film) -> NDArray[np.float64] | None:
        """The food an intake leaves on the teeth, as concentrations to add to the box.

        The amount per m² of the region is spread evenly through the film's
        depth above each face of the substratum the region covers.
        """
        retained = intake.retained
        if retained is None:
            return None
        column = np.zeros(grid.shape[-1])
        column[-film_layers(grid, film) :] = retained.amount_mol_per_m2 / (film.thickness_um * _UM)
        field = np.zeros((len(self.names), *grid.shape))
        faces = faces_within(grid, retained.region_um)
        field[self.names.index(retained.component)] = faces[..., None] * column
        return field
