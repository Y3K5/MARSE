"""Measures of a Stephan curve: how low the pH falls, when, for how long, and when it is back.

After sugar, the pH of plaque falls within minutes to a minimum and returns
over half an hour or more (Stephan 1944). Enamel dissolves below a critical pH
of about 5.5, so how long plaque stays below it, and how far, are the usual
measures of a sugar challenge. Each measure here takes the curve as the
straight lines between its records, and is exact for that curve: a crossing
of a level between two records is found where the line crosses it.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["CRITICAL_PH", "StephanCurve", "area_below", "back_above", "minutes_below"]

CRITICAL_PH = 5.5
"""The pH below which enamel dissolves, roughly: the level acid exposure is measured from."""


def _series(times_h: Sequence[float], ph: Sequence[float]) -> tuple[NDArray, NDArray]:
    t = np.asarray(times_h, dtype=float) * 60.0
    p = np.asarray(ph, dtype=float)
    if t.shape != p.shape or t.ndim != 1 or t.size < 2:
        raise ValueError("a curve needs matching times and pH values, at least two of each")
    if np.any(np.diff(t) <= 0):
        raise ValueError("the times of a curve must increase")
    return t, p


def minutes_below(times_h: Sequence[float], ph: Sequence[float], level: float) -> float:
    """How many minutes the curve spends below ``level``."""
    t, p = _series(times_h, ph)
    total = 0.0
    for t0, t1, p0, p1 in zip(t[:-1], t[1:], p[:-1], p[1:], strict=True):
        low, high = min(p0, p1), max(p0, p1)
        if high <= level:
            total += t1 - t0
        elif low < level:
            total += (t1 - t0) * (level - low) / (high - low)
    return total


def area_below(times_h: Sequence[float], ph: Sequence[float], level: float) -> float:
    """The area between ``level`` and the curve where the curve is below it, in pH x minutes."""
    t, p = _series(times_h, ph)
    total = 0.0
    for t0, t1, p0, p1 in zip(t[:-1], t[1:], p[:-1], p[1:], strict=True):
        d0, d1 = level - p0, level - p1  # depths below the level
        if d0 >= 0 and d1 >= 0:
            total += 0.5 * (d0 + d1) * (t1 - t0)
        elif d0 > 0 or d1 > 0:  # the line crosses the level: a triangle below it
            deep = max(d0, d1)
            total += 0.5 * deep * (t1 - t0) * deep / abs(d1 - d0)
    return total


def back_above(times_h: Sequence[float], ph: Sequence[float], level: float) -> float | None:
    """When, in minutes, the curve first returns to ``level`` after its minimum; None if never."""
    t, p = _series(times_h, ph)
    lowest = int(np.argmin(p))
    for i in range(lowest, t.size - 1):
        if p[i] < level <= p[i + 1]:
            return float(t[i] + (t[i + 1] - t[i]) * (level - p[i]) / (p[i + 1] - p[i]))
    return None


@dataclass(frozen=True, slots=True)
class StephanCurve:
    """A curve's start, its minimum and when, its exposure below the critical pH, and its return."""

    start: float
    minimum: float
    minimum_at_min: float
    minutes_below_critical: float
    area_below_critical: float
    back_above_6_min: float | None

    @classmethod
    def of(cls, times_h: Sequence[float], ph: Sequence[float]) -> StephanCurve:
        t, p = _series(times_h, ph)
        lowest = int(np.argmin(p))
        return cls(
            start=float(p[0]),
            minimum=float(p[lowest]),
            minimum_at_min=float(t[lowest]),
            minutes_below_critical=minutes_below(times_h, ph, CRITICAL_PH),
            area_below_critical=area_below(times_h, ph, CRITICAL_PH),
            back_above_6_min=back_above(times_h, ph, 6.0),
        )
