"""The mouth's liquid: the flow that fills it and the swallows that empty it.

Dawes's (1983) model of how the mouth clears sugar:

- the liquid in the mouth grows from its resting volume, RESID, at the
  salivary flow, until it reaches VMAX, when the person swallows
  (Lagerlof and Dawes 1984 measured the two volumes);
- a swallow is an incomplete syphon: it leaves RESID of the same liquid, so a
  swallow changes the volume at once but no concentration;
- the flow is the unstimulated flow plus a part that the taste of a stimulus,
  such as sugar, adds: Q = Q_u + Q_s c / (K + c). Chewing adds a flow of its
  own, Q_c, while it lasts (Dawes and Macpherson 1992);
- the glands secrete resting saliva at the unstimulated flow, and saliva
  closer to stimulated saliva as the flow rises; bicarbonate, for one, rises
  steeply with flow (Bardow et al. 2000).

RESID and VMAX count the liquid on every surface of the mouth, the film over
the modelled plaque included (Collins and Dawes 1987). That film is the top of
the box, so the pool the box exchanges with is the mouth's liquid less the
film. Its composition is solved with the box (:mod:`marse.core.reservoir`).
The mouth runs ahead of the box to the end of each span, or to the next
swallow, with the stimulus it holds at the start, diluted by the saliva it
secretes and added to by what is eaten or drunk (:mod:`marse.oral.diet`):
that sets its volume over the span. It expects the plaque to go on taking up,
or giving back, the stimulus at the rate it did over the span before; a
change in that rate reaches the flow from the next span on, and the swallows
it moves are booked exactly wherever they fall.

A rinse is held without swallowing, whatever the volume, and expelled at its
end down to RESID; like a swallow, that changes no concentration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from marse.oral.diet import NOTHING, Inflow
from marse.schemas.domain import Film, Mouth

__all__ = ["OralFluid", "Stretch"]

_ML = 1e-6  # m3
_UM = 1e-6  # m
_RESOLUTION_S = 0.25  # the run-ahead's step: the flow changes over seconds


@dataclass(frozen=True, slots=True)
class Stretch:
    """What the mouth's volume does over one span: samples of it, and how the span ends."""

    times_s: NDArray[np.float64]  # from the span's start
    volumes_m3: NDArray[np.float64]
    flows_m3_per_s: NDArray[np.float64]  # the volume's rate of growth at each sample
    seconds: float  # how long the span lasted: to its end, or to a swallow
    swallowed: bool


class OralFluid:
    """The volume of liquid in the mouth, the flow that fills it, and the swallows."""

    def __init__(self, mouth: Mouth, film: Film, names: tuple[str, ...]) -> None:
        self.area_m2 = mouth.plaque_area_cm2 * 1e-4
        self.film_m3 = self.area_m2 * film.thickness_um * _UM
        self.resting_m3 = mouth.resting_volume_ml * _ML
        self.full_m3 = mouth.swallow_volume_ml * _ML
        self.unstimulated_m3_per_s = mouth.unstimulated_flow_ml_per_min * _ML / 60.0
        self.stimulated_m3_per_s = mouth.stimulated_flow_ml_per_min * _ML / 60.0
        self.chewing_m3_per_s = mouth.chewing_flow_ml_per_min * _ML / 60.0
        self.stimulus = None if mouth.stimulus is None else names.index(mouth.stimulus)
        self.half = mouth.stimulus_half_mol_per_m3 or 0.0

        def full(composition: dict[str, float] | None) -> NDArray[np.float64] | None:
            if composition is None:
                return None
            return np.array([composition.get(n, 0.0) for n in names])

        resting = full(mouth.saliva_mol_per_m3)
        assert resting is not None
        self.saliva = resting
        self.stimulated_saliva = full(mouth.stimulated_saliva_mol_per_m3)
        self.volume_m3 = self.resting_m3
        self.swallows = 0

    # -- geometry -------------------------------------------------------------------------

    def thickness_um(self, volume_m3: float) -> float:
        """The pool's volume per unit area of the modelled plaque, in um."""
        return (volume_m3 - self.film_m3) / self.area_m2 / _UM

    @property
    def reference_um(self) -> float:
        """The pool's thickness at the resting volume: H_ref, which scales its unknowns."""
        return self.thickness_um(self.resting_m3)

    def growth_um_per_h(self, flow_m3_per_s: float) -> float:
        return flow_m3_per_s / self.area_m2 / _UM * 3600.0

    # -- flow and secretion ---------------------------------------------------------------

    def flow_m3_per_s(self, stimulus_mol_per_m3: float, chewing: bool = False) -> float:
        """The salivary flow at a concentration of the stimulus in the pool, and while chewing."""
        flow = self.unstimulated_m3_per_s
        if chewing:
            flow += self.chewing_m3_per_s
        if self.stimulus is None or self.stimulated_m3_per_s == 0.0:
            return flow
        c = max(stimulus_mol_per_m3, 0.0)
        return flow + self.stimulated_m3_per_s * c / (self.half + c)

    def secreted(self, flow_m3_per_s: float) -> NDArray[np.float64]:
        """The composition of the saliva secreted at a flow: resting, towards stimulated.

        It reaches stimulated saliva at the most the stimulus adds to the flow,
        or, in a mouth without one, at the chewing flow.
        """
        stimulated = self.stimulated_m3_per_s or self.chewing_m3_per_s
        if self.stimulated_saliva is None or stimulated == 0.0:
            return self.saliva
        share = (flow_m3_per_s - self.unstimulated_m3_per_s) / stimulated
        share = min(max(share, 0.0), 1.0)
        return (1.0 - share) * self.saliva + share * self.stimulated_saliva

    # -- a span ---------------------------------------------------------------------------

    def run_ahead(
        self,
        stimulus_mol: float,
        most_s: float,
        inflow: Inflow = NOTHING,
        returned_mol_per_s: float = 0.0,
    ) -> Stretch:
        """The volume from now until ``most_s`` seconds, or until the mouth is full.

        Classical Runge-Kutta in steps of a quarter of a second, the flow
        following the stimulus as secretion dilutes it and an intake adds to
        it; a drink adds its own flow. ``returned_mol_per_s`` is how fast the
        plaque gave the stimulus back over the span before, if it did, which
        the mouth expects it to go on doing. A step that would overfill the
        mouth is shortened to end at the swallow, and the volume there is set
        to exactly VMAX. A rinse is held: the mouth does not swallow it.
        """
        tasted = inflow.stimulus_mol_per_s(self.stimulus) + returned_mol_per_s
        drink = inflow.liquid_m3_per_s

        def flow(t: float, volume: float) -> float:
            held = stimulus_mol + tasted * t
            return self.flow_m3_per_s(held / (volume - self.film_m3), inflow.chewing) + drink

        volume, t = self.volume_m3, 0.0
        times, volumes, flows = [0.0], [volume], [flow(t, volume)]
        swallowed = False
        while t < most_s * (1.0 - 1e-12):
            dt = min(_RESOLUTION_S, most_s - t)
            rate = flow(t, volume)
            to_full = (self.full_m3 - volume) / rate
            # Full within this step, or within rounding of its end: swallow, landing no later
            # than the span's end, so a swallow on a span's last instant is never put off.
            if not inflow.held and to_full <= dt * (1.0 + 1e-9):
                dt, swallowed = min(max(to_full, 0.0), most_s - t), True
            k1 = rate
            k2 = flow(t + 0.5 * dt, volume + 0.5 * dt * k1)
            k3 = flow(t + 0.5 * dt, volume + 0.5 * dt * k2)
            k4 = flow(t + dt, volume + dt * k3)
            volume = volume + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            if swallowed:
                volume = self.full_m3
            times.append(t)
            volumes.append(volume)
            flows.append(flow(t, volume))
            if swallowed:
                break
        return Stretch(np.array(times), np.array(volumes), np.array(flows), t, swallowed)

    def end(self, stretch: Stretch) -> float:
        """Move to the end of a span; after a swallow, the share of the pool that stays."""
        if not stretch.swallowed:
            self.volume_m3 = float(stretch.volumes_m3[-1])
            return 1.0
        kept = self.thickness_um(self.resting_m3) / self.thickness_um(self.full_m3)
        self.volume_m3 = self.resting_m3
        self.swallows += 1
        return kept

    def take(self, volume_m3: float) -> None:
        """Take a rinse into the mouth: its volume adds to the liquid at once."""
        self.volume_m3 += volume_m3

    def expel(self) -> float:
        """Expel everything above the resting volume; the share of the pool that stays."""
        kept = self.thickness_um(self.resting_m3) / self.thickness_um(self.volume_m3)
        self.volume_m3 = self.resting_m3
        return kept
