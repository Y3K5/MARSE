"""Seeded random streams owned by the core.

Every source of randomness in MARSE draws from a stream derived here, so that a
run is reproducible from its seed alone.

The streams are keyed by **name**, not by the order providers were registered.
This matters more than it looks: with order-derived streams, adding a
diagnostic provider shifts every subsequent stream and silently changes the
trajectory of an otherwise identical run. Name-keyed streams are stable under
adding, removing and reordering providers (docs/theory.md, section 9.4).
"""

from __future__ import annotations

import hashlib

import numpy as np

__all__ = ["SeedRegistry"]


def _key_for(name: str) -> int:
    """A stable 64-bit key from a provider name.

    ``hash()`` is salted per process and would break reproducibility across
    runs, so a cryptographic digest is used instead.
    """
    digest = hashlib.blake2b(name.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big")


class SeedRegistry:
    """Hands out an independent generator per named consumer."""

    def __init__(self, seed: int) -> None:
        if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
            raise ValueError("seed must be a non-negative integer")
        self._seed = seed
        self._issued: dict[str, np.random.Generator] = {}

    @property
    def seed(self) -> int:
        return self._seed

    @property
    def issued_names(self) -> tuple[str, ...]:
        """Names that have requested a stream, sorted for a stable manifest."""
        return tuple(sorted(self._issued))

    def stream(self, name: str) -> np.random.Generator:
        """Return the generator for ``name``, creating it on first request.

        The same name always yields the same stream within a run, and the same
        stream across runs with the same seed.
        """
        if not name:
            raise ValueError("stream name must not be empty")
        if name not in self._issued:
            sequence = np.random.SeedSequence(entropy=self._seed, spawn_key=(_key_for(name),))
            self._issued[name] = np.random.default_rng(sequence)
        return self._issued[name]
