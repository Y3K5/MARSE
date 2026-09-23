"""Deterministic simulation kernel.

The core owns everything that is not biology: canonical simulation state, the
clock and run loop, orchestration of providers, spatial bookkeeping, random
seed streams, checkpoints, serialization and provenance. It never encodes a
biological rule, so alternative scientific models can be compared on the same
infrastructure.

Planned modules: ``simulation`` (run loop, clocks, checkpoints), ``state``
(canonical state), ``config`` (validated configuration) and ``provenance``
(run manifests). See docs/architecture.md.
"""
