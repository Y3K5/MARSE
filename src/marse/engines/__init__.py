"""Simulation engines selected by experiment configuration.

Five engines are registered:

- ``legacy_batch``, manifest kind ``batch``
- ``steady_biofilm``, manifest kind ``biofilm_profile``
- ``well_mixed``, manifest kind ``well_mixed``
- ``reactive_transport``, manifest kind ``reactive_transport``
- ``legacy_ecosystem_v1``, manifest kind ``ecosystem`` (unverified; see
  docs/validation.md)

A file may set a top-level ``engine`` to one of those names. The field is not
part of the scientific configuration and is not written into the manifest.
"""

from marse.engines.registry import registry

__all__ = ["registry"]
