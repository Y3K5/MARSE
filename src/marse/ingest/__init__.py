"""Read a table, check it, and write a feed an experiment file can use.

An importer does not run a simulation and does not decide that a table is a
measurement of a habitat. Abundance and transcript counts stay in the feed.
Only an environment row, or one metabolite already in mol per m3, is copied
into an experiment fragment.
"""

from marse.ingest.errors import IngestError
from marse.ingest.feed import Feed, apply_fragment
from marse.ingest.importers import (
    Importer,
    importer_for,
    register,
    registered_kinds,
    unregister,
)

__all__ = [
    "Feed",
    "Importer",
    "IngestError",
    "apply_fragment",
    "importer_for",
    "register",
    "registered_kinds",
    "unregister",
]
