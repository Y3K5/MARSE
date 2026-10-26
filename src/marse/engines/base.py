"""The engine protocol.

An engine validates a configuration, runs it, and replays a manifest of that
kind. :mod:`marse.engines.registry` selects one from the experiment file, so
the command line does not branch on file layout itself.

``engine``, when a file sets it, names the engine and is removed before the
scientific schema is read. It is not part of the configuration the manifest
records. Replay selects the engine from the manifest's ``kind``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from marse.core.provenance import Manifest


@dataclass(frozen=True, slots=True)
class Completed:
    """A finished run, plus any files already written for it."""

    result: Any
    written: tuple[Path, ...] = ()
    stored_frames: int | None = None
    frame_every: int | None = None


class Engine(Protocol):
    """One simulation engine, registered under :attr:`name`."""

    name: str
    kinds: tuple[str, ...]
    unverified: bool
    writes_during_replay: bool

    def load(self, path: Path) -> Any:
        """Read a file that does not name an engine."""

    def validate(self, document: dict[str, Any]) -> Any:
        """Read a document whose top-level ``engine`` key was already removed."""

    def default_output_dir(self, plan: Any, source: Path) -> Path:
        """Where a run writes when the caller does not pass ``--output``."""

    def run(
        self,
        plan: Any,
        output_dir: Path,
        *,
        frame_interval_h: float | None = None,
    ) -> Completed:
        """Run ``plan`` and write its outputs under ``output_dir``."""

    def replay(self, manifest: Manifest, output_dir: Path | None) -> Completed:
        """Re-run the manifest.

        When :attr:`writes_during_replay` is false, ``output_dir`` is ignored
        and the caller writes files only after the results match.
        """

    def write_outputs(self, result: Any, output_dir: Path) -> tuple[Path, ...]:
        """Write files after a replay that matched. Unused when writing already happened."""
