"""Run manifests: everything needed to replay a simulation, and nothing more.

A manifest is meant to be attached to a paper, pasted into an issue and
archived, so it must identify the *run* without identifying the person or
machine that produced it. The rules in docs/architecture.md are enforced here
rather than left to good intentions:

* no usernames, hostnames, IP addresses or environment variables;
* no absolute paths — the experiment file is recorded by name only;
* platform limited to the operating-system family;
* timestamps in UTC.

MARSE makes no network requests and collects no telemetry. Manifests are JSON
with sorted keys so that the same run produces byte-identical output.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from marse import __version__
from marse.core.config import experiment_from_dict

__all__ = ["Manifest", "RunConfig", "config_checksum"]

MANIFEST_VERSION = 2
"""Schema version of the manifest itself, so old manifests stay readable.

Version 2 added ``kind``. Version 1 manifests predate the ecosystem engine, so
their kind is recovered from the configuration when they are read.
"""

_READABLE_VERSIONS = (1, 2)

_CORE_KINDS = ("batch", "biofilm_profile")


class RunConfig(Protocol):
    """What a configuration must offer to be recorded in, and rebuilt from, a manifest."""

    experiment_id: str
    seed: int

    @property
    def kind(self) -> str: ...

    def to_dict(self) -> dict[str, Any]: ...


def _canonical(data: Any) -> str:
    """Deterministic JSON: sorted keys, compact separators, no NaN."""
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)


def config_checksum(config: RunConfig) -> str:
    """SHA-256 over the canonical form of the configuration."""
    return hashlib.sha256(_canonical(config.to_dict()).encode("utf-8")).hexdigest()


def _load_config(kind: str, raw: dict[str, Any]) -> RunConfig:
    if kind in _CORE_KINDS:
        return experiment_from_dict(raw)
    if kind == "ecosystem":
        # Imported here rather than at module level: the core does not depend on
        # the ecosystem engine, and only an ecosystem manifest needs it.
        from marse.ecosystem.model import ecosystem_from_dict

        return ecosystem_from_dict(raw)
    if kind == "well_mixed":
        from marse.schemas.experiment import experiment_from_dict as network_experiment

        return network_experiment(raw)
    raise ValueError(f"manifest kind {kind!r} is not supported by MARSE {__version__}")


def _environment_record() -> dict[str, str]:
    """Only what affects numerical results. Deliberately not the machine's identity."""
    return {
        "python": platform.python_version(),
        "implementation": sys.implementation.name,
        "os_family": platform.system(),  # "Linux", "Darwin", "Windows" - not the hostname
        "numpy": np.__version__,
        "marse": __version__,
    }


@dataclass(frozen=True, slots=True)
class Manifest:
    """The record of one run, sufficient to reproduce it."""

    run_id: str
    experiment_id: str
    config: dict[str, Any]
    config_sha256: str
    seed: int
    models: dict[str, str]
    random_streams: tuple[str, ...]
    environment: dict[str, str]
    started_utc: str
    finished_utc: str
    steps: int
    outputs: dict[str, Any]
    kind: str
    manifest_version: int = MANIFEST_VERSION

    @classmethod
    def build(
        cls,
        *,
        config: RunConfig,
        models: dict[str, str],
        random_streams: tuple[str, ...],
        started: datetime,
        finished: datetime,
        steps: int,
        outputs: dict[str, Any],
    ) -> Manifest:
        checksum = config_checksum(config)
        return cls(
            # Derived from the configuration, so the same experiment yields the
            # same run_id: identifying, but not a record of when or where.
            run_id=f"MARSE-{checksum[:12]}",
            experiment_id=config.experiment_id,
            config=config.to_dict(),
            config_sha256=checksum,
            seed=config.seed,
            models=dict(sorted(models.items())),
            random_streams=tuple(sorted(random_streams)),
            environment=_environment_record(),
            started_utc=started.astimezone(UTC).isoformat(timespec="seconds"),
            finished_utc=finished.astimezone(UTC).isoformat(timespec="seconds"),
            steps=steps,
            outputs=outputs,
            kind=config.kind,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "manifest_version": self.manifest_version,
            "kind": self.kind,
            "run_id": self.run_id,
            "experiment_id": self.experiment_id,
            "seed": self.seed,
            "config_sha256": self.config_sha256,
            "models": self.models,
            "random_streams": list(self.random_streams),
            "environment": self.environment,
            "started_utc": self.started_utc,
            "finished_utc": self.finished_utc,
            "steps": self.steps,
            "outputs": self.outputs,
            "config": self.config,
        }

    def write(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return destination

    @classmethod
    def read(cls, path: str | Path) -> Manifest:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        version = raw.get("manifest_version")
        if version not in _READABLE_VERSIONS:
            raise ValueError(
                f"manifest version {version!r} is not supported by MARSE {__version__} "
                f"(readable: {', '.join(map(str, _READABLE_VERSIONS))})"
            )
        if version == 1:
            # Version 1 had no kind field and only the two core kinds existed.
            kind = "biofilm_profile" if raw["config"].get("biofilm") is not None else "batch"
        else:
            kind = raw["kind"]
        return cls(
            run_id=raw["run_id"],
            experiment_id=raw["experiment_id"],
            config=raw["config"],
            config_sha256=raw["config_sha256"],
            seed=raw["seed"],
            models=raw["models"],
            random_streams=tuple(raw["random_streams"]),
            environment=raw["environment"],
            started_utc=raw["started_utc"],
            finished_utc=raw["finished_utc"],
            steps=raw["steps"],
            outputs=raw["outputs"],
            kind=kind,
            manifest_version=version,
        )

    def experiment(self) -> RunConfig:
        """Rebuild the configuration, verifying it still matches its checksum.

        A mismatch means the manifest was edited after the run, so replaying it
        would not reproduce the recorded results. The rebuilt configuration must
        also be of the kind the manifest declares.
        """
        config = _load_config(self.kind, self.config)
        if config.kind != self.kind:
            raise ValueError(
                f"manifest declares a {self.kind!r} run but its configuration "
                f"describes a {config.kind!r} run"
            )
        actual = config_checksum(config)
        if actual != self.config_sha256:
            raise ValueError(
                "manifest config does not match its checksum: the manifest has been "
                f"modified (recorded {self.config_sha256[:12]}, computed {actual[:12]})"
            )
        return config
