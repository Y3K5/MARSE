"""Registry of simulation engines.

``marse run`` asks :meth:`Registry.load` which engine a file uses. An explicit
top-level ``engine`` field wins. Otherwise the file is recognised the way
``marse run`` always has: a schema-version-2 document with a ``domain`` is
reactive transport, one without is well mixed, a biofilm block is the steady
one-dimensional solver, and anything else is the legacy batch kernel.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from marse.core.config import ConfigError
from marse.engines.base import Engine

__all__ = ["Registry", "registry"]


class Registry:
    """Engines keyed by name and by the manifest kind each one owns."""

    def __init__(self) -> None:
        self._by_name: dict[str, Engine] = {}
        self._by_kind: dict[str, Engine] = {}

    def register(self, engine: Engine) -> None:
        if engine.name in self._by_name:
            raise ValueError(f"engine {engine.name!r} is already registered")
        self._by_name[engine.name] = engine
        for kind in engine.kinds:
            if kind in self._by_kind:
                raise ValueError(f"manifest kind {kind!r} is already handled")
            self._by_kind[kind] = engine

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._by_name))

    def get(self, name: str) -> Engine:
        try:
            return self._by_name[name]
        except KeyError:
            known = ", ".join(self.names())
            raise ConfigError(
                f"engine {name!r} is not registered; known engines: {known}"
            ) from None

    def for_kind(self, kind: str) -> Engine:
        try:
            return self._by_kind[kind]
        except KeyError:
            raise ValueError(f"manifest kind {kind!r} has no registered engine") from None

    def load(self, path: Path) -> tuple[Engine, Any]:
        """Choose an engine for ``path`` and return it with a validated plan."""
        raw = _peek(path)
        if isinstance(raw, dict) and "engine" in raw:
            name = raw["engine"]
            if not isinstance(name, str) or not name.strip():
                raise ConfigError("engine: must be the name of a registered engine")
            engine = self.get(name.strip())
            document = {key: value for key, value in raw.items() if key != "engine"}
            return engine, engine.validate(document)
        engine = self.infer(raw)
        return engine, engine.load(path)

    def infer(self, raw: Any) -> Engine:
        """The engine ``marse run`` used before files could name one."""
        if isinstance(raw, dict) and "schema_version" in raw:
            if "domain" in raw:
                return self.get("reactive_transport")
            return self.get("well_mixed")
        if isinstance(raw, dict) and raw.get("biofilm") is not None:
            return self.get("steady_biofilm")
        return self.get("legacy_batch")


def _peek(path: Path) -> Any:
    """The parsed JSON object, or ``None`` when the regular loader should report the problem."""
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None


def _builtins() -> Registry:
    from marse.engines.legacy_batch import legacy_batch
    from marse.engines.legacy_ecosystem import legacy_ecosystem_v1
    from marse.engines.reactive_transport import reactive_transport
    from marse.engines.steady_biofilm import steady_biofilm
    from marse.engines.well_mixed import well_mixed

    engines = Registry()
    for engine in (
        legacy_batch,
        steady_biofilm,
        well_mixed,
        reactive_transport,
        legacy_ecosystem_v1,
    ):
        engines.register(engine)
    return engines


registry = _builtins()
