"""Strict reading of configuration data: every key known, every number exact.

Schema version 2 refuses what version 1 let through. An unknown key is an
error that names it and, when a known key is close, the key that was probably
meant. A key whose unit differs from the declared one is pointed at the right
name, because the unit is part of the name. A number is read as the exact
decimal written in the file, so balances computed from it are exact.

Each object in the schema is described by a mapping from field name to
:class:`Field`. The same mappings drive the reader, a test that every numeric
field names its unit, and a test that docs/networks.md documents exactly these
fields.
"""

from __future__ import annotations

import difflib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Literal

from marse.core.config import ConfigError

UNIT_SUFFIXES = {
    "_mol_per_mol": "mol of one component per mol of another",
    "_mol_per_m3": "mol per cubic metre, which is mmol per litre",
    "_per_h": "per hour",
    "_h": "hours",
}
"""Every numeric field name ends in one of these suffixes, which names its unit.

A name is matched against the longest suffix first, so ``maximum_per_h`` is
per hour, not hours.
"""

DIMENSIONLESS = {
    "schema_version": "the version number of the configuration format",
    "charge": "the electric charge of one formula unit, in elementary charges",
    "seed": "the seed of the random streams, an identifier rather than a quantity",
    "relative_tolerance": "a fraction: the accepted local error relative to each concentration",
}
"""Numeric fields that are labels or pure numbers, so they carry no unit suffix."""

Kind = Literal[
    "text", "name", "names", "choice", "integer", "number", "numbers", "object", "objects"
]
NUMERIC_KINDS = frozenset({"integer", "number", "numbers"})

_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,63}")


@dataclass(frozen=True, slots=True)
class Field:
    """How one key of a configuration object is read.

    ``text`` is any string, ``name`` an identifier, ``names`` a list of
    distinct identifiers, ``choice`` one of ``choices``, ``integer`` a whole
    number, ``number`` an exact decimal, ``numbers`` an object from names to
    exact decimals, ``object`` an object read by the caller, and ``objects``
    a list of objects read by the caller.
    """

    kind: Kind
    required: bool = True
    choices: tuple[str, ...] = ()


def _describe(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return f"the text '{value}'"
    if isinstance(value, list):
        return "a list"
    if isinstance(value, dict):
        return "an object"
    return type(value).__name__


def unit_suffix(name: str) -> str | None:
    """The unit suffix a field name ends in, longest first, or None."""
    return next((s for s in sorted(UNIT_SUFFIXES, key=len, reverse=True) if name.endswith(s)), None)


def _hint(key: str, known: Mapping[str, Field]) -> str:
    for name in known:
        suffix = unit_suffix(name)
        base = name[: -len(suffix)] if suffix else None
        if base and (key == base or key.startswith(base + "_")):
            unit = UNIT_SUFFIXES[suffix]
            return f"; the unit is part of the name, and this field is '{name}', in {unit}"
    close = difflib.get_close_matches(key, list(known), n=1, cutoff=0.75)
    return f"; did you mean '{close[0]}'?" if close else ""


def read_object(raw: Any, where: str, schema: Mapping[str, Field]) -> dict[str, Any]:
    """Check an object's keys against ``schema`` and convert every value present."""
    if not isinstance(raw, dict):
        raise ConfigError(f"{where}: expected an object, got {_describe(raw)}")
    unknown = [key for key in raw if key not in schema]
    if unknown:
        details = ", ".join(f"'{key}'{_hint(str(key), schema)}" for key in unknown)
        raise ConfigError(f"{where}: unknown field {details}")
    missing = [name for name, field in schema.items() if field.required and name not in raw]
    if missing:
        listed = ", ".join(f"'{m}'" for m in missing)
        raise ConfigError(f"{where}: missing required field {listed}")
    return {key: _convert(raw[key], f"{where}.{key}", schema[key]) for key in raw}


def _convert(value: Any, where: str, field: Field) -> Any:
    match field.kind:
        case "text":
            if not isinstance(value, str):
                raise ConfigError(f"{where}: expected text, got {_describe(value)}")
            return value
        case "name":
            return name(value, where)
        case "names":
            if not isinstance(value, list):
                raise ConfigError(f"{where}: expected a list of names, got {_describe(value)}")
            names = tuple(name(item, f"{where}[{i}]") for i, item in enumerate(value))
            repeated = sorted({n for n in names if names.count(n) > 1})
            if repeated:
                raise ConfigError(f"{where}: {', '.join(repeated)} listed more than once")
            return names
        case "choice":
            if value not in field.choices:
                options = ", ".join(f"'{c}'" for c in field.choices)
                raise ConfigError(f"{where}: expected one of {options}, got {_describe(value)}")
            return value
        case "integer":
            if isinstance(value, bool) or not isinstance(value, int):
                raise ConfigError(f"{where}: expected a whole number, got {_describe(value)}")
            return value
        case "number":
            return exact(value, where)
        case "numbers":
            if not isinstance(value, dict):
                raise ConfigError(f"{where}: expected an object of numbers, got {_describe(value)}")
            return {name(k, f"{where} key"): exact(v, f"{where}.{k}") for k, v in value.items()}
        case "object":
            if not isinstance(value, dict):
                raise ConfigError(f"{where}: expected an object, got {_describe(value)}")
            return value
        case "objects":
            if not isinstance(value, list):
                raise ConfigError(f"{where}: expected a list of objects, got {_describe(value)}")
            return value
    raise AssertionError(f"unhandled field kind {field.kind!r}")  # pragma: no cover


def exact(value: Any, where: str) -> Fraction:
    """A JSON number as the exact decimal it was written as.

    ``Fraction(repr(x))`` recovers the shortest decimal that reads back as the
    same float, which is the decimal in the file: 0.1 becomes exactly 1/10,
    not the binary value nearest to it.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{where}: expected a number, got {_describe(value)}")
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ConfigError(f"{where}: must be a finite number, got {value}")
        return Fraction(repr(value))
    return Fraction(value)


def name(value: Any, where: str) -> str:
    """An identifier: it becomes a JSON key and, in outputs, part of a file name."""
    if not isinstance(value, str) or not _NAME.fullmatch(value):
        raise ConfigError(
            f"{where}: {_describe(value)} is not a valid name; use letters, digits and "
            "underscores, starting with a letter, at most 64 characters"
        )
    return value


def plain(value: Fraction) -> int | float:
    """An exact number back as JSON: whole numbers as integers, others as floats.

    The float's shortest decimal is the one :func:`exact` read, so reading the
    output again gives back exactly ``value``.
    """
    return value.numerator if value.denominator == 1 else float(value)


def _refuse_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    keys = [key for key, _ in pairs]
    repeated = sorted({key for key in keys if keys.count(key) > 1})
    if repeated:
        listed = ", ".join(f"'{k}'" for k in repeated)
        raise ConfigError(
            f"key {listed} appears twice in one object; JSON would keep only the last, "
            "silently dropping the first"
        )
    return dict(pairs)


def load_json(path: str | Path) -> Any:
    """A configuration file's JSON, refusing duplicate keys and naming the file on error."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    try:
        return json.loads(text, object_pairs_hook=_refuse_duplicate_keys)
    except json.JSONDecodeError as error:
        raise ConfigError(
            f"{path.name}: invalid JSON ({error.msg} at line {error.lineno})"
        ) from error
