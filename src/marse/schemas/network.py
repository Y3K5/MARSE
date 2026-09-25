"""Reaction networks: components with a composition, processes that conserve it.

This is the first part of MARSE's configuration schema version 2
(docs/networks.md). A network lists the *components* a simulation tracks,
each with the chemical formula that fixes its carbon, nitrogen and electron
content, and the *processes* that turn components into one another, each a
row of a stoichiometric (Gujer) matrix.

Every process is checked as it is read: along its row, carbon, nitrogen and
electrons must balance exactly (docs/theory.md, section 3.6). A process that
would create or destroy any of them is refused with a :class:`ContinuityError`
naming the process and the quantity. A network MARSE accepts therefore cannot
make matter from nothing, whichever engine later runs it.

Rows are rarely written out in full. A growth process gives its yield and lists
in ``balanced_by`` the components whose coefficients the balances determine,
typically the electron acceptor, carbon dioxide and the nitrogen source, and
MARSE solves for them exactly. Water and protons are never listed: they close
the oxygen, hydrogen and charge balances implicitly.
"""

from __future__ import annotations

import difflib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from marse.core.config import ConfigError
from marse.schemas._reading import Field, plain, read_object
from marse.schemas.formula import QUANTITIES, QUANTITY_UNITS, Formula, parse_formula

__all__ = [
    "SCHEMA",
    "SCHEMA_VERSION",
    "Component",
    "ContinuityError",
    "Growth",
    "Network",
    "Process",
    "load_network",
    "network_from_dict",
]

SCHEMA_VERSION = 2
PHASES = ("dissolved", "particulate")
KINDS = ("growth", "reaction")

NETWORK_FIELDS = {
    "schema_version": Field("integer"),
    "description": Field("text", required=False),
    "components": Field("objects"),
    "processes": Field("objects"),
}
COMPONENT_FIELDS = {
    "name": Field("name"),
    "phase": Field("choice", choices=PHASES),
    "formula": Field("text"),
    "charge": Field("number", required=False),
}
GROWTH_FIELDS = {
    "name": Field("name"),
    "kind": Field("choice", choices=KINDS),
    "biomass": Field("name"),
    "substrate": Field("name"),
    "yield_mol_per_mol": Field("number"),
    "products_mol_per_mol": Field("numbers", required=False),
    "balanced_by": Field("names", required=False),
}
REACTION_FIELDS = {
    "name": Field("name"),
    "kind": Field("choice", choices=KINDS),
    "stoichiometry_mol_per_mol": Field("numbers"),
    "balanced_by": Field("names", required=False),
}
SCHEMA = {
    "network": NETWORK_FIELDS,
    "component": COMPONENT_FIELDS,
    "growth process": GROWTH_FIELDS,
    "reaction process": REACTION_FIELDS,
}
"""Every object of the network format and its fields, as docs/networks.md lists them."""

type Composition = Mapping[str, tuple[Fraction, Fraction, Fraction]]


class ContinuityError(ConfigError):
    """A process would create or destroy carbon, nitrogen or electrons.

    ``imbalance`` maps each quantity to the exact amount, per unit of the
    process, by which its products exceed its reactants; a balanced quantity
    maps to zero.
    """

    def __init__(
        self, message: str, process: str = "", imbalance: Mapping[str, Fraction] | None = None
    ) -> None:
        super().__init__(message)
        self.process = process
        self.imbalance = dict(imbalance or {})


@dataclass(frozen=True, slots=True)
class Component:
    """Something a simulation tracks, counted in mol of its formula unit.

    Biomass written per carbon atom, such as CH1.8O0.5N0.2, is therefore
    counted in C-mol. ``dissolved`` components are carried by the liquid;
    ``particulate`` ones, such as biomass, move only with the biofilm.
    """

    name: str
    phase: Literal["dissolved", "particulate"]
    formula: Formula

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "phase": self.phase,
            "formula": self.formula.text,
            "charge": plain(self.formula.charge),
        }


@dataclass(frozen=True, slots=True)
class Growth:
    """What a growth process states: who grows on what, with which yields."""

    biomass: str
    substrate: str
    yield_mol_per_mol: Fraction
    products_mol_per_mol: Mapping[str, Fraction]


@dataclass(frozen=True, slots=True)
class Process:
    """One row of the stoichiometric matrix, proven on loading to balance.

    ``coefficients`` is the complete row, per unit of the process: negative
    for what it consumes and positive for what it produces. Components absent
    from it take no part. ``given`` holds the coefficients the configuration
    stated, or implied through its yields, and ``balanced_by`` the components
    whose coefficients were solved from the balances. ``water`` and
    ``protons`` are the implicit amounts that close the oxygen, hydrogen and
    charge balances.

    A growth process is normalised per mol of biomass formed, so the rate of
    the process is the growth rate of that biomass.
    """

    name: str
    given: Mapping[str, Fraction]
    balanced_by: tuple[str, ...]
    coefficients: Mapping[str, Fraction]
    water: Fraction
    protons: Fraction
    growth: Growth | None = None

    @property
    def kind(self) -> str:
        return "reaction" if self.growth is None else "growth"

    def coefficient(self, component: str) -> Fraction:
        return self.coefficients.get(component, Fraction(0))

    def equation(self) -> str:
        """The row as a balanced equation, water and protons included."""
        terms = [*self.coefficients.items(), ("H2O", self.water), ("H+", self.protons)]
        left = " + ".join(_term(-amount, name) for name, amount in terms if amount < 0)
        right = " + ".join(_term(amount, name) for name, amount in terms if amount > 0)
        return f"{left or 'nothing'} -> {right or 'nothing'}"

    def to_dict(self) -> dict[str, Any]:
        """The process as it was configured, every default written out."""
        if self.growth is None:
            return {
                "name": self.name,
                "kind": "reaction",
                "stoichiometry_mol_per_mol": {n: plain(c) for n, c in self.given.items()},
                "balanced_by": list(self.balanced_by),
            }
        growth = self.growth
        return {
            "name": self.name,
            "kind": "growth",
            "biomass": growth.biomass,
            "substrate": growth.substrate,
            "yield_mol_per_mol": plain(growth.yield_mol_per_mol),
            "products_mol_per_mol": {n: plain(a) for n, a in growth.products_mol_per_mol.items()},
            "balanced_by": list(self.balanced_by),
        }


@dataclass(frozen=True, slots=True)
class Network:
    """Components and the processes between them: a checked Gujer matrix."""

    components: tuple[Component, ...]
    processes: tuple[Process, ...]
    description: str = ""

    @property
    def component_names(self) -> tuple[str, ...]:
        return tuple(c.name for c in self.components)

    def component(self, name: str) -> Component:
        for component in self.components:
            if component.name == name:
                return component
        raise KeyError(name)

    def process(self, name: str) -> Process:
        for process in self.processes:
            if process.name == name:
                return process
        raise KeyError(name)

    def stoichiometric_matrix(self) -> NDArray[np.float64]:
        """Processes by components: every coefficient, rounded once to float."""
        names = self.component_names
        rows = [[float(p.coefficient(n)) for n in names] for p in self.processes]
        return np.array(rows, dtype=float).reshape(len(self.processes), len(names))

    def composition_matrix(self) -> NDArray[np.float64]:
        """Components by quantities: carbon, nitrogen and electrons per mol."""
        rows = [[float(q) for q in c.formula.composition] for c in self.components]
        return np.array(rows, dtype=float).reshape(len(self.components), len(QUANTITIES))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "description": self.description,
            "components": [c.to_dict() for c in self.components],
            "processes": [p.to_dict() for p in self.processes],
        }


def _decimal(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{float(value):.6g}"


def _term(amount: Fraction, name: str) -> str:
    return name if amount == 1 else f"{_decimal(amount)} {name}"


def _join(parts: Sequence[str]) -> str:
    return parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]


def _label(section: str, index: int, raw: Any) -> str:
    name = raw.get("name") if isinstance(raw, dict) else None
    return f"{section}[{index}] '{name}'" if isinstance(name, str) else f"{section}[{index}]"


def _require_unique(names: Sequence[str], section: str) -> None:
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise ConfigError(f"{section}: names must be unique; {', '.join(repeated)} repeated")


def _component(raw: Any, index: int) -> Component:
    where = _label("components", index, raw)
    values = read_object(raw, where, COMPONENT_FIELDS)
    formula = parse_formula(values["formula"], values.get("charge", 0), f"{where}.formula")
    if not any(formula.composition):
        raise ConfigError(
            f"{where}: {formula.label()} holds no carbon, nitrogen or electrons, so no "
            "balance could constrain it. Water and protons are never listed as components: "
            "they close the oxygen, hydrogen and charge balances implicitly"
        )
    return Component(values["name"], values["phase"], formula)


def _net(row: Mapping[str, Fraction], composition: Composition) -> dict[str, Fraction]:
    """Per quantity, how much more the products hold than the reactants."""
    return {
        quantity: sum((c * composition[name][k] for name, c in row.items()), Fraction(0))
        for k, quantity in enumerate(QUANTITIES)
    }


def _changes(net: Mapping[str, Fraction]) -> str:
    """``creates 3 mol C and destroys 2 mol e-``, for an error message."""
    created = [f"{_decimal(a)} {QUANTITY_UNITS[q]}" for q, a in net.items() if a > 0]
    destroyed = [f"{_decimal(-a)} {QUANTITY_UNITS[q]}" for q, a in net.items() if a < 0]
    parts = [f"creates {_join(created)}"] if created else []
    parts += [f"destroys {_join(destroyed)}"] if destroyed else []
    return " and ".join(parts)


_CARRIERS = {
    "carbon": "carbon dioxide",
    "nitrogen": "the nitrogen source, such as ammonium",
    "electrons": "the electron acceptor, such as oxygen, or a reduced product",
}


def _solve(
    given: Mapping[str, Fraction],
    closers: tuple[str, ...],
    composition: Composition,
    where: str,
) -> dict[str, Fraction]:
    """The closers' coefficients that make every quantity balance, found exactly.

    Gauss-Jordan elimination over fractions on the system A x = b, where A holds
    the closers' compositions and b what the given coefficients leave
    unbalanced. A closer the balances cannot pin down is refused as ambiguous;
    a system with no solution is left for the caller to report.
    """
    size = len(QUANTITIES)
    left = _net(given, composition)
    rows = [
        [composition[name][k] for name in closers] + [-left[quantity]]
        for k, quantity in enumerate(QUANTITIES)
    ]
    pivots: list[int] = []
    for column in range(len(closers)):
        row = len(pivots)
        found = next((i for i in range(row, size) if rows[i][column] != 0), None)
        if found is None:
            partners = [closers[pivots[i]] for i in range(row) if rows[i][column] != 0]
            same = _join([f"'{p}'" for p in partners]) if partners else "the others"
            raise ConfigError(
                f"{where}.balanced_by: '{closers[column]}' holds carbon, nitrogen and "
                f"electrons in the same proportions as {same}, so the balances cannot decide "
                "between them; keep one in balanced_by and give the other a coefficient"
            )
        rows[row], rows[found] = rows[found], rows[row]
        pivot = rows[row][column]
        rows[row] = [x / pivot for x in rows[row]]
        for i in range(size):
            if i != row and rows[i][column] != 0:
                factor = rows[i][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[row], strict=True)]
        pivots.append(column)
    return {closers[column]: rows[i][-1] for i, column in enumerate(pivots)}


def _balance(
    process: str,
    per: str,
    where: str,
    given: Mapping[str, Fraction],
    closers: tuple[str, ...],
    composition: Composition,
) -> dict[str, Fraction]:
    solved = _solve(given, closers, composition, where) if closers else {}
    net = _net({**given, **solved}, composition)
    if not any(net.values()):
        return solved
    unbalanced = [q for q, amount in net.items() if amount]
    if not closers:
        raise ContinuityError(
            f"{where}: {_changes(net)} per {per}. A process may only rearrange carbon, "
            "nitrogen and electrons, never create or destroy them. List in balanced_by the "
            "components whose coefficients the balances should determine (docs/networks.md)",
            process,
            net,
        )
    carried = {q for k, q in enumerate(QUANTITIES) if any(composition[n][k] for n in closers)}
    missing = [q for q in unbalanced if q not in carried]
    if missing:
        leftover = _net(given, composition)
        amounts = _join([f"{_decimal(abs(leftover[q]))} {QUANTITY_UNITS[q]}" for q in missing])
        raise ContinuityError(
            f"{where}: the given coefficients leave {amounts} per {per} unbalanced, and "
            f"nothing in balanced_by carries {_join(missing)}; add a component that does: "
            f"{_join([_CARRIERS[q] for q in missing])}",
            process,
            net,
        )
    raise ContinuityError(
        f"{where}: no choice of coefficients for {_join([f"'{n}'" for n in closers])} "
        f"balances carbon, nitrogen and electrons together ({_join(unbalanced)} would stay "
        "unbalanced); change what balanced_by lists",
        process,
        net,
    )


def _known(names: Sequence[str], components: Mapping[str, Component], where: str) -> None:
    for name in names:
        if name not in components:
            close = difflib.get_close_matches(name, list(components), n=1, cutoff=0.75)
            hint = f"; did you mean '{close[0]}'?" if close else ""
            raise ConfigError(f"{where}: '{name}' is not a component of this network{hint}")


def _growth_row(
    values: Mapping[str, Any], where: str, components: Mapping[str, Component]
) -> tuple[dict[str, Fraction], Growth]:
    biomass, substrate = values["biomass"], values["substrate"]
    products: dict[str, Fraction] = values.get("products_mol_per_mol", {})
    _known([biomass], components, f"{where}.biomass")
    _known([substrate], components, f"{where}.substrate")
    _known(list(products), components, f"{where}.products_mol_per_mol")
    if components[biomass].phase != "particulate":
        raise ConfigError(
            f"{where}.biomass: '{biomass}' is dissolved; what grows must be a particulate component"
        )
    if substrate == biomass:
        raise ConfigError(f"{where}: '{biomass}' cannot be both the biomass and its substrate")
    yield_ = values["yield_mol_per_mol"]
    if yield_ <= 0:
        raise ConfigError(f"{where}.yield_mol_per_mol: must be positive, got {_decimal(yield_)}")
    for product, amount in products.items():
        if product in (biomass, substrate):
            role = "biomass" if product == biomass else "substrate"
            raise ConfigError(
                f"{where}.products_mol_per_mol: '{product}' is this process's {role}, not a product"
            )
        if amount <= 0:
            raise ConfigError(
                f"{where}.products_mol_per_mol.{product}: must be positive, got "
                f"{_decimal(amount)}; a component a process consumes needs a reaction process"
            )
    given = {biomass: Fraction(1), substrate: -1 / yield_}
    given |= {product: amount / yield_ for product, amount in products.items()}
    return given, Growth(biomass, substrate, yield_, dict(products))


def _reaction_row(
    values: Mapping[str, Any], where: str, components: Mapping[str, Component]
) -> dict[str, Fraction]:
    given: dict[str, Fraction] = values["stoichiometry_mol_per_mol"]
    field = f"{where}.stoichiometry_mol_per_mol"
    if not given:
        raise ConfigError(f"{field}: is empty; a process must change at least one component")
    _known(list(given), components, field)
    zero = [name for name, c in given.items() if c == 0]
    if zero:
        raise ConfigError(
            f"{field}: {_join(zero)} given a coefficient of zero; leave out what a process "
            "does not change"
        )
    return dict(given)


def _process(raw: Any, index: int, components: Mapping[str, Component]) -> Process:
    where = _label("processes", index, raw)
    if not isinstance(raw, dict):
        raise ConfigError(f"{where}: expected an object")
    if raw.get("kind") not in KINDS:
        got = "missing" if "kind" not in raw else f"{raw['kind']!r}"
        raise ConfigError(f"{where}.kind: expected 'growth' or 'reaction', got {got}")
    growth_kind = raw["kind"] == "growth"
    values = read_object(raw, where, GROWTH_FIELDS if growth_kind else REACTION_FIELDS)
    closers: tuple[str, ...] = values.get("balanced_by", ())
    _known(closers, components, f"{where}.balanced_by")
    if growth_kind:
        given, growth = _growth_row(values, where, components)
        per = f"mol of {growth.biomass} formed"
    else:
        given, growth = _reaction_row(values, where, components), None
        per = "unit of reaction"
    overlap = [name for name in closers if name in given]
    if overlap:
        raise ConfigError(
            f"{where}.balanced_by: {_join([f"'{n}'" for n in overlap])} already given a "
            "coefficient here; list only components the balances should determine"
        )
    if len(closers) > len(QUANTITIES):
        raise ConfigError(
            f"{where}.balanced_by: lists {len(closers)} components, but only "
            f"{len(QUANTITIES)} balances ({_join(list(QUANTITIES))}) can determine them"
        )
    composition = {name: c.formula.composition for name, c in components.items()}
    row = given | _balance(values["name"], per, where, given, closers, composition)
    coefficients = {name: row[name] for name in components if row.get(name, 0) != 0}
    formulas = {name: components[name].formula for name in coefficients}
    water = -sum((c * formulas[n].oxygen for n, c in coefficients.items()), Fraction(0))
    protons = -sum((c * formulas[n].charge for n, c in coefficients.items()), Fraction(0))
    hydrogen = sum((c * formulas[n].hydrogen for n, c in coefficients.items()), Fraction(0))
    if hydrogen + 2 * water + protons != 0:  # implied by the three balances (theory.md 3.6)
        raise ArithmeticError(f"{where}: hydrogen does not close; this is a bug in MARSE")
    return Process(
        name=values["name"],
        given=given,
        balanced_by=closers,
        coefficients=coefficients,
        water=water,
        protons=protons,
        growth=growth,
    )


def network_from_dict(raw: Any) -> Network:
    """Read a network and prove every process balances.

    Raises :class:`~marse.core.config.ConfigError` for a malformed network and
    :class:`ContinuityError`, a subclass, for a process that creates or
    destroys carbon, nitrogen or electrons.
    """
    if isinstance(raw, dict) and "schema_version" not in raw:
        raise ConfigError(
            "network: no schema_version, so this is not a version 2 configuration. Version 1 "
            "ecosystem configurations still run with 'marse ecosystem' (docs/networks.md)"
        )
    values = read_object(raw, "network", NETWORK_FIELDS)
    if values["schema_version"] != SCHEMA_VERSION:
        raise ConfigError(
            f"network.schema_version: {values['schema_version']} is not supported; "
            f"this version of MARSE reads version {SCHEMA_VERSION}"
        )
    components = tuple(_component(item, i) for i, item in enumerate(values["components"]))
    if not components:
        raise ConfigError("network.components: at least one component is required")
    _require_unique([c.name for c in components], "network.components")
    by_name = {c.name: c for c in components}
    processes = tuple(_process(item, i, by_name) for i, item in enumerate(values["processes"]))
    _require_unique([p.name for p in processes], "network.processes")
    return Network(components, processes, values.get("description", ""))


def _refuse_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    keys = [key for key, _ in pairs]
    repeated = sorted({key for key in keys if keys.count(key) > 1})
    if repeated:
        raise ConfigError(
            f"key {_join([f"'{k}'" for k in repeated])} appears twice in one object; "
            "JSON would keep only the last, silently dropping the first"
        )
    return dict(pairs)


def load_network(path: str | Path) -> Network:
    """Read a network from a JSON file; see :func:`network_from_dict`."""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    try:
        raw = json.loads(text, object_pairs_hook=_refuse_duplicate_keys)
    except json.JSONDecodeError as error:
        raise ConfigError(
            f"{path.name}: invalid JSON ({error.msg} at line {error.lineno})"
        ) from error
    return network_from_dict(raw)
