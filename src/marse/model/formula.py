"""Chemical formulas and the quantities every process must conserve.

Each component of a reaction network (docs/networks.md) declares its chemical
formula, and the formula fixes its content of the conserved quantities
(docs/theory.md, section 3.6):

- ``carbon``: mol C per mol of the component;
- ``nitrogen``: mol N per mol;
- ``electrons``: mol e- per mol, the degree of reduction. These are the
  electrons the component releases when it is oxidised completely to CO2, H2O,
  NH4+ and phosphate, so glucose holds 24, carbon dioxide 0, and oxygen, which
  accepts electrons, holds -4;
- ``phosphorus``, ``potassium``, ``chlorine`` and ``sodium``: mol of each
  element per mol, balanced in every network whose components contain it.

Counts are exact fractions of the decimal text as written, so CH1.8O0.5N0.2
contains exactly 9/5 hydrogen. Balances computed from them are therefore
exactly zero or visibly not, never "approximately zero".
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction

from marse.core.config import ConfigError

__all__ = [
    "ATOMIC_WEIGHTS_G_PER_MOL",
    "ELEMENT_QUANTITIES",
    "QUANTITIES",
    "QUANTITY_UNITS",
    "Formula",
    "parse_formula",
]

QUANTITIES = ("carbon", "nitrogen", "electrons")
"""The quantities every process conserves, in the order of composition vectors."""

ELEMENT_QUANTITIES = {"phosphorus": "P", "potassium": "K", "chlorine": "Cl", "sodium": "Na"}
"""Elements balanced as well, in every network whose components contain them, in this order."""

QUANTITY_UNITS = {
    "carbon": "mol C",
    "nitrogen": "mol N",
    "electrons": "mol e-",
    "phosphorus": "mol P",
    "potassium": "mol K",
    "chlorine": "mol Cl",
    "sodium": "mol Na",
}

ATOMIC_WEIGHTS_G_PER_MOL = {
    "C": 12.011,
    "H": 1.008,
    "Cl": 35.45,
    "K": 39.098,
    "N": 14.007,
    "Na": 22.990,
    "O": 15.999,
    "P": 30.974,
}
"""IUPAC abridged standard atomic weights of the elements MARSE balances."""

_ELEMENTS = ("C", "H", "Cl", "K", "N", "Na", "O", "P")  # Hill order: C, H, then alphabetical
_VALENCE = {"P": 5, "K": 1, "Cl": -1, "Na": 1}  # in the reference state of the degree of reduction
_TOKEN = re.compile(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)?")
_FORMULA = re.compile(r"(?:[A-Z][a-z]?(?:\d+(?:\.\d+)?)?)+")


@dataclass(frozen=True, slots=True)
class Formula:
    """A parsed chemical formula: exact element counts and a charge.

    Build one with :func:`parse_formula`, which validates the text.
    """

    text: str
    charge: Fraction
    counts: tuple[tuple[str, Fraction], ...]

    def count(self, element: str) -> Fraction:
        """How many atoms of ``element`` one formula unit contains."""
        return next((n for symbol, n in self.counts if symbol == element), Fraction(0))

    @property
    def carbon(self) -> Fraction:
        return self.count("C")

    @property
    def hydrogen(self) -> Fraction:
        return self.count("H")

    @property
    def nitrogen(self) -> Fraction:
        return self.count("N")

    @property
    def oxygen(self) -> Fraction:
        return self.count("O")

    @property
    def electrons(self) -> Fraction:
        """The degree of reduction: 4C + H - 2O - 3N + 5P + K - Cl + Na - charge.

        It counts the electrons released on complete oxidation to CO2, H2O,
        NH4+, phosphate, K+, Cl-, Na+ and H+, the reference state in which it is
        zero. Protonation leaves it unchanged: an added proton brings one
        hydrogen and one charge.
        """
        valences = sum((v * self.count(e) for e, v in _VALENCE.items()), Fraction(0))
        return (
            4 * self.carbon
            + self.hydrogen
            - 2 * self.oxygen
            - 3 * self.nitrogen
            + valences
            - self.charge
        )

    @property
    def composition(self) -> tuple[Fraction, Fraction, Fraction]:
        """Carbon, nitrogen and electrons per formula unit, as in ``QUANTITIES``."""
        return self.carbon, self.nitrogen, self.electrons

    def content(self, quantity: str) -> Fraction:
        """How much of a conserved quantity one formula unit holds, by its name."""
        if quantity in ELEMENT_QUANTITIES:
            return self.count(ELEMENT_QUANTITIES[quantity])
        return self.composition[QUANTITIES.index(quantity)]

    @property
    def molar_mass_g_per_mol(self) -> float:
        return sum(float(n) * ATOMIC_WEIGHTS_G_PER_MOL[symbol] for symbol, n in self.counts)

    def label(self) -> str:
        """The formula with its charge, for reports: ``NH4 (+1)``."""
        if self.charge == 0:
            return self.text
        return f"{self.text} ({'+' if self.charge > 0 else '-'}{_number(abs(self.charge))})"


def _number(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{float(value):g}"


def parse_formula(text: str, charge: Fraction | int = 0, where: str = "formula") -> Formula:
    """Read a formula such as ``C6H12O6`` or ``CH1.8O0.5N0.2``.

    Element symbols are followed by an optional count, which may be a decimal
    for a lumped composition such as biomass. A symbol may repeat (``CH3COOH``)
    and its counts are summed. The elements are C, H, N and O, and P, K, Cl and
    Na, which are balanced too; others, such as sulphur and calcium, are refused
    until they are balanced (docs/roadmap.md).

    A formula whose counts and charge are all whole numbers must have an even
    number of electrons. An odd count almost always means an ion written
    without its charge, such as ammonium as ``NH4``, which would give it an
    electron it does not have. Radicals (NO, for example) are not supported
    yet.
    """
    if not isinstance(text, str) or not text:
        raise ConfigError(f"{where}: expected a chemical formula such as C6H12O6")
    if "(" in text or ")" in text:
        raise ConfigError(
            f"{where}: '{text}' uses parentheses, which are not supported; "
            "write each element once with its total count, e.g. C16H32O2"
        )
    if "+" in text or "-" in text:
        raise ConfigError(
            f"{where}: '{text}' includes a charge sign; write the formula without it "
            "and give the charge in the separate 'charge' field, e.g. NH4 with charge 1"
        )
    if not _FORMULA.fullmatch(text):
        raise ConfigError(
            f"{where}: '{text}' is not a chemical formula; expected element symbols, "
            "each followed by an optional count, such as C6H12O6 or CH1.8O0.5N0.2"
        )
    totals: dict[str, Fraction] = {}
    for symbol, digits in _TOKEN.findall(text):
        if symbol not in _ELEMENTS:
            raise ConfigError(
                f"{where}: '{text}' contains {symbol}; formulas may use only C, H, N, O, P, K, "
                "Cl and Na until other elements are balanced as well (docs/roadmap.md)"
            )
        count = Fraction(digits) if digits else Fraction(1)
        if count <= 0:
            raise ConfigError(f"{where}: '{text}' gives {symbol} a count of zero")
        totals[symbol] = totals.get(symbol, Fraction(0)) + count
    charge = Fraction(charge)
    formula = Formula(
        text=text,
        charge=charge,
        counts=tuple((symbol, totals[symbol]) for symbol in _ELEMENTS if symbol in totals),
    )
    whole = charge.denominator == 1 and all(n.denominator == 1 for _, n in formula.counts)
    odd = formula.hydrogen + formula.nitrogen - charge + sum(formula.count(e) for e in _VALENCE)
    if whole and odd % 2:
        example = "ammonium is NH4 with charge 1, acetate C2H3O2 with charge -1"
        raise ConfigError(
            f"{where}: {text} with charge {_number(charge)} would have an odd number of "
            f"electrons. If it is an ion, give its charge ({example}); radicals are not "
            "supported yet"
        )
    return formula
