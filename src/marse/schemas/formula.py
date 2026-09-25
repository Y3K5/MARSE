"""Chemical formulas and the three quantities every process must conserve.

Each component of a reaction network (docs/networks.md) declares its chemical
formula, and the formula fixes its content of the three conserved quantities
(docs/theory.md, section 3.6):

- ``carbon``: mol C per mol of the component;
- ``nitrogen``: mol N per mol;
- ``electrons``: mol e- per mol, the degree of reduction. These are the
  electrons the component releases when it is oxidised completely to CO2, H2O
  and NH4+, so glucose holds 24, carbon dioxide 0, and oxygen, which accepts
  electrons, holds -4.

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
    "QUANTITIES",
    "QUANTITY_UNITS",
    "Formula",
    "parse_formula",
]

QUANTITIES = ("carbon", "nitrogen", "electrons")
"""The quantities every process conserves, in the order of composition vectors."""

QUANTITY_UNITS = {"carbon": "mol C", "nitrogen": "mol N", "electrons": "mol e-"}

ATOMIC_WEIGHTS_G_PER_MOL = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999}
"""IUPAC abridged standard atomic weights of the elements MARSE balances."""

_ELEMENTS = ("C", "H", "N", "O")  # Hill order: carbon, hydrogen, then alphabetical
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
        """The degree of reduction: 4C + H - 2O - 3N - charge.

        It counts the electrons released on complete oxidation to CO2, H2O,
        NH4+ and H+, the reference state in which it is zero.
        """
        return 4 * self.carbon + self.hydrogen - 2 * self.oxygen - 3 * self.nitrogen - self.charge

    @property
    def composition(self) -> tuple[Fraction, Fraction, Fraction]:
        """Carbon, nitrogen and electrons per formula unit, as in ``QUANTITIES``."""
        return self.carbon, self.nitrogen, self.electrons

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
    and its counts are summed. Only C, H, O and N are accepted until sulphur and
    phosphorus are balanced too (docs/roadmap.md, Stage 3).

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
                f"{where}: '{text}' contains {symbol}; formulas may use only C, H, O and N "
                "until sulphur and phosphorus are balanced as well (docs/roadmap.md, Stage 3)"
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
    if whole and (formula.hydrogen + formula.nitrogen - charge) % 2:
        example = "ammonium is NH4 with charge 1, acetate C2H3O2 with charge -1"
        raise ConfigError(
            f"{where}: {text} with charge {_number(charge)} would have an odd number of "
            f"electrons. If it is an ion, give its charge ({example}); radicals are not "
            "supported yet"
        )
    return formula
