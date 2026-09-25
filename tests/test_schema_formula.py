"""Chemical formulas: exact counts, and the three quantities every process conserves.

The electron content (the degree of reduction) is checked against textbook
values, and through it the chemical oxygen demand that wastewater and biofilm
models use as their electron unit: one mol of electrons is 8 g of oxygen.
"""

from fractions import Fraction

import pytest

from marse.core.config import ConfigError
from marse.schemas.formula import parse_formula


def test_counts_are_the_exact_decimals_written():
    biomass = parse_formula("CH1.8O0.5N0.2")
    assert biomass.counts == (
        ("C", 1),
        ("H", Fraction(9, 5)),
        ("N", Fraction(1, 5)),
        ("O", Fraction(1, 2)),
    )
    assert biomass.electrons == Fraction(21, 5)  # exactly 4.2, not a float near it


def test_a_repeated_symbol_is_summed():
    assert parse_formula("CH3COOH").counts == parse_formula("C2H4O2").counts


@pytest.mark.parametrize(
    ("formula", "charge", "electrons"),
    [
        ("C6H12O6", 0, 24),  # glucose, 4 per carbon
        ("CH1.8O0.5N0.2", 0, Fraction(21, 5)),  # average biomass, 4.2 per C-mol
        ("C5H7O2N", 0, 20),  # biomass as Hoover and Porges wrote it
        ("O2", 0, -4),  # the acceptor takes up four electrons
        ("CO2", 0, 0),  # carbon dioxide, water and ammonium are the reference state
        ("H2O", 0, 0),
        ("NH4", 1, 0),
        ("NO3", -1, -8),  # nitrate
        ("N2", 0, -6),
        ("C2H3O2", -1, 8),  # acetate
        ("C3H5O3", -1, 12),  # lactate
        ("C4H4O4", -2, 14),  # succinate
        ("C2H6O", 0, 12),  # ethanol
        ("CH4", 0, 8),
        ("H2", 0, 2),
    ],
)
def test_the_degree_of_reduction_matches_textbook_values(formula, charge, electrons):
    assert parse_formula(formula, charge).electrons == electrons


def test_protonation_does_not_change_the_electron_content():
    assert parse_formula("C2H4O2").electrons == parse_formula("C2H3O2", -1).electrons
    assert parse_formula("NH3").electrons == parse_formula("NH4", 1).electrons


def test_molar_masses_from_the_standard_atomic_weights():
    assert parse_formula("C6H12O6").molar_mass_g_per_mol == pytest.approx(180.156, abs=1e-9)
    assert parse_formula("CH1.8O0.5N0.2").molar_mass_g_per_mol == pytest.approx(24.6263, abs=1e-9)


@pytest.mark.parametrize(("formula", "cod_g_per_g"), [("C6H12O6", "1.07"), ("C5H7O2N", "1.42")])
def test_chemical_oxygen_demand_reproduces_the_textbook_factors(formula, cod_g_per_g):
    # One mol of electrons reduces a quarter mol of O2, which is 8 g. Textbooks
    # quote 1.07 g COD per g of glucose and 1.42 g per g of cells, worked out
    # with integer atomic masses: 192/180 and 160/113.
    parsed = parse_formula(formula)
    integer_mass = sum(n * {"C": 12, "H": 1, "N": 14, "O": 16}[s] for s, n in parsed.counts)
    assert round(8 * parsed.electrons / integer_mass, 2) == Fraction(cod_g_per_g)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "expected a chemical formula"),
        ("ch4", "is not a chemical formula"),
        ("C6H12O6!", "is not a chemical formula"),
        ("C1.", "is not a chemical formula"),
        ("Ca(OH)2", "parentheses"),
        ("NH4+", "give the charge in the separate 'charge' field"),
        ("C0H4", "a count of zero"),
        ("CH4S", "contains S; formulas may use only C, H, O and N"),
        ("FeCl3", "contains Fe"),
    ],
)
def test_malformed_formulas_are_refused(text, message):
    with pytest.raises(ConfigError, match=message.replace("(", r"\(")):
        parse_formula(text)


@pytest.mark.parametrize(("text", "charge"), [("NH4", 0), ("C2H3O2", 0), ("NO", 0)])
def test_an_odd_number_of_electrons_is_refused(text, charge):
    # Ammonium and acetate written without their charge, and a radical.
    with pytest.raises(ConfigError, match="odd number of electrons"):
        parse_formula(text, charge)


def test_lumped_formulas_have_no_parity_to_check():
    assert parse_formula("CH1.7O0.4N0.3").electrons == 4


def test_a_charged_formula_is_labelled_with_its_charge():
    assert parse_formula("NH4", 1).label() == "NH4 (+1)"
    assert parse_formula("C4H4O4", -2).label() == "C4H4O4 (-2)"
    assert parse_formula("CO2").label() == "CO2"
