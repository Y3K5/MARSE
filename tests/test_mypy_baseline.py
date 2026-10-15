"""The mypy baseline compares counts. It does not run mypy; CI does that."""

from collections import Counter

import pytest

from tools import mypy_baseline as baseline


def test_a_clean_module_that_grows_an_error_fails():
    problems = baseline.compare(Counter({"src/marse/cli/run.py": 1}), {})
    assert problems == ["src/marse/cli/run.py: 1 error(s), and it is not in the baseline"]


def test_a_module_may_not_exceed_or_silently_under_cut_its_allowance():
    allowed = {"src/marse/core/config.py": 4}
    assert baseline.compare(Counter({"src/marse/core/config.py": 5}), allowed)
    under = baseline.compare(Counter({"src/marse/core/config.py": 3}), allowed)
    assert under
    assert "lower it" in under[0]


def test_an_exact_baseline_matches():
    allowed = {"src/marse/core/config.py": 2}
    assert baseline.compare(Counter({"src/marse/core/config.py": 2}), allowed) == []


def test_the_baseline_file_is_path_then_count():
    assert baseline.parse_baseline("# comment\n\nsrc/marse/core/config.py 3\n") == {
        "src/marse/core/config.py": 3
    }
    with pytest.raises(ValueError, match="expected"):
        baseline.parse_baseline("src/marse/core/config.py\n")


def test_error_lines_are_counted_and_notes_are_not(tmp_path):
    report = "\n".join(
        [
            "src/marse/core/config.py:10: error: bad [operator]",
            "src/marse/core/config.py:10: note: left operand",
            "src/marse/core/config.py:11: error: also bad [arg-type]",
            "src/marse/cli/run.py:1: error: no [attr-defined]",
        ]
    )
    assert baseline.errors_by_module(report, tmp_path) == {
        "src/marse/core/config.py": 2,
        "src/marse/cli/run.py": 1,
    }
