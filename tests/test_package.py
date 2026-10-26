import importlib
import importlib.metadata
import re
import subprocess
import sys

import pytest

import marse
from marse import cli

SUBPACKAGES = [
    "core",
    "model",
    "processes",
    "environments",
    "space",
    "spatial",
    "validation",
    "ecosystem",
    "analysis",
    "evidence",
    "experimental",
    "engines",
    "cli",
    "io",
    "bench",
    "ingest",
    "structure",
]


def test_version_is_pep440_and_matches_installed_metadata():
    assert re.fullmatch(r"\d+\.\d+\.\d+((a|b|rc)\d+)?(\.post\d+)?(\.dev\d+)?", marse.__version__)
    assert importlib.metadata.version("marse") == marse.__version__


@pytest.mark.parametrize("name", SUBPACKAGES)
def test_architecture_subpackages_are_importable_and_documented(name):
    module = importlib.import_module(f"marse.{name}")
    assert module.__doc__
    assert module.__doc__.strip()


@pytest.mark.parametrize(
    "name",
    ["schemas", "microbes", "chemistry", "biofilm", "oral"],
)
def test_moved_science_packages_warn(name: str):
    """The pre-consolidation package names still import, and say where they went."""
    sys.modules.pop(f"marse.{name}", None)
    with pytest.warns(DeprecationWarning, match=f"marse.{name} has moved to "):
        module = importlib.import_module(f"marse.{name}")
    assert module.__doc__


@pytest.mark.parametrize("name", ["adaptation", "interventions"])
def test_removed_placeholder_packages_warn(name: str):
    """The empty packages stay importable, and say they are not implemented."""
    sys.modules.pop(f"marse.{name}", None)
    with pytest.warns(DeprecationWarning, match=f"marse.{name} is not implemented"):
        module = importlib.import_module(f"marse.{name}")
    assert module.__doc__


def test_cli_reports_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"marse {marse.__version__}"


def test_python_dash_m_entry_point():
    result = subprocess.run(
        [sys.executable, "-m", "marse", "--version"], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == f"marse {marse.__version__}"
