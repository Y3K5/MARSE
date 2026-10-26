"""Modules moved into subpackages keep working at their old paths, for one release.

Each old top-level path (``marse.niche`` and so on) is now an alias, installed
by :mod:`marse._compat`, that warns and re-exports the same objects. Nothing
inside MARSE may use an alias: the suite runs with warnings as errors, and the
last test imports every module in a fresh interpreter to prove it.
"""

import importlib
import pkgutil
import subprocess
import sys
import warnings

import pytest

import marse
from marse._compat import MOVED


@pytest.mark.parametrize(("old", "new"), MOVED.items())
def test_an_old_path_warns_and_is_the_same_module_underneath(old: str, new: str):
    names = old.split(".")
    # Load parent aliases first. Their own warnings are checked as separate cases.
    for index in range(2, len(names)):
        parent = ".".join(names[:index])
        if parent not in sys.modules:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                importlib.import_module(parent)
    sys.modules.pop(old, None)  # make the alias run again, so its warning is observable
    with pytest.warns(DeprecationWarning, match=f"{old} has moved to {new}"):
        alias = importlib.import_module(old)
    module = importlib.import_module(new)
    assert alias.__all__ == module.__all__
    for name in module.__all__:
        assert getattr(alias, name) is getattr(module, name), name


def test_nothing_inside_marse_imports_an_old_path():
    """Import every module except the aliases, in a fresh interpreter, warnings as errors."""
    modules = sorted(
        info.name
        for info in pkgutil.walk_packages(marse.__path__, prefix="marse.")
        if info.name not in MOVED
    )
    assert "marse.analysis.ensemble" in modules  # the walk found the new layout
    code = f"import importlib\nfor name in {modules!r}:\n    importlib.import_module(name)\n"
    completed = subprocess.run(
        [sys.executable, "-W", "error::DeprecationWarning", "-c", code],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr[-2000:]
