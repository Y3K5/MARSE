"""Modules moved into subpackages keep working at their old paths, for one release.

Each old top-level path (``marse.niche`` and so on) is now an alias that warns
and re-exports the same objects. Nothing inside MARSE may use an alias: the
suite runs with warnings as errors, and the last test imports every module in
a fresh interpreter to prove it.
"""

import importlib
import pkgutil
import subprocess
import sys

import pytest

import marse

MOVED = {
    "marse.niche": "marse.microbes.niche",
    "marse.genotype": "marse.microbes.genotype",
    "marse.additives": "marse.microbes.additives",
    "marse.science": "marse.evidence.science",
    "marse.calibration": "marse.analysis.calibration",
    "marse.uncertainty": "marse.analysis.uncertainty",
    "marse.ensemble": "marse.analysis.ensemble",
    "marse.immune": "marse.experimental.host.immune",
    "marse.actions": "marse.experimental.host.actions",
}


@pytest.mark.parametrize(("old", "new"), MOVED.items())
def test_an_old_path_warns_and_is_the_same_module_underneath(old: str, new: str):
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
