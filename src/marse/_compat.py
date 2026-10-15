"""Deprecated import paths, kept working for one release.

Top-level modules that moved into subpackages (``marse.niche`` and the rest)
and the empty ``adaptation`` and ``interventions`` packages are no longer
files in the tree. Importing an old name warns and, for a moved module,
re-exports the new module's public names.

Nothing inside MARSE may import an old path. The layout test loads every
current module with deprecation warnings treated as errors.
"""

from __future__ import annotations

import importlib
import sys
import warnings
from collections.abc import Sequence
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from types import ModuleType

MOVED: dict[str, str] = {
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

PLACEHOLDERS: dict[str, str] = {
    "marse.adaptation": (
        "marse.adaptation is not implemented. The empty package is deprecated and "
        "will be removed; genotype mapping lives in marse.microbes.genotype"
    ),
    "marse.interventions": (
        "marse.interventions is not implemented. The empty package is deprecated and "
        "will be removed; host actions live in marse.experimental.host.actions and "
        "dose responses in marse.microbes.additives"
    ),
}


class _AliasLoader(Loader):
    """Load a moved module's public names under the old module name."""

    def create_module(self, spec: ModuleSpec) -> None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        target_name = MOVED[module.__name__]
        warnings.warn(
            f"{module.__name__} has moved to {target_name}; "
            "the old import path will be removed in the next release",
            DeprecationWarning,
            stacklevel=2,
        )
        target = importlib.import_module(target_name)
        exported = list(vars(target)["__all__"])
        vars(module)["__all__"] = exported
        module.__doc__ = f"Deprecated alias of :mod:`{target_name}`, kept for one release."
        for name in exported:
            setattr(module, name, getattr(target, name))


class _PlaceholderLoader(Loader):
    """Import an empty package that has been removed from the tree."""

    def create_module(self, spec: ModuleSpec) -> None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        message = PLACEHOLDERS[module.__name__]
        module.__doc__ = message
        vars(module)["__all__"] = []
        warnings.warn(message, DeprecationWarning, stacklevel=2)


class _CompatFinder(MetaPathFinder):
    """Resolve the deprecated names in :data:`MOVED` and :data:`PLACEHOLDERS`."""

    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None:
        if fullname in MOVED:
            return ModuleSpec(fullname, _AliasLoader(), origin="marse._compat")
        if fullname in PLACEHOLDERS:
            spec = ModuleSpec(
                fullname, _PlaceholderLoader(), origin="marse._compat", is_package=True
            )
            spec.submodule_search_locations = []
            return spec
        return None


def install() -> None:
    """Register the compatibility finder. Safe to call more than once."""
    if any(isinstance(finder, _CompatFinder) for finder in sys.meta_path):
        return
    sys.meta_path.insert(0, _CompatFinder())
