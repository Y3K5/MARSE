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

# Old import path -> module that now holds it. Packages are listed again in PACKAGES
# so ``import marse.schemas.network`` still resolves after the directory is gone.
MOVED: dict[str, str] = {
    "marse.niche": "marse.processes.niche",
    "marse.genotype": "marse.processes.genotype",
    "marse.additives": "marse.processes.additives",
    "marse.science": "marse.evidence.science",
    "marse.calibration": "marse.analysis.calibration",
    "marse.uncertainty": "marse.analysis.uncertainty",
    "marse.ensemble": "marse.analysis.ensemble",
    "marse.immune": "marse.processes.host.immune",
    "marse.actions": "marse.processes.host.actions",
    "marse.schemas": "marse.model",
    "marse.schemas._reading": "marse.model._reading",
    "marse.schemas.domain": "marse.model.domain",
    "marse.schemas.experiment": "marse.model.experiment",
    "marse.schemas.formula": "marse.model.formula",
    "marse.schemas.network": "marse.model.network",
    "marse.microbes": "marse.processes",
    "marse.microbes.adhesion": "marse.processes.adhesion",
    "marse.microbes.additives": "marse.processes.additives",
    "marse.microbes.cardinal": "marse.processes.cardinal",
    "marse.microbes.genotype": "marse.processes.genotype",
    "marse.microbes.growth": "marse.processes.growth",
    "marse.microbes.kinetics": "marse.processes.kinetics",
    "marse.microbes.niche": "marse.processes.niche",
    "marse.chemistry": "marse.processes.acid_base",
    "marse.chemistry.acid_base": "marse.processes.acid_base",
    "marse.biofilm": "marse.processes",
    "marse.biofilm.biomass": "marse.processes.biomass",
    "marse.biofilm.spreading": "marse.processes.spreading",
    "marse.oral": "marse.environments.oral",
    "marse.oral.diet": "marse.environments.oral.diet",
    "marse.oral.film": "marse.environments.oral.film",
    "marse.oral.mouth": "marse.environments.oral.mouth",
    "marse.oral.stephan": "marse.environments.oral.stephan",
    "marse.experimental.host": "marse.processes.host",
    "marse.experimental.host.actions": "marse.processes.host.actions",
    "marse.experimental.host.immune": "marse.processes.host.immune",
    "marse.spatial.air": "marse.space.air",
    "marse.spatial.colloids": "marse.space.colloids",
    "marse.spatial.column": "marse.space.column",
    "marse.spatial.grid": "marse.space.grid",
    "marse.spatial.solutes": "marse.space.solutes",
    "marse.spatial.spreading": "marse.space.spreading",
    "marse.spatial.surface": "marse.space.surface",
    "marse.spatial.transport": "marse.space.transport",
}

PACKAGES = frozenset(
    {
        "marse.schemas",
        "marse.microbes",
        "marse.chemistry",
        "marse.biofilm",
        "marse.oral",
        "marse.experimental.host",
    }
)

PLACEHOLDERS: dict[str, str] = {
    "marse.adaptation": (
        "marse.adaptation is not implemented. The empty package is deprecated and "
        "will be removed; genotype mapping lives in marse.processes.genotype"
    ),
    "marse.interventions": (
        "marse.interventions is not implemented. The empty package is deprecated and "
        "will be removed; host actions live in marse.processes.host.actions and "
        "dose responses in marse.processes.additives"
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
        if "__all__" in vars(target):
            exported = list(vars(target)["__all__"])
        else:
            exported = [name for name in vars(target) if not name.startswith("_")]
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
            spec = ModuleSpec(fullname, _AliasLoader(), origin="marse._compat")
            if fullname in PACKAGES:
                spec.submodule_search_locations = []
            return spec
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
