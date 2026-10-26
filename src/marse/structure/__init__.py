"""Score a predicted structure. A score is not a MARSE parameter.

ESMFold and AlphaFold adapters call a function the caller supplies. This
package does not download a model and does not contact a server. Tests pass
a stand-in. pLDDT and PAE are reported with the bins AlphaFold publishes.
Nothing here sets a rate, a yield, or a biomass formula.
"""

from marse.structure.adapters import FoldAdapter, alphafold_adapter, esmfold_adapter
from marse.structure.record import StructureError, StructureRecord

__all__ = [
    "FoldAdapter",
    "StructureError",
    "StructureRecord",
    "alphafold_adapter",
    "esmfold_adapter",
]
