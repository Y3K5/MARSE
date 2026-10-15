"""Writers shared by every engine: frame stores and the run's data files.

VTK export still lives in :mod:`marse.spatial.vtk`. This package is the single
home for the frame store that used to be split across ``core`` and ``ecosystem``.
"""

from marse.io.frames import FORMAT, FORMAT_VERSION, FrameStore, FrameStoreError

__all__ = ["FORMAT", "FORMAT_VERSION", "FrameStore", "FrameStoreError"]
