"""Recorded frames, streamed to disk as standard NumPy files.

The implementation lives in :mod:`marse.io.frames`. This module remains the
historical import path and re-exports the same classes.
"""

from marse.io.frames import FORMAT, FORMAT_VERSION, FrameStore, FrameStoreError

__all__ = ["FORMAT", "FORMAT_VERSION", "FrameStore", "FrameStoreError"]
