"""VTK image files, so a run can be explored in ParaView and other 3-D viewers.

Each recorded frame becomes one ``.vti`` file (VTK XML image data): a
regular grid whose voxels carry one value per component. A ``.pvd`` file lists
the frames with their times, so a viewer plays them as a time series. The
files are written with the standard library and NumPy only, in the documented
format (Kitware, "VTK File Formats", XML image data with raw appended data).

Voxels are cells: the extent counts voxel corners, the spacing is the voxel
edge in micrometres, and height is the third axis with the substratum at
z = 0. One- and two-dimensional grids are written as a single column or slab
of voxels.

:func:`read_vti` reads the files this module writes. It is used by the tests,
and it lets a user check a file without a viewer.
"""

from __future__ import annotations

import re
from pathlib import Path
from xml.sax.saxutils import quoteattr

import numpy as np
from numpy.typing import NDArray

from marse.spatial.grid import Grid

__all__ = ["read_vti", "write_pvd", "write_vti"]

_TYPES = {np.dtype("<f4"): "Float32", np.dtype("<f8"): "Float64"}


def _as_3d(grid: Grid) -> tuple[int, int, int]:
    """Voxel counts along x, y and height; a missing lateral axis has one voxel."""
    lateral = [grid.shape[a] for a in grid.lateral_axes]
    lateral += [1] * (2 - len(lateral))
    return lateral[0], lateral[1], grid.shape[-1]


def write_vti(
    path: str | Path,
    grid: Grid,
    fields: dict[str, NDArray[np.float64]],
    *,
    dtype: str = "<f4",
) -> Path:
    """Write one frame: each field has the grid's shape. Returns the path written."""
    out = Path(path)
    nx, ny, nz = _as_3d(grid)
    kind = np.dtype(dtype)
    if kind not in _TYPES:
        raise ValueError("VTK frames are written as little-endian float32 or float64")
    arrays, blocks, offset = [], [], 0
    for name, values in fields.items():
        data = np.asarray(values, dtype=float)
        if data.shape != grid.shape:
            raise ValueError(f"field {name!r} has shape {data.shape}, not the grid's {grid.shape}")
        # VTK orders cells with x varying fastest.
        raw = data.reshape(nx, ny, nz).transpose(2, 1, 0).astype(kind).tobytes()
        block = np.uint64(len(raw)).tobytes() + raw
        arrays.append(
            f'        <DataArray type="{_TYPES[kind]}" Name={quoteattr(name)} '
            f'format="appended" offset="{offset}"/>'
        )
        blocks.append(block)
        offset += len(block)
    h = grid.voxel_um
    extent = f"0 {nx} 0 {ny} 0 {nz}"
    header = "\n".join(
        [
            '<?xml version="1.0"?>',
            '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" '
            'header_type="UInt64">',
            f'  <ImageData WholeExtent="{extent}" Origin="0 0 0" Spacing="{h!r} {h!r} {h!r}">',
            f'    <Piece Extent="{extent}">',
            "      <CellData>",
            *arrays,
            "      </CellData>",
            "    </Piece>",
            "  </ImageData>",
            '  <AppendedData encoding="raw">',
            "   _",
        ]
    )
    with out.open("wb") as handle:
        handle.write(header.encode("ascii"))
        for block in blocks:
            handle.write(block)
        handle.write(b"\n  </AppendedData>\n</VTKFile>\n")
    return out


def write_pvd(path: str | Path, frames: list[tuple[float, str]]) -> Path:
    """Write the time-series index: (time in hours, .vti file name relative to the .pvd)."""
    out = Path(path)
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        "  <Collection>",
        *(
            f'    <DataSet timestep="{time_h!r}" group="" part="0" file={quoteattr(name)}/>'
            for time_h, name in frames
        ),
        "  </Collection>",
        "</VTKFile>",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="ascii")
    return out


def read_vti(path: str | Path) -> tuple[dict[str, NDArray[np.float64]], tuple[int, ...], float]:
    """Read a file from :func:`write_vti`: fields, (x, y, height) voxel counts, spacing."""
    content = Path(path).read_bytes()
    marker = content.index(b"_", content.index(b"<AppendedData"))
    header = content[:marker].decode("ascii")
    extent = re.search(r'WholeExtent="([^"]+)"', header)
    spacing = re.search(r'Spacing="([^"]+)"', header)
    if extent is None or spacing is None:
        raise ValueError(f"{path} is not a VTK image file written by MARSE")
    bounds = [int(v) for v in extent.group(1).split()]
    nx, ny, nz = bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]
    data = content[marker + 1 :]
    fields = {}
    for kind, name, offset in re.findall(
        r'<DataArray type="(Float32|Float64)" Name="([^"]+)" format="appended" offset="(\d+)"/>',
        header,
    ):
        start = int(offset)
        size = int(np.frombuffer(data[start : start + 8], dtype="<u8")[0])
        values = np.frombuffer(
            data[start + 8 : start + 8 + size], dtype="<f4" if kind == "Float32" else "<f8"
        )
        fields[name] = values.reshape(nz, ny, nx).transpose(2, 1, 0).astype(float)
    return fields, (nx, ny, nz), float(spacing.group(1).split()[0])
