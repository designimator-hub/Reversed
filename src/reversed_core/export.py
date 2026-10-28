"""Result export: raw densities (.npy), legacy VTK for ParaView, binary STL surface.

The STL surface is the 0.5 iso-surface of the physical density, extracted with
marching cubes [L] via scikit-image [V]. The grid is padded with zeros so the surface
is closed where material touches the domain boundary.

References
----------
[L] Lorensen, W.E., Cline, H.E. (1987). Marching cubes: a high resolution 3D surface
    construction algorithm. Computer Graphics (SIGGRAPH '87) 21(4), 163-169.
[V] van der Walt, S. et al. (2014). scikit-image: image processing in Python.
    PeerJ 2, e453.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def write_vtk(path: str | Path, x_phys: FloatArray, h: float = 1.0) -> None:
    """Legacy VTK STRUCTURED_POINTS with cell data 'density'."""
    nx, ny, nz = x_phys.shape
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nreversed density\nASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx + 1} {ny + 1} {nz + 1}\n")
        f.write("ORIGIN 0 0 0\n")
        f.write(f"SPACING {h} {h} {h}\n")
        f.write(f"CELL_DATA {nx * ny * nz}\nSCALARS density float 1\nLOOKUP_TABLE default\n")
        # VTK expects x fastest, then y, then z.
        np.savetxt(f, x_phys.transpose(2, 1, 0).ravel(), fmt="%.6g")


def density_to_mesh(
    x_phys: FloatArray, h: float = 1.0, level: float = 0.5
) -> tuple[FloatArray, NDArray[np.int64]]:
    """Vertices (in physical units) and triangle indices of the iso-surface."""
    from skimage.measure import marching_cubes

    padded = np.pad(x_phys, 1, mode="constant", constant_values=0.0)
    verts, faces, _, _ = marching_cubes(padded, level=level, spacing=(h, h, h))
    # Element values sit at voxel centres: shift so centre of voxel (0,0,0) is at h/2.
    verts = verts - h + 0.5 * h
    return verts, faces.astype(np.int64)


def write_stl(path: str | Path, verts: FloatArray, faces: NDArray[np.int64]) -> int:
    """Binary STL. Returns the number of triangles written."""
    tri = verts[faces]  # (n, 3, 3)
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = np.divide(normals, norms, out=np.zeros_like(normals), where=norms > 0)
    with open(path, "wb") as f:
        f.write(b"reversed-core binary STL".ljust(80, b" "))
        f.write(struct.pack("<I", len(faces)))
        rec = np.zeros(len(faces), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec["n"] = normals
        rec["v"] = tri
        f.write(rec.tobytes())
    return len(faces)
