"""Export checks: STL file structure, closed surface, and correct physical placement."""

import numpy as np
import pytest

from reversed_core.export import density_to_mesh, write_stl, write_vtk


def test_stl_cube_is_closed_and_placed(tmp_path):
    pytest.importorskip("skimage")
    x = np.zeros((6, 6, 6))
    x[1:5, 1:5, 1:5] = 1.0  # solid block from 1 to 5 voxel edges, h = 2
    verts, faces = density_to_mesh(x, h=2.0)
    # Iso 0.5 between centres of solid and void voxels lands on the voxel faces.
    assert verts.min(axis=0) == pytest.approx([2.0] * 3)
    assert verts.max(axis=0) == pytest.approx([10.0] * 3)
    # Closed 2-manifold: every edge is shared by exactly two triangles.
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(counts == 2)
    n = write_stl(tmp_path / "c.stl", verts, faces)
    assert (tmp_path / "c.stl").stat().st_size == 84 + 50 * n


def test_vtk_header(tmp_path):
    write_vtk(tmp_path / "d.vtk", np.ones((2, 3, 4)), h=0.5)
    text = (tmp_path / "d.vtk").read_text().splitlines()
    assert "DIMENSIONS 3 4 5" in text
    assert "CELL_DATA 24" in text
    assert len(text) == 10 + 24
