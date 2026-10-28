"""Multigrid verification: interpolation exactness, symmetry of the V-cycle (required
for CG, Briggs, Henson & McCormick 2000), and agreement with the direct solver."""

import numpy as np
import pytest

from reversed_core.element import stiffness_matrix
from reversed_core.fe import Assembler, LinearSolver
from reversed_core.mesh import VoxelGrid
from reversed_core.multigrid import MultigridPreconditioner, interp_1d, prolongation


def test_interp_1d_reproduces_linear_functions():
    P = interp_1d(8, 2)
    xc = np.linspace(0, 1, 5)
    assert np.allclose(P @ (3 * xc + 1), 3 * np.linspace(0, 1, 9) + 1)
    assert np.allclose(np.asarray(P.sum(axis=1)).ravel(), 1.0)  # partition of unity


@pytest.mark.parametrize("shape", [(8, 4, 2), (6, 3, 1), (4, 4, 4)])
def test_prolongation_exact_for_affine_displacement(shape):
    g = VoxelGrid(*shape)
    P, gc = prolongation(g)
    A = np.array([[0.1, 0.2, -0.3], [0.0, 0.5, 0.1], [0.2, 0.0, 0.4]])
    b = np.array([1.0, -2.0, 0.5])

    def field(grid, scale):
        X = grid.node_coords / grid.h * scale  # node coords in fine-element units
        return (X @ A.T + b).ravel()

    scale_c = np.array(shape) / np.array(gc.shape)
    Xc = np.stack(np.meshgrid(*(np.arange(n + 1) for n in gc.shape), indexing="ij"), -1)
    uc = ((Xc.reshape(-1, 3) * scale_c) @ A.T + b).ravel()
    assert np.allclose(P @ uc, field(g, 1.0))


def _system(shape=(16, 8, 4)):
    g = VoxelGrid(*shape)
    rng = np.random.default_rng(9)
    x = np.where(rng.uniform(size=g.n_elements) > 0.4, 1.0, 1e-3)
    K = Assembler(g, stiffness_matrix()).assemble(x)
    fixed_nodes = g.nodes_in_box(np.zeros(3), np.array([0, shape[1], shape[2]]))
    fixed = (3 * fixed_nodes[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(g.n_dofs), fixed)
    F = rng.normal(size=g.n_dofs)
    return g, K, free, F


def test_vcycle_is_symmetric():
    g, K, free, _ = _system()
    mg = MultigridPreconditioner(g, free, min_coarse_dofs=100)
    mg.setup(K[free][:, free].tocsr())
    assert mg.n_levels >= 3
    rng = np.random.default_rng(10)
    u, v = rng.normal(size=(2, free.size))
    assert u @ mg.apply(v) == pytest.approx(v @ mg.apply(u), rel=1e-10)
    assert u @ mg.apply(u) > 0


def test_multigrid_matches_direct_with_high_contrast():
    g, K, free, F = _system()
    Ud, _ = LinearSolver(free, "direct").solve(K, F)
    s = LinearSolver(free, "mg", grid=g, rtol=1e-10)
    s.mg = MultigridPreconditioner(g, free, min_coarse_dofs=100)
    Um, stats = s.solve(K, F)
    assert np.allclose(Um, Ud, rtol=1e-7, atol=1e-9 * abs(Ud).max())
    assert stats[0].residual < 1e-9


def test_auto_selects_multigrid_for_large_systems():
    g = VoxelGrid(40, 20, 12)
    free = np.arange(30_000)
    assert LinearSolver(free, "auto", grid=g).method == "mg"
    assert LinearSolver(free[:1000], "auto", grid=g).method == "direct"
