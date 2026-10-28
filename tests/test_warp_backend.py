"""Warp matrix-free backend checked against the assembled CPU reference.

Runs on Warp's CPU device when no CUDA GPU is present, so the kernels are verified in
CI; on a CUDA machine the same tests exercise the GPU."""

import numpy as np
import pytest

wp = pytest.importorskip("warp")

from reversed_core import Box, Fixture, Load, LoadCase, Problem, Settings  # noqa: E402
from reversed_core.element import stiffness_matrix  # noqa: E402
from reversed_core.fe import Assembler, LinearSolver  # noqa: E402
from reversed_core.mesh import VoxelGrid  # noqa: E402
from reversed_core.simp import ComplianceModel  # noqa: E402
from reversed_core.warp_backend import WarpSolver  # noqa: E402


@pytest.fixture(scope="module")
def system():
    g = VoxelGrid(10, 6, 4)
    KE = stiffness_matrix()
    rng = np.random.default_rng(11)
    Ee = rng.uniform(0.01, 1.0, g.n_elements)
    fn = g.nodes_in_box(np.zeros(3), np.array([0, 6, 4]))
    fixed = (3 * fn[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(g.n_dofs), fixed)
    F = np.zeros(g.n_dofs)
    F[free] = rng.normal(size=free.size)
    return g, KE, Ee, fixed, free, F


def test_matvec_matches_assembled_matrix(system):
    g, KE, Ee, fixed, free, _ = system
    ws = WarpSolver(g, KE, free)
    ws.set_moduli(Ee)
    x = np.random.default_rng(12).normal(size=g.n_dofs)
    xm = x.copy()
    xm[fixed] = 0.0
    ref = Assembler(g, KE).assemble(Ee) @ xm
    ref[fixed] = x[fixed]  # masked operator acts as identity on fixed DOFs
    assert np.abs(ws.apply_K(x) - ref).max() <= 1e-13 * np.abs(ref).max()


@pytest.mark.parametrize("dtype,tol", [("float64", 1e-7), ("float32", 5e-4)])
def test_solve_matches_direct(system, dtype, tol):
    g, KE, Ee, _, free, F = system
    Ud, _ = LinearSolver(free, "direct").solve(Assembler(g, KE).assemble(Ee), F)
    ws = WarpSolver(g, KE, free, dtype=getattr(wp, dtype),
                    rtol=1e-9 if dtype == "float64" else 1e-6)
    Uw, _ = ws.solve(Ee, F)
    assert np.linalg.norm(Uw - Ud) / np.linalg.norm(Ud) < tol


def test_compliance_and_gradient_match_reference():
    def problem(solver):
        return Problem(
            nelx=8, nely=4, nelz=2, h=1.0,
            fixtures=(Fixture(Box((0, 0, 0), (0, 4, 2))),),
            load_cases=(LoadCase((Load(Box((8, 0, 0), (8, 0, 2)), (0, -1, 0)),)),),
            settings=Settings(volume_fraction=0.4, solver=solver),
        )

    x = np.random.default_rng(13).uniform(0.2, 0.9, 64)
    c_ref, dc_ref, _ = ComplianceModel(problem("direct")).evaluate(x)
    c_w, dc_w, _ = ComplianceModel(problem("warp")).evaluate(x)
    assert c_w == pytest.approx(c_ref, rel=1e-7)
    assert np.allclose(dc_w, dc_ref, rtol=1e-6, atol=1e-9 * abs(dc_ref).max())
