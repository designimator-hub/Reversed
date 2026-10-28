"""Finite-element verification against Timoshenko beam theory and solver cross-checks.

Timoshenko cantilever tip deflection under end load P:
    delta = P L^3 / (3 E I) + P L / (kappa G A),  kappa = 5/6 for a rectangle.
Reference: Timoshenko, S.P., Gere, J.M. (1972). Mechanics of Materials. Van Nostrand.
Cowper, G.R. (1966). The shear coefficient in Timoshenko's beam theory.
J. Appl. Mech. 33, 335-340.
"""

import numpy as np
import pytest

from reversed_core import Box, Fixture, Load, LoadCase, Material, Problem, Settings
from reversed_core.element import stiffness_matrix
from reversed_core.fe import Assembler, solve_free
from reversed_core.mesh import VoxelGrid
from reversed_core.simp import ComplianceModel

E, L = 1000.0, 10.0


def beam_problem(n: int, solver: str = "direct") -> Problem:
    return Problem(
        nelx=10 * n, nely=n, nelz=n, h=1.0 / n,
        fixtures=(Fixture(Box((0, 0, 0), (0, 1, 1))),),
        load_cases=(LoadCase((Load(Box((L, 0, 0), (L, 1, 1)), (0, -1, 0)),)),),
        material=Material(E=E, nu=0.0),
        settings=Settings(volume_fraction=1.0, r_min=1.0, solver=solver),
    )


def tip_deflection(n: int) -> float:
    p = beam_problem(n)
    c, _, _ = ComplianceModel(p).evaluate(np.ones(p.grid.n_elements))
    return c  # P = 1, so compliance = load-averaged tip deflection


def timoshenko() -> float:
    inertia, area, shear_mod = 1.0 / 12.0, 1.0, E / 2.0
    return L**3 / (3 * E * inertia) + L / (5.0 / 6.0 * shear_mod * area)


@pytest.mark.slow
def test_cantilever_converges_to_timoshenko():
    d2, d4, d8 = tip_deflection(2), tip_deflection(4), tip_deflection(8)
    exact = timoshenko()
    e2, e4, e8 = (abs(exact - d) for d in (d2, d4, d8))
    # Trilinear elements converge at second order in h for this quantity.
    assert e2 / e4 > 3.0 and e4 / e8 > 3.0
    # Richardson extrapolation with p = 2 must recover the theory within 0.5 %.
    extrap = d8 + (d8 - d4) / 3.0
    assert extrap == pytest.approx(exact, rel=5e-3)
    assert d8 == pytest.approx(exact, rel=1.5e-2)


def test_direct_and_cg_agree():
    grid = VoxelGrid(8, 4, 3)
    KE = stiffness_matrix(1.0, 0.3)
    rng = np.random.default_rng(1)
    K = Assembler(grid, KE).assemble(rng.uniform(0.1, 1.0, grid.n_elements))
    fixed = 3 * grid.nodes_in_box(np.zeros(3), np.array([0, 4, 3])).repeat(3) + np.tile(
        np.arange(3), grid.nodes_in_box(np.zeros(3), np.array([0, 4, 3])).size)
    free = np.setdiff1d(np.arange(grid.n_dofs), fixed)
    F = rng.normal(size=(grid.n_dofs, 2))
    Ud, _ = solve_free(K, F, free, "direct")
    Uc, stats = solve_free(K, F, free, "cg", rtol=1e-11)
    assert np.allclose(Ud, Uc, rtol=1e-6, atol=1e-8)
    assert all(s.residual < 1e-10 for s in stats)
    assert np.all(Ud[fixed] == 0)


def test_equilibrium_reactions_balance_load():
    """Sum of reactions at the support equals minus the applied load."""
    p = beam_problem(2)
    m = ComplianceModel(p)
    m.evaluate(np.ones(p.grid.n_elements))
    K = m.assembler.assemble(np.full(p.grid.n_elements, E))
    R = K @ m.U[:, 0] - m.F[:, 0]
    assert R[1::3].sum() == pytest.approx(1.0, rel=1e-9)  # applied Fy = -1
    assert abs(R[0::3].sum()) < 1e-9 and abs(R[2::3].sum()) < 1e-9
