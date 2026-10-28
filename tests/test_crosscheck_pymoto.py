"""Cross-check against an independent implementation: pyMOTO (TU Delft).

pyMOTO: Delissen, A. et al., pyMOTO - Modular Topology Optimization framework in
Python, https://github.com/aatmdelissen/pyMOTO. Its node and element numbering
differ from ours (x fastest, ours z fastest), so agreement also verifies our
numbering, boundary-condition selection and load distribution.
"""

import warnings

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from reversed_core import Box, Fixture, Load, LoadCase, Material, Problem, Settings
from reversed_core.element import stiffness_matrix
from reversed_core.simp import ComplianceModel

pym = pytest.importorskip("pymoto")


def test_element_stiffness_spectrum_matches():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dom = pym.VoxelDomain(2, 2, 2, unitx=0.5, unity=0.5, unitz=0.5)
        Kp = pym.AssembleStiffness(dom, e_modulus=2.0, poisson_ratio=0.3).stiffness_element
    Ko = stiffness_matrix(2.0, 0.3, 0.5)
    assert np.allclose(np.sort(np.linalg.eigvalsh(Kp)), np.sort(np.linalg.eigvalsh(Ko)),
                       rtol=1e-12, atol=1e-12)


def test_heterogeneous_cantilever_compliance_matches():
    nx, ny, nz = 12, 4, 3
    rng = np.random.default_rng(0)
    x_ours = rng.uniform(0.1, 1.0, nx * ny * nz)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dom = pym.VoxelDomain(nx, ny, nz)
        jj, kk = np.meshgrid(np.arange(ny + 1), np.arange(nz + 1), indexing="ij")
        face = dom.get_nodenumber(np.zeros(jj.size, int), jj.ravel(), kk.ravel())
        bc = np.concatenate([3 * face + d for d in range(3)])
        line = dom.get_nodenumber(np.full(nz + 1, nx), np.zeros(nz + 1, int), np.arange(nz + 1))
        f = np.zeros(3 * dom.nnodes)
        f[3 * line + 1] = -1.0 / (nz + 1)
        x_pym = x_ours.reshape(nx, ny, nz).transpose(2, 1, 0).ravel()
        K = pym.AssembleStiffness(dom, bc=bc, e_modulus=1.0, poisson_ratio=0.3)(
            pym.Signal("x", state=x_pym)).state
        c_pym = f @ spla.spsolve(K.tocsc(), f)

    p = Problem(
        nelx=nx, nely=ny, nelz=nz, h=1.0,
        fixtures=(Fixture(Box((0, 0, 0), (0, ny, nz))),),
        load_cases=(LoadCase((Load(Box((nx, 0, 0), (nx, 0, nz)), (0, -1, 0)),)),),
        material=Material(E=1.0, nu=0.3, E_min_ratio=0.0),
        settings=Settings(r_min=1.0, penal=1.0),  # identity filter, E_e = x_e
    )
    c_ours, _, _ = ComplianceModel(p).evaluate(x_ours)
    assert c_ours == pytest.approx(c_pym, rel=1e-10)
