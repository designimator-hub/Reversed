"""Verification of the hex8 element against properties any correct element must have.

References: Zienkiewicz, Taylor & Zhu (2013), FEM: Its Basis and Fundamentals, on
rigid-body modes and the patch test (Irons & Razzaque 1972).
"""

import numpy as np
import pytest

from reversed_core.element import NODE_SIGNS, elasticity_matrix, stiffness_matrix


@pytest.fixture(scope="module")
def KE():
    return stiffness_matrix(1.0, 0.3, 1.0)


def test_symmetric(KE):
    assert np.allclose(KE, KE.T, atol=1e-14)


def test_six_rigid_body_modes(KE):
    """A free 3D solid element has exactly 6 zero-energy modes (3 translations, 3 rotations)."""
    w = np.linalg.eigvalsh(KE)
    assert np.sum(np.abs(w) < 1e-10 * w.max()) == 6
    assert np.all(w > -1e-10 * w.max())  # positive semi-definite


def test_rigid_motions_in_null_space(KE):
    X = (NODE_SIGNS + 1.0) / 2.0  # unit cube node coordinates
    modes = []
    for d in range(3):
        t = np.zeros((8, 3))
        t[:, d] = 1.0
        modes.append(t.ravel())
    for axis in np.eye(3):
        modes.append(np.cross(axis, X - 0.5).ravel())
    for m in modes:
        assert np.linalg.norm(KE @ m) < 1e-12 * np.linalg.norm(m)


@pytest.mark.parametrize("nu", [0.0, 0.3, 0.45])
def test_patch_test_uniform_strain_energy(nu):
    """Displacement u = eps x gives uniform strain; energy must equal 1/2 eps:D:eps V exactly."""
    h = 2.5
    K = stiffness_matrix(7.0, nu, h)
    D = elasticity_matrix(7.0, nu)
    rng = np.random.default_rng(0)
    eps_t = rng.normal(size=(3, 3))
    eps_t = 0.5 * (eps_t + eps_t.T)  # symmetric strain tensor
    X = (NODE_SIGNS + 1.0) / 2.0 * h
    u = (X @ eps_t.T).ravel()
    voigt = np.array([eps_t[0, 0], eps_t[1, 1], eps_t[2, 2],
                      2 * eps_t[1, 2], 2 * eps_t[0, 2], 2 * eps_t[0, 1]])
    expected = 0.5 * voigt @ D @ voigt * h**3
    assert 0.5 * u @ K @ u == pytest.approx(expected, rel=1e-12)


def test_scales_linearly_with_size():
    assert np.allclose(stiffness_matrix(1.0, 0.3, 3.0), 3.0 * stiffness_matrix(1.0, 0.3, 1.0))


def test_invalid_poisson_ratio():
    with pytest.raises(ValueError):
        elasticity_matrix(1.0, 0.5)
