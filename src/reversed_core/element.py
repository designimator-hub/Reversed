"""Trilinear 8-node hexahedral (hex8) element for 3D linear elasticity.

The element is a cube of side ``h`` integrated with 2x2x2 Gauss-Legendre quadrature.
For an undistorted brick the integrand of B^T D B is at most quadratic in each
natural coordinate, so 2-point quadrature integrates the stiffness exactly.

References
----------
[Z] Zienkiewicz, O.C., Taylor, R.L., Zhu, J.Z. (2013). The Finite Element Method:
    Its Basis and Fundamentals, 7th ed. Butterworth-Heinemann. Chapters on shape
    functions, isoparametric mapping and Gauss quadrature.
[H] Hughes, T.J.R. (2000). The Finite Element Method: Linear Static and Dynamic
    Finite Element Analysis. Dover. Ch. 2 (elasticity in Voigt form) and Ch. 3
    (isoparametric hexahedra).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# Local node order, natural coordinates (xi, eta, zeta) in {-1, +1}. [H]
NODE_SIGNS: NDArray[np.float64] = np.array(
    [
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1],
    ],
    dtype=np.float64,
)

_GAUSS_1D = np.array([-1.0, 1.0]) / np.sqrt(3.0)


def elasticity_matrix(E: float, nu: float) -> NDArray[np.float64]:
    """Isotropic 6x6 constitutive matrix in Voigt order (xx, yy, zz, yz, xz, xy).

    Shear components use engineering strain (gamma = 2 * epsilon). [H]
    """
    if not (-1.0 < nu < 0.5):
        raise ValueError(f"Poisson's ratio must lie in (-1, 0.5), got {nu}")
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2.0 * mu
    D[np.arange(3, 6), np.arange(3, 6)] = mu
    return D


def shape_gradients_natural(xi: NDArray[np.float64]) -> NDArray[np.float64]:
    """dN_a/dxi_i for the 8 trilinear shape functions N_a = prod_j (1 + s_aj xi_j) / 8.

    Returns an (8, 3) array. [Z]
    """
    factors = 1.0 + NODE_SIGNS * xi  # (8, 3)
    grads = np.empty((8, 3))
    for i in range(3):
        others = [j for j in range(3) if j != i]
        grads[:, i] = NODE_SIGNS[:, i] * factors[:, others[0]] * factors[:, others[1]] / 8.0
    return grads


def strain_displacement(xi: NDArray[np.float64], h: float) -> NDArray[np.float64]:
    """6x24 strain-displacement matrix B at natural point xi for a cube of side h.

    The map x = h/2 (xi + 1) gives a constant Jacobian J = (h/2) I. [H]
    """
    dN = shape_gradients_natural(xi) * (2.0 / h)  # dN/dx, (8, 3)
    B = np.zeros((6, 24))
    for a in range(8):
        dx, dy, dz = dN[a]
        c = 3 * a
        B[0, c] = dx
        B[1, c + 1] = dy
        B[2, c + 2] = dz
        B[3, c + 1], B[3, c + 2] = dz, dy  # gamma_yz
        B[4, c], B[4, c + 2] = dz, dx  # gamma_xz
        B[5, c], B[5, c + 1] = dy, dx  # gamma_xy
    return B


def stiffness_matrix(E: float = 1.0, nu: float = 0.3, h: float = 1.0) -> NDArray[np.float64]:
    """24x24 element stiffness K_e = sum_g B^T D B det(J) w_g. [Z]

    DOF order is (u_x, u_y, u_z) per node, nodes in NODE_SIGNS order.
    K_e scales linearly with h for a cube (B ~ 1/h, det J ~ h^3).
    """
    D = elasticity_matrix(E, nu)
    det_j = (h / 2.0) ** 3
    K = np.zeros((24, 24))
    for gx in _GAUSS_1D:
        for gy in _GAUSS_1D:
            for gz in _GAUSS_1D:
                B = strain_displacement(np.array([gx, gy, gz]), h)
                K += B.T @ D @ B * det_j  # Gauss weights are 1 for 2-point rule
    return 0.5 * (K + K.T)


def element_centroid_B(h: float) -> NDArray[np.float64]:
    """B evaluated at the element centre, used for element stress recovery."""
    return strain_displacement(np.zeros(3), h)
