"""Geometric multigrid V-cycle as a preconditioner for conjugate gradients (MG-PCG).

Hierarchy: each level halves the element count in every direction that is even
(semi-coarsening where a direction is odd or already 1). Prolongation P is trilinear
interpolation of nodal displacements. Coarse operators are Galerkin products
A_c = P^T A P restricted to free DOFs, which stay symmetric positive definite and
handle the high stiffness contrast of SIMP designs better than rediscretization [AAL].
Smoother: damped Jacobi, the same number of sweeps before and after coarse correction,
so the V-cycle is a symmetric operator and a valid CG preconditioner [BHM, TOS].
Coarsest level: sparse LU.

References
----------
[AAL] Amir, O., Aage, N., Lazarov, B.S. (2014). On multigrid-CG for efficient topology
      optimization. Struct. Multidisc. Optim. 49, 815-829.
[BHM] Briggs, W.L., Henson, V.E., McCormick, S.F. (2000). A Multigrid Tutorial, 2nd ed.
      SIAM.
[TOS] Trottenberg, U., Oosterlee, C.W., Schuller, A. (2001). Multigrid. Academic Press.
[WDW] Wu, J., Dick, C., Westermann, R. (2016). A system for high-resolution topology
      optimization. IEEE Trans. Vis. Comput. Graph. 22(3), 1195-1208.
      (the same scheme run matrix-free on the GPU; planned for the Warp backend)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from .mesh import VoxelGrid

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def interp_1d(n_fine: int, factor: int) -> sp.csr_matrix:
    """(n_fine + 1) x (n_fine // factor + 1) linear interpolation of node values."""
    if factor == 1:
        return sp.identity(n_fine + 1, format="csr")
    n_coarse = n_fine // 2
    rows, cols, vals = [], [], []
    for i in range(n_fine + 1):
        if i % 2 == 0:
            rows.append(i)
            cols.append(i // 2)
            vals.append(1.0)
        else:
            rows += [i, i]
            cols += [i // 2, i // 2 + 1]
            vals += [0.5, 0.5]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n_fine + 1, n_coarse + 1))


def prolongation(grid: VoxelGrid) -> tuple[sp.csr_matrix, VoxelGrid]:
    """DOF prolongation from the next coarser grid, and that grid.

    Node ordering (i slowest, k fastest) matches kron(Px, kron(Py, Pz)); each node's three
    DOFs are interpolated independently, hence the trailing kron with I_3.
    """
    factors = [2 if n % 2 == 0 and n >= 2 else 1 for n in grid.shape]
    if factors == [1, 1, 1]:
        raise ValueError("grid cannot be coarsened further")
    Px, Py, Pz = (interp_1d(n, f) for n, f in zip(grid.shape, factors, strict=True))
    P_nodes = sp.kron(Px, sp.kron(Py, Pz, format="csr"), format="csr")
    P = sp.kron(P_nodes, sp.identity(3, format="csr"), format="csr")
    nx, ny, nz = (n // f for n, f in zip(grid.shape, factors, strict=True))
    # Only the node layout of the coarse grid is used (Galerkin operators need no h).
    coarse = VoxelGrid(nx, ny, nz, h=grid.h * max(factors))
    return P, coarse


@dataclass
class _Level:
    A: sp.csr_matrix
    inv_diag: FloatArray
    P: sp.csr_matrix | None  # prolongation from the next coarser level (free DOFs)


class MultigridPreconditioner:
    """Build once per grid and free-DOF set; call ``setup`` with each new K_ff."""

    def __init__(self, grid: VoxelGrid, free: IntArray, min_coarse_dofs: int = 3000,
                 max_levels: int = 8, sweeps: int = 2, omega: float = 0.6) -> None:
        self.sweeps = sweeps
        self.omega = omega
        self.P_free: list[sp.csr_matrix] = []
        g = grid
        rows = free
        for _ in range(max_levels - 1):
            if 3 * g.n_nodes <= min_coarse_dofs:
                break
            try:
                P, gc = prolongation(g)
            except ValueError:
                break
            Pf = P[rows]
            keep = np.flatnonzero(np.asarray(abs(Pf).sum(axis=0)).ravel() > 0)
            self.P_free.append(Pf[:, keep].tocsr())
            rows = keep  # coarse "free" DOFs are those reachable from fine free DOFs
            g = gc
        self.levels: list[_Level] = []
        self._coarse_lu: spla.SuperLU | None = None

    @property
    def n_levels(self) -> int:
        return len(self.P_free) + 1

    def setup(self, A: sp.csr_matrix) -> None:
        self.levels = []
        Ak = A.tocsr()
        for P in self.P_free:
            self.levels.append(_Level(Ak, 1.0 / Ak.diagonal(), P))
            Ak = (P.T @ Ak @ P).tocsr()
        self.levels.append(_Level(Ak, 1.0 / Ak.diagonal(), None))
        self._coarse_lu = spla.splu(Ak.tocsc(), permc_spec="MMD_AT_PLUS_A")

    def _vcycle(self, k: int, b: FloatArray) -> FloatArray:
        lv = self.levels[k]
        if lv.P is None:
            assert self._coarse_lu is not None
            return self._coarse_lu.solve(b)
        x = self.omega * lv.inv_diag * b
        for _ in range(self.sweeps - 1):
            x += self.omega * lv.inv_diag * (b - lv.A @ x)
        r = b - lv.A @ x
        x += lv.P @ self._vcycle(k + 1, lv.P.T @ r)
        for _ in range(self.sweeps):
            x += self.omega * lv.inv_diag * (b - lv.A @ x)
        return x

    def apply(self, r: FloatArray) -> FloatArray:
        return self._vcycle(0, r)

    def as_operator(self) -> spla.LinearOperator:
        n = self.levels[0].A.shape[0]
        return spla.LinearOperator((n, n), matvec=self.apply)
