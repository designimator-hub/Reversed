"""Finite-element assembly and linear solves on a voxel grid.

Global stiffness K = sum_e E_e * K0_e, with K0 the unit-modulus hex8 stiffness.
Only the free DOFs are solved: K_ff u_f = f_f (homogeneous Dirichlet on fixed DOFs).

References
----------
[A] Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., Sigmund, O. (2011).
    Efficient topology optimization in MATLAB using 88 lines of code.
    Struct. Multidisc. Optim. 43, 1-16. (vectorized sparse assembly with iK, jK)
[S] Saad, Y. (2003). Iterative Methods for Sparse Linear Systems, 2nd ed. SIAM.
    (preconditioned conjugate gradient)
[AAL] Amir, O., Aage, N., Lazarov, B.S. (2014). On multigrid-CG for efficient topology
    optimization. Struct. Multidisc. Optim. 49, 815-829.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from .mesh import VoxelGrid

FloatArray = NDArray[np.float64]
SolverName = Literal["auto", "direct", "cg", "mg"]

# Above this many free DOFs the "auto" solver switches from sparse LU to multigrid-PCG.
AUTO_DIRECT_LIMIT = 25_000  # measured crossover, docs/solver.md


@dataclass
class Assembler:
    """Holds the sparsity pattern so each iteration only rebuilds the values. [A]"""

    grid: VoxelGrid
    KE: FloatArray
    iK: NDArray[np.int64] = field(init=False)
    jK: NDArray[np.int64] = field(init=False)

    def __post_init__(self) -> None:
        edof = self.grid.edof
        self.iK = np.repeat(edof, 24, axis=1).ravel()
        self.jK = np.tile(edof, (1, 24)).ravel()

    def assemble(self, Ee: FloatArray) -> sp.csr_matrix:
        if Ee.shape != (self.grid.n_elements,):
            raise ValueError("one modulus per element expected")
        vals = (self.KE.ravel()[None, :] * Ee[:, None]).ravel()
        n = self.grid.n_dofs
        K = sp.coo_matrix((vals, (self.iK, self.jK)), shape=(n, n)).tocsr()
        return K


@dataclass
class SolveStats:
    method: str
    iterations: int = 0
    residual: float = 0.0


def solve_free(
    K: sp.csr_matrix,
    F: FloatArray,
    free: NDArray[np.int64],
    solver: SolverName = "auto",
    rtol: float = 1e-8,
    maxiter: int = 5000,
    x0: FloatArray | None = None,
    grid: VoxelGrid | None = None,
) -> tuple[FloatArray, list[SolveStats]]:
    """Solve K u = F with u = 0 on non-free DOFs. F may hold several load cases as columns.

    Returns U with the same shape as F. One-shot wrapper around LinearSolver.
    """
    return LinearSolver(free, solver, grid=grid, rtol=rtol, maxiter=maxiter).solve(K, F, x0)


class LinearSolver:
    """Reusable solver for a fixed set of free DOFs.

    Methods: "direct" (sparse LU), "cg" (Jacobi-PCG [S]), "mg" (multigrid-PCG [AAL], needs
    the grid), "auto" (direct for small systems, multigrid above AUTO_DIRECT_LIMIT).
    Iterative methods warm-start from the previous solution, which cuts iterations
    sharply once the design settles.
    """

    def __init__(self, free: NDArray[np.int64], method: SolverName = "auto",
                 grid: VoxelGrid | None = None, rtol: float = 1e-8, maxiter: int = 5000) -> None:
        if method == "auto":
            method = "direct" if free.size <= AUTO_DIRECT_LIMIT or grid is None else "mg"
        if method not in ("direct", "cg", "mg"):
            raise ValueError(f"unknown solver {method!r}")
        if method == "mg" and grid is None:
            raise ValueError("multigrid solver needs the grid")
        self.free = free
        self.method = method
        self.rtol = rtol
        self.maxiter = maxiter
        self.mg = None
        if method == "mg":
            from .multigrid import MultigridPreconditioner

            assert grid is not None
            self.mg = MultigridPreconditioner(grid, free)
        self._last: FloatArray | None = None

    def solve(self, K: sp.csr_matrix, F: FloatArray,
              x0: FloatArray | None = None) -> tuple[FloatArray, list[SolveStats]]:
        free = self.free
        F2 = F[:, None] if F.ndim == 1 else F
        Kff = K[free][:, free].tocsr()
        U = np.zeros_like(F2)
        stats: list[SolveStats] = []

        if self.method == "direct":
            lu = spla.splu(Kff.tocsc(), permc_spec="MMD_AT_PLUS_A")
            for c in range(F2.shape[1]):
                U[free, c] = lu.solve(F2[free, c])
                stats.append(SolveStats("direct"))
            return (U[:, 0] if F.ndim == 1 else U), stats

        if self.mg is not None:
            self.mg.setup(Kff)
            M = self.mg.as_operator()
        else:
            inv_diag = 1.0 / Kff.diagonal()
            M = spla.LinearOperator(Kff.shape, matvec=lambda r: inv_diag * r)
        start = x0 if x0 is not None else self._last
        if start is not None and start.ndim == 1:
            start = start[:, None]
        for c in range(F2.shape[1]):
            it = [0]

            def cb(_xk: FloatArray, it: list[int] = it) -> None:
                it[0] += 1

            guess = None
            if start is not None and start.shape == F2.shape:
                guess = start[free, c]
            u, info = spla.cg(Kff, F2[free, c], x0=guess, rtol=self.rtol,
                              maxiter=self.maxiter, M=M, callback=cb)
            if info != 0:
                raise RuntimeError(f"{self.method} did not converge (info={info}), case {c}")
            U[free, c] = u
            res = float(np.linalg.norm(Kff @ u - F2[free, c]) /
                        max(float(np.linalg.norm(F2[free, c])), 1e-300))
            stats.append(SolveStats(self.method, it[0], res))
        self._last = U.copy()
        return (U[:, 0] if F.ndim == 1 else U), stats


def element_energies(U: FloatArray, edof: NDArray[np.int64], KE: FloatArray) -> FloatArray:
    """u_e^T K0 u_e for every element (unit modulus). U is (n_dofs,) or (n_dofs, n_cases).

    Returns (n_elements,) or (n_elements, n_cases).
    """
    if U.ndim == 1:
        Ue = U[edof]
        return np.einsum("ei,ij,ej->e", Ue, KE, Ue)
    Ue = U[edof]  # (ne, 24, nc)
    return np.einsum("eic,ij,ejc->ec", Ue, KE, Ue)
