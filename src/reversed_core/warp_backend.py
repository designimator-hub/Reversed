"""Matrix-free finite-element solve with NVIDIA Warp (GPU when available, CPU otherwise).

The global stiffness matrix is never stored. Each product y = K x loops over elements:
    y[edof_e] += E_e * K0 @ x[edof_e]
with atomic scatter, so memory is O(n_elements) instead of O(nnz) [WDW]. Fixed DOFs are
handled by masking: K_masked = M K M + (I - M), with M the 0/1 free-DOF mask, which keeps
the operator symmetric positive definite and leaves u = 0 on fixed DOFs when f = 0 there.

Precision: float64 for validation. Consumer GPUs run float64 at 1/32 to 1/64 of the
float32 rate, so float32 operators are supported for speed; CG in float32 reaches
relative residuals near 1e-6, which is adequate for design iterations [WDW].

Preconditioner: Jacobi (diagonal) for now. A matrix-free geometric multigrid V-cycle,
following [WDW] and [AAL], is the next step and will reuse multigrid.py's hierarchy.

References
----------
[WDW] Wu, J., Dick, C., Westermann, R. (2016). A system for high-resolution topology
      optimization. IEEE Trans. Vis. Comput. Graph. 22(3), 1195-1208.
[AAL] Amir, O., Aage, N., Lazarov, B.S. (2014). Struct. Multidisc. Optim. 49, 815-829.
[NW]  NVIDIA Warp, https://github.com/NVIDIA/warp (Apache-2.0).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import warp as wp
import warp.optim.linear as wpl
from numpy.typing import NDArray

from .fe import SolveStats
from .mesh import VoxelGrid

FloatArray = NDArray[np.float64]


@wp.kernel
def _element_matvec(
    edof: wp.array2d(dtype=wp.int32),  # type: ignore[valid-type]
    Ee: wp.array(dtype=Any),  # type: ignore[valid-type]
    KE: wp.array2d(dtype=Any),  # type: ignore[valid-type]
    mask: wp.array(dtype=Any),  # type: ignore[valid-type]
    x: wp.array(dtype=Any),  # type: ignore[valid-type]
    y: wp.array(dtype=Any),  # type: ignore[valid-type]
):
    e = wp.tid()
    s = Ee[e]
    for a in range(24):
        acc = s * x.dtype(0.0)
        for b in range(24):
            db = edof[e, b]
            acc += KE[a, b] * x[db] * mask[db]
        da = edof[e, a]
        wp.atomic_add(y, da, s * acc * mask[da])


@wp.kernel
def _finish_matvec(
    Kx: wp.array(dtype=Any),  # type: ignore[valid-type]
    mask: wp.array(dtype=Any),  # type: ignore[valid-type]
    x: wp.array(dtype=Any),  # type: ignore[valid-type]
    y: wp.array(dtype=Any),  # type: ignore[valid-type]
    z: wp.array(dtype=Any),  # type: ignore[valid-type]
    alpha: Any,
    beta: Any,
):
    i = wp.tid()
    one = x.dtype(1.0)
    z[i] = alpha * (Kx[i] + (one - mask[i]) * x[i]) + beta * y[i]


@wp.kernel
def _element_diagonal(
    edof: wp.array2d(dtype=wp.int32),  # type: ignore[valid-type]
    Ee: wp.array(dtype=Any),  # type: ignore[valid-type]
    KE: wp.array2d(dtype=Any),  # type: ignore[valid-type]
    d: wp.array(dtype=Any),  # type: ignore[valid-type]
):
    e = wp.tid()
    for a in range(24):
        wp.atomic_add(d, edof[e, a], Ee[e] * KE[a, a])


@wp.kernel
def _jacobi_apply(
    inv_d: wp.array(dtype=Any),  # type: ignore[valid-type]
    x: wp.array(dtype=Any),  # type: ignore[valid-type]
    y: wp.array(dtype=Any),  # type: ignore[valid-type]
    z: wp.array(dtype=Any),  # type: ignore[valid-type]
    alpha: Any,
    beta: Any,
):
    i = wp.tid()
    z[i] = alpha * inv_d[i] * x[i] + beta * y[i]


class WarpSolver:
    """Matrix-free Jacobi-PCG on a voxel grid. ``solve(Ee, F)`` mirrors LinearSolver."""

    def __init__(self, grid: VoxelGrid, KE: FloatArray, free: NDArray[np.int64],
                 device: str | None = None, dtype: type = wp.float64,
                 rtol: float = 1e-8, maxiter: int = 20000) -> None:
        wp.init()
        self.device = wp.get_device(device)
        self.dtype = dtype
        self.np_dtype = np.float64 if dtype == wp.float64 else np.float32
        self.n = grid.n_dofs
        self.ne = grid.n_elements
        self.rtol = rtol
        self.maxiter = maxiter
        mask = np.zeros(self.n, dtype=self.np_dtype)
        mask[free] = 1.0
        self.free = free
        self.edof: Any = wp.array(grid.edof.astype(np.int32), dtype=wp.int32, device=self.device)
        self.KE: Any = wp.array(KE.astype(self.np_dtype), dtype=dtype, device=self.device)
        self.mask: Any = wp.array(mask, dtype=dtype, device=self.device)
        self.Ee = wp.zeros(self.ne, dtype=dtype, device=self.device)
        self._Kx = wp.zeros(self.n, dtype=dtype, device=self.device)
        self._inv_d = wp.zeros(self.n, dtype=dtype, device=self.device)
        self._last: FloatArray | None = None

    def set_moduli(self, Ee: FloatArray) -> None:
        self.Ee.assign(wp.array(Ee.astype(self.np_dtype), dtype=self.dtype, device=self.device))
        d = wp.zeros(self.n, dtype=self.dtype, device=self.device)
        wp.launch(_element_diagonal, dim=self.ne, inputs=[self.edof, self.Ee, self.KE, d],
                  device=self.device)
        dn = d.numpy()
        m = self.mask.numpy()
        inv = np.where(m > 0, 1.0 / np.where(dn > 0, dn, 1.0), 1.0).astype(self.np_dtype)
        self._inv_d.assign(wp.array(inv, dtype=self.dtype, device=self.device))

    def matvec(self, x: wp.array, y: wp.array, z: wp.array, alpha: Any, beta: Any) -> None:
        self._Kx.zero_()
        wp.launch(_element_matvec, dim=self.ne,
                  inputs=[self.edof, self.Ee, self.KE, self.mask, x, self._Kx],
                  device=self.device)
        wp.launch(_finish_matvec, dim=self.n,
                  inputs=[self._Kx, self.mask, x, y, z, self.dtype(alpha), self.dtype(beta)],
                  device=self.device)

    def _precond(self, x: wp.array, y: wp.array, z: wp.array, alpha: Any, beta: Any) -> None:
        wp.launch(_jacobi_apply, dim=self.n,
                  inputs=[self._inv_d, x, y, z, self.dtype(alpha), self.dtype(beta)],
                  device=self.device)

    def apply_K(self, x: FloatArray) -> FloatArray:
        """K_masked @ x on the host, for testing."""
        xa: Any = wp.array(x.astype(self.np_dtype), dtype=self.dtype, device=self.device)
        z = wp.zeros_like(xa)
        self.matvec(xa, z, z, 1.0, 0.0)
        return z.numpy().astype(np.float64)

    def solve(self, Ee: FloatArray, F: FloatArray) -> tuple[FloatArray, list[SolveStats]]:
        self.set_moduli(Ee)
        A = wpl.LinearOperator((self.n, self.n), self.dtype, self.device, self.matvec)
        M = wpl.LinearOperator((self.n, self.n), self.dtype, self.device, self._precond)
        F2 = F[:, None] if F.ndim == 1 else F
        U = np.zeros_like(F2, dtype=np.float64)
        stats: list[SolveStats] = []
        for c in range(F2.shape[1]):
            f = np.zeros(self.n, dtype=self.np_dtype)
            f[self.free] = F2[self.free, c]
            b: Any = wp.array(f, dtype=self.dtype, device=self.device)
            x0 = np.zeros(self.n, dtype=self.np_dtype)
            if self._last is not None and self._last.shape == F2.shape:
                x0[self.free] = self._last[self.free, c]
            x: Any = wp.array(x0, dtype=self.dtype, device=self.device)
            out: Any = wpl.cg(A, b, x, tol=self.rtol, maxiter=self.maxiter, M=M,
                              check_every=10, use_cuda_graph=False)
            iters, err = out[0], out[1]
            U[:, c] = x.numpy().astype(np.float64)
            stats.append(SolveStats("warp-jacobi-cg", int(iters), float(err)))
        self._last = U.copy()
        return (U[:, 0] if F.ndim == 1 else U), stats
