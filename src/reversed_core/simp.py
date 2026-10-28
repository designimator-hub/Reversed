"""SIMP compliance minimization on a voxel grid.

    min_x  c = sum_k w_k F_k^T U_k      s.t.  K(x_bar) U_k = F_k,  mean(x_bar) <= f
    E_e = E_min + x_bar_e^p (E - E_min)                                   [Be, Si]
    dc/dx_bar_e = -sum_k w_k p x_bar_e^(p-1) (E - E_min) u_ke^T K0 u_ke   [Si]
    x_bar = Heaviside(Filter(x))  (projection optional)                    [B, W]

Beta continuation for the projection doubles beta every 50 iterations or when the
design change drops below the tolerance, up to beta_max, following top88 [A].
Move limit under projection: the steepest slope of the projection is
    d x_bar / d x~ |_(x~ = eta = 0.5) = beta / (2 tanh(beta / 2))  ~  beta / 2,
so a fixed move limit lets one step flip whole members once beta is large, which
disconnects load paths. The move limit is therefore scaled by min(1, 4 / beta), which
caps the largest projected change per step at about twice the nominal move limit.
This is a stabilizing heuristic, verified by tests/test_optimize.py.
With MMA [Sv87] the objective is scaled by its first value and the constraint is
written g = mean(x_bar) / f - 1 <= 0, so both are of order one.

References
----------
[Be] Bendsoe, M.P. (1989). Optimal shape design as a material distribution problem.
     Struct. Optim. 1, 193-202.
[Si] Sigmund, O. (2001). A 99 line topology optimization code written in Matlab.
     Struct. Multidisc. Optim. 21, 120-127.
[A]  Andreassen, E. et al. (2011). Struct. Multidisc. Optim. 43, 1-16.
[LT] Liu, K., Tovar, A. (2014). An efficient 3D topology optimization code written in
     Matlab. Struct. Multidisc. Optim. 50, 1175-1196.
[B]  Bourdin, B. (2001). Int. J. Numer. Methods Eng. 50, 2143-2158.
[W]  Wang, F., Lazarov, B.S., Sigmund, O. (2011). Struct. Multidisc. Optim. 43, 767-784.
[Sv87] Svanberg, K. (1987). Int. J. Numer. Methods Eng. 24, 359-373.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from .element import stiffness_matrix
from .fe import Assembler, LinearSolver, element_energies
from .filters import DensityFilter, heaviside, heaviside_grad
from .optimizers import MMA, oc_update
from .problem import Problem

FloatArray = NDArray[np.float64]


@dataclass
class IterationRecord:
    iteration: int
    compliance: float
    volume: float
    change: float
    beta: float
    seconds: float


@dataclass
class Result:
    x: FloatArray  # design variables, shape (nelx, nely, nelz)
    x_phys: FloatArray  # physical densities after filter and projection
    compliance: float
    history: list[IterationRecord] = field(default_factory=list)
    converged: bool = False

    @property
    def iterations(self) -> int:
        return len(self.history)


class ComplianceModel:
    """Maps design variables to compliance and its gradient. Shared by the optimizer
    and by the finite-difference tests so both exercise the same code path."""

    def __init__(self, problem: Problem) -> None:
        self.problem = problem
        self.grid = problem.grid
        m = problem.material
        self.E = m.E
        self.E_min = m.E * m.E_min_ratio
        self.KE = stiffness_matrix(1.0, m.nu, problem.h)
        self.assembler = Assembler(self.grid, self.KE)
        self.filter = DensityFilter(self.grid, problem.settings.r_min)
        fixed = problem.fixed_dofs()
        self.free = np.setdiff1d(np.arange(self.grid.n_dofs), fixed)
        self.F = problem.force_matrix()
        self.weights = np.array([lc.weight for lc in problem.load_cases])
        self.warp = None
        if problem.settings.solver.startswith("warp"):
            import warp as wp

            from .warp_backend import WarpSolver  # optional dependency: pip install .[gpu]

            dtype = wp.float32 if problem.settings.solver == "warp32" else wp.float64
            self.warp = WarpSolver(self.grid, self.KE, self.free, dtype=dtype,
                                   rtol=1e-6 if dtype == wp.float32 else 1e-8)
        else:
            self.solver = LinearSolver(self.free, problem.settings.solver,  # type: ignore[arg-type]
                                       grid=self.grid)
        self.solid, self.void = problem.passive_elements()
        self.beta = 1.0
        self.U: FloatArray | None = None

    def physical(self, x: FloatArray) -> FloatArray:
        s = self.problem.settings
        xt = self.filter.apply(x)
        xp = heaviside(xt, self.beta, s.eta) if s.projection else xt
        xp = xp.copy()
        xp[self.solid] = 1.0
        xp[self.void] = 0.0
        return xp

    def physical_grad(self, x: FloatArray, g_phys: FloatArray) -> FloatArray:
        """Back-propagate d f / d x_phys to d f / d x through projection and filter."""
        s = self.problem.settings
        g = g_phys.copy()
        g[self.solid] = 0.0
        g[self.void] = 0.0
        if s.projection:
            g = g * heaviside_grad(self.filter.apply(x), self.beta, s.eta)
        return self.filter.backprop(g)

    def evaluate(self, x: FloatArray) -> tuple[float, FloatArray, FloatArray]:
        """Returns compliance, dc/dx and x_phys."""
        p = self.problem.settings.penal
        xp = self.physical(x)
        Ee = self.E_min + xp**p * (self.E - self.E_min)
        if self.warp is not None:
            U, self.last_stats = self.warp.solve(Ee, self.F)
        else:
            U, self.last_stats = self.solver.solve(self.assembler.assemble(Ee), self.F)
        self.U = U
        ce = element_energies(U, self.grid.edof, self.KE)  # (ne, n_cases)
        c = float(np.sum(self.weights * np.sum(Ee[:, None] * ce, axis=0)))
        dc_phys = -p * xp ** (p - 1) * (self.E - self.E_min) * (ce @ self.weights)
        return c, self.physical_grad(x, dc_phys), xp


def optimize(
    problem: Problem,
    callback: Callable[[IterationRecord, FloatArray], None] | None = None,
    x0: FloatArray | None = None,
) -> Result:
    """Run SIMP. ``x0`` optionally sets the starting design (flat, one value per element),
    which allows warm starts and multi-start searches."""
    s = problem.settings
    model = ComplianceModel(problem)
    ne = model.grid.n_elements
    if x0 is not None:
        x = np.clip(np.asarray(x0, dtype=np.float64).ravel(), 0.0, 1.0)
        if x.size != ne:
            raise ValueError(f"x0 has {x.size} values, grid has {ne} elements")
    else:
        x = np.full(ne, s.volume_fraction, dtype=np.float64)
    x[model.solid] = 1.0
    x[model.void] = 0.0

    dv = model.physical_grad(x, np.full(ne, 1.0 / ne))  # linear filter: constant if no projection
    history: list[IterationRecord] = []
    loop_beta = 0
    converged = False
    c = np.inf
    xp = model.physical(x)

    if s.optimizer not in ("mma", "oc"):
        raise ValueError(f"unknown optimizer {s.optimizer!r}")
    mma = MMA(ne, move=s.move) if s.optimizer == "mma" else None
    c0: float | None = None

    for it in range(1, s.max_iter + 1):
        t0 = time.perf_counter()
        loop_beta += 1
        c, dc, xp = model.evaluate(x)
        move = s.move * min(1.0, 4.0 / model.beta) if s.projection else s.move
        if s.projection:
            dv = model.physical_grad(x, np.full(ne, 1.0 / ne))
        if mma is not None:
            mma.move = move
            if c0 is None:
                c0 = c
            g = float(xp.mean()) / s.volume_fraction - 1.0
            x_new = mma.update(x, dc / c0, g, dv / s.volume_fraction,
                               fixed_solid=model.solid, fixed_void=model.void)
        else:
            x_new = oc_update(
                x, dc, dv,
                volume_of=lambda xn: float(model.physical(xn).mean()),
                target=s.volume_fraction,
                move=move,
                fixed_solid=model.solid,
                fixed_void=model.void,
            )
        change = float(np.max(np.abs(x_new - x)))
        x = x_new
        rec = IterationRecord(it, c, float(model.physical(x).mean()), change, model.beta,
                              time.perf_counter() - t0)
        history.append(rec)
        if callback is not None:
            callback(rec, x)

        if s.projection and model.beta < s.beta_max and (loop_beta >= 50 or change <= s.tol):
            model.beta = min(2.0 * model.beta, s.beta_max)
            loop_beta = 0
            if mma is not None:
                mma.reset()  # the objective landscape changed; restart asymptotes
            continue
        if change <= s.tol and (not s.projection or model.beta >= s.beta_max):
            converged = True
            break

    c_final, _, xp = model.evaluate(x)
    shape = model.grid.shape
    return Result(x.reshape(shape), xp.reshape(shape), c_final, history, converged)
