"""Design updates for one volume constraint: optimality criteria (OC) and MMA.

OC [BS, A]:
    x_new = clip( x * sqrt(-dc / (lambda dv)),  x - m, x + m ) clipped to [0, 1]
lambda is found by bisection so the projected volume matches the target. With a
Heaviside projection the volume is evaluated on the projected field, as in the
'ft = 3' variant of top88 [A]. OC is a heuristic fixed-point scheme and can fall into
period-2 oscillation when a projection steepens the response; MMA is the robust choice
there.

MMA, the method of moving asymptotes [Sv87, Sv07]: each function is replaced by a
convex separable approximation
    f~(x) = r + sum_j p_j / (U_j - x_j) + q_j / (x_j - L_j)
and the subproblem is solved through its dual. With a single constraint the dual is a
one-dimensional concave maximization in lambda >= 0, solved here by bisection on the
dual gradient (the approximated constraint value).

References
----------
[BS]   Bendsoe, M.P., Sigmund, O. (2003). Topology Optimization: Theory, Methods and
       Applications. Springer. (optimality criteria method)
[A]    Andreassen, E. et al. (2011). Efficient topology optimization in MATLAB using
       88 lines of code. Struct. Multidisc. Optim. 43, 1-16.
[Sv87] Svanberg, K. (1987). The method of moving asymptotes: a new method for
       structural optimization. Int. J. Numer. Methods Eng. 24, 359-373.
[Sv07] Svanberg, K. (2007). MMA and GCMMA, two methods for nonlinear optimization.
       Technical note, KTH Royal Institute of Technology, Stockholm.
       (asymptote update rule, p/q coefficients with the 1.001 / 0.001 split)
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def oc_update(
    x: FloatArray,
    dc: FloatArray,
    dv: FloatArray,
    volume_of: Callable[[FloatArray], float],
    target: float,
    move: float = 0.2,
    fixed_solid: NDArray[np.int64] | None = None,
    fixed_void: NDArray[np.int64] | None = None,
    tol: float = 1e-4,
    max_bisect: int = 200,
) -> FloatArray:
    """One OC step. ``volume_of`` maps design variables to the constrained volume fraction."""
    # Compliance sensitivities are non-positive in theory; clip round-off. [BS]
    ratio = np.maximum(-dc, 0.0) / np.maximum(dv, 1e-30)
    lo_x = np.maximum(0.0, x - move)
    hi_x = np.minimum(1.0, x + move)

    def candidate(lmid: float) -> FloatArray:
        xn = np.clip(x * np.sqrt(ratio / lmid), lo_x, hi_x)
        if fixed_solid is not None:
            xn[fixed_solid] = 1.0
        if fixed_void is not None:
            xn[fixed_void] = 0.0
        return xn

    l1, l2 = 0.0, 1e9
    xn = x
    for _ in range(max_bisect):
        if (l2 - l1) / (l1 + l2) <= tol:
            break
        lmid = 0.5 * (l1 + l2)
        xn = candidate(lmid)
        if volume_of(xn) > target:
            l1 = lmid
        else:
            l2 = lmid
    return xn


class MMA:
    """Method of moving asymptotes for  min f(x)  s.t.  g(x) <= 0,  xmin <= x <= xmax.

    Stateful: keeps the two previous iterates and the asymptotes. Call ``reset`` when the
    problem changes abruptly (for example after a projection beta step).
    """

    def __init__(self, n: int, move: float = 0.2, xmin: float = 0.0, xmax: float = 1.0,
                 asyinit: float = 0.5, asyincr: float = 1.2, asydecr: float = 0.7) -> None:
        self.n = n
        self.move = move
        self.xmin: FloatArray = np.full(n, xmin, dtype=np.float64)
        self.xmax: FloatArray = np.full(n, xmax, dtype=np.float64)
        self.asyinit, self.asyincr, self.asydecr = asyinit, asyincr, asydecr
        self.reset()

    def reset(self) -> None:
        self.k = 0
        self.xold1: FloatArray | None = None
        self.xold2: FloatArray | None = None
        self.low: FloatArray = self.xmin.copy()
        self.upp: FloatArray = self.xmax.copy()

    def _asymptotes(self, x: FloatArray) -> None:
        span = self.xmax - self.xmin
        if self.k <= 2 or self.xold1 is None or self.xold2 is None:
            self.low = x - self.asyinit * span
            self.upp = x + self.asyinit * span
        else:
            # Oscillation pulls the asymptotes in; steady progress pushes them out. [Sv07]
            sign = (x - self.xold1) * (self.xold1 - self.xold2)
            gamma = np.where(sign > 0, self.asyincr, np.where(sign < 0, self.asydecr, 1.0))
            self.low = x - gamma * (self.xold1 - self.low)
            self.upp = x + gamma * (self.upp - self.xold1)
            self.low = np.clip(self.low, x - 10.0 * span, x - 0.01 * span)
            self.upp = np.clip(self.upp, x + 0.01 * span, x + 10.0 * span)

    def _pq(self, x: FloatArray, df: FloatArray) -> tuple[FloatArray, FloatArray]:
        span = np.maximum(self.xmax - self.xmin, 1e-12)
        dpos = np.maximum(df, 0.0)
        dneg = np.maximum(-df, 0.0)
        eps = 1e-5 / span
        p = (self.upp - x) ** 2 * (1.001 * dpos + 0.001 * dneg + eps)
        q = (x - self.low) ** 2 * (0.001 * dpos + 1.001 * dneg + eps)
        return p, q

    def update(self, x: FloatArray, df: FloatArray, g: float, dg: FloatArray,
               fixed_solid: NDArray[np.int64] | None = None,
               fixed_void: NDArray[np.int64] | None = None) -> FloatArray:
        self.k += 1
        self._asymptotes(x)
        span = self.xmax - self.xmin
        alpha = np.maximum(np.maximum(self.xmin, self.low + 0.1 * (x - self.low)),
                           x - self.move * span)
        beta = np.minimum(np.minimum(self.xmax, self.upp - 0.1 * (self.upp - x)),
                          x + self.move * span)
        p0, q0 = self._pq(x, df)
        p1, q1 = self._pq(x, dg)
        r1 = g - np.sum(p1 / (self.upp - x) + q1 / (x - self.low))

        def x_of(lam: float) -> FloatArray:
            sp_ = np.sqrt(p0 + lam * p1)
            sq_ = np.sqrt(q0 + lam * q1)
            return np.clip((sp_ * self.low + sq_ * self.upp) / (sp_ + sq_), alpha, beta)

        def g_approx(xx: FloatArray) -> float:
            return float(r1 + np.sum(p1 / (self.upp - xx) + q1 / (xx - self.low)))

        # Dual gradient g~(x(lambda)) is non-increasing in lambda; bisect for its root.
        lam_lo, lam_hi = 0.0, 1.0
        if g_approx(x_of(0.0)) <= 0.0:
            x_new = x_of(0.0)
        else:
            while g_approx(x_of(lam_hi)) > 0.0 and lam_hi < 1e12:
                lam_hi *= 10.0
            for _ in range(100):
                lam = 0.5 * (lam_lo + lam_hi)
                if g_approx(x_of(lam)) > 0.0:
                    lam_lo = lam
                else:
                    lam_hi = lam
                if lam_hi - lam_lo <= 1e-12 * max(1.0, lam_hi):
                    break
            x_new = x_of(lam_hi)

        if fixed_solid is not None:
            x_new[fixed_solid] = 1.0
        if fixed_void is not None:
            x_new[fixed_void] = 0.0
        self.xold2 = self.xold1
        self.xold1 = x.copy()
        return x_new
