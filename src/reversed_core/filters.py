"""Density filter and Heaviside projection.

Density filter [B]:
    x~_e = sum_i w_ei x_i / sum_i w_ei,   w_ei = max(0, r_min - |c_e - c_i|)
Smoothed Heaviside projection [W]:
    x_bar = (tanh(beta eta) + tanh(beta (x~ - eta))) / (tanh(beta eta) + tanh(beta (1 - eta)))

References
----------
[B] Bourdin, B. (2001). Filters in topology optimization.
    Int. J. Numer. Methods Eng. 50, 2143-2158.
[W] Wang, F., Lazarov, B.S., Sigmund, O. (2011). On projection methods, convergence and
    robust formulations in topology optimization. Struct. Multidisc. Optim. 43, 767-784.
"""

from __future__ import annotations

import itertools

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray

from .mesh import VoxelGrid

FloatArray = NDArray[np.float64]


class DensityFilter:
    """Linear cone filter on a structured grid. r_min is in element units (voxel edges)."""

    def __init__(self, grid: VoxelGrid, r_min: float) -> None:
        if r_min <= 0:
            raise ValueError("filter radius must be positive")
        self.grid = grid
        self.r_min = r_min
        self.H = self._build(grid, r_min)
        self.Hs = np.asarray(self.H.sum(axis=1)).ravel()

    @staticmethod
    def _build(grid: VoxelGrid, r_min: float) -> sp.csr_matrix:
        nx, ny, nz = grid.shape
        reach = int(np.ceil(r_min)) - 1
        gi, gj, gk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
        i, j, k = gi.ravel(), gj.ravel(), gk.ravel()
        e = (i * ny + j) * nz + k
        rows, cols, vals = [], [], []
        rng = range(-reach, reach + 1)
        for di, dj, dk in itertools.product(rng, rng, rng):
            w = r_min - np.sqrt(di * di + dj * dj + dk * dk)
            if w <= 0:
                continue
            ii, jj, kk = i + di, j + dj, k + dk
            ok = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny) & (kk >= 0) & (kk < nz)
            rows.append(e[ok])
            cols.append(((ii * ny + jj) * nz + kk)[ok])
            vals.append(np.full(int(ok.sum()), w))
        n = grid.n_elements
        return sp.coo_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        ).tocsr()

    def apply(self, x: FloatArray) -> FloatArray:
        return (self.H @ x) / self.Hs

    def backprop(self, g: FloatArray) -> FloatArray:
        """Chain rule: d f / d x = H^T (g / Hs) where g = d f / d x~."""
        return self.H.T @ (g / self.Hs)


def heaviside(x: FloatArray, beta: float, eta: float = 0.5) -> FloatArray:
    """Smoothed Heaviside projection [W]. Maps 0 -> 0 and 1 -> 1 for any beta > 0."""
    tb = np.tanh(beta * eta)
    return (tb + np.tanh(beta * (x - eta))) / (tb + np.tanh(beta * (1.0 - eta)))


def heaviside_grad(x: FloatArray, beta: float, eta: float = 0.5) -> FloatArray:
    """d x_bar / d x~ for the projection above."""
    tb = np.tanh(beta * eta)
    return beta * (1.0 - np.tanh(beta * (x - eta)) ** 2) / (tb + np.tanh(beta * (1.0 - eta)))
