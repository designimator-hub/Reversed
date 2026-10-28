"""Structured voxel grid: node/element numbering and region selection.

Elements are indexed (i, j, k) over (nelx, nely, nelz) in C order:
    e = (i * nely + j) * nelz + k
Nodes are indexed (i, j, k) over (nelx+1, nely+1, nelz+1) in C order:
    n = (i * (nely + 1) + j) * (nelz + 1) + k
Degree of freedom d (0=x, 1=y, 2=z) of node n is 3 * n + d.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import NDArray

from .element import NODE_SIGNS

IntArray = NDArray[np.int64]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class VoxelGrid:
    nelx: int
    nely: int
    nelz: int
    h: float = 1.0  # voxel edge length

    def __post_init__(self) -> None:
        if min(self.nelx, self.nely, self.nelz) < 1:
            raise ValueError("grid needs at least one element in each direction")
        if self.h <= 0:
            raise ValueError("voxel size must be positive")

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.nelx, self.nely, self.nelz)

    @property
    def node_shape(self) -> tuple[int, int, int]:
        return (self.nelx + 1, self.nely + 1, self.nelz + 1)

    @property
    def n_elements(self) -> int:
        return self.nelx * self.nely * self.nelz

    @property
    def n_nodes(self) -> int:
        return (self.nelx + 1) * (self.nely + 1) * (self.nelz + 1)

    @property
    def n_dofs(self) -> int:
        return 3 * self.n_nodes

    @property
    def extent(self) -> FloatArray:
        return np.array(self.shape, dtype=np.float64) * self.h

    def node_id(self, i: IntArray, j: IntArray, k: IntArray) -> IntArray:
        return (i * (self.nely + 1) + j) * (self.nelz + 1) + k

    @cached_property
    def edof(self) -> IntArray:
        """(n_elements, 24) global DOF indices per element, in element.NODE_SIGNS order."""
        gi, gj, gk = np.meshgrid(
            np.arange(self.nelx), np.arange(self.nely), np.arange(self.nelz), indexing="ij"
        )
        i, j, k = gi.ravel(), gj.ravel(), gk.ravel()
        offs = (NODE_SIGNS > 0).astype(np.int64)
        nodes = np.stack(
            [self.node_id(i + o[0], j + o[1], k + o[2]) for o in offs], axis=1
        )  # (ne, 8)
        dofs = 3 * nodes[:, :, None] + np.arange(3)[None, None, :]
        return dofs.reshape(self.n_elements, 24)

    @cached_property
    def node_coords(self) -> FloatArray:
        g = np.meshgrid(*(np.arange(n) for n in self.node_shape), indexing="ij")
        return np.stack([a.ravel() for a in g], axis=1).astype(np.float64) * self.h

    @cached_property
    def element_centers(self) -> FloatArray:
        g = np.meshgrid(*(np.arange(n) for n in self.shape), indexing="ij")
        return (np.stack([a.ravel() for a in g], axis=1) + 0.5).astype(np.float64) * self.h

    def nodes_in_box(self, lo: FloatArray, hi: FloatArray) -> IntArray:
        """Indices of nodes whose coordinates lie in the closed box [lo, hi].

        A tolerance of 1e-6 h absorbs round-off so faces at exact grid lines are caught.
        """
        tol = 1e-6 * self.h
        c = self.node_coords
        mask = np.all((c >= np.asarray(lo) - tol) & (c <= np.asarray(hi) + tol), axis=1)
        return np.flatnonzero(mask)

    def elements_in_box(self, lo: FloatArray, hi: FloatArray) -> IntArray:
        """Indices of elements whose centres lie in the closed box [lo, hi]."""
        c = self.element_centers
        mask = np.all((c >= np.asarray(lo)) & (c <= np.asarray(hi)), axis=1)
        return np.flatnonzero(mask)
