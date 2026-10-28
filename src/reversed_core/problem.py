"""Job specification: the complete, serializable input for one optimization run.

Every region is an axis-aligned box in physical coordinates, [lo, hi], with the domain
spanning [0, nel * h] in each direction. A job spec is plain JSON so any run can be
repeated, compared, or sent to another compute node (docs/architecture.md).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .mesh import VoxelGrid

Vec3 = tuple[float, float, float]


@dataclass(frozen=True)
class Box:
    lo: Vec3
    hi: Vec3

    def __post_init__(self) -> None:
        if len(self.lo) != 3 or len(self.hi) != 3:
            raise ValueError("box corners need three coordinates")
        if any(a > b for a, b in zip(self.lo, self.hi, strict=True)):
            raise ValueError(f"box lo must not exceed hi: {self.lo} > {self.hi}")


@dataclass(frozen=True)
class Fixture:
    box: Box
    dofs: str = "xyz"  # which displacement components are held at zero


@dataclass(frozen=True)
class Load:
    box: Box
    force: Vec3  # total force, shared equally among the nodes in the box


@dataclass(frozen=True)
class LoadCase:
    loads: tuple[Load, ...]
    weight: float = 1.0


@dataclass(frozen=True)
class Material:
    E: float = 1.0
    nu: float = 0.3
    E_min_ratio: float = 1e-9  # void stiffness as a fraction of E, avoids a singular K


@dataclass(frozen=True)
class Settings:
    volume_fraction: float = 0.3
    penal: float = 3.0
    r_min: float = 1.5  # filter radius in voxel edges
    projection: bool = False
    beta_max: float = 32.0
    eta: float = 0.5
    move: float = 0.2
    max_iter: int = 200
    tol: float = 0.01
    solver: str = "auto"
    optimizer: str = "mma"  # "mma" or "oc"


@dataclass(frozen=True)
class Problem:
    nelx: int
    nely: int
    nelz: int
    h: float
    fixtures: tuple[Fixture, ...]
    load_cases: tuple[LoadCase, ...]
    material: Material = field(default_factory=Material)
    settings: Settings = field(default_factory=Settings)
    keep: tuple[Box, ...] = ()  # forced solid
    avoid: tuple[Box, ...] = ()  # forced void

    @property
    def grid(self) -> VoxelGrid:
        return VoxelGrid(self.nelx, self.nely, self.nelz, self.h)

    # --- derived arrays -------------------------------------------------------------

    def fixed_dofs(self) -> NDArray[np.int64]:
        g = self.grid
        out: list[NDArray[np.int64]] = []
        for fx in self.fixtures:
            nodes = g.nodes_in_box(np.array(fx.box.lo), np.array(fx.box.hi))
            if nodes.size == 0:
                raise ValueError(f"fixture box {fx.box} selects no nodes")
            for c in fx.dofs:
                if c not in "xyz":
                    raise ValueError(f"unknown dof {c!r}; use x, y, z")
                out.append(3 * nodes + "xyz".index(c))
        if not out:
            raise ValueError("problem has no fixtures; stiffness matrix would be singular")
        return np.unique(np.concatenate(out))

    def force_matrix(self) -> NDArray[np.float64]:
        g = self.grid
        F = np.zeros((g.n_dofs, len(self.load_cases)))
        for c, lc in enumerate(self.load_cases):
            for ld in lc.loads:
                nodes = g.nodes_in_box(np.array(ld.box.lo), np.array(ld.box.hi))
                if nodes.size == 0:
                    raise ValueError(f"load box {ld.box} selects no nodes")
                for d in range(3):
                    np.add.at(F[:, c], 3 * nodes + d, ld.force[d] / nodes.size)
        return F

    def passive_elements(self) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        g = self.grid
        solid = [g.elements_in_box(np.array(b.lo), np.array(b.hi)) for b in self.keep]
        void = [g.elements_in_box(np.array(b.lo), np.array(b.hi)) for b in self.avoid]
        s = np.unique(np.concatenate(solid)) if solid else np.zeros(0, np.int64)
        v = np.unique(np.concatenate(void)) if void else np.zeros(0, np.int64)
        if np.intersect1d(s, v).size:
            raise ValueError("keep and avoid regions overlap")
        if s.size / g.n_elements > self.settings.volume_fraction:
            raise ValueError("keep regions alone exceed the volume fraction")
        return s, v

    # --- JSON round trip --------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict[str, Any]) -> Problem:
        def box(b: dict[str, Any]) -> Box:
            return Box(tuple(b["lo"]), tuple(b["hi"]))

        return Problem(
            nelx=int(d["nelx"]),
            nely=int(d["nely"]),
            nelz=int(d["nelz"]),
            h=float(d.get("h", 1.0)),
            fixtures=tuple(
                Fixture(box(f["box"]), f.get("dofs", "xyz")) for f in d["fixtures"]
            ),
            load_cases=tuple(
                LoadCase(
                    tuple(Load(box(ld["box"]), tuple(ld["force"])) for ld in lc["loads"]),
                    float(lc.get("weight", 1.0)),
                )
                for lc in d["load_cases"]
            ),
            material=Material(**d.get("material", {})),
            settings=Settings(**d.get("settings", {})),
            keep=tuple(box(b) for b in d.get("keep", [])),
            avoid=tuple(box(b) for b in d.get("avoid", [])),
        )
