"""Run this on the machine with the GPU:  python scripts/gpu_check.py

It (1) reports the devices Warp sees, (2) verifies the GPU matrix-free solve against the
CPU direct solver on a small problem, and (3) times one finite-element solve at several
sizes on CPU multigrid and on the GPU in float32 and float64.
Paste the printed output back into the conversation or an issue.
"""

from __future__ import annotations

import time

import numpy as np
import warp as wp

from reversed_core.element import stiffness_matrix
from reversed_core.fe import Assembler, LinearSolver
from reversed_core.mesh import VoxelGrid
from reversed_core.warp_backend import WarpSolver


def cantilever_system(nx: int, ny: int, nz: int, seed: int = 0):
    g = VoxelGrid(nx, ny, nz)
    KE = stiffness_matrix()
    rng = np.random.default_rng(seed)
    Ee = 1e-9 + rng.uniform(0.3, 1.0, g.n_elements) ** 3
    fn = g.nodes_in_box(np.zeros(3), np.array([0, ny, nz]))
    fixed = (3 * fn[:, None] + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(g.n_dofs), fixed)
    F = np.zeros(g.n_dofs)
    tip = g.nodes_in_box(np.array([nx, 0, 0]), np.array([nx, 0, nz]))
    F[3 * tip + 1] = -1.0 / tip.size
    return g, KE, Ee, free, F


def main() -> None:
    wp.init()
    cuda = [d for d in wp.get_devices() if d.is_cuda]
    print("Warp devices:", [str(d) for d in wp.get_devices()])
    if not cuda:
        print("No CUDA device found. Check the NVIDIA driver; GPU timings skipped.")
    dev = str(cuda[0]) if cuda else "cpu"

    g, KE, Ee, free, F = cantilever_system(20, 8, 4)
    Ud, _ = LinearSolver(free, "direct").solve(Assembler(g, KE).assemble(Ee), F)
    for dtype in (wp.float64, wp.float32):
        Uw, st = WarpSolver(g, KE, free, device=dev, dtype=dtype,
                            rtol=1e-9 if dtype == wp.float64 else 1e-6).solve(Ee, F)
        err = np.linalg.norm(Uw - Ud) / np.linalg.norm(Ud)
        ok = err < (1e-7 if dtype == wp.float64 else 5e-4)
        print(f"verify {dev} {dtype.__name__}: rel. error {err:.2e}  {'PASS' if ok else 'FAIL'}")

    print(f"\n{'grid':>14} {'elements':>9} {'CPU mg [s]':>11} {dev + ' f32 [s]':>16} "
          f"{dev + ' f64 [s]':>16}")
    for nx, ny, nz in [(60, 20, 4), (120, 40, 8), (160, 80, 40), (240, 120, 60)]:
        g, KE, Ee, free, F = cantilever_system(nx, ny, nz)
        row = [f"{nx}x{ny}x{nz}", f"{g.n_elements}"]
        if g.n_elements <= 300_000:
            t = time.perf_counter()
            LinearSolver(free, "mg", grid=g).solve(Assembler(g, KE).assemble(Ee), F)
            row.append(f"{time.perf_counter() - t:.2f}")
        else:
            row.append("skipped")
        for dtype in (wp.float32, wp.float64):
            if not cuda and g.n_elements > 10_000:  # CPU fallback is slow; keep the dry run short
                row.append("skipped")
                continue
            ws = WarpSolver(g, KE, free, device=dev, dtype=dtype,
                            rtol=1e-6 if dtype == wp.float32 else 1e-8, maxiter=50_000)
            t = time.perf_counter()
            _, st = ws.solve(Ee, F)
            wp.synchronize()
            row.append(f"{time.perf_counter() - t:.2f} ({st[0].iterations} it)")
        print(f"{row[0]:>14} {row[1]:>9} {row[2]:>11} {row[3]:>16} {row[4]:>16}")


if __name__ == "__main__":
    main()
