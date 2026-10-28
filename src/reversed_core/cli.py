"""Command line: run a JSON job spec and write results.

    reversed-solve examples/cantilever.json --out results/cantilever
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .export import density_to_mesh, write_stl, write_vtk
from .problem import Problem
from .simp import IterationRecord, optimize


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Voxel SIMP topology optimization")
    ap.add_argument("job", type=Path, help="JSON job specification")
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--no-stl", action="store_true", help="skip surface export")
    args = ap.parse_args(argv)

    problem = Problem.from_dict(json.loads(args.job.read_text()))
    args.out.mkdir(parents=True, exist_ok=True)

    def report(rec: IterationRecord, _x: np.ndarray) -> None:
        print(f"it {rec.iteration:4d}  c {rec.compliance:12.5g}  vol {rec.volume:.4f}  "
              f"change {rec.change:.4f}  beta {rec.beta:5.1f}  {rec.seconds:6.2f}s")

    result = optimize(problem, callback=report)
    np.save(args.out / "x_phys.npy", result.x_phys)
    write_vtk(args.out / "x_phys.vtk", result.x_phys, problem.h)
    summary = {
        "compliance": result.compliance,
        "iterations": result.iterations,
        "converged": result.converged,
        "volume": float(result.x_phys.mean()),
        "history": [r.__dict__ for r in result.history],
    }
    if not args.no_stl:
        verts, faces = density_to_mesh(result.x_phys, problem.h)
        summary["triangles"] = write_stl(args.out / "surface.stl", verts, faces)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"done: compliance {result.compliance:.6g}, {result.iterations} iterations, "
          f"converged={result.converged}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
