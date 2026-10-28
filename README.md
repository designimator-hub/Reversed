# Reversed

A generative design engine: voxel topology optimization now, with bio-inspired lattices and multiscale structure to follow.

- Architecture: [docs/architecture.md](docs/architecture.md). The interface is hosted; computation runs on your own GPU.
- Methods survey with references: [docs/research/bio-inspired-methods.md](docs/research/bio-inspired-methods.md)
- Patent landscape (preliminary): [docs/research/patent-landscape.md](docs/research/patent-landscape.md)
- Solver details, verification and benchmarks: [docs/solver.md](docs/solver.md)

## Quick start

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest -q                                   # full verification suite
reversed-solve examples/cantilever.json --out results/cantilever
```

Outputs: `x_phys.npy` (densities), `x_phys.vtk` (ParaView), `surface.stl` (Blender, Fusion, slicers) and `summary.json` (iteration history).

## Job specification

A job is one JSON file: grid size, voxel size, fixtures, load cases, keep/avoid regions, material and settings. Regions are boxes in physical coordinates. See `examples/` and `src/reversed_core/problem.py`.

## Checks

Every push runs ruff, mypy and pytest on Python 3.10, 3.12 and 3.13, plus one end-to-end example (`.github/workflows/ci.yml`). The suite tests the code against beam theory, finite differences, analytic optima and an independent code (pyMOTO).
