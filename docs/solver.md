# Solver core: what it does, how it is checked

Package: `src/reversed_core`. This is the CPU reference path. Every later backend (GPU through NVIDIA Warp, matrix-free multigrid) must reproduce its answers within tolerance before it is trusted.

## Pipeline per iteration

| Step | Module | Method | Reference |
|---|---|---|---|
| Element stiffness | `element.py` | Trilinear hex8, 2x2x2 Gauss (exact for bricks) | Zienkiewicz, Taylor & Zhu (2013); Hughes (2000) |
| Assembly | `fe.py` | Vectorized COO with fixed sparsity pattern | Andreassen et al. (2011) |
| Linear solve | `fe.py`, `multigrid.py` | Sparse LU below 25k free DOFs, multigrid-PCG above | Amir, Aage & Lazarov (2014) |
| Linear solve, GPU | `warp_backend.py` | Matrix-free element-by-element K·u with Jacobi-PCG (float32 or float64) | Wu, Dick & Westermann (2016); NVIDIA Warp |
| Interpolation | `simp.py` | SIMP, E = E_min + x^p (E − E_min) | Bendsøe (1989); Sigmund (2001) |
| Filter | `filters.py` | Linear cone density filter | Bourdin (2001) |
| Projection (optional) | `filters.py` | Smoothed Heaviside with beta continuation | Wang, Lazarov & Sigmund (2011) |
| Update | `optimizers.py` | MMA (default) or OC | Svanberg (1987, 2007); Bendsøe & Sigmund (2003) |
| Export | `export.py` | .npy, legacy VTK, STL by marching cubes | Lorensen & Cline (1987) |

The equations are written out in each module's docstring next to the code that implements them.

## Verification (all in `tests/`, run by CI on every push)

| Check | What it proves | Result |
|---|---|---|
| Element: symmetry, exactly 6 zero eigenvalues, rigid motions in null space | Element is a valid 3D solid element | pass |
| Element: patch test with random uniform strain, ν = 0, 0.3, 0.45 | Strain energy exact to 1e-12 | pass |
| Cantilever vs Timoshenko beam theory, 3 mesh levels | Error falls ~4x per halving of h (second order); Richardson extrapolation within 0.5 % | measured: 11.2 %, 3.1 %, 0.78 % error; extrapolated 0.02 % |
| Reactions balance applied load | Global equilibrium | pass to 1e-9 |
| Direct vs Jacobi-CG vs multigrid-PCG | Solvers agree | pass to 1e-7 relative |
| Multigrid: prolongation exact for affine fields, V-cycle symmetric | Valid CG preconditioner | pass |
| Warp matrix-free K·u vs assembled matrix | GPU kernel is the same operator | agree to 1e-13 (measured 6e-16) |
| Warp solve vs direct, float64 and float32 | GPU path gives the same displacements | measured 2e-10 (f64), 8e-5 (f32) |
| Warp compliance and gradient vs CPU reference | Optimizer sees identical numbers on GPU | rel. 1e-7 |
| **Independent code: pyMOTO (TU Delft)** | Element spectrum and heterogeneous cantilever compliance | agree to 1e-10 (measured 13 digits) |
| Filter: preserves constants, adjoint consistency, weights by formula | Filter and its chain rule are correct | pass |
| Heaviside: endpoints, derivative vs finite differences | Projection gradient correct | pass |
| Compliance and volume gradients vs central finite differences, with and without projection, 1 and 2 load cases, with passive elements | Optimizer follows the true gradient | pass, rel. 1e-5 |
| OC and MMA on an analytic problem with closed-form optimum | Update schemes converge to the right point | pass, rel. 2e-3 |
| OC vs MMA on the same cantilever | Two independent schemes agree | within 2 % |
| Volume constraint, monotone compliance, symmetric design for a symmetric problem, passive regions | Optimization behaves as theory says | pass |
| Projection run is crisp and stays connected | Regression test for the bug below | pass |
| STL closed 2-manifold and placed correctly; VTK header | Export is usable downstream | pass |

### Bug found by the checks during development

Without the beta-scaled move limit, a beta step in the projection run let one design update flip whole members. The load path disconnected and compliance rose about 40-fold (16 801 against 383) while the design looked crisp. The fix scales the move limit by min(1, 4/β). The reasoning is in `simp.py`, and `test_projection_gives_crisp_connected_design` guards against it coming back.

## Benchmark run: Liu and Tovar cantilever (60x20x4, f = 0.3, p = 3, r_min = 1.5)

| Optimizer | Iterations to converge (max change ≤ 0.01) | Final compliance |
|---|---|---|
| MMA | 256 | 94.84 |
| OC | 265 | 96.62 |

MMA reaches a 1.8 % stiffer design here. Both use the same filter, solver and stopping rule.

## Local optima and round-off sensitivity (found by CI, 2026-10-28)

On the 24x8x2 cantilever, the GitHub runners (Python 3.10 and 3.12) gave OC compliance 358.8, while this container gave 382.1. Repeating the run with the starting design changed by 1e-10 gave OC results of 359.2, 350.8, 350.8, 350.7 and 350.8. MMA gave 383.32 in every case.

- OC's trajectory is chaotic here. Round-off picks the local optimum, and on this problem some of those optima are 9 % stiffer than MMA's.
- MMA is reproducible, but it can settle in a worse local optimum. On 60x20x4, MMA was the better of the two by 1.8 %.
- Consequence for the product: one run is one sample from a landscape with many optima. That makes the case for multi-start or replica-exchange outcome generation (methods doc section 6.5). `optimize(..., x0=...)` now accepts a starting design to support it.

## Measured performance (4-core cloud container, CPU only)

| Problem | Free DOFs | Sparse LU | Multigrid-PCG |
|---|---|---|---|
| 60x20x4 uniform density | 17 640 | 1.6 s (COLAMD), 1.3 s (MMD) | 0.21 s, 11 iterations |
| 60x20x4 optimized design | 17 640 | 1.0 s | 0.97 s, 56 iterations |
| 120x40x8 optimized design | 132 840 | 56.9 s | 6.3 s, 53 iterations |

Direct solves grow roughly with the square of problem size in 3D. Multigrid grows close to linearly. The GPU path should add a further 10–50x on a desktop card, according to Wu, Dick and Westermann (2016). That number still has to be measured on your machine.

## Not yet verified against a published number

The Liu and Tovar (2014) and Ferrari and Sigmund (2020) papers print benchmark results, but this session's network blocked arxiv.org and the publisher sites, so their exact compliance values have not been compared. To close this gap, either allow `arxiv.org` in the environment's network settings, or run `top3d.m` in MATLAB or Octave on the same problem and compare against `examples/cantilever.json`.

## Known limits (honest list)

- The Warp backend has been verified only on Warp's CPU device in this session (no GPU here). `scripts/gpu_check.py` verifies and times it on a CUDA machine.
- The GPU path uses Jacobi preconditioning; matrix-free multigrid on the GPU is the next step. On the CPU, multigrid needed 11 CG iterations where Jacobi needed 398 for the uniform 60x20x4 case, so a large cut in GPU iterations is expected but not yet measured.
- Assembled matrices: the index arrays, values and CSR matrix take about 20 KB per element. That allows roughly 700 000 elements in 16 GB on CPU. The matrix-free GPU version removes this limit.
- Only cubic voxels, linear elasticity, compliance objective, one volume constraint.
- No stress constraint, overhang filter or symmetry planes yet (methods doc sections 1.5 and 10).

## References

- Amir, O., Aage, N., Lazarov, B.S. (2014). On multigrid-CG for efficient topology optimization. Struct. Multidisc. Optim. 49, 815–829.
- Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Struct. Multidisc. Optim. 43, 1–16.
- Bendsøe, M.P. (1989). Optimal shape design as a material distribution problem. Struct. Optim. 1, 193–202.
- Bendsøe, M.P., Sigmund, O. (2003). Topology Optimization: Theory, Methods and Applications. Springer.
- Bourdin, B. (2001). Filters in topology optimization. Int. J. Numer. Methods Eng. 50, 2143–2158.
- Briggs, W.L., Henson, V.E., McCormick, S.F. (2000). A Multigrid Tutorial, 2nd ed. SIAM.
- Cowper, G.R. (1966). The shear coefficient in Timoshenko's beam theory. J. Appl. Mech. 33, 335–340.
- Ferrari, F., Sigmund, O. (2020). A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D. Struct. Multidisc. Optim. 62, 2211–2228.
- Hughes, T.J.R. (2000). The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Dover.
- Liu, K., Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. Struct. Multidisc. Optim. 50, 1175–1196.
- Lorensen, W.E., Cline, H.E. (1987). Marching cubes: a high resolution 3D surface construction algorithm. Computer Graphics 21(4), 163–169.
- NVIDIA Warp, Python framework for GPU simulation. https://github.com/NVIDIA/warp
- pyMOTO, Modular Topology Optimization framework in Python. https://github.com/aatmdelissen/pyMOTO
- Sigmund, O. (2001). A 99 line topology optimization code written in Matlab. Struct. Multidisc. Optim. 21, 120–127.
- Svanberg, K. (1987). The method of moving asymptotes: a new method for structural optimization. Int. J. Numer. Methods Eng. 24, 359–373.
- Svanberg, K. (2007). MMA and GCMMA, two methods for nonlinear optimization. Technical note, KTH, Stockholm.
- Timoshenko, S.P., Gere, J.M. (1972). Mechanics of Materials. Van Nostrand Reinhold.
- Wang, F., Lazarov, B.S., Sigmund, O. (2011). On projection methods, convergence and robust formulations in topology optimization. Struct. Multidisc. Optim. 43, 767–784.
- Wu, J., Dick, C., Westermann, R. (2016). A system for high-resolution topology optimization. IEEE Trans. Vis. Comput. Graph. 22(3), 1195–1208.
- Zienkiewicz, O.C., Taylor, R.L., Zhu, J.Z. (2013). The Finite Element Method: Its Basis and Fundamentals, 7th ed. Butterworth-Heinemann.
