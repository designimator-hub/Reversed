"""End-to-end optimization checks: constraint satisfaction, monotone progress,
symmetry, passive regions, projection crispness, and the OC update in isolation."""

import json

import numpy as np
import pytest

from reversed_core import Box, Fixture, Load, LoadCase, Problem, Settings, optimize
from reversed_core.optimizers import MMA, oc_update


def cantilever(nelx=24, nely=8, nelz=4, **settings) -> Problem:
    """Liu & Tovar (2014) style cantilever: clamped at x = 0, line load at the free end."""
    s = dict(volume_fraction=0.3, max_iter=60)
    s.update(settings)
    return Problem(
        nelx=nelx, nely=nely, nelz=nelz, h=1.0,
        fixtures=(Fixture(Box((0, 0, 0), (0, nely, nelz))),),
        load_cases=(LoadCase((Load(Box((nelx, 0, 0), (nelx, 0, nelz)), (0, -1, 0)),)),),
        settings=Settings(**s),
    )


def test_oc_meets_volume_and_move_limit():
    rng = np.random.default_rng(7)
    x = rng.uniform(0.2, 0.6, 500)
    dc = -rng.uniform(0.1, 10, 500)
    dv = np.full(500, 1 / 500)
    xn = oc_update(x, dc, dv, lambda v: float(v.mean()), 0.4, move=0.1)
    assert xn.mean() == pytest.approx(0.4, abs=1e-3)
    assert np.all(np.abs(xn - x) <= 0.1 + 1e-12)
    assert xn.min() >= 0 and xn.max() <= 1


@pytest.mark.parametrize("method", ["oc", "mma"])
def test_update_solves_analytic_problem(method):
    """min sum a_j / x_j  s.t. mean(x) <= v has the closed-form optimum
    x_j = n v sqrt(a_j) / sum_k sqrt(a_k)  (Lagrange stationarity: a_j / x_j^2 = lambda / n).
    This is the compliance-like model problem behind both OC and MMA (Svanberg 1987)."""
    rng = np.random.default_rng(8)
    n, v = 40, 0.4
    a = rng.uniform(0.5, 2.0, n)
    exact = n * v * np.sqrt(a) / np.sqrt(a).sum()
    x = np.full(n, v)
    opt = MMA(n, move=0.2, xmin=1e-3)
    for _ in range(200):
        df = -a / x**2
        dv = np.full(n, 1.0 / n)
        if method == "oc":
            x = oc_update(x, df, dv, lambda z: float(z.mean()), v, move=0.2)
        else:
            x = opt.update(x, df / abs(df).max(), float(x.mean()) / v - 1.0, dv / v)
    assert np.allclose(x, exact, rtol=2e-3)


def test_mma_and_oc_reach_comparable_local_optima():
    """Compliance minimization is non-convex, so two schemes may stop in different local
    optima. Measured on this problem (2026-10-28): MMA 383.3 on every platform; OC between
    350.7 and 382.1 depending on round-off (see test below and docs/solver.md)."""
    r_oc = optimize(cantilever(24, 8, 2, optimizer="oc", max_iter=300))
    r_mma = optimize(cantilever(24, 8, 2, optimizer="mma", max_iter=300))
    assert r_oc.converged and r_mma.converged
    assert r_mma.compliance == pytest.approx(r_oc.compliance, rel=0.12)
    for r in (r_oc, r_mma):
        assert r.compliance < 0.2 * r.history[0].compliance


def test_mma_is_stable_under_tiny_perturbations():
    """MMA's result must not depend on round-off: a 1e-10 change to the start must not move
    the converged compliance. (OC fails this on the same problem, measured spread 9 %.)"""
    p = cantilever(24, 8, 2, optimizer="mma", max_iter=300)
    ne = p.grid.n_elements
    base = optimize(p)
    nudged = optimize(p, x0=0.3 + 1e-10 * np.random.default_rng(14).standard_normal(ne))
    assert nudged.compliance == pytest.approx(base.compliance, rel=1e-4)


@pytest.fixture(scope="module")
def cantilever_result():
    return optimize(cantilever())


def test_volume_constraint_satisfied(cantilever_result):
    assert cantilever_result.x_phys.mean() == pytest.approx(0.3, abs=2e-3)


def test_compliance_drops_substantially(cantilever_result):
    h = cantilever_result.history
    assert cantilever_result.compliance < 0.25 * h[0].compliance
    # After the first 10 iterations OC with filtering should be close to monotone.
    tail = np.array([r.compliance for r in h[10:]])
    assert np.all(np.diff(tail) <= 1e-3 * tail[:-1])


def test_symmetric_problem_gives_symmetric_design():
    """Load and supports symmetric about z mid-plane -> design symmetric about it."""
    r = optimize(cantilever(16, 6, 4, max_iter=25))
    assert np.allclose(r.x_phys, r.x_phys[:, :, ::-1], atol=1e-6)


def test_passive_regions_respected():
    p = cantilever(16, 6, 2, max_iter=10)
    p = Problem(**{**p.__dict__,
                   "keep": (Box((14, 0, 0), (16, 2, 2)),),
                   "avoid": (Box((6, 2, 0), (10, 4, 2)),)})
    r = optimize(p)
    g = p.grid
    s, v = p.passive_elements()
    assert np.all(r.x_phys.ravel()[s] == 1.0)
    assert np.all(r.x_phys.ravel()[v] == 0.0)
    assert s.size == 2 * 2 * 2 and v.size == 4 * 2 * 2
    assert g.n_elements == 16 * 6 * 2


@pytest.mark.slow
def test_projection_gives_crisp_connected_design():
    """Regression test: before the beta-scaled move limit, a beta step disconnected the
    load path and compliance rose about 40-fold (found during development, 2026-10-28)."""
    r_plain = optimize(cantilever(24, 8, 2, max_iter=300))
    r = optimize(cantilever(24, 8, 2, projection=True, beta_max=32, max_iter=500))
    gray = np.mean((r.x_phys > 0.1) & (r.x_phys < 0.9))
    assert r.converged
    assert gray < 0.05
    assert r.x_phys.mean() == pytest.approx(0.3, abs=5e-3)
    # Crisp 0/1 designs have no penalized gray material, so they should not be worse.
    assert r.compliance < 1.05 * r_plain.compliance


def test_job_spec_round_trip(tmp_path):
    p = cantilever()
    path = tmp_path / "job.json"
    path.write_text(json.dumps(p.to_dict()))
    assert Problem.from_dict(json.loads(path.read_text())) == p


def test_bad_specs_rejected():
    p = cantilever()
    with pytest.raises(ValueError, match="selects no nodes"):
        Problem(**{**p.__dict__, "fixtures": (Fixture(Box((50, 50, 50), (60, 60, 60))),)}
                ).fixed_dofs()
    with pytest.raises(ValueError, match="exceed the volume fraction"):
        Problem(**{**p.__dict__, "keep": (Box((0, 0, 0), (24, 8, 4)),)}).passive_elements()
    with pytest.raises(ValueError):
        Box((1, 0, 0), (0, 1, 1))
