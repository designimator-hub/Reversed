"""Adjoint sensitivities checked against central finite differences.

A mismatch here means the optimizer is following the wrong gradient, so this test
guards the whole method. Reference: Bendsoe & Sigmund (2003), Topology Optimization,
on compliance sensitivity (the problem is self-adjoint).
"""

import numpy as np
import pytest

from reversed_core import Box, Fixture, Load, LoadCase, Problem, Settings
from reversed_core.simp import ComplianceModel


def small_problem(projection: bool, two_cases: bool) -> Problem:
    cases = [LoadCase((Load(Box((4, 0, 0), (4, 0, 2)), (0, -1, 0)),))]
    if two_cases:
        cases.append(LoadCase((Load(Box((4, 2, 0), (4, 2, 2)), (0.5, 0, 0.3)),), weight=0.7))
    return Problem(
        nelx=4, nely=2, nelz=2, h=1.0,
        fixtures=(Fixture(Box((0, 0, 0), (0, 2, 2))),),
        load_cases=tuple(cases),
        keep=(Box((3.5, 0, 0), (4, 1, 1)),),
        settings=Settings(volume_fraction=0.5, r_min=1.6, projection=projection),
    )


@pytest.mark.parametrize("projection", [False, True])
@pytest.mark.parametrize("two_cases", [False, True])
def test_compliance_gradient_matches_finite_difference(projection, two_cases):
    m = ComplianceModel(small_problem(projection, two_cases))
    m.beta = 4.0
    rng = np.random.default_rng(5)
    x = rng.uniform(0.2, 0.8, m.grid.n_elements)
    _, dc, _ = m.evaluate(x)
    eps = 1e-6
    for e in range(m.grid.n_elements):
        xp, xm = x.copy(), x.copy()
        xp[e] += eps
        xm[e] -= eps
        fd = (m.evaluate(xp)[0] - m.evaluate(xm)[0]) / (2 * eps)
        assert dc[e] == pytest.approx(fd, rel=1e-5, abs=1e-8 * abs(dc).max())


def test_volume_gradient_matches_finite_difference():
    m = ComplianceModel(small_problem(True, False))
    m.beta = 4.0
    x = np.random.default_rng(6).uniform(0.2, 0.8, m.grid.n_elements)
    ne = m.grid.n_elements
    dv = m.physical_grad(x, np.full(ne, 1.0 / ne))
    eps = 1e-6
    for e in range(ne):
        xp, xm = x.copy(), x.copy()
        xp[e] += eps
        xm[e] -= eps
        fd = (m.physical(xp).mean() - m.physical(xm).mean()) / (2 * eps)
        assert dv[e] == pytest.approx(fd, rel=1e-6, abs=1e-12)
