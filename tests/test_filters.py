"""Filter and projection verification (Bourdin 2001; Wang, Lazarov & Sigmund 2011)."""

import numpy as np
import pytest

from reversed_core.filters import DensityFilter, heaviside, heaviside_grad
from reversed_core.mesh import VoxelGrid


@pytest.fixture(scope="module")
def filt():
    return DensityFilter(VoxelGrid(7, 5, 4), r_min=2.3)


def test_preserves_constant_field(filt):
    x = np.full(filt.grid.n_elements, 0.37)
    assert np.allclose(filt.apply(x), 0.37)


def test_backprop_is_adjoint_of_apply(filt):
    rng = np.random.default_rng(2)
    x, y = rng.normal(size=(2, filt.grid.n_elements))
    assert filt.apply(x) @ y == pytest.approx(x @ filt.backprop(y), rel=1e-12)


def test_output_within_input_bounds(filt):
    x = np.random.default_rng(3).uniform(size=filt.grid.n_elements)
    xt = filt.apply(x)
    assert xt.min() >= x.min() - 1e-12 and xt.max() <= x.max() + 1e-12


def test_radius_one_is_identity():
    f = DensityFilter(VoxelGrid(4, 3, 2), r_min=1.0)
    x = np.random.default_rng(4).uniform(size=24)
    assert np.allclose(f.apply(x), x)


def test_weights_match_formula():
    g = VoxelGrid(5, 5, 5)
    f = DensityFilter(g, r_min=1.8)
    centre = (2 * 5 + 2) * 5 + 2
    row = f.H.getrow(centre).toarray().ravel()
    # face neighbour (distance 1) -> 0.8, edge (sqrt 2) -> 0.3858, corner (sqrt 3) -> 0.0679
    assert row[centre] == pytest.approx(1.8)
    assert row[centre + 1] == pytest.approx(0.8)
    assert row[centre + 5 + 1] == pytest.approx(1.8 - np.sqrt(2))
    assert row[centre + 25 + 5 + 1] == pytest.approx(1.8 - np.sqrt(3))
    assert np.count_nonzero(row) == 27


@pytest.mark.parametrize("beta", [1.0, 8.0, 64.0])
def test_heaviside_endpoints_and_gradient(beta):
    assert heaviside(np.array([0.0, 1.0]), beta) == pytest.approx([0.0, 1.0])
    x = np.linspace(0.05, 0.95, 19)
    eps = 1e-6
    fd = (heaviside(x + eps, beta) - heaviside(x - eps, beta)) / (2 * eps)
    assert np.allclose(heaviside_grad(x, beta), fd, rtol=1e-5, atol=1e-8)
