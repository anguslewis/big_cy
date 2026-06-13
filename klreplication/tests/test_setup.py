"""Validation — init_setup pieces (shock grid + Smolyak wiring) vs benchmark."""
import numpy as np
import pytest
import torch

from klrep.params import load_param_file
from klrep.setup.shock_grid import build_shock_grid
from klrep.setup.smolyak_setup import build_smolyak
from klrep.grid.interp_fit import basis_condition_number


@pytest.fixture(scope="module")
def bm():
    return load_param_file(1)


def test_shock_grid_shapes(bm):
    sg = build_shock_grid(bm)
    # benchmark n_quad_shocks = [3,3,5,5] -> n_gh = 225, n_quad = 227.
    assert sg.n_gh == 3 * 3 * 5 * 5 == 225
    assert sg.n_quad == 227
    assert sg.shock_grid.shape == (227, 4)
    assert sg.quad_weight_vec.shape == (226,)  # n_gh + 1
    assert sg.no_shock_idx == 225


def test_product_weights_sum_to_one(bm):
    sg = build_shock_grid(bm)
    w = sg.quad_weight_vec[: sg.n_gh]
    assert float(w.sum()) == pytest.approx(1.0, abs=1e-12)
    assert float(sg.quad_weight_vec[sg.n_gh]) == 0.0  # no-shock node weight


def test_product_nodes_first_two_moments(bm):
    sg = build_shock_grid(bm)
    shock = sg.shock_grid[: sg.n_gh]            # (n_gh, 4)
    w = sg.quad_weight_vec[: sg.n_gh]            # (n_gh,)
    # Weighted mean ~ 0.
    mean = (w.unsqueeze(1) * shock).sum(0)
    assert torch.allclose(mean, torch.zeros(4, dtype=torch.float64), atol=1e-12)
    # Weighted second-moment matrix == U Uᵀ (the Fortran shock=U z convention).
    M = torch.einsum("k,ki,kj->ij", w, shock, shock)
    UUt = sg.chol_upper @ sg.chol_upper.T
    assert torch.allclose(M, UUt, atol=1e-10)
    # z (dim 0) is independent: its variance is exactly sig_z^2.
    assert float(M[0, 0]) == pytest.approx(bm.sig_z ** 2, rel=1e-9)


def test_disaster_and_no_shock_nodes(bm):
    sg = build_shock_grid(bm)
    # No-shock node is all zeros.
    assert torch.allclose(sg.shock_grid[sg.no_shock_idx],
                          torch.zeros(4, dtype=torch.float64))
    # Disaster node: only a z drop of -disast_shock.
    dis = sg.shock_grid[sg.n_quad - 1]
    assert float(dis[0]) == pytest.approx(-bm["disast_shock"])
    assert torch.allclose(dis[1:], torch.zeros(3, dtype=torch.float64))
    # dz_vec / dz_vec_adj.
    assert torch.allclose(sg.dz_vec, sg.shock_grid[:, 0])
    assert float(sg.dz_vec_adj[sg.n_quad - 1]) == 0.0


def test_smolyak_setup_benchmark(bm):
    s = build_smolyak(bm)
    # 7 active dims: k, tht, zf, wh, wf, p, omg (ih/if inactive).
    assert s.active_dims == [0, 1, 2, 3, 4, 5, 8]
    assert s.mus_redux == [3, 3, 3, 3, 3, 3, 3]
    assert s.mu == 3
    # Square, invertible, well-conditioned basis on the grid.
    assert s.smol_polynom.shape == (s.n_states, s.n_states)
    assert s.smol_grid.shape == (s.n_states, 7)
    cond = basis_condition_number(s.smol_polynom)
    assert np.isfinite(cond) and cond < 1e8
    # Grid points lie in [-1, 1].
    assert s.smol_grid.min() >= -1.0 - 1e-12
    assert s.smol_grid.max() <= 1.0 + 1e-12
