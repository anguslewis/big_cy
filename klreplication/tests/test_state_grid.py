"""Validation — grid_setup (state grid + exogenous AR(1) transitions)."""
import pytest
import torch

from klrep.params import (load_param_file, IDX_K, IDX_THH, IDX_ZF, IDX_WH,
                          IDX_WF, IDX_DIS, IDX_IH, IDX_IF, IDX_OMG)
from klrep.solve.steady_state import calc_steady
from klrep.setup.smolyak_setup import build_smolyak
from klrep.setup.shock_grid import build_shock_grid
from klrep.setup.state_grid import build_state_grid


@pytest.fixture(scope="module")
def built():
    p = load_param_file(1)
    ss = calc_steady(p)
    smol = build_smolyak(p)
    sg = build_shock_grid(p)
    stg = build_state_grid(p, smol, ss, sg)
    return p, ss, smol, sg, stg


def test_shapes(built):
    p, ss, smol, sg, stg = built
    assert stg.state_grid.shape == (smol.n_states, 9)
    assert stg.wealth_share_grid.shape == (smol.n_states, 2)
    assert stg.next_zf_mat.shape == (smol.n_states, sg.n_quad)
    assert stg.next_dis_mat.shape == (smol.n_states, sg.n_quad)
    assert stg.next_omg_mat.shape == (smol.n_states, sg.n_quad)


def test_inactive_dims_pinned_to_center(built):
    p, ss, smol, sg, stg = built
    # ih, if are inactive in the benchmark -> constant at their grid mean.
    assert torch.allclose(stg.state_grid[:, IDX_IH],
                          torch.full((smol.n_states,), p["ih_grid_mean"], dtype=torch.float64))
    assert torch.allclose(stg.state_grid[:, IDX_IF],
                          torch.full((smol.n_states,), p["if_grid_mean"], dtype=torch.float64))


def test_active_dim_centers_and_ranges(built):
    p, ss, smol, sg, stg = built
    # Capital centered at k_grid_mean, within +/- k_grid_dev of it.
    k = stg.state_grid[:, IDX_K]
    assert k.min() >= ss.k_grid_mean - ss.k_grid_dev - 1e-9
    assert k.max() <= ss.k_grid_mean + ss.k_grid_dev + 1e-9
    # zf centered at 0; theta around its target.
    assert abs(float(stg.state_grid[:, IDX_ZF].mean())) < ss.k_grid_mean  # finite, near 0-centered
    th = stg.state_grid[:, IDX_THH]
    assert th.min() >= p["tht_trgt_h"] - p["theta_h_dev"] - 1e-9
    assert th.max() <= p["tht_trgt_h"] + p["theta_h_dev"] + 1e-9


def test_wealth_share_sums_to_one(built):
    p, ss, smol, sg, stg = built
    s = stg.wealth_share_grid.sum(dim=1)
    assert torch.allclose(s, torch.ones_like(s), atol=1e-12)


def test_exogenous_transition_formula(built):
    p, ss, smol, sg, stg = built
    # next_zf = (1-rho)*0 + rho*zf_current + zf_shock, checked at the no-shock node.
    nsi = sg.no_shock_idx
    zf_cur = stg.state_grid[:, IDX_ZF]
    expect_zf = p.rho_zf * zf_cur  # mean 0, no shock
    assert torch.allclose(stg.next_zf_mat[:, nsi], expect_zf, atol=1e-10)
    # next omega at no-shock node = (1-rho)*omg_mean + rho*omg_current.
    omg_cur = stg.state_grid[:, IDX_OMG]
    expect_omg = (1 - p.omg_rho) * p["omg_mean"] + p.omg_rho * omg_cur
    assert torch.allclose(stg.next_omg_mat[:, nsi], expect_omg, atol=1e-10)
    # next log-p at no-shock node = (1-rho)*disast_p_in + rho*dis_current.
    dis_cur = stg.state_grid[:, IDX_DIS]
    expect_dis = (1 - p.disast_rho) * p.disast_p_in + p.disast_rho * dis_cur
    assert torch.allclose(stg.next_dis_mat[:, nsi], expect_dis, atol=1e-10)
