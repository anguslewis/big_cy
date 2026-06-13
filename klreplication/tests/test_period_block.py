"""Validation — current-period production + price of capital reproduce the
steady state (the static algebra of the equilibrium step, checked at SS)."""
import torch
import pytest

from klrep.params import load_param_file
from klrep.solve.steady_state import calc_steady
from klrep.model.period_block import compute_current_period, compute_q_new


@pytest.fixture(scope="module")
def setup():
    p = load_param_file(1)
    ss = calc_steady(p)
    return p, ss


def test_current_period_reproduces_ss(setup):
    p, ss = setup
    f64 = torch.float64
    one = lambda x: torch.tensor([x], dtype=f64)
    k = one(ss.k_ss)
    s = one(ss.s_ss)
    q = one(ss.q_ss)
    l = torch.tensor([[1.0, 1.0]], dtype=f64)
    zf_state = one(0.0)
    wshare = torch.tensor([[p["tht_trgt_h"], 1 - p["tht_trgt_h"]]], dtype=f64)
    varsigma = torch.tensor([p["varsigma_h"], p["varsigma_f"]], dtype=f64)

    cp = compute_current_period(
        k, s, l, q, zf_state, wshare,
        aalpha=p["aalpha"], sigma=p["sigma"], ddelta=p["delta"],
        zeta=p.zeta, varsigma=varsigma,
    )
    # Capital allocation, wages reproduce SS.
    assert torch.allclose(cp.kappa[0], torch.tensor(ss.kappa_ss, dtype=f64), rtol=1e-6)
    assert torch.allclose(cp.w_choice[0], torch.tensor(ss.w_ss, dtype=f64), rtol=1e-6)
    # Profit per unit capital equals pi_ss in both countries.
    assert torch.allclose(cp.pi_current[0],
                          torch.tensor([ss.pi_ss, ss.pi_ss], dtype=f64), rtol=1e-6)
    # Aggregate wealth = (pi_ss + (1-delta) q_ss) * k_ss.
    expect_w = (ss.pi_ss + (1 - p["delta"]) * ss.q_ss) * ss.k_ss
    assert float(cp.aggr_wealth[0]) == pytest.approx(expect_w, rel=1e-6)


def test_q_new_reproduces_q_ss(setup):
    p, ss = setup
    f64 = torch.float64
    one = lambda x: torch.tensor([x], dtype=f64)
    q_new = compute_q_new(one(ss.s_ss), one(ss.k_ss), one(ss.k_ss),
                          inv_share_h=p["inv_share_h"], sigma=p["sigma"],
                          chiX=p["chiX"])
    assert float(q_new[0]) == pytest.approx(ss.q_ss, rel=1e-9)


def test_capital_return_is_one_plus_rk_ss(setup):
    # rk = ((1-delta) q_nxt + pi_nxt) / q_current; at SS q_nxt=q_current=q_ss,
    # pi_nxt = pi_ss  ->  rk = 1 + rk_ss.
    p, ss = setup
    rk = ((1 - p["delta"]) * ss.q_ss + ss.pi_ss) / ss.q_ss
    assert rk == pytest.approx(1.0 + ss.rk_ss, rel=1e-9)
