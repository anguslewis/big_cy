"""Tests for the bond-ladder port (calc_bond_prices): the SDF pricing block matches
equilibrium_step's M_vec at converged-ish policies; the ladder produces finite
positive prices; and the bond-dependent moments 8/9/10 compute on synthetic series."""
import math

import torch

from klrep.params import load_param_file
from klrep.solve.solve_model import init_solver, solve_step, _interp_input
from klrep.grid.interp_fit import interp_factor, interp_solve
from klrep.grid.smolyak_polynomial import smolyak_polynomial
from klrep.model.equilibrium_step import equilibrium_step
from klrep.model.bond_ladder import (compute_bond_columns, compute_pricing_block,
                                     STORED_MATS, N_BOND)
from klrep.post.moments import (compute_table2_moments_per_sim, _ols,
                                compute_extended_tables)


def _coarse_solved(n_iters=120):
    p = load_param_file(1)
    p.mus_dims = [1, 1, 1, 1, 1, 1, 0, 0, 1]
    st = init_solver(p)
    diff = 1.0
    for it in range(1, n_iters + 1):
        st, diff = solve_step(st, it)
    return p, st, diff


def test_pricing_block_matches_equilibrium_step():
    """compute_pricing_block reproduces equilibrium_step's SDF (same formula, one
    evaluation at the current policies — agreement scales with convergence)."""
    p, st, diff = _coarse_solved(120)
    const = st.const
    S, Q = const.state_grid.shape[0], const.n_quad
    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, _interp_input(st))
    Bn = smolyak_polynomial(st.next_state.reshape(S * Q, -1), const.elem, const.mu)
    nxt_interp = (Bn @ coeffs).reshape(S, Q, 12)
    nxt_mat2 = torch.stack([st.k_next_mat, const.next_zf_mat, const.next_omg_mat], dim=2)
    out = equilibrium_step(const, st.g, nxt_interp, nxt_mat2, st.big_w, 9999)
    pb = compute_pricing_block(const, st.g, st.v_mat, st.mc_mat, st.next_state)
    assert pb.M_vec.shape == (S, Q, 2)
    rel = ((pb.M_vec - out.M_vec).abs() / out.M_vec.abs().clamp_min(1e-12)).max()
    # same formula; residual is the single extra damped-Euler step in equilibrium_step.
    assert float(rel) < 0.05, f"M_vec rel diff {float(rel):.2e} at solve diff {diff:.2e}"


def test_bond_columns_finite_positive():
    p, st, _ = _coarse_solved(120)
    bg = compute_bond_columns(st.const, st.g, st.v_mat, st.mc_mat, st.next_state)
    S = st.const.state_grid.shape[0]
    assert set(bg) == {"q1_h", "q1_f", "q1_hw", "q19_h", "q19_f", "q20_h", "q20_f"}
    for k, v in bg.items():
        assert v.shape == (S,), k
        assert torch.isfinite(v).all(), k
        assert (v > 0).all(), k          # bond prices are positive
    assert STORED_MATS == [1, 2, 3, 4, 19, 20] and N_BOND == 20


def test_ols_recovers_known_coeffs():
    # y = 2 + 3*x1 - 1*x2 exactly -> OLS recovers (2,3,-1), residual ~0.
    torch.manual_seed(0)
    n_sims, L = 3, 50
    x1 = torch.randn(n_sims, L, dtype=torch.float64)
    x2 = torch.randn(n_sims, L, dtype=torch.float64)
    y = 2.0 + 3.0 * x1 - 1.0 * x2
    beta, resid = _ols([x1, x2], y)
    assert torch.allclose(beta, torch.tensor([2.0, 3.0, -1.0], dtype=torch.float64).expand(n_sims, 3), atol=1e-8)
    assert float(resid.abs().max()) < 1e-8


def test_full_15_moments_synthetic():
    # Provide all series (12-moment + bond-dependent) -> all 15 moments finite.
    n_sims, T = 2, 40
    g = torch.Generator().manual_seed(1)
    pos = lambda: 1.0 + 0.05 * torch.rand(n_sims, T, dtype=torch.float64, generator=g)
    s = {k: pos() for k in ["yh", "yf", "s", "ch", "cf", "inv", "lh", "lf",
                            "infl_h", "infl_f", "qx", "P_Phh", "P_Phf", "q", "pi"]}
    s["rfh"] = 0.005 + 0.001 * torch.randn(n_sims, T, dtype=torch.float64, generator=g)
    s["nfa_rel"] = -0.8 + 0.1 * torch.randn(n_sims, T, dtype=torch.float64, generator=g)
    s["nfa_rel_growth"] = 0.01 * torch.randn(n_sims, T, dtype=torch.float64, generator=g)
    for k in ["rA", "exc_retA", "exc_rf", "uip_pvt", "fama_yield_pvt", "y_growth",
              "div_price_smoothed_1"]:
        s[k] = torch.randn(n_sims, T, dtype=torch.float64, generator=g)
    moms = compute_table2_moments_per_sim(s, bg_yss=0.127)
    assert set(moms.keys()) == set(range(1, 16))
    for k, v in moms.items():
        assert torch.isfinite(v).all(), f"moment {k} non-finite"


def test_extended_tables_synthetic():
    # Series needed by Tables 3/4/5/10 -> all 17 moments finite.
    n_sims, T = 3, 40
    g = torch.Generator().manual_seed(2)
    pos = lambda: 1.0 + 0.05 * torch.rand(n_sims, T, dtype=torch.float64, generator=g)
    s = {k: pos() for k in ["yh", "yf", "ch", "cf", "qx", "h_sav", "h_ksav",
                            "h_kap", "h_bh_sav"]}
    for k in ["exc_retA", "exc_rf", "rfh", "uip_pvt", "y_growth", "nfa_rel_growth",
              "E_change"]:
        s[k] = torch.randn(n_sims, T, dtype=torch.float64, generator=g)
    ext = compute_extended_tables(s)
    assert set(ext) == {"t3_1", "t3_2", "t3_3", "t3_4", "t3_5", "t3_6",
                        "t4_1", "t4_2", "t4_3", "t4_4", "t4_5", "t4_6",
                        "t5_1", "t5_2", "t10_k", "t10_bH", "t10_bF"}
    assert all(math.isfinite(v) for v in ext.values())
