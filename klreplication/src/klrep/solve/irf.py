"""Generalized (MIT-shock) impulse responses — port of the IRF preparation in
`calc_sol` (mod_calc.f90:989-1213, incl. `calc_unexpected_transition` 1940-2174)
and the IRF simulation in `mod_results.f90` (:489-792).

A generalized IRF: from each of n_sample ergodic starting states, sit at the
stochastic steady state for two periods, apply ONE shock at t=3 (jump into the
"shock-transition" state), then propagate with the no-shock transition; difference
the shocked path against the matched no-shock path and average. The shock-transition
state (standardized next-state when a specific shock of a given size hits) is a
damped fixed point because the realized wealth share theta' depends on next-period
values interpolated at that very state (`calc_unexpected_transition`).

Shocks (sidx, irf_shock_sizes from init_setup): z (1, -5*sig_z), zf (2, -2*sig_zf),
p (3, 2*sig_dis), omg (4, 2*sig_omg), disaster (period-3 = disaster quad node).
"""
from dataclasses import dataclass

import torch

from ..config import DTYPE
from ..params import (IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF, IDX_DIS, IDX_IH,
                      IDX_IF, IDX_OMG)
from ..grid.smolyak_polynomial import smolyak_polynomial
from ..grid.interp_fit import interp_factor, interp_solve
from ..model.period_block import compute_current_period

SIDX_Z, SIDX_ZF, SIDX_P, SIDX_OMG = 0, 1, 2, 3   # 0-based shock-vector indices
SQRT_EPS = 1.49e-8


@dataclass
class _Current:
    """Iteration-invariant current-period block (depends only on the converged g)."""
    w_choice: torch.Tensor    # (S,2)
    pi_home: torch.Tensor     # (S,)
    wealth_vec: torch.Tensor  # (S,2)
    P: torch.Tensor           # (S,2)
    omg: torch.Tensor         # (S,)


def _interp_input(g, v_mat, mc_mat):
    return torch.stack([
        v_mat[:, 0], v_mat[:, 1], mc_mat[:, 0], mc_mat[:, 1],
        g.s, g.q, g.l_aggr[:, 0], g.l_aggr[:, 1], g.infl[:, 0], g.infl[:, 1],
        g.c_spending[:, 0], g.c_spending[:, 1]], dim=1)


def _current(const, g) -> _Current:
    sg = const.state_grid
    cp = compute_current_period(sg[:, IDX_K], g.s, g.l_aggr, g.q, sg[:, IDX_ZF],
                                const.wealth_share_grid, aalpha=const.aalpha,
                                sigma=const.sigma, ddelta=const.ddelta,
                                zeta=const.zeta, varsigma=const.varsigma_vec)
    # pi_home + aggregate wealth per calc_unexpected_transition (2032,2044).
    pi_home = const.aalpha * (g.l_aggr[:, 0] / cp.kappa[:, 0]) ** (1.0 - const.aalpha)
    aggr_wealth = sg[:, IDX_K] * ((1.0 - const.ddelta) * g.q + pi_home)
    wealth_vec = aggr_wealth.unsqueeze(1) * const.wealth_share_grid
    omg = torch.exp(sg[:, IDX_OMG]) + const.omg_shift
    return _Current(w_choice=cp.w_choice, pi_home=pi_home, wealth_vec=wealth_vec,
                    P=cp.p_div_ph, omg=omg)


def _unexpected_theta(const, g, cur, nxt, omg_nxt_exo, dz_temp):
    """Realized next wealth share theta' (S,) for one shocked node — port of
    calc_unexpected_transition. `nxt` (S,12) = next-period vars interpolated at the
    transition state; `omg_nxt_exo` (S,) the next (log) omega state; dz_temp the z
    innovation."""
    sigma, ddelta = const.sigma, const.ddelta
    s_nxt, q_nxt = nxt[:, 4], nxt[:, 5]
    infl_nxt = nxt[:, 8:10]
    dz_exp = torch.exp(torch.tensor(dz_temp, dtype=DTYPE, device=s_nxt.device))
    cf_sp_nxt = nxt[:, 11] * dz_exp
    vs = const.varsigma_vec
    hg1 = infl_nxt[:, 0] * ((vs[0] + (1 - vs[0]) * g.s ** (sigma - 1))
                            / (vs[0] + (1 - vs[0]) * s_nxt ** (sigma - 1))) ** (1.0 / (1.0 - sigma))
    hg2 = infl_nxt[:, 1] * s_nxt / g.s * ((vs[1] * g.s ** (1 - sigma) + (1 - vs[1]))
                            / (vs[1] * s_nxt ** (1 - sigma) + (1 - vs[1]))) ** (1.0 / (1.0 - sigma))
    omg_nxt = torch.exp(omg_nxt_exo) + const.omg_shift
    rk = (1.0 - ddelta) * q_nxt / (g.q - cur.pi_home)          # realized return (2112)
    rf0 = g.nom_i[:, 0] / hg1
    rf1 = (g.nom_i[:, 0] + g.nom_i[:, 1]) / hg2
    rf_home = rf0 / (1.0 - cur.omg)
    if const.bg_yss > 0:
        seig = cf_sp_nxt * const.bg_yss * omg_nxt
        seignorage = torch.stack([seig, -seig], dim=1)
    else:
        seignorage = torch.zeros((s_nxt.shape[0], 2), dtype=DTYPE, device=s_nxt.device)
    msc = const.min_cons_sav
    nxt_wealth = torch.empty((s_nxt.shape[0], 2), dtype=DTYPE, device=s_nxt.device)
    for iii in range(2):
        wl = cur.wealth_vec[:, iii]
        cs = torch.minimum(torch.maximum(g.c_spending[:, iii], torch.full_like(wl, msc)),
                           wl + cur.w_choice[:, iii] * g.l_aggr[:, iii] - msc)
        savings = wl + cur.w_choice[:, iii] * g.l_aggr[:, iii] - cs
        r_alpha = ((g.share[:, iii] - g.bF_share[:, iii]) * rf_home
                   + g.bF_share[:, iii] * rf1 + (1.0 - g.share[:, iii]) * rk)
        nxt_wealth[:, iii] = savings * r_alpha + seignorage[:, iii]
    return nxt_wealth[:, 0] / nxt_wealth.sum(dim=1)


def build_shock_transition(const, g, v_mat, mc_mat, next_state, k_next_new,
                           sidx, shock_size, *, max_iter=500, conv=1e-8, damp=0.2):
    """Damped fixed point for the standardized shock-transition state (S, d) when
    shock `sidx` (0..3 = z/zf/p/omg) of `shock_size` hits — port of mod_calc.f90:
    1015-1191."""
    S, Q, d = next_state.shape
    dev = next_state.device
    no_shock = const.no_shock_idx
    shock = torch.zeros(4, dtype=DTYPE, device=dev)
    shock[sidx] = shock_size
    cur = _current(const, g)
    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, _interp_input(g, v_mat, mc_mat))         # (S,12)

    zf_no = const.next_zf_mat[:, no_shock]
    omg_no = const.next_omg_mat[:, no_shock]
    dis_no = const.next_dis_mat[:, no_shock]
    dz_adj_exp = torch.exp(const.dz_vec_adj[no_shock])
    m, dv = const.grid_means, const.grid_devs
    sigma, vs = const.sigma, const.varsigma_vec
    s_col = g.s
    dz_z_exp = torch.exp(shock[SIDX_Z])

    state = next_state[:, no_shock, :].clone()
    for it in range(1, max_iter + 1):
        nxt = smolyak_polynomial(state, const.elem, const.mu) @ coeffs   # (S,12)
        omg_nxt_exo = omg_no + shock[SIDX_OMG]
        theta_nxt = _unexpected_theta(const, g, cur, nxt, omg_nxt_exo, float(shock[SIDX_Z]))

        cols = {}
        cols[IDX_K] = (k_next_new / dz_z_exp / dz_adj_exp - m[IDX_K]) / dv[IDX_K]
        cols[IDX_THH] = (theta_nxt - m[IDX_THH]) / dv[IDX_THH]
        cols[IDX_ZF] = (zf_no + shock[SIDX_ZF] - m[IDX_ZF]) / dv[IDX_ZF]
        A1 = vs[0] + (1 - vs[0]) * s_col ** (sigma - 1.0)
        cols[IDX_WH] = (cur.w_choice[:, 0] / dz_z_exp / A1 ** (1.0 / (sigma - 1.0))
                        - m[IDX_WH]) / dv[IDX_WH]
        A2 = vs[1] * s_col ** (1.0 - sigma) + (1 - vs[1])
        cols[IDX_WF] = (cur.w_choice[:, 1] / dz_z_exp * s_col / A2 ** (1.0 / (sigma - 1.0))
                        - m[IDX_WF]) / dv[IDX_WF]
        cols[IDX_DIS] = (torch.zeros(S, dtype=DTYPE, device=dev) if dv[IDX_DIS] < SQRT_EPS
                         else (dis_no + shock[SIDX_P] - m[IDX_DIS]) / dv[IDX_DIS])
        cols[IDX_IH] = ((g.nom_i[:, 0] - 1.0) - m[IDX_IH]) / dv[IDX_IH]
        cols[IDX_IF] = ((g.nom_i[:, 0] + g.nom_i[:, 1] - 1.0) - m[IDX_IF]) / dv[IDX_IF]
        cols[IDX_OMG] = (torch.zeros(S, dtype=DTYPE, device=dev) if dv[IDX_OMG] < SQRT_EPS
                         else (omg_no + shock[SIDX_OMG] - m[IDX_OMG]) / dv[IDX_OMG])

        new = torch.stack([cols[dd] for dd in const.active_dims], dim=1).clamp(-1.0, 1.0)
        diff = float((new - state).abs().max()) if it >= 10 else 1.0
        state = state + damp * (new - state)
        if diff < conv:
            break
    return state
