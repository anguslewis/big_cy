"""The equilibrium step — port of `calc_equilibrium_and_update` (mod_calc.f90:
2176-3099), batched tensor-native over ALL grid points at once.

The 10 sub-steps (see PLAN sect 3.1 / 13): (1) current production, (2) next-period
returns, (3) price of capital, (4) home bond clearing (Brent on nom_i + inner
Brent on capital share), (5) foreign bond clearing (Brent on the spread + inner
Brent on the foreign-bond share), (6) value functions + theta transition + next
aggregate capital, (7) terms of trade (damped fixed point), (8) consumption-
savings (damped Euler) which also fills the SDF M_vec_mat, (9) labor (damped wage
Phillips curve), (10) inflation (Taylor-rule inversion).

KEY (matches Fortran): steps 6/8/9 use the *incoming* portfolio shares
(share/bF_share held fixed within the iteration); the shares freshly solved in
steps 4/5 become next iteration's guess (returned). nom_i, by contrast, is solved
in steps 4/5 and used immediately in 6/8.

This is the tightly-coupled core (one module, not function-per-file) because the
sub-steps share scratch — mirroring the monolithic Fortran routine.
"""
import math
from dataclasses import dataclass

import torch

from ..config import DTYPE
from ..params import IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF, IDX_DIS, IDX_OMG
from .calc_bundle import calc_bundle  # noqa: F401  (kept for parity / future use)
from .portfolio import portfolio_return, portfolio_foc
from .period_block import compute_current_period, compute_q_new
from ..solve.brent import batched_brent, batched_bracket_expand

SQRT_EPS = math.sqrt(2.220446049250313e-16)
_CONV = 1e-8


@dataclass
class Guesses:
    c_spending: torch.Tensor   # (S,2)
    s: torch.Tensor            # (S,)
    l_aggr: torch.Tensor       # (S,2)
    q: torch.Tensor            # (S,)
    infl: torch.Tensor         # (S,2)
    nom_i: torch.Tensor        # (S,2)  [home level, spread i*-i]
    share: torch.Tensor        # (S,2)  capital-vs-bond bond share per agent
    bF_share: torch.Tensor     # (S,2)  foreign-bond share per agent


@dataclass
class StepOut:
    g: Guesses                 # updated guesses
    v: torch.Tensor            # (S,2) new value
    mc: torch.Tensor           # (S,2) new marginal value
    w_choice: torch.Tensor     # (S,2) production wage
    theta_nxt: torch.Tensor    # (S,Q) next wealth share (home)
    k_next_new: torch.Tensor   # (S,) next aggregate capital
    M_vec: torch.Tensor        # (S,Q,2) SDF


def _solve_share(savings, share_guess, bF_share, rf_home, rf_foreign, rf_raw0, rk,
                 v, mc, seign, big_w, c_cost_growth, gma, ies, q_l_ss, tot_wealth_ss,
                 dz_exp, *, low_fixed, high_fixed, is_capital):
    """Batched portfolio-share solve (over states S). For the capital share
    (is_capital=True) r1=rf_home,r2=rk and `share_guess` is the variable, holding
    bF_share fixed. For the foreign-bond share (is_capital=False) r1=rf_foreign,
    r2=rf_home and bF_share is the variable, holding share_guess (capital) fixed.
    Returns the solved share (S,). Mirrors calc_portfolio_share /
    calc_bond_portfolio_share + the leverage bounds in the excess routines."""
    dz_exp = dz_exp.unsqueeze(0)  # (1,Q)

    if is_capital:
        # bounds: temp = (-rk - bF*(rf_raw1 - rf_raw0)) / (rf_home - rk + eps)
        temp = (-rk - bF_share.unsqueeze(1) * (rf_foreign - rf_raw0)) / (rf_home - rk + SQRT_EPS)
    else:
        # foreign: temp = (-(1-share_cap)*rk - share_cap*rf_raw0) / (rf_foreign - rf_home)
        temp_r = rf_foreign - rf_home
        temp = (-(1.0 - share_guess.unsqueeze(1)) * rk
                - share_guess.unsqueeze(1) * rf_raw0) / temp_r
        temp = torch.where(temp_r == 0, torch.full_like(temp, -1000.0), temp)

    low_masked = torch.where(temp >= 0, torch.full_like(temp, -1e4), temp)
    share_low = torch.maximum(low_masked.max(dim=1).values + SQRT_EPS,
                              torch.full_like(savings, low_fixed))
    high_masked = torch.where(temp <= 0, torch.full_like(temp, 1e4), temp)
    share_high = torch.minimum(high_masked.min(dim=1).values - SQRT_EPS,
                               torch.full_like(savings, high_fixed))

    def foc(x):  # x: (S,)
        if is_capital:
            r_alpha = portfolio_return(x.unsqueeze(1), bF_share.unsqueeze(1),
                                       rf_home, rf_foreign, rk)
            r1, r2 = rf_home, rk
        else:
            r_alpha = portfolio_return(share_guess.unsqueeze(1), x.unsqueeze(1),
                                       rf_home, rf_foreign, rk)
            r1, r2 = rf_foreign, rf_home
        nps = (savings.unsqueeze(1) * r_alpha + dz_exp * q_l_ss + seign) / tot_wealth_ss
        return portfolio_foc(nps, v, mc, big_w, r1, r2, c_cost_growth, gma, ies)

    foc_low = foc(share_low)
    foc_high = foc(share_high)
    root = batched_brent(foc, share_low, share_high)
    # Constraint binding -> clamp to the boundary (signs depend on which share).
    if is_capital:
        share = torch.where(foc_low <= 0, share_low,
                            torch.where(foc_high >= 0, share_high, root))
    else:
        share = torch.where(foc_low <= 0, share_low,
                            torch.where(foc_high >= 0, share_high, root))
    return share


def equilibrium_step(const, g: Guesses, nxt_interp, nxt_mat2, big_w, outer_iter,
                     *, bond_supply_shock=0.0, bond_supply_shock_nxt=0.0):
    """One equilibrium update over all grid points. nxt_interp (S,Q,12),
    nxt_mat2 (S,Q,3)=[k_next,next_zf,next_omg], big_w (S,Q). Returns StepOut."""
    f64 = DTYPE
    S, Q = big_w.shape
    n_I = 2
    aalpha, sigma, ddelta = const.aalpha, const.sigma, const.ddelta
    zeta = const.zeta
    dz = const.dz_vec                      # (Q,)
    dz_exp = torch.exp(dz)                 # (Q,)
    sg = const.state_grid
    k_aggr = sg[:, IDX_K]
    zf_state = sg[:, IDX_ZF]
    theta = sg[:, IDX_THH]
    w_current = torch.stack([sg[:, IDX_WH], sg[:, IDX_WF]], dim=1)   # (S,2)
    omg = torch.exp(sg[:, IDX_OMG]) + const.omg_shift - const.b_lmbd * bond_supply_shock  # (S,)
    nom_i_in = g.nom_i.clone()             # incoming nom_i for the Taylor rule

    # bbeta adjustment for stationarity (per state, per agent).
    bbeta_adj = torch.stack([-const.bbeta_coeff * (theta - const.tht_h_grid_mean),
                             +const.bbeta_coeff * (theta - const.tht_h_grid_mean)], dim=1)
    bbeta = const.bbeta_vec.unsqueeze(0) + bbeta_adj  # (S,2)

    # ---- STEP 1: current production ----
    cp = compute_current_period(k_aggr, g.s, g.l_aggr, g.q, zf_state,
                                const.wealth_share_grid, aalpha=aalpha, sigma=sigma,
                                ddelta=ddelta, zeta=zeta, varsigma=const.varsigma_vec)
    kappa, pi_current, y_current = cp.kappa, cp.pi_current, cp.y_current
    w_choice = cp.w_choice
    P = cp.p_div_ph                        # (S,2)
    aggr_wealth, wealth_vec = cp.aggr_wealth, cp.wealth_vec

    # ---- STEP 2: next period ----
    v_temp = nxt_interp[:, :, 0:2]                       # (S,Q,2)
    mc_temp = nxt_interp[:, :, 2:4]
    s_nxt = nxt_interp[:, :, 4]
    q_nxt = nxt_interp[:, :, 5]
    l_nxt = nxt_interp[:, :, 6:8]
    infl_nxt = nxt_interp[:, :, 8:10]
    ch_sp_nxt = nxt_interp[:, :, 10] * dz_exp
    cf_sp_nxt = nxt_interp[:, :, 11] * dz_exp
    k_nxt = nxt_mat2[:, :, 0]
    next_k = k_nxt[:, 0]
    zf_nxt = torch.stack([torch.ones_like(s_nxt), zeta * torch.exp(nxt_mat2[:, :, 1])], dim=2)
    omg_nxt = torch.exp(nxt_mat2[:, :, 2]) + const.omg_shift - const.b_lmbd * bond_supply_shock_nxt
    vs = const.varsigma_vec
    s_col = g.s.unsqueeze(1)               # (S,1)
    hg1 = infl_nxt[:, :, 0] * ((vs[0] + (1 - vs[0]) * s_col ** (sigma - 1))
                               / (vs[0] + (1 - vs[0]) * s_nxt ** (sigma - 1))) ** (1.0 / (1.0 - sigma))
    hg2 = infl_nxt[:, :, 1] * s_nxt / s_col * ((vs[1] * s_col ** (1 - sigma) + (1 - vs[1]))
                               / (vs[1] * s_nxt ** (1 - sigma) + (1 - vs[1]))) ** (1.0 / (1.0 - sigma))
    homegood_infl_nxt = torch.stack([hg1, hg2], dim=2)   # (S,Q,2)
    ratio_nxt = s_nxt ** (1.0 / (1.0 - aalpha)) * zf_nxt[:, :, 0] / zf_nxt[:, :, 1] * l_nxt[:, :, 0] / l_nxt[:, :, 1]
    kappa_nxt_h = k_nxt * ratio_nxt / (1.0 + ratio_nxt)
    kappa_nxt = torch.stack([kappa_nxt_h, k_nxt - kappa_nxt_h], dim=2)
    pi_nxt = aalpha * torch.exp((torch.log(zf_nxt[:, :, 0]) + dz) * (1.0 - aalpha)) \
        * l_nxt[:, :, 0] ** (1.0 - aalpha) * kappa_nxt[:, :, 0] ** (aalpha - 1.0)
    rk_vec = ((1.0 - ddelta) * q_nxt + pi_nxt) / g.q.unsqueeze(1)
    rk_vec = rk_vec.clone()
    rk_vec[:, Q - 1] = rk_vec[:, Q - 1] * dz_exp[Q - 1]   # disaster-node capital scaling
    Pn1 = (vs[0] + (1 - vs[0]) * s_nxt ** (sigma - 1)) ** (1.0 / (1.0 - sigma))
    Pn2 = s_nxt ** (-1) * (vs[1] * s_nxt ** (1 - sigma) + (1 - vs[1])) ** (1.0 / (1.0 - sigma))
    P_nxt = torch.stack([Pn1, Pn2], dim=2)               # (S,Q,2)

    # seigniorage transfer
    if const.bg_yss > 0:
        seig = cf_sp_nxt * (const.bg_yss + bond_supply_shock_nxt) * omg_nxt   # (S,Q)
        seignorage = torch.stack([seig, -seig], dim=2)
    else:
        seignorage = torch.zeros((S, Q, n_I), dtype=f64, device=big_w.device)

    # ---- STEP 3: price of capital ----
    q_new = compute_q_new(g.s, next_k, k_aggr, inv_share_h=const.inv_share_h,
                          sigma=sigma, chiX=const.chiX)

    # savings per agent (independent of nom_i; clamp c_spend as in Fortran)
    msc = const.min_cons_sav
    c_spend = torch.minimum(torch.maximum(g.c_spending, torch.full_like(g.c_spending, msc)),
                            wealth_vec + w_choice * g.l_aggr - msc)
    savings = wealth_vec + w_choice * g.l_aggr - c_spend       # (S,2)

    nom_i = g.nom_i.clone()
    share_solved = g.share.clone()
    bF_solved = g.bF_share.clone()

    def rf_from_nom(nom_home, spread):
        rf0 = nom_home.unsqueeze(1) / homegood_infl_nxt[:, :, 0]
        rf1 = (nom_home + spread).unsqueeze(1) / homegood_infl_nxt[:, :, 1]
        return rf0, rf1

    # ---- STEP 4: home bond market clearing (Brent on nom_i home) ----
    if outer_iter > 10:
        def excess_home(nom_home):
            rf0, rf1 = rf_from_nom(nom_home, nom_i[:, 1])
            rf_home = rf0 / (1.0 - omg).unsqueeze(1)
            total = torch.zeros_like(nom_home)
            for iii in range(n_I):
                ccg = P_nxt[:, :, iii] / P[:, iii].unsqueeze(1)
                sh = _solve_share(savings[:, iii], g.share[:, iii], g.bF_share[:, iii],
                                  rf_home, rf1, rf0, rk_vec, v_temp[:, :, iii],
                                  mc_temp[:, :, iii], seignorage[:, :, iii], big_w, ccg,
                                  const.gma_vec[iii], const.ies_vec[iii],
                                  const.q_l_ss_vec[iii], const.tot_wealth_ss_vec[iii],
                                  dz_exp, low_fixed=const.low_guess_fixed,
                                  high_fixed=const.high_guess_fixed, is_capital=True)
                total = total + savings[:, iii] * sh
            return total
        a, b, _, _ = batched_bracket_expand(excess_home, nom_i[:, 0], step=1e-2)
        nom_home = batched_brent(excess_home, a, b)
        nom_i[:, 0] = nom_home
        # recover the equilibrium capital shares at the solved nom_i
        rf0, rf1 = rf_from_nom(nom_i[:, 0], nom_i[:, 1])
        rf_home = rf0 / (1.0 - omg).unsqueeze(1)
        for iii in range(n_I):
            ccg = P_nxt[:, :, iii] / P[:, iii].unsqueeze(1)
            share_solved[:, iii] = _solve_share(
                savings[:, iii], g.share[:, iii], g.bF_share[:, iii], rf_home, rf1, rf0,
                rk_vec, v_temp[:, :, iii], mc_temp[:, :, iii], seignorage[:, :, iii],
                big_w, ccg, const.gma_vec[iii], const.ies_vec[iii],
                const.q_l_ss_vec[iii], const.tot_wealth_ss_vec[iii], dz_exp,
                low_fixed=const.low_guess_fixed, high_fixed=const.high_guess_fixed,
                is_capital=True)

    # rf at the (possibly updated) nom_i for steps 6/8
    rf0, rf1 = rf_from_nom(nom_i[:, 0], nom_i[:, 1])
    rf_home = rf0 / (1.0 - omg).unsqueeze(1)

    # ---- STEP 5: foreign bond market clearing (Brent on the spread) ----
    if const.foreign_trading == 1 and outer_iter > 100:
        def excess_foreign(spread):
            rf0s = nom_i[:, 0].unsqueeze(1) / homegood_infl_nxt[:, :, 0]
            rf1s = (nom_i[:, 0] + spread).unsqueeze(1) / homegood_infl_nxt[:, :, 1]
            rf_home_s = rf0s / (1.0 - omg).unsqueeze(1)
            total = torch.zeros_like(spread)
            for iii in range(n_I):
                ccg = P_nxt[:, :, iii] / P[:, iii].unsqueeze(1)
                bf = _solve_share(savings[:, iii], share_solved[:, iii], g.bF_share[:, iii],
                                  rf_home_s, rf1s, rf0s, rk_vec, v_temp[:, :, iii],
                                  mc_temp[:, :, iii], seignorage[:, :, iii], big_w, ccg,
                                  const.gma_vec[iii], const.ies_vec[iii],
                                  const.q_l_ss_vec[iii], const.tot_wealth_ss_vec[iii],
                                  dz_exp, low_fixed=const.safe_low_guess_fixed,
                                  high_fixed=const.safe_high_guess_fixed, is_capital=False)
                total = total - (share_solved[:, iii] - bf) * savings[:, iii]
            return total
        a, b, _, _ = batched_bracket_expand(excess_foreign, nom_i[:, 1], step=1e-3)
        spread = batched_brent(excess_foreign, a, b)
        nom_i[:, 1] = spread
        rf0, rf1 = rf_from_nom(nom_i[:, 0], nom_i[:, 1])
        rf_home = rf0 / (1.0 - omg).unsqueeze(1)
        for iii in range(n_I):
            ccg = P_nxt[:, :, iii] / P[:, iii].unsqueeze(1)
            bF_solved[:, iii] = _solve_share(
                savings[:, iii], share_solved[:, iii], g.bF_share[:, iii], rf_home, rf1,
                rf0, rk_vec, v_temp[:, :, iii], mc_temp[:, :, iii], seignorage[:, :, iii],
                big_w, ccg, const.gma_vec[iii], const.ies_vec[iii],
                const.q_l_ss_vec[iii], const.tot_wealth_ss_vec[iii], dz_exp,
                low_fixed=const.safe_low_guess_fixed, high_fixed=const.safe_high_guess_fixed,
                is_capital=False)

    # ---- STEP 6: value functions, theta transition, next aggregate capital ----
    # NOTE: uses the INCOMING shares g.share / g.bF_share (held fixed this iteration).
    v_new = torch.empty((S, n_I), dtype=f64, device=big_w.device)
    c_vec = torch.empty((S, n_I), dtype=f64, device=big_w.device)
    savings6 = torch.empty((S, n_I), dtype=f64, device=big_w.device)
    nxt_wealth = torch.empty((S, Q, n_I), dtype=f64, device=big_w.device)
    for iii in range(n_I):
        ies = const.ies_vec[iii]; gma = const.gma_vec[iii]; chi0 = const.chi0_vec[iii]
        c_cost = P[:, iii]
        wealth = wealth_vec[:, iii]
        labor = g.l_aggr[:, iii]
        cs = torch.minimum(torch.maximum(g.c_spending[:, iii], torch.full_like(wealth, msc)),
                           wealth + w_choice[:, iii] * labor - msc)
        cons = cs / c_cost
        c_vec[:, iii] = cons
        r_alpha = (g.share[:, iii] - g.bF_share[:, iii]).unsqueeze(1) * rf_home \
            + g.bF_share[:, iii].unsqueeze(1) * rf1 + (1.0 - g.share[:, iii]).unsqueeze(1) * rk_vec
        labor_part = 1.0 + (1.0 / ies - 1.0) * chi0 * const.chi / (1.0 + const.chi) \
            * labor ** ((1.0 + const.chi) / const.chi)
        sav = wealth + w_choice[:, iii] * labor - cons * c_cost
        savings6[:, iii] = sav
        nps = (sav.unsqueeze(1) * r_alpha + dz_exp.unsqueeze(0) * const.q_l_ss_vec[iii]
               + seignorage[:, :, iii]) / const.tot_wealth_ss_vec[iii]
        EV = ((big_w * (v_temp[:, :, iii] * nps) ** (1.0 - gma)).sum(dim=1)) ** (1.0 / (1.0 - gma))
        nxt_wealth[:, :, iii] = sav.unsqueeze(1) * r_alpha + seignorage[:, :, iii]
        objf = (const.v_normalization_vec[iii] * labor_part ** (1.0 / ies) * cons ** ((ies - 1.0) / ies)
                + bbeta[:, iii] * EV ** ((ies - 1.0) / ies)) ** (ies / (ies - 1.0))
        v_new[:, iii] = objf * ((wealth + const.q_l_ss_vec[iii]) / const.tot_wealth_ss_vec[iii]) ** (-1)

    theta_nxt = nxt_wealth[:, :, 0] / nxt_wealth.sum(dim=2)
    k_next_new = (savings6 * (1.0 - g.share) / g.q.unsqueeze(1)).sum(dim=1)

    # ---- STEP 7: terms of trade (damped fixed point, weight 0.1) ----
    investment = next_k - k_aggr * (1.0 - ddelta)
    s_new = g.s.clone()
    for _ in range(2000):
        inv_h = investment * g.q / (1.0 + s_new ** (sigma - 1.0) * (1 - const.inv_share_h) / const.inv_share_h)
        inv_f = inv_h * s_new ** sigma * (1 - const.inv_share_h) / const.inv_share_h
        ch1 = g.c_spending[:, 0] / ((1 - vs[0]) / vs[0] * s_new ** (sigma - 1.0) + 1.0)
        ch2 = g.c_spending[:, 1] / ((1 - vs[1]) / vs[1] * s_new ** (sigma - 1.0) + 1.0)
        ch_sum = ch1 + ch2
        csp_minus_ch = (g.c_spending[:, 0] - ch1) + (g.c_spending[:, 1] - ch2)
        s_update = (csp_minus_ch / ch_sum * (y_current[:, 0] - inv_h) / (y_current[:, 1] - inv_f)) ** (-1)
        diff = (s_update - s_new).abs()
        s_new = s_new + 0.1 * (s_update - s_new)
        if bool((diff < _CONV).all()):
            break

    # ---- STEP 8: consumption-savings (damped Euler, weight 0.5) ----
    c_spending_new = torch.empty((S, n_I), dtype=f64, device=big_w.device)
    mc_new = torch.empty((S, n_I), dtype=f64, device=big_w.device)
    M_vec = torch.empty((S, Q, n_I), dtype=f64, device=big_w.device)
    l_temp = g.l_aggr.clone()
    for iii in range(n_I):
        ies = const.ies_vec[iii]; gma = const.gma_vec[iii]; chi0 = const.chi0_vec[iii]
        c_cost = P[:, iii]; c_cost_nxt = P_nxt[:, :, iii]
        wealth = wealth_vec[:, iii]
        labor = g.l_aggr[:, iii].clone()
        cs = torch.minimum(torch.maximum(g.c_spending[:, iii], torch.full_like(wealth, msc)),
                           wealth + w_choice[:, iii] * labor - msc)
        cons = cs / c_cost
        r_alpha = (g.share[:, iii] - g.bF_share[:, iii]).unsqueeze(1) * rf_home \
            + g.bF_share[:, iii].unsqueeze(1) * rf1 + (1.0 - g.share[:, iii]).unsqueeze(1) * rk_vec
        for _ in range(1000):
            labor_part = 1.0 + (1.0 / ies - 1.0) * chi0 * const.chi / (1.0 + const.chi) \
                * labor ** ((1.0 + const.chi) / const.chi)
            sav = wealth + w_choice[:, iii] * labor - cons * c_cost
            nps = (sav.unsqueeze(1) * r_alpha + dz_exp.unsqueeze(0) * const.q_l_ss_vec[iii]
                   + seignorage[:, :, iii]) / const.tot_wealth_ss_vec[iii]
            util_c_deriv = (labor_part / cons) ** (1.0 / ies)
            mc_cur = util_c_deriv * ((const.q_l_ss_vec[iii] + wealth) / const.tot_wealth_ss_vec[iii]) ** (1.0 / ies)
            EV = ((big_w * (v_temp[:, :, iii] * nps) ** (1.0 - gma)).sum(dim=1)) ** (1.0 / (1.0 - gma))
            M = bbeta[:, iii].unsqueeze(1) * mc_temp[:, :, iii] / util_c_deriv.unsqueeze(1) \
                * (v_temp[:, :, iii] / EV.unsqueeze(1)) ** (1.0 / ies - gma) * nps ** (-gma)
            temp = (M * r_alpha * c_cost.unsqueeze(1) / c_cost_nxt * big_w).sum(dim=1)
            cons_update = cons / temp ** ies
            cons_update = torch.minimum(torch.maximum(cons_update, torch.full_like(cons, msc)),
                                        (wealth + w_choice[:, iii] * labor - msc) / c_cost)
            diff = (cons - cons_update).abs()
            cons = cons + 0.5 * (cons_update - cons)
            if const.phi_w == 0:
                labor_new = (w_choice[:, iii] / c_cost
                             * (1.0 + (1.0 / ies - 1.0) * chi0 * const.chi / (1.0 + const.chi)
                                * labor ** ((1.0 + const.chi) / const.chi))
                             / (1.0 / ies * cons * chi0)) ** const.chi
                labor = labor + 0.2 * (labor_new - labor)
            if bool((diff < _CONV).all()):
                break
        c_spending_new[:, iii] = cons * c_cost
        mc_new[:, iii] = mc_cur
        M_vec[:, :, iii] = M
        l_temp[:, iii] = labor

    # ---- STEP 9: labor / wage Phillips curve (damped, weight 0.005) ----
    if const.phi_w > 0:
        l_aggr_new = g.l_aggr.clone()
        zf2 = zeta * torch.exp(zf_state)
        zf_vec = torch.stack([torch.ones_like(zf2), zf2], dim=1)
        s_vec = torch.stack([torch.ones_like(g.s), g.s], dim=1)
        for iii in range(n_I):
            ies = const.ies_vec[iii]; chi0 = const.chi0_vec[iii]
            w_temp = w_choice[:, iii].clone()
            w_next_choice_i = ((1 - aalpha) * (nxt_interp[:, :, 6 + iii] * 0 + 1))  # placeholder fixed below
            # next-period production wage for agent iii
            ynext_i = torch.exp((torch.log(zf_nxt[:, :, iii]) + dz) * (1.0 - aalpha)) \
                * l_nxt[:, :, iii] * (kappa_nxt[:, :, iii] / l_nxt[:, :, iii]) ** aalpha
            if iii == 0:
                w_next_choice_i = (1 - aalpha) * ynext_i / l_nxt[:, :, 0]
            else:
                w_next_choice_i = s_nxt ** (-1) * (1 - aalpha) * ynext_i / l_nxt[:, :, 1]
            for it in range(100):
                w_cn = w_temp.unsqueeze(1).expand(S, Q).clone()
                w_cn[:, Q - 1] = w_temp * dz_exp[Q - 1]
                infl_term = (w_temp / w_current[:, iii] * g.infl[:, iii] / P[:, iii] - 1.0) \
                    * w_temp / w_current[:, iii] * g.infl[:, iii] / P[:, iii]
                cont = (M_vec[:, :, iii] * big_w * homegood_infl_nxt[:, :, iii] / infl_nxt[:, :, iii]
                        * (w_next_choice_i / w_cn * homegood_infl_nxt[:, :, iii] - 1.0)
                        * (w_next_choice_i / w_cn) ** 2 * homegood_infl_nxt[:, :, iii]
                        * l_nxt[:, :, iii] / l_aggr_new[:, iii].unsqueeze(1)).sum(dim=1)
                disutil = -(1.0 / ies) * c_vec[:, iii] * chi0 * (l_aggr_new[:, iii] ** (1.0 / const.chi)) \
                    / (1.0 + (1.0 / ies - 1.0) * chi0 * const.chi / (1.0 + const.chi)
                       * l_aggr_new[:, iii] ** ((1.0 + const.chi) / const.chi))
                w_temp_new = const.vareps_w / (1.0 - const.vareps_w) * P[:, iii] * disutil \
                    + w_temp * const.phi_w / (1.0 - const.vareps_w) * (infl_term - cont)
                w_diff = (w_temp_new - w_temp).abs()
                w_temp = w_temp + 0.005 * (w_temp_new - w_temp)
                l_aggr_new[:, iii] = kappa[:, iii] * ((zf_vec[:, iii] ** (1.0 - aalpha) * (1.0 - aalpha)
                                                       / w_temp / s_vec[:, iii]) ** (1.0 / aalpha))
                if it >= 10 and bool((w_diff < _CONV).all()):
                    break
    else:
        l_aggr_new = l_temp
    w_choice_new = w_choice.clone()

    # ---- STEP 10: inflation (Taylor-rule inversion) ----
    ih_last = sg[:, 6]   # IDX_IH
    if_last = sg[:, 7]   # IDX_IF
    infl_new0 = ((nom_i_in[:, 0] / (1.0 + ih_last) ** const.rho_i) ** (1.0 / (1.0 - const.rho_i))
                 * math.exp(-const.tayl_ic_h) / (y_current[:, 0] ** const.phi_yh)) ** (1.0 / const.phi_h)
    infl_new1 = (((nom_i_in[:, 0] + nom_i_in[:, 1]) / (1.0 + if_last) ** const.rho_i) ** (1.0 / (1.0 - const.rho_i))
                 * math.exp(-const.tayl_ic_f) / (y_current[:, 1] ** const.phi_yf)) ** (1.0 / const.phi_f)
    infl_new = torch.stack([infl_new0, infl_new1], dim=1)

    # ---- bond-pricing fallback (outer_iter <= 10) + share/nom_i finalization ----
    if outer_iter <= 10:
        share_out = torch.zeros((S, n_I), dtype=f64, device=big_w.device)
        nh = torch.stack([
            (big_w * M_vec[:, :, 0] / (1.0 - omg).unsqueeze(1) * (1.0 / homegood_infl_nxt[:, :, 0])
             * P[:, 0].unsqueeze(1) / P_nxt[:, :, 0]).sum(dim=1),
            (big_w * M_vec[:, :, 1] / (1.0 - omg).unsqueeze(1) * (1.0 / homegood_infl_nxt[:, :, 0])
             * P[:, 1].unsqueeze(1) / P_nxt[:, :, 1]).sum(dim=1)], dim=1)
        nom_i[:, 0] = 1.0 / nh.max(dim=1).values
    else:
        share_out = share_solved

    if const.foreign_trading == 1 and outer_iter > 100:
        bF_out = bF_solved
    else:
        bF_out = torch.zeros((S, n_I), dtype=f64, device=big_w.device)
        nf = torch.stack([
            (big_w * M_vec[:, :, 0] * (1.0 / homegood_infl_nxt[:, :, 1])
             * P[:, 0].unsqueeze(1) / P_nxt[:, :, 0]).sum(dim=1),
            (big_w * M_vec[:, :, 1] * (1.0 / homegood_infl_nxt[:, :, 1])
             * P[:, 1].unsqueeze(1) / P_nxt[:, :, 1]).sum(dim=1)], dim=1)
        nom_i[:, 1] = 1.0 / nf.max(dim=1).values - nom_i[:, 0]

    g_new = Guesses(c_spending=c_spending_new, s=s_new, l_aggr=l_aggr_new, q=q_new,
                    infl=infl_new, nom_i=nom_i, share=share_out, bF_share=bF_out)
    return StepOut(g=g_new, v=v_new, mc=mc_new, w_choice=w_choice_new,
                   theta_nxt=theta_nxt, k_next_new=k_next_new, M_vec=M_vec)
