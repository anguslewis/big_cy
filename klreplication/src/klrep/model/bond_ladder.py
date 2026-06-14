"""Recursive nominal-bond pricing — port of `calc_bond_prices` (mod_calc.f90:
3101-3501) and its maturity recursion (the bbb=1..n_bond loop, 716-800).

Post-solve: given the converged policies, price zero-coupon nominal bonds of every
maturity 1..n_bond in three "currencies" (h = home private, f = foreign, hw = home
omega-adjusted / convenience). Recursion: an m-maturity bond today is a claim that
next period is worth the (m-1)-maturity bond, so

    q_bond[:,:,1] = 1                      (maturity-0 payoff)
    for m = 1..n_bond:
        coeffs       = fit Smolyak to q_bond[:,:,m]        (the 3-currency prices)
        nxt_bond     = basis(next_state) @ coeffs          (price next period, per node)
        q_bond[:,:,m+1] = price_bonds(nxt_bond)            (one calc_bond_prices step)

`price_bonds` discounts by the SDF M of the *marginal investor* (the agent who
values the bond most — `maxval` over agents), converting the nominal payoff to real
home-good units via home/foreign-good inflation and the basket-cost ratio
P/P_nxt; the hw bond additionally divides by (1-omg).

The SDF M_vec here reproduces `equilibrium_step`'s STEP-6/8 kernel evaluated ONCE
at the converged policies (no Euler iteration) — validated against it in tests.

Only the maturities the Table-2 bond moments need are returned (1-period yields,
and the 19/20-period prices for the levered-equity 20-period bond return). All
columns are bond prices (ratios) — NOT de-trended downstream.
"""
from dataclasses import dataclass

import torch

from ..config import DTYPE
from ..params import IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF, IDX_OMG
from ..grid.smolyak_polynomial import smolyak_polynomial
from ..grid.interp_fit import interp_factor, interp_solve
from ..solve.model_const import big_weight
from .period_block import compute_current_period

N_BOND = 20
# Maturities KL stores (bbb<=4 or bbb>n_bond-2): [1,2,3,4,19,20]. We expose the
# ones the Table-2 bond moments use.
STORED_MATS = [1, 2, 3, 4, 19, 20]


@dataclass
class PricingBlock:
    M_vec: torch.Tensor       # (S,Q,2) SDF per agent
    hg_infl: torch.Tensor     # (S,Q,2) home/foreign good gross inflation next period
    P: torch.Tensor           # (S,2) basket cost in own good, this period
    P_nxt: torch.Tensor       # (S,Q,2) basket cost next period
    omg: torch.Tensor         # (S,) convenience wedge
    next_state: torch.Tensor  # (S,Q,d) standardized next states (active dims)


def _interp_input(const, g, v_mat, mc_mat):
    return torch.stack([
        v_mat[:, 0], v_mat[:, 1], mc_mat[:, 0], mc_mat[:, 1],
        g.s, g.q, g.l_aggr[:, 0], g.l_aggr[:, 1], g.infl[:, 0], g.infl[:, 1],
        g.c_spending[:, 0], g.c_spending[:, 1]], dim=1)


def compute_pricing_block(const, g, v_mat, mc_mat, next_state) -> PricingBlock:
    """Compute the SDF + price-conversion tensors at the converged policies."""
    f64 = DTYPE
    S, Q, d = next_state.shape
    n_I = 2
    aalpha, sigma, ddelta, zeta = const.aalpha, const.sigma, const.ddelta, const.zeta
    dz = const.dz_vec
    dz_exp = torch.exp(dz)
    sg = const.state_grid
    k_aggr = sg[:, IDX_K]
    zf_state = sg[:, IDX_ZF]
    theta = sg[:, IDX_THH]
    omg = torch.exp(sg[:, IDX_OMG]) + const.omg_shift                       # (S,)

    bbeta_adj = torch.stack([-const.bbeta_coeff * (theta - const.tht_h_grid_mean),
                             +const.bbeta_coeff * (theta - const.tht_h_grid_mean)], dim=1)
    bbeta = const.bbeta_vec.unsqueeze(0) + bbeta_adj                        # (S,2)

    # STEP 1 current production
    cp = compute_current_period(k_aggr, g.s, g.l_aggr, g.q, zf_state,
                                const.wealth_share_grid, aalpha=aalpha, sigma=sigma,
                                ddelta=ddelta, zeta=zeta, varsigma=const.varsigma_vec)
    P, w_choice, wealth_vec = cp.p_div_ph, cp.w_choice, cp.wealth_vec

    # interpolate next-period objects at next_state
    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, _interp_input(const, g, v_mat, mc_mat))       # (S,12)
    B = smolyak_polynomial(next_state.reshape(S * Q, d), const.elem, const.mu)
    nxt = (B @ coeffs).reshape(S, Q, 12)
    v_temp = nxt[:, :, 0:2]
    mc_temp = nxt[:, :, 2:4]
    s_nxt = nxt[:, :, 4]
    q_nxt = nxt[:, :, 5]
    l_nxt = nxt[:, :, 6:8]
    infl_nxt = nxt[:, :, 8:10]
    cf_sp_nxt = nxt[:, :, 11] * dz_exp

    # next non-interp states: recompute k_next_mat from converged g
    msc = const.min_cons_sav
    c_spend = torch.minimum(torch.maximum(g.c_spending, torch.full_like(g.c_spending, msc)),
                            wealth_vec + w_choice * g.l_aggr - msc)
    savings = wealth_vec + w_choice * g.l_aggr - c_spend                    # (S,2)
    k_next_new = (savings * (1.0 - g.share) / g.q.unsqueeze(1)).sum(dim=1)  # (S,)
    dz_adj_exp = torch.exp(const.dz_vec_adj)
    k_nxt = k_next_new.unsqueeze(1) / dz_adj_exp.unsqueeze(0) * dz_exp.unsqueeze(0)  # (S,Q)

    zf_nxt0 = torch.ones_like(s_nxt)
    zf_nxt1 = zeta * torch.exp(const.next_zf_mat)
    vs = const.varsigma_vec
    s_col = g.s.unsqueeze(1)
    hg1 = infl_nxt[:, :, 0] * ((vs[0] + (1 - vs[0]) * s_col ** (sigma - 1))
                               / (vs[0] + (1 - vs[0]) * s_nxt ** (sigma - 1))) ** (1.0 / (1.0 - sigma))
    hg2 = infl_nxt[:, :, 1] * s_nxt / s_col * ((vs[1] * s_col ** (1 - sigma) + (1 - vs[1]))
                               / (vs[1] * s_nxt ** (1 - sigma) + (1 - vs[1]))) ** (1.0 / (1.0 - sigma))
    hg_infl = torch.stack([hg1, hg2], dim=2)                                # (S,Q,2)

    ratio_nxt = s_nxt ** (1.0 / (1.0 - aalpha)) * zf_nxt0 / zf_nxt1 * l_nxt[:, :, 0] / l_nxt[:, :, 1]
    kappa_nxt_h = k_nxt * ratio_nxt / (1.0 + ratio_nxt)
    pi_nxt = aalpha * torch.exp((torch.log(zf_nxt0) + dz) * (1.0 - aalpha)) \
        * l_nxt[:, :, 0] ** (1.0 - aalpha) * kappa_nxt_h ** (aalpha - 1.0)
    rk_vec = ((1.0 - ddelta) * q_nxt + pi_nxt) / g.q.unsqueeze(1)
    rk_vec = rk_vec.clone()
    rk_vec[:, Q - 1] = rk_vec[:, Q - 1] * dz_exp[Q - 1]

    Pn1 = (vs[0] + (1 - vs[0]) * s_nxt ** (sigma - 1)) ** (1.0 / (1.0 - sigma))
    Pn2 = s_nxt ** (-1) * (vs[1] * s_nxt ** (1 - sigma) + (1 - vs[1])) ** (1.0 / (1.0 - sigma))
    P_nxt = torch.stack([Pn1, Pn2], dim=2)                                  # (S,Q,2)

    omg_nxt = torch.exp(const.next_omg_mat) + const.omg_shift
    if const.bg_yss > 0:
        seig = cf_sp_nxt * const.bg_yss * omg_nxt
        seignorage = torch.stack([seig, -seig], dim=2)
    else:
        seignorage = torch.zeros((S, Q, n_I), dtype=f64, device=next_state.device)

    rf0 = g.nom_i[:, 0].unsqueeze(1) / hg_infl[:, :, 0]
    rf1 = (g.nom_i[:, 0] + g.nom_i[:, 1]).unsqueeze(1) / hg_infl[:, :, 1]
    rf_home = rf0 / (1.0 - omg).unsqueeze(1)

    big_w = big_weight(const)
    M_vec = torch.empty((S, Q, n_I), dtype=f64, device=next_state.device)
    for iii in range(n_I):
        ies = const.ies_vec[iii]; gma = const.gma_vec[iii]; chi0 = const.chi0_vec[iii]
        cons = c_spend[:, iii] / P[:, iii]
        labor = g.l_aggr[:, iii]
        r_alpha = (g.share[:, iii] - g.bF_share[:, iii]).unsqueeze(1) * rf_home \
            + g.bF_share[:, iii].unsqueeze(1) * rf1 + (1.0 - g.share[:, iii]).unsqueeze(1) * rk_vec
        nps = (savings[:, iii].unsqueeze(1) * r_alpha
               + dz_exp.unsqueeze(0) * const.q_l_ss_vec[iii] + seignorage[:, :, iii]) \
            / const.tot_wealth_ss_vec[iii]
        EV = ((big_w * (v_temp[:, :, iii] * nps) ** (1.0 - gma)).sum(dim=1)) ** (1.0 / (1.0 - gma))
        labor_part = 1.0 + (1.0 / ies - 1.0) * chi0 * const.chi / (1.0 + const.chi) \
            * labor ** ((1.0 + const.chi) / const.chi)
        util_c_deriv = (labor_part / cons) ** (1.0 / ies)
        M_vec[:, :, iii] = bbeta[:, iii].unsqueeze(1) * mc_temp[:, :, iii] / util_c_deriv.unsqueeze(1) \
            * (v_temp[:, :, iii] / EV.unsqueeze(1)) ** (1.0 / ies - gma) * nps ** (-gma)

    return PricingBlock(M_vec=M_vec, hg_infl=hg_infl, P=P, P_nxt=P_nxt, omg=omg,
                        next_state=next_state)


def _price_one_step(nxt_bp, pb: PricingBlock, big_w):
    """One calc_bond_prices step: nxt_bp (S,Q,3) [h,f,hw] -> q_bond (S,3), each the
    marginal-investor (maxval over agent) valuation."""
    S, Q, _ = nxt_bp.shape
    M, hg, P, Pn, omg = pb.M_vec, pb.hg_infl, pb.P, pb.P_nxt, pb.omg
    out = torch.empty((S, 3), dtype=nxt_bp.dtype, device=nxt_bp.device)
    # currency h: home-good inflation (col 0), no wedge
    vh = torch.stack([(big_w * M[:, :, i] * (nxt_bp[:, :, 0] / hg[:, :, 0])
                       * P[:, i].unsqueeze(1) / Pn[:, :, i]).sum(dim=1) for i in range(2)], dim=1)
    out[:, 0] = vh.max(dim=1).values
    # currency f: foreign-good inflation (col 1)
    vf = torch.stack([(big_w * M[:, :, i] * (nxt_bp[:, :, 1] / hg[:, :, 1])
                       * P[:, i].unsqueeze(1) / Pn[:, :, i]).sum(dim=1) for i in range(2)], dim=1)
    out[:, 1] = vf.max(dim=1).values
    # currency hw: home-good inflation, /(1-omg)
    vhw = torch.stack([(big_w * M[:, :, i] * (nxt_bp[:, :, 2] / hg[:, :, 0]) / (1.0 - omg).unsqueeze(1)
                        * P[:, i].unsqueeze(1) / Pn[:, :, i]).sum(dim=1) for i in range(2)], dim=1)
    out[:, 2] = vhw.max(dim=1).values
    return out


def compute_bond_columns(const, g, v_mat, mc_mat, next_state):
    """Run the maturity recursion and return grid-level bond-price columns needed by
    the Table-2 bond moments: dict name -> (S,) tensor. Currencies h/f/hw."""
    pb = compute_pricing_block(const, g, v_mat, mc_mat, next_state)
    S, Q, d = next_state.shape
    lu = interp_factor(const.smol_polynom)
    B = smolyak_polynomial(next_state.reshape(S * Q, d), const.elem, const.mu)  # (S*Q, S)
    big_w = big_weight(const)

    q_bond = torch.ones((S, 3), dtype=DTYPE, device=next_state.device)      # maturity 0
    stored = {}
    for m in range(1, N_BOND + 1):
        coeffs = interp_solve(lu, q_bond)                                   # (S,3)
        nxt_bp = (B @ coeffs).reshape(S, Q, 3)
        q_bond = _price_one_step(nxt_bp, pb, big_w)
        if m in STORED_MATS:
            stored[m] = q_bond.clone()

    cols = {
        "q1_h": stored[1][:, 0], "q1_f": stored[1][:, 1], "q1_hw": stored[1][:, 2],
        "q19_h": stored[19][:, 0], "q19_f": stored[19][:, 1],
        "q20_h": stored[20][:, 0], "q20_f": stored[20][:, 1],
    }
    return cols
