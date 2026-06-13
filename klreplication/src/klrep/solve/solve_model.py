"""Outer backward-iteration solver — port of the iteration loop in
`calc_sol` (mod_calc.f90:204-709). Dampened fixed point on the policy/value/price
arrays over the Smolyak grid, with staged activation (bond clearing after iter 10,
foreign trading after iter 100) and per-object damping weights (PLAN sect 13.2).

Tensor-native: every per-grid-point quantity carries a leading n_states batch dim;
the equilibrium step is evaluated for all states at once.
"""
import math
from dataclasses import dataclass

import torch

from ..config import DTYPE
from ..params import (Params, IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF, IDX_DIS,
                      IDX_IH, IDX_IF, IDX_OMG)
from ..solve.steady_state import calc_steady, SteadyState
from ..setup.smolyak_setup import build_smolyak
from ..setup.shock_grid import build_shock_grid
from ..setup.state_grid import build_state_grid
from ..solve.model_const import build_model_const, big_weight, ModelConst
from ..grid.smolyak_polynomial import smolyak_polynomial
from ..grid.interp_fit import interp_factor, interp_solve
from ..model.equilibrium_step import equilibrium_step, Guesses

_CONV = 1e-8


@dataclass
class SolverState:
    const: ModelConst
    g: Guesses
    v_mat: torch.Tensor      # (S,2)
    mc_mat: torch.Tensor     # (S,2)
    next_state: torch.Tensor  # (S, Q, d) standardized next states (active dims)
    big_w: torch.Tensor      # (S, Q)
    k_next_mat: torch.Tensor  # (S, Q)
    w_choice: torch.Tensor   # (S, 2) production wage (recomputed; stored for diff)


def _standardized_next(const, k_next, theta_nxt, w_choice, nom_i, s_in,
                       *, theta_is_zero=False):
    """Build the (S, Q, d) standardized next-state for the active dims, clamped."""
    S = const.state_grid.shape[0]
    Q = const.n_quad
    dz_adj_exp = torch.exp(const.dz_vec_adj).unsqueeze(0)   # (1,Q)
    m, dv = const.grid_means, const.grid_devs
    sigma = const.sigma
    vs = const.varsigma_vec
    s_col = s_in.unsqueeze(1)                                # (S,1)

    cols = {}
    cols[IDX_K] = (k_next.unsqueeze(1) / dz_adj_exp - m[IDX_K]) / dv[IDX_K]
    if theta_is_zero:
        cols[IDX_THH] = torch.zeros((S, Q), dtype=DTYPE, device=const.state_grid.device)
    else:
        cols[IDX_THH] = (theta_nxt - m[IDX_THH]) / dv[IDX_THH]
    cols[IDX_ZF] = (const.next_zf_mat - m[IDX_ZF]) / dv[IDX_ZF]
    A1 = (vs[0] + (1 - vs[0]) * s_col ** (sigma - 1.0))
    cols[IDX_WH] = (w_choice[:, 0].unsqueeze(1) / dz_adj_exp / A1 ** (1.0 / (sigma - 1.0))
                    - m[IDX_WH]) / dv[IDX_WH]
    A2 = (vs[1] * s_col ** (1.0 - sigma) + (1 - vs[1]))
    cols[IDX_WF] = (w_choice[:, 1].unsqueeze(1) / dz_adj_exp * s_col / A2 ** (1.0 / (sigma - 1.0))
                    - m[IDX_WF]) / dv[IDX_WF]
    cols[IDX_DIS] = (const.next_dis_mat - m[IDX_DIS]) / dv[IDX_DIS]
    cols[IDX_IH] = ((nom_i[:, 0] - 1.0).unsqueeze(1).expand(S, Q) - m[IDX_IH]) / dv[IDX_IH]
    cols[IDX_IF] = ((nom_i[:, 0] + nom_i[:, 1] - 1.0).unsqueeze(1).expand(S, Q) - m[IDX_IF]) / dv[IDX_IF]
    cols[IDX_OMG] = (const.next_omg_mat - m[IDX_OMG]) / dv[IDX_OMG]

    stacked = torch.stack([cols[d] for d in const.active_dims], dim=2)  # (S,Q,d)
    return stacked.clamp(-1.0, 1.0)


def init_solver(p: Params, ss: SteadyState = None, *, device=None) -> SolverState:
    if ss is None:
        ss = calc_steady(p)
    smol = build_smolyak(p, device=device)
    sg = build_shock_grid(p, device=device)
    stg = build_state_grid(p, smol, ss, sg, device=device)
    const = build_model_const(p, ss, sg, stg, smol, device=device)
    S, Q = const.state_grid.shape[0], const.n_quad
    f64 = DTYPE
    state_grid = const.state_grid

    # initial guesses (calc_sol:304-352)
    l_aggr = torch.ones((S, 2), dtype=f64, device=device)
    q = torch.full((S,), ss.q_ss, dtype=f64, device=device)
    nom_i0 = (1.0 + ss.rf_ss) / (1.0 + torch.exp(state_grid[:, IDX_OMG]))
    nom_i = torch.stack([nom_i0, (1.0 + ss.rf_ss) - nom_i0], dim=1)
    s = torch.full((S,), ss.s_ss, dtype=f64, device=device)
    infl = torch.ones((S, 2), dtype=f64, device=device)
    share = torch.zeros((S, 2), dtype=f64, device=device)
    bF_share = torch.zeros((S, 2), dtype=f64, device=device)

    zf2 = const.zeta * torch.exp(state_grid[:, IDX_ZF])
    zf_tmp = torch.stack([torch.ones_like(zf2), zf2], dim=1)
    aalpha = const.aalpha
    ratio = (s ** (1.0 / (1.0 - aalpha)) * zf_tmp[:, 0] / zf_tmp[:, 1]
             * l_aggr[:, 0] / l_aggr[:, 1])
    kappa_h = state_grid[:, IDX_K] * ratio / (1.0 + ratio)
    kappa = torch.stack([kappa_h, state_grid[:, IDX_K] - kappa_h], dim=1)
    w_choice = torch.stack([
        (1 - aalpha) * kappa[:, 0] ** aalpha * 1.0 ** (-aalpha) * zf_tmp[:, 0] ** (1 - aalpha),
        (1 - aalpha) * kappa[:, 1] ** aalpha * 1.0 ** (-aalpha) * ss.s_ss ** (-1) * zf_tmp[:, 1] ** (1 - aalpha),
    ], dim=1)
    c_spending = torch.stack([
        w_choice[:, 0] * 1.0 + state_grid[:, IDX_K] * const.wealth_share_grid[:, 0] * (1.0 / const.bbeta_vec[0] - 1.0),
        w_choice[:, 1] * 1.0 + state_grid[:, IDX_K] * const.wealth_share_grid[:, 1] * (1.0 / const.bbeta_vec[1] - 1.0),
    ], dim=1)

    g = Guesses(c_spending=c_spending, s=s, l_aggr=l_aggr, q=q, infl=infl,
                nom_i=nom_i, share=share, bF_share=bF_share)
    v_mat = const.tot_wealth_ss_vec.new_tensor([0.0])  # placeholder
    v_mat = torch.stack([torch.full((S,), 1.0, dtype=f64, device=device),
                         torch.full((S,), 1.0, dtype=f64, device=device)], dim=1)
    # v_ss == 1 by normalization; mc from SS
    import math as _m
    mc_mat = torch.stack([torch.full((S,), ss.mc_ss[0], dtype=f64, device=device),
                          torch.full((S,), ss.mc_ss[1], dtype=f64, device=device)], dim=1)

    # initial k_next_mat (calc_sol:362) and next_state (theta dim = 0 initially)
    dz_exp = torch.exp(const.dz_vec).unsqueeze(0)
    dz_adj_exp = torch.exp(const.dz_vec_adj).unsqueeze(0)
    k_next_mat = state_grid[:, IDX_K].unsqueeze(1) * dz_exp / dz_adj_exp   # (S,Q)
    next_state = _standardized_next(const, state_grid[:, IDX_K], None, w_choice, nom_i, s,
                                    theta_is_zero=True)
    big_w = big_weight(const)
    return SolverState(const=const, g=g, v_mat=v_mat, mc_mat=mc_mat,
                       next_state=next_state, big_w=big_w, k_next_mat=k_next_mat,
                       w_choice=w_choice)


def _interp_input(st: SolverState):
    g = st.g
    return torch.stack([
        st.v_mat[:, 0], st.v_mat[:, 1], st.mc_mat[:, 0], st.mc_mat[:, 1],
        g.s, g.q, g.l_aggr[:, 0], g.l_aggr[:, 1], g.infl[:, 0], g.infl[:, 1],
        g.c_spending[:, 0], g.c_spending[:, 1]], dim=1)   # (S,12)


def solve_step(st: SolverState, outer_iter: int):
    """One outer iteration (fit -> interpolate -> equilibrium_step -> damp ->
    rebuild next_state). Returns (new SolverState, diff)."""
    const = st.const
    S, Q = const.state_grid.shape[0], const.n_quad
    # fit interpolation coefficients
    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, _interp_input(st))             # (T,12)
    # evaluate basis at all Q next-states for all S grid points
    pts = st.next_state.reshape(S * Q, -1)                   # (S*Q, d)
    B = smolyak_polynomial(pts, const.elem, const.mu)        # (S*Q, T)
    nxt_interp = (B @ coeffs).reshape(S, Q, 12)
    nxt_mat2 = torch.stack([st.k_next_mat, const.next_zf_mat, const.next_omg_mat], dim=2)

    out = equilibrium_step(const, st.g, nxt_interp, nxt_mat2, st.big_w, outer_iter)
    gn, g = out.g, st.g

    # next_state rebuild from new policies (uses input s)
    next_state_new = _standardized_next(const, out.k_next_new, out.theta_nxt,
                                        out.w_choice, gn.nom_i, g.s)
    k_next_mat_new = out.k_next_new.unsqueeze(1) / torch.exp(const.dz_vec_adj).unsqueeze(0) \
        * torch.exp(const.dz_vec).unsqueeze(0)

    # convergence diff (PLAN 13.2: max of abs-log-ratios + abs diffs)
    def amaxlog(a, b):
        return (torch.log(a / b)).abs().max()
    diff = float(torch.stack([
        amaxlog(out.v, st.v_mat), amaxlog(out.mc, st.mc_mat),
        amaxlog(gn.c_spending, g.c_spending), amaxlog(k_next_mat_new, st.k_next_mat),
        amaxlog(gn.l_aggr, g.l_aggr), amaxlog(out.w_choice, st.w_choice),
        amaxlog(gn.q, g.q), amaxlog(gn.infl, g.infl),
        (gn.share - g.share).abs().max(), (gn.nom_i - g.nom_i).abs().max(),
        (next_state_new - st.next_state).abs().max(),
    ]).max())

    # damping (PLAN 13.2)
    v_mat = st.v_mat + 1.0 * (out.v - st.v_mat)
    mc_mat = st.mc_mat + 1.0 * (out.mc - st.mc_mat)
    next_state = st.next_state + 0.1 * (next_state_new - st.next_state)
    q = g.q + 0.2 * (gn.q - g.q)
    s = g.s + 0.1 * (gn.s - g.s)
    l_aggr = g.l_aggr + 0.2 * (gn.l_aggr - g.l_aggr)
    w_choice = st.w_choice * (1.0 + 0.2 * torch.log(out.w_choice / st.w_choice))  # log-damp
    c_spending = g.c_spending + 0.2 * (gn.c_spending - g.c_spending)
    nom_i = g.nom_i + 0.1 * (gn.nom_i - g.nom_i)
    infl = g.infl + 0.05 * (gn.infl - g.infl)
    k_next_mat = st.k_next_mat + 0.2 * (k_next_mat_new - st.k_next_mat)
    share = g.share + 0.1 * (gn.share - g.share)
    bF_w = 0.0 if outer_iter < 101 else 0.05
    bF_share = g.bF_share + bF_w * (gn.bF_share - g.bF_share)

    g2 = Guesses(c_spending=c_spending, s=s, l_aggr=l_aggr, q=q, infl=infl,
                 nom_i=nom_i, share=share, bF_share=bF_share)
    st2 = SolverState(const=const, g=g2, v_mat=v_mat, mc_mat=mc_mat,
                      next_state=next_state, big_w=st.big_w, k_next_mat=k_next_mat,
                      w_choice=w_choice)
    return st2, diff


def solve_model(p: Params, *, max_iter=5000, conv=_CONV, device=None, verbose=False):
    # `conv` is the OUTER convergence tolerance (the inner step fixed points stay
    # tight at _CONV=1e-8). 1e-8 is KL's; 1e-6 is ample for SS/Table-2 validation
    # and converges much sooner (the linear tail is slow).
    st = init_solver(p, device=device)
    diff = 1.0
    it = 0
    while diff > conv and it < max_iter:
        it += 1
        st, diff = solve_step(st, it)
        if verbose and (it <= 12 or it % 50 == 0):
            print(f"iter {it:5d}  diff {diff:.3e}")
    return st, diff, it
