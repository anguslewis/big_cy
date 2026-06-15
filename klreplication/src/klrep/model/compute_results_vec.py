"""Grid-level Table-2 economic variables — port of the `results_vec` assembly in
`calc_bond_prices` (mod_calc.f90:3449-3493), restricted to the deterministic
columns the Table-2 moments need (the bond-ladder columns are deferred; see
PLAN sect 13.6).

KL computes `results_mat` (one results_vec per grid point) from the CONVERGED
policies, then the MATLAB post-processing interpolates that matrix along the
simulated standardized state path. We replicate that: this module produces the
grid-level columns; `post/table2_series.py` fits Smolyak coeffs to them and
interpolates at the simulated `std_series`. Computing `inv = k_next_new -
(1-delta)*k` and `savings` on the grid (where `k_next_new = sum savings*(1-share)/q`
is exact) avoids the de-trending pitfalls of reconstructing them from consecutive
simulated economic states.

All quantities are in the solver's de-trended ("re-scaled economy") units; the
caller re-trends level columns by `* exp(cumsum(z_shock))`. Returns are scale-free
and are NOT re-trended.
"""
from dataclasses import dataclass

import torch

from ..params import IDX_K, IDX_ZF
from .period_block import compute_current_period

# Column order of the returned matrix (one (S,) column each). Matches the names
# the de-trending layer expects; the comment gives the Fortran results_vec index.
RESULT_NAMES = [
    "yh",       # results_vec(1)  y_current(1)            [level -> *z]
    "yf",       # results_vec(2)  y_current(2)            [level -> *z]
    "lh",       # results_vec(3)  l_aggr(1)               [stationary]
    "lf",       # results_vec(4)  l_aggr(2)               [stationary]
    "ch",       # results_vec(5)  c_vec(1)=c_spend/P1     [level -> *z]
    "cf",       # results_vec(6)  c_vec(2)=c_spend/P2     [level -> *z]
    "inv",      # results_vec(7)  k_nxt(1)-(1-d)k         [level -> *z]
    "inv_h",    # results_vec(38) inv_h/P1 (home-good inv) [level -> *z]
    "pi",       # results_vec(8)  pi_current(1)/P1        [ratio]
    "P_Phh",    # results_vec(10) P_div_P_h(1)            [ratio]
    "P_Phf",    # results_vec(11) P_div_P_h(2)            [ratio]
    "qx",       # results_vec(12) P_div_P_h(1)/P_div_P_h(2)[ratio]
    "s",        # results_vec(13) terms of trade          [ratio]
    "q",        # results_vec(17) q_current/P1            [ratio]
    "h_kap",    # results_vec(30) kappa_vec(1)            [level -> *z]
    "h_ksav",   # results_vec(29) savings(1)*(1-share)/P1 [level -> *z]
    "h_bh_sav", # results_vec(40) savings(1)*(share-bF)/P1[level -> *z]
    "h_sav",    # results_vec(44) savings(1)/P1           [level -> *z]
    "infl_h",   # results_vec(25) infl_vec(1)             [gross]
    "infl_f",   # results_vec(26) infl_vec(2)             [gross]
    "nom_ih",   # results_vec(27) nom_i(1)                [gross]
    "nom_if",   # results_vec(28) nom_i(1)+nom_i(2)       [gross]
]

# Names that are de-trended levels (re-trended by * exp(cumsum z_shock)).
LEVEL_NAMES = {"yh", "yf", "ch", "cf", "inv", "inv_h", "h_kap", "h_ksav",
               "h_bh_sav", "h_sav"}


@dataclass
class ResultsGrid:
    mat: torch.Tensor             # (S, len(RESULT_NAMES))
    names: list                   # RESULT_NAMES
    level_mask: torch.Tensor      # (ncols,) bool, True where a de-trended level


def compute_results_vec(const, g) -> ResultsGrid:
    """Build the grid-level Table-2 columns from the converged guesses `g`."""
    sg = const.state_grid
    k_aggr = sg[:, IDX_K]
    zf_state = sg[:, IDX_ZF]
    ddelta = const.ddelta

    # STEP 1 current production (same block the equilibrium step uses).
    cp = compute_current_period(
        k_aggr, g.s, g.l_aggr, g.q, zf_state, const.wealth_share_grid,
        aalpha=const.aalpha, sigma=const.sigma, ddelta=ddelta,
        zeta=const.zeta, varsigma=const.varsigma_vec,
    )
    P = cp.p_div_ph                                   # (S,2)  [P1, P2]
    w_choice = cp.w_choice
    wealth_vec = cp.wealth_vec

    # savings + consumption (clamp c_spend exactly as Fortran / equilibrium_step).
    msc = const.min_cons_sav
    c_spend = torch.minimum(
        torch.maximum(g.c_spending, torch.full_like(g.c_spending, msc)),
        wealth_vec + w_choice * g.l_aggr - msc,
    )
    savings = wealth_vec + w_choice * g.l_aggr - c_spend       # (S,2)
    c_vec = c_spend / P                                        # (S,2) consumption

    # next aggregate capital and investment (exact on the grid).
    k_next_new = (savings * (1.0 - g.share) / g.q.unsqueeze(1)).sum(dim=1)
    inv = k_next_new - (1.0 - ddelta) * k_aggr
    # home-good investment spending (equilibrium_step STEP 7 / results_vec 38).
    ish = const.inv_share_h
    inv_h = inv * g.q / (1.0 + g.s ** (const.sigma - 1.0) * (1.0 - ish) / ish)

    cols = {
        "yh": cp.y_current[:, 0],
        "yf": cp.y_current[:, 1],
        "lh": g.l_aggr[:, 0],
        "lf": g.l_aggr[:, 1],
        "ch": c_vec[:, 0],
        "cf": c_vec[:, 1],
        "inv": inv,
        "inv_h": inv_h / P[:, 0],
        "pi": cp.pi_current[:, 0] / P[:, 0],
        "P_Phh": P[:, 0],
        "P_Phf": P[:, 1],
        "qx": P[:, 0] / P[:, 1],
        "s": g.s,
        "q": g.q / P[:, 0],
        "h_kap": cp.kappa[:, 0],
        "h_ksav": savings[:, 0] * (1.0 - g.share[:, 0]) / P[:, 0],
        "h_bh_sav": savings[:, 0] * (g.share[:, 0] - g.bF_share[:, 0]) / P[:, 0],
        "h_sav": savings[:, 0] / P[:, 0],
        "infl_h": g.infl[:, 0],
        "infl_f": g.infl[:, 1],
        "nom_ih": g.nom_i[:, 0],
        "nom_if": g.nom_i[:, 0] + g.nom_i[:, 1],
    }
    mat = torch.stack([cols[n] for n in RESULT_NAMES], dim=1)
    level_mask = torch.tensor([n in LEVEL_NAMES for n in RESULT_NAMES],
                              dtype=torch.bool, device=mat.device)
    return ResultsGrid(mat=mat, names=list(RESULT_NAMES), level_mask=level_mask)
