"""Build the named Table-2 time series from the simulated path — port of the
relevant parts of `extract_series.m` (KL's MATLAB post-processing).

Pipeline (PLAN sect 13.6):
1. Fit Smolyak coeffs to the grid-level `results_vec` columns (compute_results_vec).
2. Interpolate them along the simulated standardized state path `std_series`
   (basis @ coeffs) — exactly what KL does (interpolate `results_mat`).
3. De-trend level columns by `* exp(cumsum(z_shock))`; ratios/returns untouched.
4. Construct the derived series the moments need: `dis`, `rfh`, `nfa`, `nfa_rel`.

Each returned series is a (n_sims, T) tensor (time along dim 1). Only the series
needed by the 12 deliverable moments {1-7, 11-15} are built here; the bond-ladder
series (rA/exc_retA, moments 8/9/10) are deferred (PLAN sect 13.6).
"""
import math

import torch

from ..config import DTYPE
from ..grid.interp_fit import interp_factor, interp_solve
from ..grid.smolyak_polynomial import smolyak_polynomial
from ..model.compute_results_vec import compute_results_vec

SQRT_EPS = math.sqrt(2.220446049250313e-16)


def _prepend_first(x):
    """MATLAB `[v(1); v]`: prepend the first element along time (dim 1)."""
    return torch.cat([x[:, :1], x], dim=1)


def build_table2_series(const, g, sim, *, disast_shock, bond_grid=None,
                        debt_to_equity=1.0, zeta=None, delta=None):
    """Interpolate the grid results at the simulated path and build named series.

    Parameters
    ----------
    const, g : solver constants and converged guesses (for compute_results_vec).
    sim : dict from simulate_ensemble — needs `std_series` (n_sims,T,d) and
          `z_shock_series` (n_sims,T).
    disast_shock : float — the disaster TFP drop (to isolate `dis_series`).
    bond_grid : optional dict name->(S,) of grid-level bond prices (from
        `bond_ladder.compute_bond_columns`). If given, the bond-dependent return
        series (yields, rq20, rfh_lev, rk, rA, exc_retA, exc_rf, uip, ...) needed
        by moments 8/9/10 are also built.
    debt_to_equity, zeta, delta : params for the levered-equity / 20-period blend
        (debt_to_equity=1.0; zeta=const.zeta; delta=const.ddelta if None).

    Returns
    -------
    dict[str -> (n_sims, T) tensor] of de-trended named series.
    """
    zeta = const.zeta if zeta is None else zeta
    delta = const.ddelta if delta is None else delta
    rg = compute_results_vec(const, g)                      # ResultsGrid

    # stack results + bond columns; interpolate all at the simulated path at once.
    names = list(rg.names)
    level = list(rg.level_mask)
    mat = rg.mat
    if bond_grid is not None:
        bnames = list(bond_grid.keys())
        names += bnames
        level += [False] * len(bnames)                     # bond prices are ratios
        mat = torch.cat([mat, torch.stack([bond_grid[n] for n in bnames], dim=1)], dim=1)

    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, mat)                          # (S, ncols)
    std = sim["std_series"]                                 # (n_sims, T, d)
    n_sims, T, d = std.shape
    B = smolyak_polynomial(std.reshape(n_sims * T, d), const.elem, const.mu)
    interp = (B @ coeffs).reshape(n_sims, T, -1)           # (n_sims, T, ncols)

    z_shock = sim["z_shock_series"]                         # (n_sims, T)
    z_series = torch.exp(torch.cumsum(z_shock, dim=1))     # (n_sims, T)
    cols = {}
    for j, name in enumerate(names):
        c = interp[:, :, j]
        cols[name] = c * z_series if level[j] else c

    # disaster component of the z-shock (extract_series.m:15-17).
    is_dis = (z_shock >= -disast_shock - SQRT_EPS) & (z_shock <= -disast_shock + SQRT_EPS)
    dis = torch.where(is_dis, z_shock, torch.zeros_like(z_shock))
    cols["dis"] = dis

    # realized home real safe rate: log(nom_ih_{t-1} / infl_h_t), prepend first.
    cols["rfh"] = _prepend_first(torch.log(cols["nom_ih"][:, :-1] / cols["infl_h"][:, 1:]))

    # net foreign assets: h_sav_t - h_kap_{t+1}*exp(-dis_{t+1})*q_t, duplicate LAST.
    nfa_base = (cols["h_sav"][:, :-1]
                - cols["h_kap"][:, 1:] * torch.exp(-cols["dis"][:, 1:]) * cols["q"][:, :-1])
    cols["nfa"] = torch.cat([nfa_base, nfa_base[:, -1:]], dim=1)
    cols["nfa_rel"] = cols["nfa"] / cols["yh"]
    # nfa growth relative to output (extract_series.m:320-321), prepend first.
    cols["nfa_rel_growth"] = _prepend_first(
        (cols["nfa"][:, 1:] - cols["nfa"][:, :-1]) / cols["yh"][:, 1:])

    if bond_grid is not None:
        _build_bond_series(cols, debt_to_equity=debt_to_equity, zeta=zeta, delta=delta)
    return cols


def _build_bond_series(cols, *, debt_to_equity, zeta, delta):
    """Add the bond-dependent return + regression series (moments 8/9/10) to `cols`,
    in place. Mirrors extract_series.m: realized bond returns, levered-equity return
    rA, excess returns, UIP residual, smoothed dividend-price, 4q output growth."""
    d2e = debt_to_equity
    infl_h, infl_f = cols["infl_h"], cols["infl_f"]
    P_hh, P_hf = cols["P_Phh"], cols["P_Phf"]
    q, pi = cols["q"], cols["pi"]

    # 1-period yields from the 1-period bond prices.
    cols["yield1_h"] = -torch.log(cols["q1_h"])
    cols["yield1_f"] = -torch.log(cols["q1_f"])
    cols["yield1_hw"] = -torch.log(cols["q1_hw"])

    # realized real-exchange-rate change term (home<-foreign basket), length T-1.
    xrate = torch.log(P_hf[:, 1:] / P_hf[:, :-1] / P_hh[:, 1:] * P_hh[:, :-1])

    # 20-period bond realized holding returns (buy 20-per at t-1 -> 19-per at t).
    rq20_h = _prepend_first(torch.log(cols["q19_h"][:, 1:] / infl_h[:, 1:] / cols["q20_h"][:, :-1]))
    rq20_f = _prepend_first(torch.log(cols["q19_f"][:, 1:] / infl_f[:, 1:] / cols["q20_f"][:, :-1]))
    rq20_fh = _prepend_first(torch.log(cols["q19_f"][:, 1:] / infl_f[:, 1:] / cols["q20_f"][:, :-1]) + xrate)

    # levered risk-free blend (1/(1+zeta) home-20 + zeta/(1+zeta) foreign-20-in-home).
    rfh_lev = torch.log(1.0 / (1.0 + zeta) * torch.exp(rq20_h)
                        + zeta / (1.0 + zeta) * torch.exp(rq20_fh))

    # realized capital return (with disaster correction), then levered equity return.
    rk = _prepend_first(torch.log(torch.exp(cols["dis"][:, 1:])
                                  * (pi[:, 1:] + (1.0 - delta) * q[:, 1:]) / q[:, :-1]))
    cols["rk"] = rk
    cols["rA"] = torch.log((1.0 + d2e) * torch.exp(rk) - d2e * torch.exp(rfh_lev))
    cols["exc_retA"] = cols["rA"] - cols["rfh"]

    # foreign safe rate in home basket, and the foreign-minus-home excess.
    rff = _prepend_first(torch.log(cols["nom_if"][:, :-1] / infl_f[:, 1:]))
    rff_h = _prepend_first(rff[:, 1:] + xrate)
    cols["exc_rf"] = rff_h - cols["rfh"]

    # expected log nominal depreciation (UIP forward term), prepend first.
    qx = cols["qx"]
    E_change = _prepend_first(-torch.log(qx[:, :-1] / qx[:, 1:] * infl_h[:, 1:] / infl_f[:, 1:]))
    cols["E_change"] = E_change
    # UIP residual on private/omega-adjusted yields (extract_series.m:469-470), prepend 0.
    uip_base = 100.0 * (cols["yield1_f"][:, :-1] - cols["yield1_hw"][:, :-1] - E_change[:, 1:])
    cols["uip_pvt"] = torch.cat([torch.zeros_like(uip_base[:, :1]), uip_base], dim=1)
    cols["fama_yield_pvt"] = 100.0 * 4.0 * (cols["yield1_f"] - cols["yield1_h"])

    # 4-quarter output growth (pct), zero-pad first 4.
    yg = 100.0 * (torch.log(cols["yh"][:, 4:]) - torch.log(cols["yh"][:, :-4]))
    cols["y_growth"] = torch.cat([torch.zeros_like(yg[:, :4]), yg], dim=1)

    # smoothed dividend-price on the levered equity claim.
    div1 = _prepend_first(pi[:, 1:] + (1.0 - delta) * q[:, 1:]
                          - d2e / (1.0 + d2e) * q[:, :-1] * torch.exp(rfh_lev[:, 1:])
                          - 1.0 / (1.0 + d2e) * q[:, 1:])
    qE1 = q / (1.0 + d2e)
    # trailing 4-quarter sum of dividends (movsum(div,[3,0])).
    cum = torch.cumsum(div1, dim=1)
    ms = cum.clone()
    ms[:, 4:] = cum[:, 4:] - cum[:, :-4]
    dp = ms / qE1
    dp0 = dp[:, :1].clone()
    dp[:, :3] = dp0 * 4.0
    cols["div_price_smoothed_1"] = dp
