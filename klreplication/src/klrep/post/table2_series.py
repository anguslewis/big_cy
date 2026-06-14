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


def build_table2_series(const, g, sim, *, disast_shock):
    """Interpolate the grid results at the simulated path and build named series.

    Parameters
    ----------
    const, g : solver constants and converged guesses (for compute_results_vec).
    sim : dict from simulate_ensemble — needs `std_series` (n_sims,T,d) and
          `z_shock_series` (n_sims,T).
    disast_shock : float — the disaster TFP drop (to isolate `dis_series`).

    Returns
    -------
    dict[str -> (n_sims, T) tensor] of de-trended named series.
    """
    rg = compute_results_vec(const, g)                      # ResultsGrid
    lu = interp_factor(const.smol_polynom)
    coeffs = interp_solve(lu, rg.mat)                       # (S, ncols)

    std = sim["std_series"]                                 # (n_sims, T, d)
    n_sims, T, d = std.shape
    B = smolyak_polynomial(std.reshape(n_sims * T, d), const.elem, const.mu)
    interp = (B @ coeffs).reshape(n_sims, T, -1)           # (n_sims, T, ncols)

    # de-trend: level columns * exp(cumsum z_shock); ratios/returns untouched.
    z_shock = sim["z_shock_series"]                         # (n_sims, T)
    z_series = torch.exp(torch.cumsum(z_shock, dim=1))     # (n_sims, T)
    cols = {}
    for j, name in enumerate(rg.names):
        c = interp[:, :, j]
        cols[name] = c * z_series if rg.level_mask[j] else c

    # disaster component of the z-shock (extract_series.m:15-17).
    is_dis = (z_shock >= -disast_shock - SQRT_EPS) & (z_shock <= -disast_shock + SQRT_EPS)
    dis = torch.where(is_dis, z_shock, torch.zeros_like(z_shock))
    cols["dis"] = dis

    # realized home real safe rate: log(nom_ih_{t-1} / infl_h_t), prepend first.
    rfh = torch.log(cols["nom_ih"][:, :-1] / cols["infl_h"][:, 1:])
    cols["rfh"] = _prepend_first(rfh)

    # net foreign assets: h_sav_t - h_kap_{t+1}*exp(-dis_{t+1})*q_t, duplicate LAST.
    nfa_base = (cols["h_sav"][:, :-1]
                - cols["h_kap"][:, 1:] * torch.exp(-cols["dis"][:, 1:]) * cols["q"][:, :-1])
    cols["nfa"] = torch.cat([nfa_base, nfa_base[:, -1:]], dim=1)
    cols["nfa_rel"] = cols["nfa"] / cols["yh"]

    return cols
