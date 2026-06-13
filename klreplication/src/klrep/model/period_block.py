"""Current-period production/prices, price of capital, and next-period returns —
ports of STEP 1, STEP 3, and STEP 2 of `calc_equilibrium_and_update`
(mod_calc.f90:2269-2396, 2398-2406).

These are the static (non-root-finding) algebra blocks of the equilibrium step,
batched tensor-native: leading dim S = grid points, quadrature dim Q last where
present. They are validated against the deterministic steady state (a no-shock
fixed point reproduces SS capital allocation, wages, q, and the SS capital
return 1+rk_ss).
"""
import math
from dataclasses import dataclass

import torch

SQRT_EPS = math.sqrt(2.220446049250313e-16)


@dataclass
class CurrentPeriod:
    kappa: torch.Tensor        # (S, 2) capital allocation [home, foreign]
    pi_current: torch.Tensor   # (S, 2) profit per unit capital, by country
    y_current: torch.Tensor    # (S, 2) output
    w_choice: torch.Tensor     # (S, 2) real wages (home-good units)
    p_div_ph: torch.Tensor     # (S, 2) P/P_h, P*/P*_f (bundle cost in home good)
    aggr_wealth: torch.Tensor  # (S,)
    wealth_vec: torch.Tensor   # (S, 2) wealth by agent


def compute_current_period(k_aggr, s, l_aggr, q_current, zf_state, wealth_share,
                           *, aalpha, sigma, ddelta, zeta, varsigma):
    """STEP 1 (mod_calc.f90:2273-2302). All inputs batched on leading dim S;
    l_aggr, wealth_share, varsigma are (S,2)/(2,); returns CurrentPeriod."""
    zf2 = zeta * torch.exp(zf_state)              # (S,)
    ratio = s ** (1.0 / (1.0 - aalpha)) * (1.0 / zf2) * (l_aggr[:, 0] / l_aggr[:, 1])
    kappa_h = k_aggr * ratio / (1.0 + ratio)
    kappa = torch.stack([kappa_h, k_aggr - kappa_h], dim=1)  # (S,2)

    zf_vec = torch.stack([torch.ones_like(zf2), zf2], dim=1)  # (S,2)
    s_vec = torch.stack([torch.ones_like(s), s], dim=1)       # (S,2)
    pi_current = (s_vec ** (-1)) * aalpha * (zf_vec * l_aggr / kappa) ** (1.0 - aalpha)
    y_current = zf_vec ** (1.0 - aalpha) * kappa ** aalpha * l_aggr ** (1.0 - aalpha)
    w_choice = (s_vec ** (-1)) * (1.0 - aalpha) * y_current / l_aggr

    vs = varsigma if torch.is_tensor(varsigma) else torch.tensor(varsigma, dtype=s.dtype)
    p1 = (vs[0] + (1.0 - vs[0]) * s ** (sigma - 1.0)) ** (1.0 / (1.0 - sigma))
    p2 = s ** (-1) * (vs[1] * s ** (1.0 - sigma) + (1.0 - vs[1])) ** (1.0 / (1.0 - sigma))
    p_div_ph = torch.stack([p1, p2], dim=1)

    aggr_wealth = (kappa[:, 0] * ((1 - ddelta) * q_current + pi_current[:, 0])
                   + kappa[:, 1] * ((1 - ddelta) * q_current + pi_current[:, 1]))
    wealth_vec = aggr_wealth.unsqueeze(1) * wealth_share
    return CurrentPeriod(kappa, pi_current, y_current, w_choice, p_div_ph,
                         aggr_wealth, wealth_vec)


def compute_q_new(s, next_k, k_aggr, *, inv_share_h, sigma, chiX):
    """STEP 3 price of capital (mod_calc.f90:2403-2405)."""
    base = inv_share_h ** (1.0 / (1.0 - sigma)) * (
        1.0 + (1.0 - inv_share_h) / inv_share_h * s ** (sigma - 1.0)
    ) ** (1.0 / (1.0 - sigma))
    return base * (next_k / k_aggr) ** chiX
