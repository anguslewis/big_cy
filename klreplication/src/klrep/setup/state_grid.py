"""State grid + exogenous transitions — port of `mod_param.f90 :: grid_setup`.

Maps the Smolyak collocation grid (in [-1,1], active dims only) to economic state
values via per-dimension center/half-width, with inactive dims pinned to their
center. Then builds the next-period AR(1) values of the exogenous states
(z_F, p, omega) at every quadrature node from every current grid state.

Centers/half-widths: capital & wages come from the steady state (calc_steady);
the rest from parameters. Runs after calc_steady (main.f90 ordering:
init_setup -> calc_steady -> grid_setup).

Layout note: tensor-native (batch dim = states first). state_grid is
(n_states, 9) and next_*_mat are (n_states, n_quad) — the TRANSPOSE of the
Fortran (9, n_states) / (n_quad, n_states), so the leading dim batches over grid
points for the vectorized solver.
"""
from dataclasses import dataclass

import numpy as np
import torch

from ..config import DTYPE
from ..params import (Params, SQRT_EPS, IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF,
                      IDX_DIS, IDX_IH, IDX_IF, IDX_OMG)
from ..solve.steady_state import SteadyState
from .smolyak_setup import SmolyakSetup
from .shock_grid import ShockGrid


@dataclass
class StateGrid:
    state_grid: torch.Tensor       # (n_states, 9) economic state values
    wealth_share_grid: torch.Tensor  # (n_states, 2) [theta, 1-theta]
    next_zf_mat: torch.Tensor      # (n_states, n_quad) next z_F per node
    next_dis_mat: torch.Tensor     # (n_states, n_quad) next log-p per node
    next_omg_mat: torch.Tensor     # (n_states, n_quad) next omega per node
    grid_means: torch.Tensor       # (9,) per-dim centers
    grid_devs: torch.Tensor        # (9,) per-dim half-widths


def build_state_grid(p: Params, smol: SmolyakSetup, ss: SteadyState,
                     sg: ShockGrid, *, dtype=DTYPE, device=None) -> StateGrid:
    # Per-dimension centers and half-widths (9 slots; order = idx_*).
    means = np.zeros(9)
    devs = np.zeros(9)
    means[IDX_K] = ss.k_grid_mean
    devs[IDX_K] = ss.k_grid_dev
    means[IDX_THH] = p["tht_trgt_h"]
    devs[IDX_THH] = p["theta_h_dev"]
    means[IDX_ZF] = 0.0
    devs[IDX_ZF] = max(SQRT_EPS, 3.0 * p.zf_std)
    means[IDX_WH] = ss.wh_grid_mean
    devs[IDX_WH] = ss.wh_grid_dev
    means[IDX_WF] = ss.wf_grid_mean
    devs[IDX_WF] = ss.wf_grid_dev
    means[IDX_DIS] = p.disast_p_in
    devs[IDX_DIS] = 3.0 * p.disast_std_in
    means[IDX_IH] = p["ih_grid_mean"]
    devs[IDX_IH] = p["ih_grid_dev"]
    means[IDX_IF] = p["if_grid_mean"]
    devs[IDX_IF] = p["if_grid_dev"]
    means[IDX_OMG] = p["omg_mean"]
    devs[IDX_OMG] = max(SQRT_EPS, 3.0 * p.omg_std)

    # Map smol_grid columns (active dims, in order) to the 9 state dims.
    n_states = smol.n_states
    state_grid = np.empty((n_states, 9))
    col = 0
    for d in range(9):
        if p.mus_dims[d] > 0:
            state_grid[:, d] = smol.smol_grid[:, col] * devs[d] + means[d]
            col += 1
        else:
            state_grid[:, d] = means[d]

    wealth_share = np.stack(
        [state_grid[:, IDX_THH], 1.0 - state_grid[:, IDX_THH]], axis=1
    )

    # Exogenous AR(1) next-period values at each quad node (shock cols: z=0, zf=1,
    # p=2, omg=3). next = (1-rho)*mean + rho*current + shock.
    shock = sg.shock_grid.cpu().numpy()  # (n_quad, 4)
    zf_cur = state_grid[:, IDX_ZF][:, None]     # (n_states, 1)
    dis_cur = state_grid[:, IDX_DIS][:, None]
    omg_cur = state_grid[:, IDX_OMG][:, None]
    next_zf = (1.0 - p.rho_zf) * means[IDX_ZF] + p.rho_zf * zf_cur + shock[:, 1][None, :]
    next_dis = (1.0 - p.disast_rho) * p.disast_p_in + p.disast_rho * dis_cur + shock[:, 2][None, :]
    next_omg = (1.0 - p.omg_rho) * means[IDX_OMG] + p.omg_rho * omg_cur + shock[:, 3][None, :]

    t = lambda a: torch.tensor(a, dtype=dtype, device=device)
    return StateGrid(
        state_grid=t(state_grid),
        wealth_share_grid=t(wealth_share),
        next_zf_mat=t(next_zf),
        next_dis_mat=t(next_dis),
        next_omg_mat=t(next_omg),
        grid_means=t(means),
        grid_devs=t(devs),
    )
