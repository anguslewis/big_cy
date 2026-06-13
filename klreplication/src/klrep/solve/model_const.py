"""ModelConst — all the static tensors/scalars the equilibrium step + outer loop
need, assembled once from Params + SteadyState + ShockGrid + StateGrid +
SmolyakSetup. Keeps the per-step functions free of 30-argument signatures.
"""
from dataclasses import dataclass
from typing import List

import torch

from ..config import DTYPE
from ..params import Params
from ..solve.steady_state import SteadyState
from ..setup.shock_grid import ShockGrid
from ..setup.state_grid import StateGrid
from ..setup.smolyak_setup import SmolyakSetup


@dataclass
class ModelConst:
    # scalar params
    aalpha: float; sigma: float; ddelta: float; inv_share_h: float; zeta: float
    chi: float; chiX: float; vareps_w: float; phi_w: float
    phi_h: float; phi_f: float; phi_yh: float; phi_yf: float
    tayl_ic_h: float; tayl_ic_f: float; rho_i: float
    bg_yss: float; b_lmbd: float; omg_shift: float; bbeta_coeff: float
    tht_h_grid_mean: float
    # per-agent tensors (2,)
    bbeta_vec: torch.Tensor; ies_vec: torch.Tensor; gma_vec: torch.Tensor
    varsigma_vec: torch.Tensor; chi0_vec: torch.Tensor
    v_normalization_vec: torch.Tensor
    q_l_ss_vec: torch.Tensor; tot_wealth_ss_vec: torch.Tensor
    # shock grid
    dz_vec: torch.Tensor; dz_vec_adj: torch.Tensor
    quad_weight_vec: torch.Tensor       # (n_gh+1,)
    n_quad: int; no_shock_idx: int
    # state grid
    state_grid: torch.Tensor            # (S,9)
    wealth_share_grid: torch.Tensor     # (S,2)
    next_zf_mat: torch.Tensor; next_dis_mat: torch.Tensor; next_omg_mat: torch.Tensor
    grid_means: torch.Tensor; grid_devs: torch.Tensor  # (9,)
    # smolyak
    elem: object; mu: int; active_dims: List[int]
    smol_polynom: torch.Tensor          # (S,S)
    # numerical constants
    min_cons_sav: float = 1e-8
    low_guess_fixed: float = -10.0
    high_guess_fixed: float = 5.0
    safe_low_guess_fixed: float = -10.0
    safe_high_guess_fixed: float = 10.0
    foreign_trading: int = 1


def build_model_const(p: Params, ss: SteadyState, sg: ShockGrid, stg: StateGrid,
                      smol: SmolyakSetup, *, dtype=DTYPE, device=None) -> ModelConst:
    t = lambda lst: torch.tensor(lst, dtype=dtype, device=device)
    return ModelConst(
        aalpha=p["aalpha"], sigma=p["sigma"], ddelta=p["delta"],
        inv_share_h=p["inv_share_h"], zeta=p.zeta, chi=p["chi"], chiX=p["chiX"],
        vareps_w=p["vareps_w"], phi_w=p["varphi_w"],
        phi_h=p["phi_pi_h"], phi_f=p["phi_pi_f"], phi_yh=p["phi_y_h"], phi_yf=p["phi_y_f"],
        tayl_ic_h=p["tayl_ic_h"], tayl_ic_f=p["tayl_ic_f"], rho_i=p["rho_i"],
        bg_yss=p["bg_yss"], b_lmbd=p["b_lmbd"], omg_shift=p["omg_shift"],
        bbeta_coeff=0.001, tht_h_grid_mean=p["tht_trgt_h"],
        bbeta_vec=t([p["bbeta_h"], p["bbeta_f"]]),
        ies_vec=t([p["ies_h"], p["ies_f"]]),
        gma_vec=t([p["gma_h"], p["gma_f"]]),
        varsigma_vec=t([p["varsigma_h"], p["varsigma_f"]]),
        chi0_vec=t(ss.chi0_vec), v_normalization_vec=t(ss.v_normalization_vec),
        q_l_ss_vec=t(ss.q_l_ss), tot_wealth_ss_vec=t(ss.tot_wealth_ss),
        dz_vec=sg.dz_vec.to(dtype=dtype, device=device),
        dz_vec_adj=sg.dz_vec_adj.to(dtype=dtype, device=device),
        quad_weight_vec=sg.quad_weight_vec.to(dtype=dtype, device=device),
        n_quad=sg.n_quad, no_shock_idx=sg.no_shock_idx,
        state_grid=stg.state_grid, wealth_share_grid=stg.wealth_share_grid,
        next_zf_mat=stg.next_zf_mat, next_dis_mat=stg.next_dis_mat,
        next_omg_mat=stg.next_omg_mat,
        grid_means=stg.grid_means, grid_devs=stg.grid_devs,
        elem=smol.elem, mu=smol.mu, active_dims=smol.active_dims,
        smol_polynom=smol.smol_polynom,
        foreign_trading=int(p["foreign_trading"]),
    )


def big_weight(const: ModelConst) -> torch.Tensor:
    """Per-state quadrature weights (S, n_quad) including the state-dependent
    disaster mass: [:, :n_quad-1] = quad_weight*(1-p_dis), [:, -1] = p_dis."""
    from ..params import IDX_DIS
    p_dis = torch.exp(const.state_grid[:, IDX_DIS])          # (S,)
    S = p_dis.shape[0]
    w = torch.empty((S, const.n_quad), dtype=p_dis.dtype, device=p_dis.device)
    w[:, : const.n_quad - 1] = const.quad_weight_vec.unsqueeze(0) * (1.0 - p_dis).unsqueeze(1)
    w[:, const.n_quad - 1] = p_dis
    return w
