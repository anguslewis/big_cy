"""Calibration parameters: load a param_file_*.csv into a typed Params object.

Faithful to `mod_param.f90 :: init_setup`'s read of `param_input_vec(1..65)` and
the scalar quantities it derives immediately after (disaster lognormal moments,
the AR(1) innovation SDs, the per-dimension Smolyak levels with the rho_i
adjustment, and the per-shock quadrature node counts). The covariance/Cholesky,
shock grid, and Smolyak grid construction live in the (separate) init_setup port;
this module is just the calibration inputs.

Path-touching module: imports project_strings for `raw` (the inputs live in
…/big_cy/data/raw/klreplication/params/, copied from the KL package — see
MANIFEST.md). Override the directory with `params_dir=` for tests.
"""
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np

# --- project_strings import (per coding conventions) ----------------------
work_dir = os.environ.get("work_dir", str(Path.home() / "code"))
sys.path.insert(0, work_dir + "/big_cy")
from project_strings import raw  # noqa: E402

SQRT_EPS = math.sqrt(np.finfo(np.float64).eps)  # matches base_lib sqrt_eps

# 65 column names, in the exact order of param_input_vec (CSV header order).
PARAM_NAMES: List[str] = [
    "zeta", "bbeta_h", "bbeta_f", "gma_h", "gma_f", "ies_h", "ies_f", "chi",
    "sigma", "varsigma_h", "varsigma_f", "l_target", "delta", "aalpha", "chiX",
    "inv_share_h", "varphi_w", "vareps_w", "phi_pi_h", "phi_pi_f", "phi_y_h",
    "phi_y_f", "tayl_ic_h", "tayl_ic_f", "rho_i", "b_lmbd", "bg_yss", "sig_z",
    "std_zf", "rho_zf", "disast_p", "disast_shock", "p_rho", "p_std", "omg_mean",
    "omg_std", "omg_shift", "omg_rho", "corr_omg_dis", "tht_trgt_h",
    "theta_h_dev", "k_dev_param", "k_grid_adj", "w_dev_param", "w_grid_adj",
    "ih_grid_mean", "if_grid_mean", "ih_grid_dev", "if_grid_dev", "mu_k",
    "mu_tht", "mu_zf", "mu_wh", "mu_wf", "mu_p", "mu_ih", "mu_if", "mu_omg",
    "n_quad_z", "n_quad_zf", "n_quad_p", "n_quad_omg", "foreign_trading",
    "run_bg", "run_samp",
]

# State-dimension indices (0-based) into the 9-slot mus vector; mirrors the
# Fortran idx_* (1-based): k, tht_h, zf, wh, wf, dis(p), ih, if, omg.
IDX_K, IDX_THH, IDX_ZF, IDX_WH, IDX_WF, IDX_DIS, IDX_IH, IDX_IF, IDX_OMG = range(9)


@dataclass
class Params:
    """Calibration for one KL specification (one param_file_*.csv)."""

    raw: Dict[str, float]  # the 65 values keyed by PARAM_NAMES

    # Derived scalars (computed in __post_init__, matching init_setup).
    disast_p_in: float = field(init=False)     # log-mean input (CSV 'disast_p')
    disast_rho: float = field(init=False)       # CSV 'p_rho'
    disast_std_in: float = field(init=False)    # CSV 'p_std', floored at sqrt_eps
    disast_p: float = field(init=False)         # lognormal mean
    disast_std: float = field(init=False)       # lognormal SD
    sig_dis: float = field(init=False)          # AR(1) innovation SD of log p
    sig_z: float = field(init=False)
    zf_std: float = field(init=False)
    rho_zf: float = field(init=False)
    sig_zf: float = field(init=False)
    omg_std: float = field(init=False)
    omg_rho: float = field(init=False)
    sig_omg: float = field(init=False)
    corr_omg_dis: float = field(init=False)
    mus_dims: List[int] = field(init=False)     # length-9 per-dim Smolyak levels
    n_quad_shocks: List[int] = field(init=False)  # [z, zf, p, omg] node counts

    def __post_init__(self):
        r = self.raw
        # Disaster process (mod_param.f90:196-201,246).
        self.disast_p_in = r["disast_p"]
        self.disast_shock = r["disast_shock"]
        self.disast_rho = r["p_rho"]
        self.disast_std_in = max(SQRT_EPS, r["p_std"])
        self.disast_p = math.exp(self.disast_p_in + self.disast_std_in ** 2 / 2)
        self.disast_std = math.sqrt(
            (math.exp(self.disast_std_in ** 2) - 1.0)
            * math.exp(2 * self.disast_p_in + self.disast_std_in ** 2)
        )
        self.sig_dis = self.disast_std_in * math.sqrt(1.0 - self.disast_rho ** 2)

        # Other shocks (mod_param.f90:191-194,204-209,247-248).
        self.sig_z = max(SQRT_EPS, r["sig_z"])
        self.zf_std = max(SQRT_EPS, r["std_zf"])
        self.rho_zf = r["rho_zf"]
        self.sig_zf = self.zf_std * math.sqrt(1.0 - self.rho_zf ** 2)
        self.omg_std = r["omg_std"]
        self.omg_rho = r["omg_rho"]
        self.sig_omg = self.omg_std * math.sqrt(1.0 - self.omg_rho ** 2)
        self.corr_omg_dis = r["corr_omg_dis"]
        # Special case (mod_param.f90:250-252): perfectly-correlated, zero-vol omg.
        if self.omg_std == 0.0 and self.corr_omg_dis == 1.0:
            self.sig_omg = 1.0

        # Per-dimension Smolyak levels (mod_param.f90:224-233,324-327).
        mus = [int(round(r[name])) for name in (
            "mu_k", "mu_tht", "mu_zf", "mu_wh", "mu_wf", "mu_p",
            "mu_ih", "mu_if", "mu_omg",
        )]
        if r["rho_i"] > SQRT_EPS:  # interest smoothing activates the lagged-rate dims
            mus[IDX_IH] = 3
            mus[IDX_IF] = 3
        self.mus_dims = mus

        # Per-shock Gauss-Hermite node counts (mod_param.f90:235-238,259-261).
        n_quad = [int(round(r[name])) for name in
                  ("n_quad_z", "n_quad_zf", "n_quad_p", "n_quad_omg")]
        if self.corr_omg_dis > 1.0 - SQRT_EPS:
            n_quad[3] = 1  # omg collapses onto the disaster node
        self.n_quad_shocks = n_quad

    # Convenient typed access to a few hot raw params.
    @property
    def zeta(self) -> float:
        return self.raw["zeta"]

    @property
    def n_active_dims(self) -> int:
        return sum(1 for m in self.mus_dims if m > 0)

    def __getitem__(self, name: str) -> float:
        return self.raw[name]


def load_param_file(idx: int, params_dir: Path = None) -> Params:
    """Load calibration `idx` (1..9). param_file_1 is the benchmark."""
    if params_dir is None:
        params_dir = Path(raw) / "klreplication" / "params"
    path = Path(params_dir) / f"param_file_{idx}.csv"
    with open(path) as fh:
        header = [h.strip() for h in fh.readline().split(",")]
        values = [float(v) for v in fh.readline().split(",")]
    if len(values) != len(PARAM_NAMES):
        raise ValueError(
            f"{path}: expected {len(PARAM_NAMES)} values, got {len(values)}"
        )
    # Trust positional order (matches the Fortran read), but verify the header.
    if header != PARAM_NAMES:
        mism = [(i, h, n) for i, (h, n) in enumerate(zip(header, PARAM_NAMES))
                if h != n]
        raise ValueError(f"{path}: header mismatch at {mism[:3]}")
    return Params(raw=dict(zip(PARAM_NAMES, values)))
