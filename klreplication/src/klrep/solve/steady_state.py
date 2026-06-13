"""Deterministic steady state — port of `mod_calc.f90 :: calc_steady`.

Solves the no-shock, no-disaster steady state to (a) center the Smolyak grid for
capital and wages and (b) calibrate the per-country labor-disutility level
`chi0_vec` and value-function normalization. Two nested scalar fixed points:
the terms of trade `s_ss` (outer, damped 0.01) and, per country, `chi0` (damped
0.5). One-time scalar computation (n_I=2), so plain Python floats — not a hot
path. Returns a SteadyState carrying everything calc_sol/grid_setup consume.

Mirrors the Fortran block at mod_calc.f90:28-200 line-for-line; also writes the
`extra_data.csv` row (chi0_h, chi0_f, sig_dis, sig_omg, dis_grid_mean,
omg_grid_mean, y_h_ss, aggr_css) that the MATLAB post-processing reads.
"""
import math
from dataclasses import dataclass
from typing import List

from ..params import Params, SQRT_EPS
from ..model.calc_bundle import calc_bundle


@dataclass
class SteadyState:
    # Core SS objects (per-country lists are [Home, Foreign]).
    s_ss: float
    k_ss: float
    q_ss: float
    rf_ss: float
    rk_ss: float
    pi_ss: float
    kappa_ss: List[float]      # [kappa_h_ss, kappa_f_ss]
    l_ss: List[float]
    w_ss: List[float]
    c_h_ss: List[float]
    c_f_ss: List[float]
    c_cost_ss: List[float]
    c_bundle_ss: List[float]
    chi0_vec: List[float]
    v_normalization_vec: List[float]
    v_ss: List[float]
    mc_ss: List[float]
    q_l_ss: List[float]
    tot_wealth_ss: List[float]
    y_h_ss: float
    y_f_ss: float
    aggr_css: float
    # Grid centers / half-widths set here (consumed by grid_setup).
    k_grid_mean: float
    wh_grid_mean: float
    wf_grid_mean: float
    k_grid_dev: float
    wh_grid_dev: float
    wf_grid_dev: float


def calc_steady(p: Params) -> SteadyState:
    sigma, aalpha, ddelta, chi = p["sigma"], p["aalpha"], p["delta"], p["chi"]
    inv_share_h = p["inv_share_h"]
    zeta = p.zeta
    bbeta = [p["bbeta_h"], p["bbeta_f"]]
    ies = [p["ies_h"], p["ies_f"]]
    varsigma = [p["varsigma_h"], p["varsigma_f"]]
    tht_trgt = [p["tht_trgt_h"], 1.0 - p["tht_trgt_h"]]

    zf_ss = [1.0, zeta]
    bbeta_avrg = sum(t * b for t, b in zip(tht_trgt, bbeta))
    rk_ss = 1.0 / bbeta_avrg - 1.0
    rf_ss = 1.0 / bbeta_avrg - 1.0
    l_ss = [p["l_target"], p["l_target"]]
    theta_ss = list(tht_trgt)

    s_ss = 1.0
    diff, it = 1.0, 0
    while diff > SQRT_EPS:
        it += 1
        q_ss = inv_share_h ** (1.0 / (1.0 - sigma)) * (
            1.0 + (1.0 - inv_share_h) / inv_share_h * s_ss ** (sigma - 1.0)
        ) ** (1.0 / (1.0 - sigma))
        pi_ss = (rk_ss + ddelta) * q_ss
        kappa_h = zf_ss[0] * l_ss[0] * (pi_ss / aalpha) ** (1.0 / (aalpha - 1.0))
        kappa_f = zf_ss[1] * l_ss[1] * (s_ss * pi_ss / aalpha) ** (1.0 / (aalpha - 1.0))
        k_ss = kappa_h + kappa_f
        inv_h = ddelta * k_ss * q_ss / (
            1.0 + (s_ss ** (sigma - 1.0)) * (1.0 - inv_share_h) / inv_share_h
        )
        inv_f = inv_h * s_ss ** sigma * (1.0 - inv_share_h) / inv_share_h
        y_h_ss = kappa_h ** aalpha * (zf_ss[0] * l_ss[0]) ** (1.0 - aalpha)
        y_f_ss = kappa_f ** aalpha * (zf_ss[1] * l_ss[1]) ** (1.0 - aalpha)
        c_f_aggr = (zf_ss[1] * l_ss[1]) ** (1.0 - aalpha) * kappa_f ** aalpha - inv_f

        w_ss = [
            (1.0 - aalpha) * kappa_h ** aalpha * l_ss[0] ** (-aalpha)
            * zf_ss[0] ** (1.0 - aalpha),
            (1.0 - aalpha) * kappa_f ** aalpha * l_ss[1] ** (-aalpha)
            * s_ss ** (-1.0) * zf_ss[1] ** (1.0 - aalpha),
        ]
        c_h_ss = [
            (theta_ss[i] * q_ss * k_ss * rf_ss + w_ss[i] * l_ss[i])
            / (1.0 + (s_ss ** (sigma - 1.0)) * (1.0 - varsigma[i]) / varsigma[i])
            for i in range(2)
        ]
        c_f_ss = [
            c_h_ss[i] * s_ss ** sigma * (1.0 - varsigma[i]) / varsigma[i]
            for i in range(2)
        ]
        denom = sum(c_h_ss[i] * (1.0 - varsigma[i]) / varsigma[i] for i in range(2))
        s_ss_new = (c_f_aggr / denom) ** (1.0 / sigma)
        diff = abs(s_ss_new - s_ss)
        s_ss = s_ss + 0.01 * (s_ss_new - s_ss)
        if it > 10000:
            raise RuntimeError("no convergence in s_ss")

    c_cost_ss = [
        varsigma[i] ** (1.0 / (1.0 - sigma)) * (
            1.0 + (1.0 - varsigma[i]) / varsigma[i] * s_ss ** (sigma - 1.0)
        ) ** (1.0 / (1.0 - sigma))
        for i in range(2)
    ]

    chi0_vec, v_norm, v_ss, mc_ss, c_bundle_ss = [], [], [], [], []
    for i in range(2):
        c_bundle = calc_bundle(c_h_ss[i], c_f_ss[i], varsigma[i], sigma)
        chi0, diff, it = 0.5, 1.0, 0
        while diff > 1e-12:
            it += 1
            chi0_new = (
                ies[i] * w_ss[i] / c_cost_ss[i] / c_bundle * l_ss[i] ** (-1.0)
                * (1.0 + (1.0 / ies[i] - 1.0) * chi0 * chi / (1.0 + chi))
            )
            diff = abs(chi0 - chi0_new)
            chi0 = chi0 + 0.5 * (chi0_new - chi0)
            if it > 1000:
                raise RuntimeError("no convergence in chi0")
        chi0_vec.append(chi0)
        vn = (1.0 - bbeta[i]) * (
            c_bundle ** (1.0 - 1.0 / ies[i])
            * (1.0 + (1.0 / ies[i] - 1.0) * chi0 / (1.0 + 1.0 / chi)) ** (1.0 / ies[i])
        ) ** (-1.0)
        v_norm.append(vn)
        v_ss.append(
            (vn / (1.0 - bbeta[i])) ** (1.0 / (1.0 - 1.0 / ies[i]))
            * c_bundle
            * (1.0 + (1.0 / ies[i] - 1.0) * chi0 / (1.0 + 1.0 / chi))
            ** ((1.0 / ies[i]) / (1.0 - 1.0 / ies[i]))
        )
        mc_ss.append(
            ((1.0 + (1.0 / ies[i] - 1.0) * chi0 / (1.0 + 1.0 / chi)) / c_bundle)
            ** (1.0 / ies[i])
        )
        c_bundle_ss.append(c_bundle)

    q_l_ss = [w_ss[i] / (1.0 - bbeta[i]) for i in range(2)]
    tot_wealth_ss = [q_l_ss[i] + theta_ss[i] * k_ss * (1.0 + rk_ss) for i in range(2)]

    k_grid_mean = p["k_grid_adj"] * k_ss
    wh_grid_mean = p["w_grid_adj"] * w_ss[0]
    wf_grid_mean = p["w_grid_adj"] * w_ss[1]
    k_grid_dev = p["k_dev_param"] * k_grid_mean
    wh_grid_dev = p["w_dev_param"] * wh_grid_mean
    wf_grid_dev = p["w_dev_param"] * wf_grid_mean

    aggr_css = sum(c_cost_ss[i] * c_bundle_ss[i] for i in range(2))

    return SteadyState(
        s_ss=s_ss, k_ss=k_ss, q_ss=q_ss, rf_ss=rf_ss, rk_ss=rk_ss, pi_ss=pi_ss,
        kappa_ss=[kappa_h, kappa_f], l_ss=l_ss, w_ss=w_ss,
        c_h_ss=c_h_ss, c_f_ss=c_f_ss, c_cost_ss=c_cost_ss, c_bundle_ss=c_bundle_ss,
        chi0_vec=chi0_vec, v_normalization_vec=v_norm, v_ss=v_ss, mc_ss=mc_ss,
        q_l_ss=q_l_ss, tot_wealth_ss=tot_wealth_ss, y_h_ss=y_h_ss, y_f_ss=y_f_ss,
        aggr_css=aggr_css,
        k_grid_mean=k_grid_mean, wh_grid_mean=wh_grid_mean, wf_grid_mean=wf_grid_mean,
        k_grid_dev=k_grid_dev, wh_grid_dev=wh_grid_dev, wf_grid_dev=wf_grid_dev,
    )
