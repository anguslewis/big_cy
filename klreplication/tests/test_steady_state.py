"""Validation ladder #4 — deterministic steady state (calc_steady).

No Fortran to compare against, so checks internal consistency and the appendix SS
relations: returns equal the inverse average discount factor; the value-function
normalization yields v_ss == 1 exactly; labor hits its target; capital/wages/
consumption are positive with the right home-bias ordering; grid centers follow
the adj/dev parameters; and the result is deterministic.
"""
import math

import pytest

from klrep.params import load_param_file
from klrep.solve.steady_state import calc_steady


@pytest.fixture(scope="module")
def bm():
    return load_param_file(1)


@pytest.fixture(scope="module")
def ss(bm):
    return calc_steady(bm)


def test_returns_equal_inverse_beta(bm, ss):
    tht = [bm["tht_trgt_h"], 1 - bm["tht_trgt_h"]]
    beta = [bm["bbeta_h"], bm["bbeta_f"]]
    beta_avg = sum(t * b for t, b in zip(tht, beta))
    assert ss.rk_ss == pytest.approx(1.0 / beta_avg - 1.0)
    assert ss.rf_ss == ss.rk_ss
    assert 0.0 < ss.rk_ss < 0.05  # plausible quarterly real rate


def test_value_normalization_gives_unit_v(ss):
    # v_normalization_vec is defined precisely so v_ss == 1.
    assert ss.v_ss[0] == pytest.approx(1.0, abs=1e-9)
    assert ss.v_ss[1] == pytest.approx(1.0, abs=1e-9)


def test_labor_at_target(bm, ss):
    assert ss.l_ss == [bm["l_target"], bm["l_target"]]


def test_positivity_and_home_bias(ss):
    assert ss.k_ss > 0 and ss.q_ss > 0 and ss.s_ss > 0
    assert all(k > 0 for k in ss.kappa_ss)
    assert all(c > 0 for c in ss.c_h_ss + ss.c_f_ss)
    assert all(x > 0 for x in ss.chi0_vec)
    # Home consumes more of the Home good; Foreign more of the Foreign good.
    assert ss.c_h_ss[0] > ss.c_f_ss[0]
    assert ss.c_f_ss[1] > ss.c_h_ss[1]
    # Foreign is larger (zeta = 1.6) -> more Foreign capital.
    assert ss.kappa_ss[1] > ss.kappa_ss[0]


def test_grid_centers_follow_params(bm, ss):
    assert ss.k_grid_mean == pytest.approx(bm["k_grid_adj"] * ss.k_ss)
    assert ss.wh_grid_mean == pytest.approx(bm["w_grid_adj"] * ss.w_ss[0])
    assert ss.k_grid_dev == pytest.approx(bm["k_dev_param"] * ss.k_grid_mean)
    assert ss.wh_grid_dev == pytest.approx(bm["w_dev_param"] * ss.wh_grid_mean)


def test_s_ss_converged_near_one(ss):
    # Terms of trade close to 1 in this near-symmetric calibration.
    assert abs(ss.s_ss - 1.0) < 0.1


def test_deterministic(bm):
    a, b = calc_steady(bm), calc_steady(bm)
    assert a.k_ss == b.k_ss and a.chi0_vec == b.chi0_vec and a.v_ss == b.v_ss


def test_all_calibrations_solve_steady_state():
    for i in range(1, 10):
        s = calc_steady(load_param_file(i))
        assert s.k_ss > 0 and all(math.isfinite(v) for v in s.v_ss)
