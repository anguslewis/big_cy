"""Validation — calibration loader vs the shipped param_file_1 (benchmark).

Reads the real copied input (…/big_cy/data/raw/klreplication/params/) and checks
raw values, the derived shock scalars (disaster lognormal moments, AR(1)
innovation SDs), and the per-dim Smolyak levels / quadrature counts against the
known benchmark calibration.
"""
import math

import pytest

from klrep.params import load_param_file


@pytest.fixture(scope="module")
def bm():
    return load_param_file(1)


def test_raw_benchmark_values(bm):
    assert bm["zeta"] == pytest.approx(1.6)
    assert bm["gma_h"] == pytest.approx(21.0)
    assert bm["gma_f"] == pytest.approx(23.8)
    assert bm["ies_h"] == pytest.approx(0.75)
    assert bm["ies_f"] == pytest.approx(0.75)
    assert bm["bbeta_h"] == pytest.approx(0.988708, abs=1e-6)
    assert bm["sigma"] == pytest.approx(1.5)
    assert bm["aalpha"] == pytest.approx(0.33)
    assert bm["delta"] == pytest.approx(0.025)
    assert bm["chi"] == pytest.approx(0.75)


def test_derived_disaster_moments(bm):
    # disast_p (lognormal mean) = exp(p_in + std_in^2 / 2).
    p_in = bm["disast_p"]            # CSV 'disast_p' is the log-mean input
    assert bm.disast_p_in == pytest.approx(-5.644891, abs=1e-6)
    assert bm.disast_std_in == pytest.approx(0.832555, abs=1e-6)
    assert bm.disast_p == pytest.approx(math.exp(p_in + bm.disast_std_in ** 2 / 2))
    assert bm.disast_p > 0
    # AR(1) innovation SD of log p.
    assert bm.sig_dis == pytest.approx(
        bm.disast_std_in * math.sqrt(1 - bm.disast_rho ** 2)
    )


def test_derived_shock_sds(bm):
    assert bm.sig_zf == pytest.approx(bm.zf_std * math.sqrt(1 - bm.rho_zf ** 2))
    assert bm.sig_omg == pytest.approx(bm.omg_std * math.sqrt(1 - bm.omg_rho ** 2))
    assert bm.sig_z == pytest.approx(bm["sig_z"])


def test_smolyak_levels_and_active_dims(bm):
    # Benchmark: 7 active dims (k, tht, zf, wh, wf, p, omg); ih/if inactive
    # because rho_i == 0.
    assert bm.mus_dims == [3, 3, 3, 3, 3, 3, 0, 0, 3]
    assert bm.n_active_dims == 7
    assert bm["rho_i"] == 0.0


def test_quadrature_counts(bm):
    # [z, zf, p, omg]; corr_omg_dis = 0.5 (< 1) so omg is not collapsed.
    assert bm.n_quad_shocks == [3, 3, 5, 5]
    assert bm.corr_omg_dis == pytest.approx(0.5)


def test_all_nine_calibrations_load():
    for i in range(1, 10):
        p = load_param_file(i)
        assert len(p.raw) == 65
        assert 1 <= p.n_active_dims <= 9
