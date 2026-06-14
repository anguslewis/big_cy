"""Unit tests for the Table-2 post-processing: compute_results_vec columns,
build_table2_series derived series, and the moment formulas (synthetic series)."""
import math

import torch

from klrep.params import load_param_file
from klrep.solve.solve_model import init_solver
from klrep.model.compute_results_vec import compute_results_vec, RESULT_NAMES, LEVEL_NAMES
from klrep.post.table2_series import build_table2_series, _prepend_first
from klrep.post.moments import compute_table2_moments, DELIVERABLE


def _coarse_params():
    p = load_param_file(1)
    p.mus_dims = [1, 1, 1, 1, 1, 1, 0, 0, 1]
    return p


def test_compute_results_vec_columns_finite():
    p = _coarse_params()
    st = init_solver(p)
    rg = compute_results_vec(st.const, st.g)
    S = st.const.state_grid.shape[0]
    assert rg.mat.shape == (S, len(RESULT_NAMES))
    assert torch.isfinite(rg.mat).all()
    assert rg.names == RESULT_NAMES
    assert rg.level_mask.sum().item() == len(LEVEL_NAMES)
    # economically sane sign checks on the initial grid
    j = {n: i for i, n in enumerate(RESULT_NAMES)}
    assert (rg.mat[:, j["yh"]] > 0).all()
    assert (rg.mat[:, j["q"]] > 0).all()
    assert (rg.mat[:, j["s"]] > 0).all()


def test_build_table2_series_shapes_derived():
    p = _coarse_params()
    st = init_solver(p)
    n_sims, T = 4, 12
    d = len(st.const.active_dims)
    # synthetic standardized path inside [-1,1]; tiny z-shocks
    std = torch.zeros((n_sims, T, d), dtype=torch.float64)
    z_shock = torch.zeros((n_sims, T), dtype=torch.float64)
    z_shock[:, 1:] = 0.001
    sim = {"std_series": std, "z_shock_series": z_shock}
    s = build_table2_series(st.const, st.g, sim, disast_shock=p["disast_shock"])
    for name in ("yh", "yf", "rfh", "nfa", "nfa_rel", "dis", "infl_h", "lh"):
        assert s[name].shape == (n_sims, T), name
        assert torch.isfinite(s[name]).all(), name
    # de-trending applied: with constant std path, level column grows by exp(cumsum z)
    z_series = torch.exp(torch.cumsum(z_shock, dim=1))
    assert torch.allclose(s["yh"][:, -1] / s["yh"][:, 0], z_series[:, -1] / z_series[:, 0])


def test_prepend_first():
    x = torch.tensor([[1.0, 2.0, 3.0]])
    assert torch.equal(_prepend_first(x), torch.tensor([[1.0, 1.0, 2.0, 3.0]]))


def test_moment_formulas_synthetic():
    # Deterministic synthetic series; check a few moments against hand computation.
    n_sims, T = 1, 6
    ones = torch.ones((n_sims, T), dtype=torch.float64)
    yf = torch.tensor([[3.0, 3.1, 2.9, 3.05, 2.95, 3.02]], dtype=torch.float64)
    s = {
        "yh": ones * 2.0, "yf": yf, "s": ones * 1.0,
        "ch": torch.tensor([[1.0, 1.01, 1.0, 1.01, 1.0, 1.01]], dtype=torch.float64),
        "cf": ones * 4.0, "inv": ones * 1.0,
        "lh": ones * 0.9, "lf": ones * 1.1,
        "infl_h": ones * 1.02, "infl_f": ones * 0.98,
        "qx": ones * 1.0,
        "rfh": ones * 0.005,
        "nfa_rel": ones * (-0.8),
    }
    moms = compute_table2_moments(s, bg_yss=0.127)
    assert set(moms.keys()) == set(DELIVERABLE)
    # mom1 = mean(yf/(s*yh))
    assert math.isclose(moms[1], float((yf / (1.0 * 2.0)).mean()), rel_tol=1e-9)
    # mom6 = 4*100*mean(rfh) = 4*100*0.005 = 2.0
    assert math.isclose(moms[6], 2.0, rel_tol=1e-9)
    # mom7 = 100*mean(nfa_rel/4) = 100*(-0.8/4) = -20.0
    assert math.isclose(moms[7], -20.0, rel_tol=1e-9)
    # mom12 = mean(lh) = 0.9 ; mom13 = mean(lf) = 1.1
    assert math.isclose(moms[12], 0.9, rel_tol=1e-9)
    assert math.isclose(moms[13], 1.1, rel_tol=1e-9)
    # mom14 = 100*mean(log(infl_h)) = 100*log(1.02)
    assert math.isclose(moms[14], 100 * math.log(1.02), rel_tol=1e-9)
    # mom11 = 100*bg_yss*mean((cf/qx)/(4*yh)) = 100*0.127*(4/(4*2))
    assert math.isclose(moms[11], 100 * 0.127 * (4.0 / (4.0 * 2.0)), rel_tol=1e-9)
    # all finite
    assert all(math.isfinite(v) for v in moms.values())
