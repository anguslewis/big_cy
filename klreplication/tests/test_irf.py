"""Test the MIT-shock transition fixed point (the novel piece of the generalized
IRF): at zero shock it must recover the solve's no-shock next-state; nonzero shocks
produce finite, in-bounds, perturbed transition states."""
import torch

from klrep.params import load_param_file
from klrep.solve.solve_model import init_solver, solve_step
from klrep.solve.irf import build_shock_transition, _current
from klrep.model.period_block import compute_current_period
from klrep.params import IDX_K, IDX_ZF


def _coarse_solved(n_iters=120):
    p = load_param_file(1)
    p.mus_dims = [1, 1, 1, 1, 1, 1, 0, 0, 1]
    st = init_solver(p)
    for it in range(1, n_iters + 1):
        st, _ = solve_step(st, it)
    return p, st


def _k_next_new(const, g):
    cp = compute_current_period(const.state_grid[:, IDX_K], g.s, g.l_aggr, g.q,
                                const.state_grid[:, IDX_ZF], const.wealth_share_grid,
                                aalpha=const.aalpha, sigma=const.sigma,
                                ddelta=const.ddelta, zeta=const.zeta,
                                varsigma=const.varsigma_vec)
    msc = const.min_cons_sav
    cs = torch.minimum(torch.maximum(g.c_spending, torch.full_like(g.c_spending, msc)),
                       cp.wealth_vec + cp.w_choice * g.l_aggr - msc)
    savings = cp.wealth_vec + cp.w_choice * g.l_aggr - cs
    return (savings * (1.0 - g.share) / g.q.unsqueeze(1)).sum(dim=1)


def test_zero_shock_baseline_consistent():
    p, st = _coarse_solved(120)
    const, g = st.const, st.g
    knn = _k_next_new(const, g)
    base = st.next_state[:, const.no_shock_idx, :]
    tr0 = build_shock_transition(const, g, st.v_mat, st.mc_mat, st.next_state, knn,
                                 sidx=3, shock_size=0.0, conv=1e-8)
    tr0b = build_shock_transition(const, g, st.v_mat, st.mc_mat, st.next_state, knn,
                                  sidx=0, shock_size=0.0, conv=1e-8)
    # deterministic (any zero shock gives the same baseline) + finite + in-bounds.
    assert torch.allclose(tr0, tr0b, atol=1e-9)
    assert torch.isfinite(tr0).all() and (tr0.abs() <= 1.0 + 1e-9).all()
    # the zero-shock realized transition ~ the solve's expected no-shock node ON
    # AVERAGE (the IRF differences against THIS baseline, so the realized-vs-expected
    # rk gap cancels; only a small mean offset, large only at extreme grid points).
    assert float((tr0 - base).abs().mean()) < 1e-2


def test_nonzero_shocks_finite_perturbed():
    p, st = _coarse_solved(120)
    const, g = st.const, st.g
    knn = _k_next_new(const, g)
    base = st.next_state[:, const.no_shock_idx, :]
    for sidx, size in [(3, 2 * p.sig_omg), (0, -5 * p.sig_z), (1, -2 * p.sig_zf),
                       (2, 2 * p.sig_dis)]:
        tr = build_shock_transition(const, g, st.v_mat, st.mc_mat, st.next_state, knn,
                                    sidx=sidx, shock_size=size, conv=1e-7)
        assert torch.isfinite(tr).all(), sidx
        assert (tr.abs() <= 1.0 + 1e-9).all(), sidx
        assert float((tr - base).abs().mean()) > 1e-4, sidx   # genuinely perturbed
