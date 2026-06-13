"""Smoke test — the outer solver runs a few iterations on a COARSE grid without
crashing and produces finite values. (Convergence/accuracy gating is separate.)"""
import math

import torch

from klrep.params import load_param_file
from klrep.solve.solve_model import init_solver, solve_step


def _coarse_params():
    p = load_param_file(1)
    # Reduce Smolyak levels to a tiny grid (level 1 -> 3 nodes/active dim).
    p.mus_dims = [1, 1, 1, 1, 1, 1, 0, 0, 1]
    return p


def test_init_solver_shapes_finite():
    p = _coarse_params()
    st = init_solver(p)
    S = st.const.state_grid.shape[0]
    Q = st.const.n_quad
    assert st.v_mat.shape == (S, 2)
    assert st.next_state.shape[0] == S and st.next_state.shape[1] == Q
    for t in [st.v_mat, st.mc_mat, st.g.q, st.g.s, st.g.c_spending, st.next_state,
              st.big_w, st.k_next_mat]:
        assert torch.isfinite(t).all()
    # big_w rows sum to 1.
    assert torch.allclose(st.big_w.sum(dim=1), torch.ones(S, dtype=torch.float64), atol=1e-10)


def test_outer_iterations_run_finite():
    p = _coarse_params()
    st = init_solver(p)
    diffs = []
    for it in range(1, 6):  # a few early (pre-bond-clearing) iterations
        st, diff = solve_step(st, it)
        diffs.append(diff)
        assert math.isfinite(diff), f"diff not finite at iter {it}"
        for t in [st.v_mat, st.mc_mat, st.g.q, st.g.s, st.g.c_spending, st.g.nom_i,
                  st.next_state, st.k_next_mat]:
            assert torch.isfinite(t).all(), f"non-finite state at iter {it}"
    # Values stay economically sane (positive q, s, consumption).
    assert (st.g.q > 0).all() and (st.g.s > 0).all()
    assert (st.g.c_spending > 0).all()
