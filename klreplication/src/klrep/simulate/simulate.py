"""Forward simulation of the solved model — port of the simulation phases of
`mod_results.f90 :: create_results` (stochastic steady state, burn-in, and the
moment ensemble).

The solved policies live on the Smolyak grid in standardized [-1,1] coords. We fit
Smolyak coefficients to (a) the per-quad-node next-state transition, (b) the 9
economic state variables, and (c) the policy variables we need downstream, then
propagate a CONTINUOUS standardized state: each period draw a quadrature node from
the (state-dependent) weight distribution and apply that node's transition.

Matches the Fortran exactly except the RNG: NAG's stream isn't reproducible here,
so simulated-moment agreement is statistical, not bit-exact (n_sims×T averages
out). Tensor-native: all sims advance together; transitions are a batched matmul.

Phases (digest fortran_mod_results.md §2):
- stochastic steady state: iterate the no-shock transition to a fixed point.
- burn-in: 2·n_burn periods, disaster ON (big_weight_vec); keep last n_burn states.
- ensemble: n_sims starts sampled from the burn-in pool; n_sim_periods each.
  no-disaster run excludes the disaster node and clips next-state to [-1,1];
  with-disaster run includes it and does NOT clip.
"""
from dataclasses import dataclass

import torch

from ..config import DTYPE
from ..params import IDX_DIS
from ..grid.smolyak_polynomial import smolyak_polynomial
from ..grid.interp_fit import interp_factor, interp_solve

# policy columns carried along the path (order is the POL_* index below)
POL_NAMES = ["s", "q", "l_h", "l_f", "c_h", "c_f", "nom_ih", "spread",
             "infl_h", "infl_f", "share_h", "share_f", "bF_h", "bF_f"]


@dataclass
class SimCoeffs:
    trans: torch.Tensor      # (Q, n_states, d) per-node next-state coeffs
    state: torch.Tensor      # (n_states, 9) economic-state coeffs
    pol: torch.Tensor        # (n_states, n_pol) policy coeffs
    elem: object
    mu: int
    d: int
    n_quad: int
    no_shock_idx: int
    quad_weight_vec: torch.Tensor  # (n_gh+1,)
    shock_grid: torch.Tensor       # (n_quad, 4)


def build_sim_coeffs(st) -> SimCoeffs:
    """Fit Smolyak coeffs to the solved transition, state grid, and policies."""
    const = st.const
    B = const.smol_polynom
    S = const.state_grid.shape[0]
    Q = const.n_quad
    d = len(const.active_dims)
    lu = interp_factor(B)

    # transition: next_state (S, Q, d) -> coeffs (S, Q, d) -> (Q, S, d)
    rhs = st.next_state.reshape(S, Q * d)
    trans = interp_solve(lu, rhs).reshape(S, Q, d).permute(1, 0, 2).contiguous()

    state_coeffs = interp_solve(lu, const.state_grid)               # (S, 9)

    g = st.g
    pol = torch.stack([
        g.s, g.q, g.l_aggr[:, 0], g.l_aggr[:, 1], g.c_spending[:, 0], g.c_spending[:, 1],
        g.nom_i[:, 0], g.nom_i[:, 1], g.infl[:, 0], g.infl[:, 1],
        g.share[:, 0], g.share[:, 1], g.bF_share[:, 0], g.bF_share[:, 1],
    ], dim=1)                                                       # (S, 14)
    pol_coeffs = interp_solve(lu, pol)

    return SimCoeffs(trans=trans, state=state_coeffs, pol=pol_coeffs,
                     elem=const.elem, mu=const.mu, d=d, n_quad=Q,
                     no_shock_idx=const.no_shock_idx,
                     quad_weight_vec=const.quad_weight_vec,
                     shock_grid=st.const.dz_vec.new_zeros(0) if False else _shock_grid(const))


def _shock_grid(const):
    # Reconstruct the (n_quad,4) shock grid from dz_vec + zeros (only z col is
    # stored on const; the other cols aren't needed for Table-2 realized vars, but
    # we keep the z column which drives de-trending). Full grid lives in ShockGrid;
    # callers needing all 4 cols pass it in. Here we expose dz (col 0).
    g = torch.zeros((const.n_quad, 4), dtype=DTYPE, device=const.dz_vec.device)
    g[:, 0] = const.dz_vec
    return g


def _basis(state, sc: SimCoeffs):
    return smolyak_polynomial(state, sc.elem, sc.mu)               # (N, n_states)


def stochastic_ss(sc: SimCoeffs, *, tol=1.49e-8, max_iter=100000):
    """Iterate the no-shock transition from the first grid corner to a fixed point.
    Returns (state_std (1,d), econ_state (9,))."""
    state = torch.zeros((1, sc.d), dtype=DTYPE, device=sc.state.device)  # grid center
    for _ in range(max_iter):
        b = _basis(state, sc)
        nxt = b @ sc.trans[sc.no_shock_idx]
        if float((state - nxt).abs().max()) < tol:
            state = nxt
            break
        state = nxt
    econ = (_basis(state, sc) @ sc.state).squeeze(0)
    return state, econ


def _draw_nodes(weights_cum, u):
    """Given per-sim cumulative weights (N, n_quad) and uniforms u (N,), return the
    first node index whose cumulative weight exceeds u (matches the Fortran walk).
    The uniforms are drawn on a CPU generator (reproducible across devices); move
    them onto the weights' device before comparing (the solve runs on GPU)."""
    u = u.to(weights_cum.device)
    q = (weights_cum <= u.unsqueeze(1)).sum(dim=1)
    return q.clamp(max=weights_cum.shape[1] - 1)


def _advance(state, sc, q_nxt):
    """next_state for each sim using its drawn node q_nxt (N,). Computes all-node
    transitions then gathers. Returns (N, d)."""
    b = _basis(state, sc)                                          # (N, S)
    nxt_all = torch.einsum("ns,qsd->nqd", b, sc.trans)            # (N, Q, d)
    idx = q_nxt.view(-1, 1, 1).expand(-1, 1, sc.d)
    return nxt_all.gather(1, idx).squeeze(1)                       # (N, d)


def _weight_matrix(econ_states, sc, *, disaster):
    """Per-sim quadrature weights (N, n_quad). no-disaster: product nodes only
    (disaster + no-shock get 0 mass). with-disaster: big_weight_vec(p_dis)."""
    N = econ_states.shape[0]
    Q = sc.n_quad
    n_gh = sc.quad_weight_vec.shape[0] - 1   # product-node count
    w = torch.zeros((N, Q), dtype=DTYPE, device=econ_states.device)
    if disaster:
        p_dis = torch.exp(econ_states[:, IDX_DIS])                # (N,)
        w[:, : Q - 1] = sc.quad_weight_vec.unsqueeze(0) * (1.0 - p_dis).unsqueeze(1)
        w[:, Q - 1] = p_dis
    else:
        # only product nodes 0..n_gh-1 carry mass (renormalized to 1 already)
        w[:, :n_gh] = sc.quad_weight_vec[:n_gh].unsqueeze(0)
    return w


def burn_in(sc: SimCoeffs, ss_state, *, n_burn=10000, seed=712):
    """2·n_burn periods with disaster ON; return the last n_burn standardized
    states (n_burn, d) as the ergodic pool."""
    g = torch.Generator(device="cpu").manual_seed(seed)
    n_periods = 2 * n_burn
    state = ss_state.clone()                                       # (1, d)
    pool = torch.empty((n_burn, sc.d), dtype=DTYPE, device=sc.state.device)
    for t in range(n_periods):
        econ = _basis(state, sc) @ sc.state                       # (1, 9)
        w = _weight_matrix(econ, sc, disaster=True)               # (1, Q)
        u = torch.rand(1, generator=g, dtype=DTYPE)
        q = _draw_nodes(w.cumsum(dim=1), u)
        if t >= n_burn:
            pool[t - n_burn] = state.squeeze(0)
        state = _advance(state, sc, q).clamp(-1.0, 1.0)
    return pool


def simulate_ensemble(sc: SimCoeffs, pool, *, n_sims=100, n_periods=400,
                      disaster=False, seed=99):
    """Simulate n_sims paths of length n_periods from random burn-in starts.
    Returns dict with state_series (n_sims,T,9), pol_series (n_sims,T,14),
    z_shock_series (n_sims,T) [the realized z-shock, for de-trending]."""
    g = torch.Generator(device="cpu").manual_seed(seed)
    n_burn = pool.shape[0]
    starts = torch.randperm(n_burn, generator=g)[:n_sims]        # CPU index tensor
    state = pool[starts.to(pool.device)].clone()                 # (n_sims, d)

    state_series = torch.empty((n_sims, n_periods, 9), dtype=DTYPE, device=sc.state.device)
    pol_series = torch.empty((n_sims, n_periods, sc.pol.shape[1]), dtype=DTYPE, device=sc.state.device)
    std_series = torch.empty((n_sims, n_periods, sc.d), dtype=DTYPE, device=sc.state.device)
    z_shock = torch.zeros((n_sims, n_periods), dtype=DTYPE, device=sc.state.device)

    u_all = torch.rand((n_periods, n_sims), generator=g, dtype=DTYPE)
    for t in range(n_periods):
        b = _basis(state, sc)
        econ = b @ sc.state                                       # (n_sims, 9)
        pol = b @ sc.pol                                          # (n_sims, n_pol)
        state_series[:, t, :] = econ
        pol_series[:, t, :] = pol
        std_series[:, t, :] = state                               # standardized state (for grid interp)
        w = _weight_matrix(econ, sc, disaster=disaster)
        q = _draw_nodes(w.cumsum(dim=1), u_all[t])
        if t < n_periods - 1:
            z_shock[:, t + 1] = sc.shock_grid[q, 0]               # realized z-shock
        state = _advance(state, sc, q)
        if not disaster:
            state = state.clamp(-1.0, 1.0)
    return {"state_series": state_series, "pol_series": pol_series,
            "std_series": std_series, "z_shock_series": z_shock}
