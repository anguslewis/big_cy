"""Sherlock/local entry — load a converged KL solution, simulate, and print the
Table-2 moment comparison to stdout (so it returns via the synced .out log; no
data-path allowlist needed).

Invoked by batch_controller.sh as:
    python -u klreplication/run_klrep_moments.py <SLURM_JOB_ID> <SLURM_ARRAY_TASK_ID>
The ARRAY TASK ID selects the spec (1 = benchmark). It loads
$big_cy/data/output/klreplication/solution_spec<N>.pt (written by run_klrep.py),
rebuilds the deterministic solver constants from param_file_N, assembles a
SolverState with the loaded converged policies, then runs the simulation +
moments pipeline.

Env knobs:
    KLREP_SPEC      override spec (else argv[2], else 1)
    KLREP_DEVICE    torch device (cpu/cuda); else auto
    KLREP_NBURN     burn-in length (default 10000)
    KLREP_NSIMS     ensemble paths (default 100)
    KLREP_NPER      periods per path (default 400)
"""
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE / "src"))
work_dir = os.environ.get("work_dir", str(REPO.parent))
sys.path.insert(0, str(REPO))
from project_strings import output                  # noqa: E402

import torch                                         # noqa: E402
from klrep.params import load_param_file             # noqa: E402
from klrep.config import get_device                  # noqa: E402
from klrep.solve.steady_state import calc_steady     # noqa: E402
from klrep.setup.smolyak_setup import build_smolyak  # noqa: E402
from klrep.setup.shock_grid import build_shock_grid  # noqa: E402
from klrep.setup.state_grid import build_state_grid  # noqa: E402
from klrep.solve.model_const import build_model_const, big_weight  # noqa: E402
from klrep.solve.solve_model import SolverState      # noqa: E402
from klrep.simulate.simulate import (build_sim_coeffs, stochastic_ss,  # noqa: E402
                                     burn_in, simulate_ensemble)
from klrep.post.table2_series import build_table2_series  # noqa: E402
from klrep.post.moments import (compute_table2_moments,  # noqa: E402
                                compute_table2_moments_per_sim, format_table2)


def rebuild_state(spec, device):
    """Rebuild a SolverState from param_file_spec + the saved solution_spec.pt."""
    params_dir = HERE / "inputs" / "params"
    p = load_param_file(spec, params_dir=params_dir)
    ss = calc_steady(p)
    smol = build_smolyak(p, device=device)
    sg = build_shock_grid(p, device=device)
    stg = build_state_grid(p, smol, ss, sg, device=device)
    const = build_model_const(p, ss, sg, stg, smol, device=device)

    sol_path = Path(output) / "klreplication" / f"solution_spec{spec}.pt"
    sol = torch.load(sol_path, map_location=device, weights_only=False)
    S = const.state_grid.shape[0]
    if sol["n_states"] != S:
        raise ValueError(f"grid mismatch: solution has {sol['n_states']} states, "
                         f"rebuilt const has {S} (mus_dims must match the solve).")
    # move loaded guesses/tensors onto the target device
    g = sol["g"]
    for fld in ("c_spending", "s", "l_aggr", "q", "infl", "nom_i", "share", "bF_share"):
        setattr(g, fld, getattr(g, fld).to(device))
    st = SolverState(
        const=const, g=g,
        v_mat=sol["v_mat"].to(device), mc_mat=sol["mc_mat"].to(device),
        next_state=sol["next_state"].to(device), big_w=big_weight(const),
        k_next_mat=const.state_grid.new_zeros(0), w_choice=const.state_grid.new_zeros(0),
    )
    return p, st, sol


def main():
    spec = int(os.environ.get("KLREP_SPEC", sys.argv[2] if len(sys.argv) > 2 else 1))
    job_id = sys.argv[1] if len(sys.argv) > 1 else "local"
    n_burn = int(os.environ.get("KLREP_NBURN", 10000))
    n_sims = int(os.environ.get("KLREP_NSIMS", 100))
    n_per = int(os.environ.get("KLREP_NPER", 400))
    device = get_device()
    print(f"=== run_klrep_moments: spec={spec} job={job_id} device={device} "
          f"n_burn={n_burn} n_sims={n_sims} n_per={n_per} ===", flush=True)

    p, st, sol = rebuild_state(spec, device)
    print(f"loaded solution_spec{spec}.pt: diff={sol['diff']:.3e} iters={sol['iters']} "
          f"n_states={sol['n_states']}", flush=True)

    t0 = time.time()
    sc = build_sim_coeffs(st)
    ss_state, econ = stochastic_ss(sc)
    print(f"stochastic-ss econ state: {[round(float(x), 5) for x in econ]}", flush=True)
    pool = burn_in(sc, ss_state, n_burn=n_burn)
    sim = simulate_ensemble(sc, pool, n_sims=n_sims, n_periods=n_per, disaster=False)
    print(f"simulation done: {time.time() - t0:.1f}s | "
          f"state_series {tuple(sim['state_series'].shape)}", flush=True)

    series = build_table2_series(st.const, st.g, sim, disast_shock=p["disast_shock"])
    per_sim = compute_table2_moments_per_sim(series, bg_yss=p["bg_yss"])
    moms = {k: float(v.mean()) for k, v in per_sim.items()}
    print(format_table2(moms, per_sim=per_sim), flush=True)
    print(f"\ntotal {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
