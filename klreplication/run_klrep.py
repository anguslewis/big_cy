"""Sherlock entry point — solve one KL calibration on the full grid (GPU) and
validate against the deterministic steady state.

Invoked by batch_controller.sh as:
    python -u klreplication/run_klrep.py <SLURM_JOB_ID> <SLURM_ARRAY_TASK_ID>
The ARRAY TASK ID selects the specification: array index N -> param_file_N
(index 1 = benchmark). Outputs go to $big_cy/data/output/klreplication/.

Local testing knobs (env vars):
    KLREP_SPEC      override the spec (else argv[2], else 1)
    KLREP_MUS       comma list of 9 Smolyak levels (coarse grid), e.g. 1,1,1,1,1,1,0,0,1
    KLREP_MAXITER   max outer iterations (default 5000)
    KLREP_DEVICE    torch device (cpu/cuda); else auto (cuda if available)
"""
import os
import sys
import time
from pathlib import Path

# --- resolve package + project paths -------------------------------------
HERE = Path(__file__).resolve().parent              # .../big_cy/klreplication
REPO = HERE.parent                                  # .../big_cy
sys.path.insert(0, str(HERE / "src"))               # klrep package
work_dir = os.environ.get("work_dir", str(REPO.parent))
sys.path.insert(0, str(REPO))                       # project_strings lives here
from project_strings import output                  # noqa: E402

import torch  # noqa: E402
from klrep.params import load_param_file            # noqa: E402
from klrep.solve.steady_state import calc_steady    # noqa: E402
from klrep.solve.solve_model import solve_model     # noqa: E402
from klrep.config import get_device                 # noqa: E402
from klrep.params import (IDX_K, IDX_THH, IDX_ZF, IDX_OMG)  # noqa: E402


def main():
    spec = int(os.environ.get("KLREP_SPEC", sys.argv[2] if len(sys.argv) > 2 else 1))
    job_id = sys.argv[1] if len(sys.argv) > 1 else "local"
    max_iter = int(os.environ.get("KLREP_MAXITER", 5000))
    device = get_device()
    print(f"=== run_klrep: spec={spec} job={job_id} device={device} max_iter={max_iter} ===",
          flush=True)

    params_dir = HERE / "inputs" / "params"   # tracked in-repo (reaches cluster via sync_git)
    p = load_param_file(spec, params_dir=params_dir)
    mus_env = os.environ.get("KLREP_MUS")
    if mus_env:
        p.mus_dims = [int(x) for x in mus_env.split(",")]
        print(f"coarse override mus_dims={p.mus_dims}", flush=True)

    ss = calc_steady(p)
    t0 = time.time()
    st, diff, it = solve_model(p, max_iter=max_iter, device=device, verbose=True)
    elapsed = time.time() - t0
    print(f"solve done: {it} iters, diff={diff:.3e}, {elapsed:.1f}s", flush=True)

    # --- outputs ---
    outdir = Path(output) / "klreplication"
    outdir.mkdir(parents=True, exist_ok=True)
    torch.save({"spec": spec, "diff": diff, "iters": it,
                "g": st.g, "v_mat": st.v_mat, "mc_mat": st.mc_mat,
                "next_state": st.next_state, "n_states": st.const.state_grid.shape[0]},
               outdir / f"solution_spec{spec}.pt")

    # central-node steady-state comparison
    sg = st.const.state_grid
    center = ((sg[:, IDX_ZF]) ** 2 + (sg[:, IDX_THH] - p["tht_trgt_h"]) ** 2).argmin()
    lines = [
        f"spec {spec}  iters {it}  diff {diff:.3e}  elapsed {elapsed:.1f}s  device {device}",
        f"n_states {st.const.state_grid.shape[0]}",
        "central node vs deterministic steady state:",
        f"  q      {float(st.g.q[center]):.5f}   q_ss   {ss.q_ss:.5f}",
        f"  s      {float(st.g.s[center]):.5f}   s_ss   {ss.s_ss:.5f}",
        f"  l_h    {float(st.g.l_aggr[center,0]):.5f}   l_ss   {ss.l_ss[0]:.5f}",
        f"  l_f    {float(st.g.l_aggr[center,1]):.5f}",
        f"  v_h    {float(st.v_mat[center,0]):.5f}   v_ss   {ss.v_ss[0]:.5f}",
        f"  nom_i  {float(st.g.nom_i[center,0]):.6f}  1+rf_ss {1+ss.rf_ss:.6f}",
        f"  share_h {float(st.g.share[center,0]):.4f}  bF_h {float(st.g.bF_share[center,0]):.4f}",
    ]
    report = "\n".join(lines)
    print(report, flush=True)
    (outdir / f"ss_check_spec{spec}.txt").write_text(report + "\n")
    print(f"wrote {outdir}/solution_spec{spec}.pt and ss_check_spec{spec}.txt", flush=True)


if __name__ == "__main__":
    main()
