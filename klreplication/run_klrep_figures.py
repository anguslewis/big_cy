"""Cluster entry — generalized-IRF figures (PDF) for the KL specifications.

For each requested spec + shock {p, omg, z, disaster}: build the MIT-shock transition,
simulate n_sample shocked + baseline IRF paths, difference them (generalized IRF =
mean_i(shocked - baseline + baseline_t1)), build the named series, and plot the KL
6-panel IRF figures (shock; r^e-r; r*-Dq-r; log q; theta; log y) as PDFs under
$big_cy/data/output/klreplication/figures/. Impact + peak values are printed to the
.out log for numeric validation. Port of mod_results IRF sim + extract_irfs +
create_figures (figs 2/3 for p/omg).

Env: KLREP_SPECS (default 1), KLREP_NSAMPLE (default 100), KLREP_NBURN (default 10000).
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
import matplotlib                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402

from klrep.config import get_device                  # noqa: E402
from klrep.params import IDX_THH, IDX_DIS, IDX_OMG   # noqa: E402
from klrep.grid.interp_fit import interp_factor, interp_solve  # noqa: E402
from klrep.model.period_block import compute_current_period    # noqa: E402
from klrep.model.bond_ladder import compute_bond_columns       # noqa: E402
from klrep.simulate.simulate import build_sim_coeffs, stochastic_ss, burn_in  # noqa: E402
from klrep.solve.irf import build_shock_transition, simulate_irf_paths  # noqa: E402
from klrep.post.table2_series import build_table2_series        # noqa: E402
from klrep.params import IDX_K, IDX_ZF                          # noqa: E402
from run_klrep_moments import rebuild_state                     # noqa: E402

# shock spec: (label, sidx or None for disaster, size-fn, is_disaster)
def shock_specs(p):
    return {
        "p":   dict(sidx=2, size=2.0 * p.sig_dis, is_disaster=False, z_at_jump=0.0),
        "omg": dict(sidx=3, size=2.0 * p.sig_omg, is_disaster=False, z_at_jump=0.0),
        "z":   dict(sidx=0, size=-5.0 * p.sig_z,  is_disaster=False, z_at_jump=-5.0 * p.sig_z),
        "dis": dict(sidx=None, size=0.0, is_disaster=True, z_at_jump=-p["disast_shock"]),
    }


def _k_next_new(const, g):
    cp = compute_current_period(const.state_grid[:, IDX_K], g.s, g.l_aggr, g.q,
                                const.state_grid[:, IDX_ZF], const.wealth_share_grid,
                                aalpha=const.aalpha, sigma=const.sigma,
                                ddelta=const.ddelta, zeta=const.zeta,
                                varsigma=const.varsigma_vec)
    msc = const.min_cons_sav
    cs = torch.minimum(torch.maximum(g.c_spending, torch.full_like(g.c_spending, msc)),
                       cp.wealth_vec + cp.w_choice * g.l_aggr - msc)
    sav = cp.wealth_vec + cp.w_choice * g.l_aggr - cs
    return (sav * (1.0 - g.share) / g.q.unsqueeze(1)).sum(dim=1)


def _irf_series(panel_S, panel_N):
    """Generalized-IRF level series (T,) from shocked/baseline panels (n_sample,T):
    mean_i(S - N + N[:,0:1])."""
    return (panel_S - panel_N + panel_N[:, :1]).mean(dim=0)


def _dev(series, *, log):
    """extract_irfs transform: 10000*( f(series)[2:] - f(series)[0] ), f=log or id."""
    f = torch.log(series) if log else series
    return (1e4 * (f[2:] - f[0])).cpu().numpy()


def compute_irf(const, g, st, sc, starts, kperm, p, shock, base_coeffs):
    """Return dict var -> IRF deviation array (198,) for one shock. The BASELINE
    (no-shock) path jumps at t=3 via the ZERO-shock transition `base_coeffs`; the
    shocked path jumps via the shock-specific transition (or the disaster node)."""
    if shock["is_disaster"]:
        coeffs = None
    else:
        tr = build_shock_transition(const, g, st.v_mat, st.mc_mat, st.next_state, kperm,
                                    shock["sidx"], shock["size"], conv=1e-8)
        coeffs = interp_solve(interp_factor(const.smol_polynom), tr)
    simS = simulate_irf_paths(sc, coeffs, starts, is_disaster=shock["is_disaster"],
                              z_at_jump=shock["z_at_jump"])
    simN = simulate_irf_paths(sc, base_coeffs, starts, is_disaster=False, z_at_jump=0.0)
    bg = compute_bond_columns(const, g, st.v_mat, st.mc_mat, st.next_state)
    S = build_table2_series(const, g, simS, disast_shock=p["disast_shock"], bond_grid=bg)
    N = build_table2_series(const, g, simN, disast_shock=p["disast_shock"], bond_grid=bg)
    # state-path panels (theta/omg/p) come from the econ state_series.
    out = {}
    out["shock_p"] = _dev(_irf_series(torch.exp(simS["state_series"][:, :, IDX_DIS]),
                                      torch.exp(simN["state_series"][:, :, IDX_DIS])), log=False)
    omgS = torch.exp(simS["state_series"][:, :, IDX_OMG]) + const.omg_shift
    omgN = torch.exp(simN["state_series"][:, :, IDX_OMG]) + const.omg_shift
    out["shock_omg"] = _dev(_irf_series(omgS, omgN), log=False)
    out["shock_z"] = _dev(_irf_series(torch.exp(torch.cumsum(simS["z_shock_series"], 1)),
                                      torch.exp(torch.cumsum(simN["z_shock_series"], 1))), log=True)
    out["excA"] = _dev(_irf_series(S["exc_retA"], N["exc_retA"]), log=False)
    out["excrf"] = _dev(_irf_series(S["exc_rf"], N["exc_rf"]), log=False)
    out["log_qx"] = _dev(_irf_series(S["qx"], N["qx"]), log=True)
    out["thth"] = _dev(_irf_series(simS["state_series"][:, :, IDX_THH],
                                   simN["state_series"][:, :, IDX_THH]), log=False)
    out["log_yh"] = _dev(_irf_series(S["yh"], N["yh"]), log=True)
    return out


PANELS6 = [("excA", "r^e - r"), ("excrf", "r* - dlog q - r"), ("log_qx", "log q"),
           ("thth", "theta"), ("log_yh", "log y")]


def main():
    specs = [int(x) for x in os.environ.get("KLREP_SPECS", "1").split(",")]
    n_sample = int(os.environ.get("KLREP_NSAMPLE", 100))
    n_burn = int(os.environ.get("KLREP_NBURN", 10000))
    device = get_device()
    figdir = Path(output) / "klreplication" / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    print(f"=== run_klrep_figures: specs={specs} n_sample={n_sample} device={device} ===", flush=True)

    for spec in specs:
        t0 = time.time()
        p, st, _ = rebuild_state(spec, device)
        sc = build_sim_coeffs(st)
        ss_state, _ = stochastic_ss(sc)
        pool = burn_in(sc, ss_state, n_burn=n_burn)
        starts = torch.cat([ss_state, pool[:n_sample - 1]], dim=0)   # mean + pool draws
        kperm = _k_next_new(st.const, st.g)
        # zero-shock baseline transition (shared across shocks).
        tr0 = build_shock_transition(st.const, st.g, st.v_mat, st.mc_mat, st.next_state,
                                     kperm, 0, 0.0, conv=1e-8)
        base_coeffs = interp_solve(interp_factor(st.const.smol_polynom), tr0)
        sspecs = shock_specs(p)
        for shock_name in ("p", "omg", "z", "dis"):
            irf = compute_irf(st.const, st.g, st, sc, starts, kperm, p, sspecs[shock_name],
                              base_coeffs)
            shock_key = {"p": "shock_p", "omg": "shock_omg", "z": "shock_z", "dis": "shock_z"}[shock_name]
            panels = [(shock_key, shock_name)] + PANELS6
            fig, axes = plt.subplots(2, 3, figsize=(12, 7))
            for ax, (key, title) in zip(axes.flat, panels):
                ax.plot(irf[key][:60], lw=1.5)
                ax.axhline(0, color="k", lw=0.5)
                ax.set_title(title); ax.set_xlabel("quarters")
            fig.suptitle(f"spec {spec} — {shock_name}-shock IRF (x1e4)")
            fig.tight_layout()
            out_pdf = figdir / f"irf_spec{spec}_{shock_name}.pdf"
            fig.savefig(out_pdf, format="pdf", dpi=300)
            plt.close(fig)
            # numeric validation: impact (h=0) + peak of the shock + key panels
            imp = {k: float(irf[k][0]) for k in [shock_key, "excA", "excrf", "log_qx", "thth", "log_yh"]}
            print(f"spec {spec} {shock_name}: impact " +
                  " ".join(f"{k}={v:.1f}" for k, v in imp.items()) +
                  f"  -> {out_pdf.name}", flush=True)
        print(f"spec {spec} figures done ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
