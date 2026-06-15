"""Cluster/local entry — compute the per-spec moment TABLES across all solved KL
specifications and print them grouped by KL table (with the spec->column map), so
the comparison returns via the synced .out log.

Loads solution_spec1..9.pt from $big_cy/data/output/klreplication, rebuilds each
spec's solver constants, simulates (no-disaster ensemble), prices the bond ladder,
and computes Table-2 (spec 1), Table-3/4/5 (specs 1,2,3), Table-9 (specs 1,5,6),
and Table-10 portfolio shares (specs 7,2,8,9,1). Tables 6(nx)/7(valuation)/8(swap)
and the Table-10 conditional-corr rows need additional machinery — flagged, not
computed here.

Env: KLREP_SPECS (comma list, default 1-9), KLREP_NBURN/NSIMS/NPER (sim sizes).
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
from klrep.config import get_device                  # noqa: E402
from klrep.post.moments import (compute_table2_moments_per_sim,  # noqa: E402
                                compute_extended_tables, KL_TARGETS)
from klrep.post.table2_series import build_table2_series  # noqa: E402
from klrep.model.bond_ladder import compute_bond_columns  # noqa: E402
from klrep.simulate.simulate import (build_sim_coeffs, stochastic_ss,  # noqa: E402
                                     burn_in, simulate_ensemble)
# rebuild_state is shared with the moments entry.
from run_klrep_moments import rebuild_state          # noqa: E402

# spec -> column maps (1-based spec numbers), from create_tables.m + main.m ix_*.
SPEC = dict(bm=1, no_omg=2, symm=3, symm_flex=4, tayl_y0=5, tayl_rho=6,
            no_omg_symm=7, nocorr_nobg=8, nocorr=9)


def spec_series(spec, device, n_burn, n_sims, n_per):
    p, st, sol = rebuild_state(spec, device)
    sc = build_sim_coeffs(st)
    ss_state, _ = stochastic_ss(sc)
    pool = burn_in(sc, ss_state, n_burn=n_burn)
    sim = simulate_ensemble(sc, pool, n_sims=n_sims, n_periods=n_per, disaster=False)
    bond_grid = compute_bond_columns(st.const, st.g, st.v_mat, st.mc_mat, st.next_state)
    s = build_table2_series(st.const, st.g, sim, disast_shock=p["disast_shock"],
                            bond_grid=bond_grid)
    return p, s


def main():
    specs = [int(x) for x in os.environ.get("KLREP_SPECS", "1,2,3,4,5,6,7,8,9").split(",")]
    n_burn = int(os.environ.get("KLREP_NBURN", 10000))
    n_sims = int(os.environ.get("KLREP_NSIMS", 100))
    n_per = int(os.environ.get("KLREP_NPER", 400))
    device = get_device()
    print(f"=== run_klrep_tables: specs={specs} device={device} "
          f"n_burn={n_burn} n_sims={n_sims} n_per={n_per} ===", flush=True)

    t0 = time.time()
    M, ext = {}, {}
    for spec in specs:
        p, s = spec_series(spec, device, n_burn, n_sims, n_per)
        per = compute_table2_moments_per_sim(s, bg_yss=p["bg_yss"])
        M[spec] = {k: float(v.mean()) for k, v in per.items()}
        ext[spec] = compute_extended_tables(s)
        print(f"spec {spec} done ({time.time() - t0:.0f}s)", flush=True)

    def row(label, vals, fmt="{:>9.3f}"):
        return f"  {label:<26} " + " ".join(fmt.format(v) for v in vals)

    # ---- Table 2 (spec 1, all 15) ----
    if 1 in M:
        print("\n=== TABLE 2 (spec 1) model vs KL ===", flush=True)
        for k in range(1, 16):
            print(f"  m{k:<2} model {M[1][k]:>9.3f}   KL {KL_TARGETS[k]:>7.2f}", flush=True)

    def table_block(name, cols, keys, labels):
        present = [c for c in cols if c in ext]
        if not present:
            return
        print(f"\n=== {name} (specs {present}) ===", flush=True)
        print(row("moment", present, "{:>9d}"), flush=True)
        for key, lab in zip(keys, labels):
            print(row(lab, [ext[c][key] for c in present]), flush=True)

    t3keys = ["t3_1", "t3_2", "t3_3", "t3_4", "t3_5", "t3_6"]
    t3lab = ["b(uip, dy)", "b(uip, r^e)", "b(dnfa, r*-r)",
             "(k-kap)/4y %", "b_H/4y %", "b_F/4y %"]
    table_block("TABLE 3 comovements", [SPEC["bm"], SPEC["no_omg"], SPEC["symm"]], t3keys, t3lab)
    table_block("TABLE 9 alt policy", [SPEC["bm"], SPEC["tayl_y0"], SPEC["tayl_rho"]], t3keys, t3lab)

    t4keys = ["t4_1", "t4_2", "t4_3", "t4_4", "t4_5", "t4_6"]
    t4lab = ["sd(4r) %", "sd(4 r^e-r) %", "sd(4 exc_rf) %", "sd(dlog q) %",
             "sd(dlog E) %", "corr(dq, dc*-dc)"]
    table_block("TABLE 4 second moments", [SPEC["bm"], SPEC["no_omg"], SPEC["symm"]], t4keys, t4lab)

    t5keys = ["t5_1", "t5_2"]
    t5lab = ["sd(dlog y) %", "sd(dlog y*) %"]
    table_block("TABLE 5 output vol", [SPEC["bm"], SPEC["no_omg"], SPEC["symm"]], t5keys, t5lab)

    t6keys = ["t6_1", "t6_2", "t6_3", "t6_4", "t6_5", "t6_6"]
    t6lab = ["sd(dnfa/y) %", "sd(nx/y) %", "sd((dnfa-nx)/y) %",
             "mean dnfa/y %", "mean nx/y %", "mean (dnfa-nx)/y %"]
    table_block("TABLE 6 nfa/nx", [SPEC["bm"], SPEC["no_omg"], SPEC["symm"]], t6keys, t6lab)

    t10keys = ["t10_k", "t10_bH", "t10_bF", "t10_4", "t10_5"]
    t10lab = ["k/a %", "b_H/a %", "b_F/a %",
              "corr(rf-sp, m)", "corr(rf-sp, m*)"]
    table_block("TABLE 10 portfolios + corr",
                [SPEC["no_omg_symm"], SPEC["no_omg"], SPEC["nocorr_nobg"],
                 SPEC["nocorr"], SPEC["bm"]], t10keys, t10lab)

    print("\nNOTE: Table 7(valuation,disaster) + Table 8(swap) + figures(IRFs/"
          "sample-path) need additional machinery (calc_valuation + disaster "
          "ensemble, swap-line MIT experiment, generalized IRFs) — in progress.",
          flush=True)
    print(f"\ntotal {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
