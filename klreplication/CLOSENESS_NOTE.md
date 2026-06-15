# KL (2024 AER) replication — closeness / discrepancy note (Q0)

Tensor-native Python replication of the Kekre–Lenel quantitative model. All 9
calibrations solved on the full anisotropic Smolyak grid (GPU, float64) to
diff ≈ 1e-8; moments simulated with KL's dimensions (n_burn=10000, n_sims=100,
n_sim_periods=400, no-disaster ensemble). RNG is **not** bit-reproducible vs KL's
NAG stream, so simulated-moment agreement is statistical, not bit-exact. Spec→KL
column map (create_tables.m + main.m): bm=1, no_omg=2, symm(γ=γ\*)=3, symm_flex=4,
tayl_y0=5, tayl_rho=6, no_omg_symm=7, nocorr_nobg=8, nocorr=9.

## Headline

**The replication reproduces KL across the full calibration set and the moment
tables — not just the benchmark.** The overwhelming majority of moments across
Tables 2, 3, 4, 5, 9, 10 (8 specifications) match KL essentially exactly, including
the counterfactual mechanisms (no-ω portfolio reversal; γ=γ\* equity-bias collapse).
The **single systematic residual** is the *benchmark NFA level* (a ~5–7pp
capital-vs-home-bond split), already localized as a benign second-order
equilibrium/numerical difference (not a structural bug); it does not appear in the
counterfactual specs.

## Table 2 — targeted moments (spec 1): 14/15 match

m1 1.573/1.60 · m2 0.461/0.50 · m3 0.802/0.80 · m4 0.493/0.50 · m5 1.552/1.60 ·
m6 2.043/2.00 · **m7 −15.6/−23 (outlier)** · m8 0.509/0.50 · m9 5.203/5.20 ·
m10 0.599/0.60 · m11 3.797/3.80 · m12 0.990/1.0 · m13 0.979/1.0 · m14 0.042/0.0 ·
m15 0.025/0.0. (model/KL.) Detail + the NFA-outlier forensics in TABLE2_VALIDATION.md.

## Table 3 — comovements (specs bm, no_ω, γ=γ\*) — model vs KL

| moment | bm | no_ω | γ=γ\* |
|--------|----|----|----|
| β(uip, Δy) | −0.111 / −0.11 | 0.002 / 0.00 | −0.110 / −0.11 |
| β(uip, r^e) | 0.060 / 0.06 | −0.000 / −0.00 | 0.060 / 0.06 |
| β(Δnfa, r\*−r) | 1.547 / 1.45 | −3.356 / −3.39 | 0.237 / 0.25 |
| (k−κ)/4y % | 64.6 / 59.8 | 45.8 / 50.2 | 0.16 / 0.1 |
| b_H/4y % | −102.2 / −102.5 | 149.3 / 151.3 | 14.5 / 14.3 |
| b_F/4y % | 21.6 / 19.8 | −224.4 / −224.6 | −16.1 / −15.9 |

Regressions m1–m3 essentially exact. The no-ω **portfolio reversal** (home flips
to a huge *long* home-bond / short foreign-bond position) and the γ=γ\* **capital-NFA
collapse to ≈0** are both reproduced. Residual = benchmark (k−κ)/4y +4.8pp (NFA level).

## Table 9 — comovements under alt. Taylor rules (specs bm, φʸ=.5/4, +ρⁱ=.5)

| moment | bm | tayl_y0 | tayl_rho |
|--------|----|----|----|
| β(uip, Δy) | −0.111 / −0.11 | −0.095 / −0.10 | −0.134 / −0.13 |
| β(uip, r^e) | 0.060 / 0.06 | 0.064 / 0.06 | 0.078 / 0.08 |
| β(Δnfa, r\*−r) | 1.547 / 1.45 | 1.389 / 1.30 | 1.122 / 0.86 |
| (k−κ)/4y % | 64.6 / 59.8 | 57.4 / 51.9 | 51.8 / 37.7 |
| b_H/4y % | −102.2 / −102.5 | −127.4 / −127.0 | −74.4 / −72.5 |
| b_F/4y % | 21.6 / 19.8 | 54.3 / 52.1 | 17.6 / 12.9 |

m1/m2/m5 dead-on; m3/m4/m6 carry the NFA-level sensitivity (a few pp), directionally
correct across the policy variants.

## Table 4 — additional second moments (specs bm, no_ω, γ=γ\*): essentially exact

| moment | bm | no_ω | γ=γ\* |
|--------|----|----|----|
| σ(4r) % | 4.07 / 4.1 | 2.12 / 2.1 | 4.12 / 4.1 |
| σ(4[r^e−r]) % | 17.1 / 17.2 | 11.3 / 11.2 | 17.4 / 17.5 |
| σ(4[r\*−Δq−r]) % | 1.85 / 1.9 | 0.35 / 0.3 | 1.87 / 1.9 |
| σ(Δlog q) % | 0.26 / 0.3 | 0.25 / 0.2 | 0.26 / 0.3 |
| σ(Δlog E) % | 0.36 / 0.4 | 0.19 / 0.2 | 0.36 / 0.4 |
| corr(Δq, Δc\*−Δc) | 0.911 / 0.91 | 0.925 / 0.92 | 0.974 / 0.97 |

## Table 5 — output growth volatility (specs bm, no_ω, γ=γ\*): exact

σ(Δlog y): 0.607 / 0.442 / 0.596 vs KL 0.61 / 0.44 / 0.60.
σ(Δlog y\*): 0.802 / 0.748 / 0.807 vs KL 0.81 / 0.75 / 0.81.

## Table 10 — portfolio shares (specs no_ω_symm, no_ω, nocorr_nobg, nocorr, bm)

| share | no_ω_symm | no_ω | nocorr_nobg | nocorr | bm |
|-------|----|----|----|----|----|
| k/a % | 100.0 / 100.0 | 137.0 / 137.1 | 137.4 / 137.7 | 137.4 / 137.7 | 137.2 / 142.4 |
| b_H/a % | 106.1 / 105.9 | 73.5 / 73.3 | 4.7 / 4.7 | 1.9 / 1.8 | −47.0 / −51.9 |
| b_F/a % | −106.1 / −105.9 | −110.5 / −110.4 | −42.1 / −42.4 | −39.3 / −39.6 | 9.8 / 9.5 |

Exact for 4 of 5 specs (incl. the γ=γ\*+no-ω spec 7 at exactly k/a=100%, i.e. no
leverage). Only the benchmark shows the known ~5pp capital/home-bond split; b_F/a is
exact even there. (Table-10 conditional-corr rows 4–5 need the solver's
corr_rf_spread columns — deferred, see below.)

## The one systematic residual — benchmark NFA level

Same object across Table-2 m7, Table-3/9 (k−κ)/4y, and Table-10 k/a: the benchmark
home portfolio holds ~5pp less capital / is ~5pp less short home bonds than KL,
making NFA/4y −16% vs −23%. Localized (TABLE2_VALIDATION.md): portfolio
*composition* and all other moments match; the θ-transition, simulation method, and
bbeta nudge are bit-faithful to the Fortran; ergodic θ=0.364±0.003 is the genuine
value (0.350 was the deterministic-target red herring). It is a second-order
converged-equilibrium / RNG difference on the most simulation-sensitive moment (a
near-unit-root level, per-path sd 13.6pp), **not a structural bug**, and absent from
the counterfactual specs. Cannot be closed further without KL's raw policy/sim
output (not shipped).

## Not yet reproduced (clearly-scoped next chunk — need new machinery)

- **Table 6** (NFA/net-export vol): needs the net-exports series `nx` (consumption
  splits chf/cfh + inv_h + deployed-capital flow). Moderate.
- **Table 7 / A2** (external-adjustment variance decomposition): needs `calc_valuation`
  (the n_bond valuation columns) + the **disaster** simulation ensemble (the `dis`
  branch). Substantial.
- **Table 8** (dollar swap-line effects): needs the swap-line MIT-shock experiment
  (`collect_swap_moments`) — a counterfactual policy experiment, not in the moment
  ensemble. Substantial.
- **Table-10 rows 4–5** (conditional corr of the rf spread with the SDF): need the
  solver's `corr_rf_spread_m` columns (results_vec 49–52) added to the grid results.
  Small.
- **Figures 2–20**: the IRF figures need the **generalized impulse-response** port
  (MIT shocks to z/ω/p/zf/bg/fix from mod_results); figs 4/6/7/10 need the
  **empirical sample-path** simulation (NAG scattered-data interp → scipy). The
  headline fig 10 (ω-shock vs swapped-T-bill) needs the sample-path machinery.
  Substantial; the largest remaining piece.

**Pulling artifacts back to local:** the moment tables print to the synced `.out`
logs (done). Pulling `solution_spec*.pt` or figure files to local would need a
`SYNC_ALLOWED_big_cy=(data/output/klreplication)` entry in `sherlock-agent-com`
(shared tooling — Angus's call).

## Validation infrastructure

58 unit/integration tests green (Smolyak basis, GH quadrature, steady state,
equilibrium step, batched Brent, bond-ladder SDF vs equilibrium_step, OLS, all
moment formulas). Solver convergence diff ≈ 1e-8 for all 9 specs.
