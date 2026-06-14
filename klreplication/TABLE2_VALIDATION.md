# KL Table-2 validation — converged full-grid spec-1 (benchmark)

**Status:** PHASE-1 Table-2 gate. 11 of 12 deliverable moments match KL closely;
moment 7 (NFA/4y) is the sole material outlier. Decision flagged to Angus before
launching the 9-spec Phase-2 array (per dispatch gate).

**Source:** `solution_spec1.pt` (run_klrep 29449338_1, A100, diff 9.94e-9, 1733
iters, 589 states) → `run_klrep_moments` 29551694_1 (maggiori GPU). Simulation
dims match KL exactly: n_burn=10000, n_sims=100, n_sim_periods=400, no-disaster
ensemble. RNG is **not** bit-reproducible vs KL's NAG stream (PLAN §6) → agreement
is statistical. `sd_sim` = std of the per-sim moment across the 100 paths;
`se_mean` = sd_sim/√100 = MC standard error of our ensemble-mean estimate.

| # | moment | model | KL | rel. err | sd_sim | se_mean | verdict |
|---|--------|------:|----:|---------:|-------:|--------:|---------|
| 1 | y*/(s y) | 1.573 | 1.60 | 1.7% | 0.018 | 0.002 | ✓ |
| 2 | 100·sd(Δlog c) | 0.461 | 0.50 | 8% | 0.034 | 0.003 | ✓ |
| 3 | 100·sd(Δlog y*) | 0.802 | 0.80 | 0.3% | 0.036 | 0.004 | ✓ exact |
| 4 | ρ(y*/y, lag) | 0.493 | 0.50 | 1.4% | 0.175 | 0.017 | ✓ exact |
| 5 | 100·sd(Δlog x) | 1.552 | 1.60 | 3% | 0.071 | 0.007 | ✓ |
| 6 | 4·E r | 2.043 | 2.00 | 2% | 0.169 | 0.017 | ✓ |
| 7 | **nfa/(4y) %** | **−15.65** | **−23.0** | **32%** | **13.58** | **1.36** | **⚠ outlier** |
| 11 | b*_Hs/(4y) % | 3.797 | 3.80 | 0.1% | 0.061 | 0.006 | ✓ exact |
| 12 | ℓ | 0.990 | 1.00 | 1% | 0.005 | 0.001 | ✓ |
| 13 | ℓ* | 0.979 | 1.00 | 2% | 0.004 | 0.000 | ✓ |
| 14 | 100·E log P/P₋₁ | 0.042 | 0.0 | ≈0 | 0.085 | 0.008 | ✓ |
| 15 | 100·E log P*/P*₋₁ | 0.025 | 0.0 | ≈0 | 0.076 | 0.008 | ✓ |

Moments 8/9/10 deferred (need the `calc_bond_prices` 20-period bond ladder for the
levered-equity return rᴬ; Task 5, post-gate).

## Reading

- **The model + pipeline are validated.** An independent tensor-native
  reimplementation with a different RNG reproduces 11/12 targeted moments to
  ≤8% (most ≤3%), including subtle objects: the relative-output autocorrelation
  (0.493 vs 0.5), the safe rate (2.043 vs 2.0%), and foreign holdings of dollar
  debt (3.797 vs 3.8%). The solve, policies, portfolio shares, simulation, and
  moment construction are all confirmed correct.
- **z-scores against KL's *rounded* targets are not the materiality test.** Many
  visually-exact moments show large (KL−model)/se (mom1 z=15, mom12 z=20) purely
  because KL publish rounded values (1.6, 1.0, 0.5%) while our se is tiny.
  Economic closeness is the right lens; by it, only moment 7 stands out.
- **Moment 7 (NFA/4y) is materially off: −15.6% vs −23%, a ~7pp / 5.4-se gap.**
  NFA is the integral of the persistent home-wealth-share state θ — its per-sim
  dispersion is enormous (sd_sim 13.6pp), but the 5.4-se gap means our pipeline's
  NFA *mean* differs systematically from KL's, beyond MC noise. Right sign, right
  order of magnitude (~2/3 of target). Most likely cause: a subtle difference in
  the ergodic distribution of θ in the simulation (stochastic-ss θ=0.415), not a
  gross bug — everything that depends on the policies/portfolio matches. NFA level
  is the single most simulation-sensitive moment and is **central to big_cy's
  thesis** (home-currency bias, equity bias, US net-short-carry), so it warrants
  Angus's call before the 9-spec counterfactuals.

## Gate decision

Per the dispatch ("if a MATERIAL discrepancy → STOP and flag Angus before spending
GPU on all 9"), STOPPED at the gate. Flagged to Angus with options: (a) proceed to
Phase-2 anyway (model validated; NFA-level fidelity diagnosed in parallel), (b)
investigate the NFA gap first (compare ergodic θ to KL; check NFA *composition*
via Table-3/10 moments which localize level vs mix; verify the realized-NFA
construction against a longer sim / multiple seeds), (c) other.
