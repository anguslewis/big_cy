# Phase-1 baseline — KL benchmark solved on Sherlock GPU (closeness note)

**Run:** `run_klrep` job `29449338_1`, spec 1 (benchmark, `param_file_1`), on a
maggiori **NVIDIA A100-SXM4-80GB** node. Completed 2026-06-13, exit 0.

## Pipeline / GPU validation (the Sherlock migration goal)
- **GPU confirmed used:** `torch 2.5.1+cu121`, `cuda.is_available=True`,
  `GPU[0]: NVIDIA A100-SXM4-80GB`, diagnostic **float64 matmul on `cuda:0` OK**,
  and the converged **solver tensors live on `cuda:0`**. Peak GPU mem 5.1 GB.
- **Full grid, full fidelity:** n_states = 589 (7 active Smolyak dims, μ=3),
  converged to **diff 9.94e-9** (< 1e-8) in **1733 dampened iterations**,
  wall **3:09:26** (~6.5 s/iter; linear tail rate ≈0.68 per 50 iters).
- End-to-end cluster path works: clone → sync_git → GPU array submit → conda env
  (`big_cy_klrep_env`, torch via pip CUDA wheel) → solve → outputs on oak
  (`$big_cy/data/output/klreplication/solution_spec1.pt`, `ss_check_spec1.txt`).

## Baseline solution — central node vs deterministic steady state
| object | model (central node) | det. SS | read |
|---|---|---|---|
| q (price of capital) | 1.00127 | 1.00509 | at SS (−0.4%) |
| s (terms of trade) | 0.99687 | 0.99177 | at SS (+0.5%) |
| ℓ_h | 0.99470 | 1.0 | at SS |
| ℓ_f | 0.97673 | 1.0 | ~at SS (−2.3%) |
| v_h (value) | 0.82647 | 1.0 | **lower** ✓ |
| 1+i (Home nominal) | 1.00697 | 1.01106 | **lower** ✓ |
| share_h (bond share of savings) | **−0.383** | — | **negative = equity bias** ✓ |
| bF_h (foreign-bond share) | +0.111 | — | small RoW-bond holding |

## Verdict: baseline ALIGNS (qualitatively + at SS level)
- **Real quantities (q, s, ℓ) sit at the deterministic steady state** (within
  ~0.5%; the small gaps are expected — the central grid node is at
  k_grid_mean=1.025·k_ss, and the *stochastic* solution differs from the
  deterministic SS).
- **v_h < 1 and the Home nominal/safe rate < 1+rf_ss** are the **correct
  disaster-risk-adjusted deviations**: the global solution prices the rare
  disaster (γ=21) that the deterministic SS ignores, so welfare and the safe rate
  fall. (See PLAN §13.3.)
- **Home is short safe bonds / leveraged into capital (share_h<0)** — the model's
  central **equity-bias** mechanism, exactly as expected.
- Consistent with the local coarse-grid validation (15-state run: v_h≈0.81,
  share_h≈−0.5, same signs) — the finer grid + tighter convergence shift
  magnitudes but not the economics.

## ⚠ What this gate does NOT yet cover (next code task)
The **quantitative KL Table-2 moment comparison** (NFA, portfolio composition,
carry, return moments, etc.) requires the **simulation + moment pipeline**
(`mod_results` port: simulate with the 3 clipping regimes, generalized IRFs,
`get_var_indices`/de-trending/moments) — not yet ported. This Phase-1 gate
validates: (1) the Sherlock GPU pipeline end-to-end, (2) the full-grid
full-fidelity solve converges, (3) the solution is economically correct at the
SS level and reproduces the KL equity-bias mechanism. The full Table-2 number-by-
number comparison is the next deliverable.

**Recommendation:** baseline is sound — pipeline works, solve converges, economics
correct. Proceed (after Angus's OK) to (a) port the simulation/moments for the
full Table-2 comparison, and/or (b) Phase-2: the 9-spec GPU array (`specs=1-9`,
`KLREP_CONV=1e-8`, 7h wall).
