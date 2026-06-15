# PLAN — Tensorized-Python replication of the Kekre–Lenel (2024 AER) quantitative model

**Agent:** big_cy-klrep · **Task:** Assignment 04 (Q0) · **Status:** drafted 2026-06-13.

This is a **language refactor**, not an improvement pass: implement *exactly* what KL
do — same equations, same algorithm, same inputs — to reproduce their point estimates
and figures. Extensions come later. Everything is **tensor-native** (PyTorch),
device-agnostic, GPU-ready, and runs **locally only**.

> **Reference files** (read these for detail; this PLAN does not duplicate them):
> - `reference/kl_model_equations.md` — equations (1)–(21), (23)–(72), (N1)–(N3); the
>   ω convenience-yield mechanism; symbol↔name map. **The equation authority.**
> - `reference/fortran_mod_calc.md` — the solver (`calc_steady`, `calc_sol`,
>   `calc_equilibrium_and_update`, the 4 Brent root-finders, damping, staged activation).
> - `reference/fortran_mod_results.md` — simulation, generalized IRFs, output files.
> - `reference/fortran_auxcodes.md` — Smolyak grid/basis, Gauss–Hermite quadrature.
> - `reference/matlab_postproc.md` — `get_var_indices` column map (151 rows), figures,
>   tables, de-trending conventions.
> - `…/big_cy/data/raw/klreplication/MANIFEST.md` — inputs + the 65-parameter layout.
> - Online appendix saved at `…/papers/Kekre_Lenel_2024_AER_Online_appendix.pdf`
>   (equations (33)–(72) are in Appendix A.6, pp. 9–13).

---

## 0. The model in one paragraph (what we are reproducing)

A two-country (Home = US, Foreign = RoW, relative size ζ\*) open-economy NK model with
Epstein–Zin recursive preferences and **bonds in the utility function**. Households value
**safe dollar bonds** through a convenience term Ω(·); the marginal convenience yield ω_t
is driven by an exogenous safety-demand shock ω^d and the (fixed-ratio) supply of
government T-bills. Foreign *also* values safe dollar bonds — the channel for RoW dollar
demand. Four aggregate shocks: global TFP z (unit-root + rare disaster φ), relative TFP
z_F, disaster probability p, and safety demand ω^d. Rotemberg sticky wages, CPI-targeting
Taylor rules, freely-deployable capital with investment adjustment costs. Solved globally
on an **anisotropic Smolyak grid** by **dampened backward iteration**, with all real
variables scaled by the TFP trend z_t (stationary "re-scaled economy") and value/consumption
additionally scaled by agent wealth μ̃ for numerical conditioning. 9 calibrations
(`param_file_1..9`), file 1 = benchmark.

---

## 1. Notation & symbol ↔ name map

The full math-symbol ↔ Fortran-name map is in `reference/kl_model_equations.md §1`. The
Python convention adds a third column. **Rule (from coding_conventions): paper symbol →
descriptive Python name, with the symbol in a comment at definition.** We keep the Fortran
names as the bridge because the Fortran is the executable ground truth.

### 1.1 Parameters (the 65 in `param_file_*.csv`; order = `param_input_vec`)

| Symbol | Meaning | Fortran | Python |
|---|---|---|---|
| ζ\* | Foreign relative size | `zeta` | `zeta` |
| β, β\* | discount factors (H/F) | `bbeta_vec` | `beta[2]` |
| γ, γ\* | risk aversion (H/F) | `gma_vec` | `gamma[2]` |
| ψ | IES (both) | `ies_vec` | `ies[2]` |
| χ^W | Rotemberg wage cost | `chi` | `chi_w` |
| σ | trade (Armington) elasticity | `sigma` | `sigma` |
| ς | home-bias in CES bundle (H/F) | `varsigma_vec` | `varsigma[2]` |
| ℓ̄ | labor target (normalization) | `l_target` | `l_target` |
| δ | depreciation | `ddelta` | `delta` |
| α | capital share | `aalpha` | `alpha` |
| χ^x | investment adj. curvature | `chiX` | `chi_x` |
| — | Home investment share | `inv_share_h` | `inv_share_h` |
| φ_w, ε_w | wage PC objects | `phi_w`,`vareps_w` | `phi_w`,`eps_w` |
| φ (Taylor) | inflation coeff (H/F) | `phi_h`,`phi_f` | `phi_pi[2]` |
| φ_y | output-gap coeff (H/F) | `phi_yh`,`phi_yf` | `phi_y[2]` |
| — | Taylor intercept (H/F) | `tayl_ic_h`,`tayl_ic_f` | `tayl_ic[2]` |
| ρ_i | rate smoothing | `rho_i` | `rho_i` |
| b̄g | T-bill supply rule | `b_lmbd`,`bg_yss` | `b_lmbd`,`bg_yss` |
| σ^z | global TFP innov. SD | `sig_z` | `sig_z` |
| σ^F, ρ^F | rel-TFP SD / AR(1) | `zf_std`,`rho_zf` | `zf_std`,`rho_zf` |
| disaster | p̄ level, jump, AR, SD | `disast_p_in`,`disast_shock`,`disast_rho`,`disast_std_in` | `disast_*` |
| ω^d | safety demand: mean/SD/shift/AR | `omg_mean`,`omg_std`,`omg_shift`,`omg_rho` | `omg_*` |
| ρ^{pω} | disaster–safety corr | `corr_omg_dis` | `corr_omg_dis` |
| θ̄ | wealth-share target (H) | `tht_trgt_vec` | `tht_trgt[2]` |
| grid devs/means | per-dim grid half-widths & centers | `*_grid_dev`,`*_grid_mean`,`k_dev_param`,`w_dev_param`,`*_grid_adj` | same |
| μ-levels | per-dim Smolyak level (9) | `vector_mus_dimensions` | `mus_dims[9]` |
| n_quad\_\* | GH nodes per shock (4) | `n_uni_quad_*` | `n_quad[4]` |
| flags | `foreign_trading`,`run_bg`,`run_samp` | same | same |

Derived in `init_setup` (port verbatim): `disast_p = exp(p_in + std²/2)`,
`disast_std = sqrt((exp(std²)−1)·exp(2 p_in+std²))`, `sig_dis = std_in·√(1−ρ²)`,
`sig_omg`, `sig_zf` analogous; the shock covariance + Cholesky (`F07FDF`='U' = `cholesky(upper)`);
grid devs `= 3·SD`; IRF shock sizes (`irf_shock_sizes`, `irf_ll_shock_sizes`).

### 1.2 State vector (Smolyak dims; `idx_*`), shocks (`sidx_*`), sizes

```
state[9] = [k (capital), θ_h (Home wealth share), z_F (rel TFP), w_h, w_f (wages),
            p (disaster prob), i_h, i_f (lagged nom rates), ω (convenience)]
```
Baseline: dims 7 (`i_h`,`i_f`) are **inactive** unless `rho_i>0` → 7 active dims.
Shocks[4] = [z, z_F, p, ω] (`sidx_z,sidx_zf,sidx_p,sidx_omg`); + disaster node + no-shock
node ⇒ `n_quad = n_GH + 2`. `n_I=2`, `n_bond=20` (bond maturity ladder used in valuation),
`n_interp = 52 + 6·12 = 124` interpolated variables.

### 1.3 Scaling devices (must be reproduced exactly — see eq file §10–§11)

- **Re-scaled economy:** all real vars ÷ z_t; the double-tilde `x̃̃` further deflates lagged
  capital/wage *states* by `exp(σ^z ε^z + φ)`. Makes the system stationary.
- **Wealth scaling:** `v̂ ≡ μ̃⁻¹ ṽ`, `ĉ ≡ μ̃⁻¹ c̃`, with
  `μ̃_t = (θ_t(π_t+(1−δ)q^k_t) k̄̃̃_{t−1} + ā)/b̄`. (N1)–(N3) give SDF, CE, and μ̃_{t+1}.
  ω enters via `(1+r)/(1−ω)` and the seigniorage term in μ̃_{t+1} / θ transition (51).

---

## 2. Module / function structure (function-per-file, tensor-native)

Layout under `~/code/big_cy/klreplication/`:

```
klreplication/
├── PLAN.md                      (this file)
├── reference/                   (the 5 digests above)
├── src/klrep/
│   ├── config.py                # device/dtype; global tensor backend selection
│   ├── params.py                # Params dataclass; load_param_file(i) -> Params
│   ├── create_param_files.py    # PORT of create_param_files.m (+ get_moments, sim_moments)
│   ├── grid/
│   │   ├── smolyak_elem_isotrop.py
│   │   ├── smolyak_elem_anisotrop.py
│   │   ├── smolyak_grid.py          # Chebyshev-extrema (Lobatto) nodes, m(i)=2^(i-1)+1
│   │   ├── smolyak_polynomial.py    # 1st-kind Chebyshev T basis (batched)
│   │   ├── interp_fit.py            # solve B·c = v  (torch.linalg.solve)
│   │   ├── interp_eval.py           # B_new @ c      (batched matmul)
│   │   └── state_grid.py            # map smol_grid -> economic state_grid (grid_setup)
│   ├── quad/
│   │   ├── gauss_hermite.py         # standard-normal GH (hermegauss / √(2π)); prob-normalized
│   │   └── shock_grid.py            # product nodes × Cholesky; + no-shock + disaster nodes
│   ├── model/
│   │   ├── util_fun.py              # Φ(ℓ), Φ'(ℓ); felicity bundle  (eq 3)
│   │   ├── production.py            # κ allocation, wages, profits  (54-57,63,70-72)
│   │   ├── prices.py               # P/P_H, P*/P*_F, q (terms of trade) (70-72)
│   │   ├── capital_price.py        # q^k  (60)
│   │   ├── returns.py              # r, r*, r^k, bond returns (67-69)
│   │   ├── sdf.py                  # m̃, m̃*, CE via (N1),(N2)
│   │   ├── portfolio_foc.py        # scaled portfolio residual (eq file; the FOC kernel)
│   │   ├── bond_clearing_home.py   # solve i_t : Home bond mkt (66) [Brent]
│   │   ├── bond_clearing_foreign.py# solve i*−i : Foreign bond mkt   [Brent]
│   │   ├── capital_share.py        # portfolio share capital vs bond [Brent]
│   │   ├── terms_of_trade.py       # solve s_t : goods mkt (35,36,44,45,58,59,61,62,65) [Brent]
│   │   ├── consumption.py          # consumption-savings via (N1) (37,40,41,46,49,50)
│   │   ├── labor_supply.py         # union FOC (52,53)
│   │   ├── inflation.py            # Taylor rules -> P/P_{-1} (10)
│   │   ├── wealth_transition.py    # θ_{t+1} (51); aggregate savings k̄̃_t (64)
│   │   └── equilibrium_step.py     # calc_equilibrium_and_update: the 10 steps, batched over grid
│   ├── solve/
│   │   ├── steady_state.py         # calc_steady
│   │   ├── solve_model.py          # calc_sol: dampened backward iteration + staged activation
│   │   ├── damping.py              # per-object relaxation weights (match Fortran schedule)
│   │   └── irf.py                  # MIT-shock generalized IRFs (+ long-lived run_bg)
│   ├── simulate/
│   │   ├── stochastic_ss.py        # no-shock fixed point
│   │   ├── simulate.py             # burn-in + ensemble; sampleK_mat path-driven sims
│   │   └── write_series.py         # emit the series the post-proc reads (or in-memory)
│   ├── post/
│   │   ├── var_indices.py          # PORT of get_var_indices.m (151-row column map)
│   │   ├── read_series.py          # assemble stacked [shocks;states;vars;value]
│   │   ├── extract_series.py / extract_irfs.py
│   │   ├── moments.py              # calc_moments + collect_moments (+ swap)
│   │   ├── figures.py              # fig 2,3,5,10,11-20 (+ recession 4,6,7) — STATA? see §6
│   │   └── tables.py               # table 1-10
│   └── run/
│       ├── solve_calibration.py    # entry: solve_calibration(i) end-to-end (≈ ./main i)
│       └── reproduce_all.py        # loop calibrations -> figures/tables (≈ runme.sh steps 4-5)
└── tests/
    ├── test_smolyak_basis.py       # vs Fortran basis on fixed elems + random pts
    ├── test_quadrature.py          # GH nodes/weights vs analytic moments
    ├── test_steady_state.py        # vs Fortran SS print (output_*.txt)
    └── test_grid_csv.py            # grid.csv / grid_locs.csv vs Fortran
```

Each non-trivial routine is **one file, one function** (per C++-derived discipline). Pure
functions take tensors + a `Params`/`Grid` struct and return tensors; no hidden globals.

---

## 3. Solve algorithm (port of `calc_sol`)

Detailed mapping in `reference/fortran_mod_calc.md`. Skeleton:

```
init_setup(param_file_i)         # params, covariance+Cholesky, GH+shock grid, Smolyak grid
calc_steady()                    # deterministic SS -> grid means (k,w_h,w_f,...), ā,b̄
grid_setup()                     # smol_grid -> state_grid; next_zf/dis/omg transition mats
# --- dampened backward iteration ---
initialize guesses (v̂,v̂*,ĉ,ĉ*,ℓ,ℓ*, i,i*,q^k,s, P/P_{-1},P*/P*_{-1}, bond shares)
outer_iter = 0; diff = ∞
while diff > 1e-8 and outer_iter < max_iter (=5000):
    fit interpolation coeffs:  solve smol_polynom · coeffs = current_values   (torch.linalg.solve)
    parallel over n_states grid points (vectorized, not OMP loop):
        evaluate Smolyak basis at the n_quad next-states (interp_eval)
        calc_equilibrium_and_update(point)   # the 10 steps below
    apply per-object damping (relaxation weights); compute diff vs previous
    STAGED ACTIVATION:  bond clearing on after iter 10;  foreign_trading on after iter 100
post: calc_bond_prices, calc_valuation (n_bond ladder), IRFs
```

**The 10 steps inside `equilibrium_step.py`** (eq groups from `kl_model_equations.md §11`):
1. current production (54-57,63,70-72) → RHS of budgets (41,50)
2. next-period production → expected capital return (69)
3. price of capital q^k (60)
4. **Home bond clearing**: solve i_t (66) holding c̃,c̃\*, Foreign-bond pos., r^k, spread fixed — Brent
5. **Foreign bond clearing**: solve i\*−i holding c̃,c̃\*, capital pos., r^k, i_t fixed — Brent
6. value functions (33,34,42,43) + θ transition (51) + aggregate savings (64); μ̃ via (N1)–(N3)
7. **terms of trade**: solve s_t so relative Home-good demand = relative supply — Brent
8. consumption-savings (37 via N1, 40,41,46,49,50): solve savings holding all but μ̃_{t+1}, c̃_t fixed
9. labor supply: union FOC (52,53)
10. inflation: Taylor rules (10)

### 3.1 Within-point solve structure — THREE NESTED LEVELS (critic B1/B2; do not under-design)

Vectorizing **across** `n_states` grid points is correct, but **within** each point the
equilibrium step is NOT closed-form — it is a nest of iterative solves. Get this structure
written down before coding `equilibrium_step.py`:

- **Level 1 (outer):** the backward iteration (§3), damped, to `diff < 1e-8`.
- **Level 2 (within-point root-finds, SEQUENTIAL):**
  - Step 4 (Home bond clearing) is **Brent-inside-Brent**: an outer Brent on `i_h` (residual
    `calc_excess_bond_nom`) that, per evaluation, runs an **inner Brent on each agent's
    capital share** (`calc_portfolio_share`). 
  - Step 5 (Foreign bond clearing) runs **after** step 4 with `i_h` held fixed: outer Brent on
    the spread `i*−i` (`calc_excess_foreign_bond_nom`) wrapping `calc_bond_portfolio_share`.
  - So steps 4→5 are ordered in the within-point dimension; each is batched over `n_states`.
- **Level 3 (within-point damped fixed-point loops):** steps 7, 8, 9 are each iterated to
  convergence with their own relaxation weights — **terms of trade s (weight 0.1), consumption
  Euler (0.5), wage Phillips curve (0.005, the slowest)**. Implement as **batched fixed-point
  iterations with a per-point convergence mask** (iterate until all points' residuals < tol).

Implementation rule: the 10 step-functions take and return a typed `EquilibriumState`
container (carrying shared scratch like the SDF matrix `M_vec_mat`, `r_alpha`, and the
θ-dependent stationarity nudge `bbeta_adj(n_I)` — critic N1/N4) so cross-step coupling is
explicit and the per-file functions never re-derive the SDF and drift out of sync.

### 3.2 Disaster / no-shock quadrature nodes (critic S2 — exact, state-dependent)

`n_quad = n_GH + 2`. The expectation weights are **state-dependent** through `p_dis`:
```
big_weight_vec[1:n_quad-1] = quad_weight_vec * (1 − p_dis)   # GH nodes + no-shock node
big_weight_vec[n_quad]     = p_dis                            # disaster node
```
The disaster node carries only `z = −disast_shock` (other shocks 0) and an extra
depreciation hit on the capital return (`r^k *= exp(dz)`). Forgetting the `(1−p_dis)` rescale
double-counts probability mass — a classic silent bug. `p_dis` varies by grid point (the
disaster-probability state), so this rescale is inside the batched expectation.

**Root-finders.** All four are **hand-rolled Brent** in the Fortran (Home rate, Foreign
spread, capital share, Foreign-bond share); a library Brent handles the ZLB/fixed-rate IRF
branch. Port as **batched Brent** (`solve/brent_batched`) solving all grid points at once
(vectorized bracketing + iteration to a tol matching Fortran), but **respect the nesting and
ordering in §3.1** — this is NOT one flat batch. The residual functions
(`calc_excess_bond_nom`, `calc_excess_foreign_bond_nom`, `portfolio_foc`,
`calc_portfolio_share`, `calc_bond_portfolio_share`) must match the scaled FOC definitions
exactly, including the `/(1−ω)` convenience wedge and natural-leverage bounds.

**Damping & staged activation are load-bearing for convergence** — replicate the exact
per-object relaxation weights and the iter-10 / iter-100 switches from `calc_sol`.

---

## 4. Calibration / "estimation" loop

There is **no outer SMM loop in the released solver** — the 9 `param_file_*.csv` are
pre-computed calibrations/counterfactuals. The mapping targets→parameters lives in
`create_param_files.m` (with `get_moments.m`, `sim_moments.m`, `additional_params.mat`).

**⚠ Dependency correction (critic B4).** `chi0_vec` (labor-disutility level) and the
steady-state terms-of-trade `s_ss` are **calibrated inside `calc_steady`** (a per-agent damped
fixed point), NOT in `create_param_files.m` — the Fortran emits `chi0(1..2)` to `extra_data.csv`
for the tables. So: **port `calc_steady` FIRST** (it owns the SS inversions), and treat
`create_param_files.py` as **best-effort, validated against the shipped CSVs** — which are the
ground-truth inputs. Do not block the solver on a perfect generator port. (This reverses the
naive ordering; see §9.)

Port plan (`create_param_files.py`):
- Reproduce each of the 9 specifications and write identical 65-column rows (validate byte-
  for-byte against the copied `param_file_*.csv` — they are our ground truth; the port is
  for reproducibility/extension, not to *replace* the inputs).
- Some parameters are set by **steady-state inversion** (e.g. β to hit an NFA/wealth-share
  target, ν̄ via `chi0_vec` to hit ℓ=1, T-bill level to hit `bg_yss`). Identify these in
  `create_param_files.m` and reproduce the inversion. `additional_params.mat` → `scipy.io.loadmat`.
- Emit `n_comp` (number of calibrations) as `create_param_files.m` does.

What each calibration is (to confirm when porting): file 1 = benchmark; others toggle
`foreign_trading`, the safety/disaster shocks, fiscal `run_bg`, sample runs `run_samp`,
output-gap Taylor terms, and `rho_i` (which activates the `i_h,i_f` grid dims). Document the
1→9 map in `create_param_files.py`'s header.

---

## 5. Vectorization / tensor-native strategy

**Where tensor-native discipline matters most** (these are the hot paths):
1. **Smolyak basis evaluation** — `B_new @ c` for all `n_states × n_quad` next-states is a
   single batched matmul. Never loop over grid points or quadrature nodes in Python.
2. **The equilibrium step over grid points** — the Fortran OMP-parallelizes over `n_states`;
   we **vectorize** it: every per-point quantity carries a leading `n_states` (and where
   relevant `n_quad`) batch dim. The 10 steps become tensor ops on `[n_states, ...]`.
3. **Batched Brent root-finders** — solve all grid points' scalar root problems at once:
   maintain `[n_states]` brackets, iterate Brent updates elementwise to a shared tol. No
   per-point scipy calls.
4. **Expectations** = weighted sums over `n_quad` → `einsum`/`tensordot` with `quad_weight_vec`.

**Discipline (from CLAUDE.md standing policy):** PyTorch, `dtype=float64` (match Fortran
`dp`), device-agnostic via `config.device` (CPU now, `cuda` later — no code change). **No
non-tensorized NumPy in hot loops; no CPU↔GPU round-trips mid-solve.** NumPy/scipy allowed
only at the edges (loading `.mat`, `create_param_files` setup, post-processing I/O).
Sanity checks after non-trivial ops (portfolio shares sum to 1; market-clearing residuals ≈ 0;
μ̃ > 0). Minimal try/except (let it crash — coding_conventions).

**Smolyak linear solve conditioning (critic B3).** The Fortran fits coefficients with NAG
`F07ABF` = LAPACK **`dgesvx`**, the *expert* driver with **equilibration + condition-number
estimation** — chosen because the Chebyshev-extrema basis can be ill-conditioned at high
anisotropic levels, and small errors compound over up to 5000 backward iterations. Port as
`torch.linalg.lu_factor` once per outer iteration + `lu_solve` for the multiple RHS, in
float64; **report the basis-matrix condition number at setup** (validation step 1) and add
equilibration (column scaling) if it is large. Do not silently use a plain `solve`.

**Why this keeps the GPU door open:** the entire solve is grid-batched tensor algebra +
batched root-finding + matmul interpolation — moves to H100 with only a `device` change.

---

## 6. Post-processing → reproduce figures/numbers

Port the MATLAB pipeline to Python (no MATLAB runtime). Detail in
`reference/matlab_postproc.md`.

- **`var_indices.py`** = exact port of `get_var_indices.m` (151 absolute row indices into the
  stacked `[shocks(4); states(9); vars(124); value(14)]` series). This is the keystone —
  get the column map exactly right. **Do not hardcode 124/9/151** (critic S7): read
  `n_interp`, `smolyak_d`, `n_active_dims` from the solver-emitted `grid.csv` and assert the
  map matches (`stopifnot`-style) at runtime.
- **Calibration-index map is already known (critic S6) — lock it now**, don't defer: benchmark
  `ix_bm`, plus `ix_no_omg`, `ix_symm`, … `ix_nocorr` (see `reference/matlab_postproc.md §0`).
  Post-processing selects **different calibration subsets per table** (`comp_idxs`: Table 3 →
  `[1,2,3]`, Table 9 → `[1,5,6]`, Table 10 → `[7,2,8,9,1]`). Encode this map in `tables.py`.
- **`read_series.py`/`extract_*`** — assemble series, apply the **de-trending convention**
  (`× exp(cumsum(z_shock))`, with `exp(dis)` return corrections and exact first-element
  padding after diffs/shifts). This is the single trickiest reproduction.
- **`moments.py`** — `calc_moments` (per-sim) → `collect_moments` (average across sims);
  `collect_swap_moments` for Table 8. Fills the "Model" columns of tables 1–10.
- **`tables.py`** — emit `table_1..10.tex` matching `Code/output/tables/`.

**Figures: STATA vs Python.** The project's `coding_conventions.md` says *figures that go in
a document must be made in Stata*. KL's figures are made in MATLAB. **Decision to confirm
with Angus (low priority, defer to step 6):** for this *replication* the figures are a
faithfulness check, not big_cy paper figures — so reproducing them in Python/matplotlib to
compare against KL's PDFs is the natural first pass. If any figure becomes a big_cy
deliverable, re-make it in Stata then. I will produce comparison plots in Python first and
flag this choice in the closeness note. (Flagged in `reference/matlab_postproc.md`.)

**Simulator details that change moments (critics S4/S5/N5):**
- **Three different state-clipping regimes** (must be carried into `simulate.py`): no-disaster
  sim clips next-state to `[−1,1]`; with-disaster sim does **not** clip; sample-path sim clips
  to `[−1+eps, 1−eps]`. Silent if wrong; changes the disaster-sim moments (Table 7).
- **Sample-path advance** (`run_samp`, needed for figs 4/6/7/10 + Table A1) uses NAG
  `e01zmf`/`e01znf` scattered-data interpolation over the active-shock node grid. Replace with
  `scipy.interpolate` (`RegularGridInterpolator` if the active-shock nodes form a regular
  product grid — confirm; else `griddata`). This is the only piece needed for the **headline
  fig 10** (ω-shock vs swapped-T-bill data), so it's not optional.
- **RNG is not bit-reproducible.** The Fortran draws via NAG `g05*` (SEED=712, GENID=3). We do
  not replicate that stream. So **simulated-moment agreement is statistical, not bit-exact**
  (n_sims=100 × 400 periods averages out). The **bit-comparable** targets are the deterministic
  objects: policy functions, steady state, and IRFs from fixed start states. Set this
  expectation in the closeness note.

**Comparison targets** (`Code/output/`): `figures/fig_2..20.pdf`, `tables/table_1..10.tex`,
`results.pdf`. Closeness note: per-figure / per-number deviation + likely cause.

---

## 7. Validation ladder (build confidence bottom-up before trusting model output)

1. **Smolyak basis** vs Fortran on a fixed `smol_elem` + random points (`test_smolyak_basis`);
   also report the basis condition number (critic B3).
2. **GH quadrature** integrates polynomials / normal moments exactly AND assert the nodes equal
   the eigenvalues of the **probabilists'** Hermite Jacobi matrix (off-diagonal `√i`, not
   `√(i/2)`) — guards the classic `hermegauss`-vs-`hermgauss` √2 footgun (critic S3). Moment
   tests alone won't catch a uniform node-scaling error.
3. **`grid.csv`/`grid_locs.csv`** vs Fortran (requires running Fortran once OR deriving from
   params; we derive from params and cross-check dims). ⚠ confirm `smolyak_d`, `n_interp`.
4. **Steady state** vs the Fortran SS print (`output_*.txt`) — needs the values; we compute
   from params and sanity-check against appendix SS relations.
5. **One full solve of benchmark (calibration 1)**; **proactively** cross-check the converged
   stochastic-steady-state values (capital, wealth share, wages, rates) and the **Table-2
   *targeted* moments** against KL's *published* numbers (critic S1) — a wrong sign in the (51)
   seigniorage term or a `q` vs `1/q` in (N3) makes the model converge to the *wrong* fixed
   point without crashing, so check early rather than waiting for a residual to misbehave.
6. **All 9 calibrations** → all figures/tables.

Because we are **not installing Fortran** (authorized last resort only), steps 1–2 are
self-validating (analytic), and 3–5 cross-check against (a) appendix relations and (b) KL's
*published* output files in `Code/output/`. If a discrepancy can't be resolved by reading,
flag to Angus before resorting to gfortran.

---

## 8. Runtime / feasibility (local-only)

KL: ~112 min for all 9 calibrations on a 32-core Xeon ⇒ ~12 min/calibration there, in
compiled Fortran. Our vectorized Python on CPU will likely be **slower per iteration** but
benefits from batching; expect the benchmark solve to take from tens of minutes to a few
hours on the Mac. **Policy:** if local runtime/RAM becomes prohibitive, **stop and flag
Angus to rethink** — do **not** move to Sherlock (CLAUDE.md standing policy). The
tensor-native design means the eventual fix is a GPU, not the cluster. We will first get
*correctness* on a coarse grid / few iterations, then scale.

---

## 9. Sequencing (next sessions)

1. Backend + Smolyak grid/basis + quadrature + their tests (validation ladder 1–2). ✔ gate.
2. `params.py` (load the shipped `param_file_*.csv` directly — the ground-truth inputs) +
   `init_setup` (covariance/Cholesky, shock+disaster nodes, Smolyak grid) + `grid_setup`.
3. `calc_steady` + steady-state checks (this owns the `chi0`/`s_ss` inversions — critic B4).
   *Then* (optional, best-effort) `create_param_files.py`, validated against the shipped CSVs.
4. `equilibrium_step` (the 10 steps, §3.1 nesting) + batched Brent + the FOC residuals — the bulk.
5. `solve_model` (damping + staged activation) → benchmark converges; proactive SS/moment
   cross-check (validation step 5) before scaling to the full grid.
6. Simulation (3 clipping regimes) + IRFs; `var_indices` (assert vs `grid.csv`) + post-processing
   (calib-subset map) + sample-path interp; figures/tables; closeness note.

**Coarse-grid timing gate (critic N2):** before committing to a full benchmark run (full
`vector_mus_dimensions` × up to 5000 iters), solve on a *reduced* anisotropic grid for a few
outer iterations to measure per-iteration wall-clock, extrapolate, and confirm local
feasibility — if prohibitive, stop and flag Angus (do not move to Sherlock).

Open items to confirm during implementation (all flagged in the digests): the handful of
`⚠ VERIFY` equation transcriptions (re-check vs appendix pp. 9–13/59–62 if a residual
won't zero); exact per-object damping weights; the `sampleK_mat` column layout; the
figures-language decision (§6).

---

## 13. Solver implementation notes (from reading `calc_sol`, mod_calc.f90:204-523)

Concrete details for porting the outer loop + equilibrium step (next session).

**Interpolated next-period variables: `n_nxt = 8 + 2·n_I = 12`.** The
`interp_input_mat` (n_states × 12) columns, in order, are:
`[v_h, v_f, mc_h, mc_f, s, q, l_aggr_h, l_aggr_f, infl_h, infl_f, c_spending_h, c_spending_f]`.
Each outer iteration: fit `smol_coeffs = solve(smol_polynom, interp_input_mat)`
(LAPACK `dgesvx` w/ equilibration → our `interp_factor`/`interp_solve`). Then per
grid point `sss`: `nxt_mat (n_quad×12) = Smolyak_basis(next_state_mat[:,:,sss]) @ smol_coeffs`,
and non-interpolated `nxt_mat_2 (n_quad×3) = [k_next, next_zf, next_omg]`.

**Outer loop:** `while diff > 1e-8 and outer_iter < max_iter(=5000)`. For each
state, call `calc_equilibrium_and_update(...)` → new policies/prices/values; then
damp + recompute `diff`. `M_vec_mat(n_quad, n_I, n_states)` carries the SDF across
steps (shared scratch — our `EquilibriumState`).

**Iteration-state initialization (mod_calc.f90:304-429):**
- `l_aggr=1`, `q=q_ss`, `s=s_ss`, `infl=1`, `q_bond=1`, `share=0`, `bF_share=0`
  (start dollar-only), `v=v_ss`, `mc=mc_ss`.
- `nom_i(1,:) = (1+rf_ss)/(1+exp(omg_state))`; `nom_i(2,:) = 1+rf_ss − nom_i(1,:)`
  (col 2 is the spread i*−i, NOT the level).
- `inv_h = δ·k·inv_share_h`, `inv_f = δ·k·(1−inv_share_h)`.
- `kappa` init via the s/zf/l ratio (lines 332-336); `w_choice` from kappa (341-344);
  `c_spending(iii) = w_choice·l_ss + k·wealth_share·(1/β−1)` (349-352).
- `k_next_mat(:,sss) = k·exp(dz_vec)/exp(dz_vec_adj)` → = k for all non-disaster
  nodes, `k·exp(−disast_shock)` at the disaster node (capital destruction).
- `next_state_mat(n_quad, n_active_dims, n_states)` = standardized next states:
  capital `(k/exp(dz_adj) − k_grid_mean)/k_grid_dev`; wealth-share dim 0 (guess);
  zf `(next_zf − mean)/dev`; wh/wf `(w_choice − mean)/dev`; dis `(next_dis − mean)/dev`
  (0 if dev<sqrt_eps); ih/if 0; omg `(next_omg − mean)/dev`. **Clamped to [−1,1].**

**The bulk still to port (exact Fortran line ranges):**
- `calc_equilibrium_and_update` — mod_calc.f90:2176-3099 (the 10 steps).
- `calc_unexpected_transition` — 1940-2174 (θ transition + seigniorage).
- Brent helpers: `calc_excess_bond_nom` 3776-3868, `calc_portfolio_share`
  3871-4082, `calc_excess_foreign_bond_nom` 4086-4175, `calc_bond_portfolio_share`
  4178-4360. `portfolio_foc`/`portfolio_return` 4388-4410 (done: small).
- `calc_bond_prices` 3101-3501, `calc_valuation` 3503-3774 (post-solve, n_bond ladder).
- ZLB/fixed-rate IRF branch: `calc_equilibrium_ifixed` 4424-5175, `update_ifixed`
  5178-5324.

These are tightly coupled and only testable once the loop converges; port together,
then gate on the proactive SS/Table-2 cross-check (§7 step 5).

### 13.7 RESUME NOTE — Table-2 moments pipeline BUILT + submitted (session 8)

**DONE + committed (4087d92, pushed origin/main, 53 tests green):**
- `model/compute_results_vec.py` — grid-level Table-2 columns from converged `g`
  (reuses `compute_current_period`): yh,yf,lh,lf,ch,cf, inv=k_next_new-(1-δ)k,
  pi, P_Phh,P_Phf,qx,s,q=q/P1, h_kap=kappa_h, h_ksav, h_sav=savings_h/P1,
  infl_h/f, nom_ih/if. `level_mask` flags the de-trended levels.
- `post/table2_series.py` — `build_table2_series(const,g,sim,disast_shock)`: fits
  Smolyak coeffs to results_grid, interpolates at `sim['std_series']`, de-trends
  levels by `*exp(cumsum z_shock)`, builds dis/rfh/nfa/nfa_rel. (Matches
  extract_series.m exactly: rfh=log(nom_ih_{t-1}/infl_h_t) prepend-first; nfa=
  h_sav_t - h_kap_{t+1}·exp(-dis_{t+1})·q_t duplicate-LAST.)
- `post/moments.py` — `compute_table2_moments` (12 moments {1-7,11-15}, avg over
  sims) + `format_table2`/KL_TARGETS. Verified against calc_moments.m lines 33-47.
- `run_klrep_moments.py` — cluster entry: rebuild const from param_file_N, load
  solution_spec<N>.pt, reassemble SolverState, simulate (n_burn=10000, n_sims=100,
  n_per=400, no-disaster), print KL Table-2 comparison to stdout (.out log).
- `batch_submit.sh` — `run_klrep_moments` GPU block (depends-afterok on run_klrep
  if chained; else loads existing solution). Committed default `to_run` set to it.
- **Local validation:** coarse smoke (full pipeline) + reload-path smoke both OK.
  Even coarse/unconverged, moms 1,3,5,11,12,13,15 ≈ KL; mom 4/6/7 are the
  grid-resolution/convergence-sensitive ones (nfa, mean rate, autocorr).

**SUBMITTED:** job **29550357_1** (specs=1, maggiori GPU) — RUNNING. When done:
sync `logs/klreplication`, read `run_klrep_moments_29550357_1.out`, compare the
12 moments to KL table_2.tex. **GATE: align → closeness note + Phase-2 (specs=1-9,
KLREP_CONV=1e-8, to_run=" run_klrep "); material discrepancy → STOP + flag Angus.**

**Bond ladder (Task 5) — DONE + committed (000eb2a, 57 tests green).**
`model/bond_ladder.py`: `compute_pricing_block` (SDF M_vec + hg_infl/P/P_nxt/omg
from converged g; reproduces equilibrium_step M_vec, validated) + `compute_bond_columns`
(recursion q_bond[:,:,1]=1; for m=1..20 fit Smolyak coeffs → interpolate at
next_state → `_price_one_step` = maxval_agent Σ big_w·M·(nxt_bp/hg_infl)·P/P_nxt,
hw adds /(1-omg); store maturities [1,2,3,4,19,20] → q1_h/f/hw, q19/q20 h/f).
`post/table2_series._build_bond_series`: yields, rq20_h/fh, rfh_lev (1/(1+ζ) blend),
rk (disaster-corrected), rA=log((1+d2e)e^rk−d2e·e^rfh_lev) [d2e=1], exc_retA,
exc_rf (rff_h−rfh), E_change, uip_pvt, fama_yield_pvt, y_growth(4q),
div_price_smoothed_1 (movsum[3,0]). `post/moments._bond_moments_per_sim`: mom8=
corr(resid_dp,resid_carry) via R1/R2 batched OLS, mom9=4·100·(mean rA−mean rfh),
mom10=R3 slope on exc_retA. run_klrep_moments prices the ladder + passes bond_grid.
Coarse smoke: all 15 moments finite; 8/9/10 in KL ballpark. Converged full-15
cluster run submitted (29555360).

### 13.8 RESUME NOTE — moment tables validated; figures + T6/7/8 remain (session 8)

**DONE + committed (d2b1ceb, 58 tests green):** all 9 specs solved (1e-8); full
Table-2 (14/15) + Tables 3/4/5/9/10 validated across specs vs KL — see
`CLOSENESS_NOTE.md`. Entries: run_klrep (solve), run_klrep_moments (Table-2 + NFA
forensics), run_klrep_tables (Tables 3/4/5/9/10, loops all specs). Code:
model/{compute_results_vec,bond_ladder}, post/{table2_series,moments}.

**REMAINING (next chunk — each needs NEW machinery; port from the read-only Fortran/
MATLAB; do in committable sub-milestones; gate each on KL .tex):**
1. **Table-10 corr rows 4-5** (SMALL): add `corr_rf_spread_m1B/m2B` (results_vec
   49,50) to a grid-results path — the pricing block already has M_vec + the rf
   spread; port the corr_rf_spread block (calc_bond_prices:3406-3443, the
   quad_weight_vec-only branch → results_vec(49)=corr(5), (50)=corr(6)). Mean over sim.
2. **Table 6** (MODERATE): net-exports `nx` series (extract_series:409-414) — needs
   chf/cfh (consumption splits via varsigma/s), inv_h (results_vec 38, ADD column),
   deployed-capital flow (h_kap, q, dis). Then table_6_moms 1-6 (nfa_rel_growth +
   nx vol/means). nfa_rel_growth already built.
3. **Table 7 / A2** (SUBSTANTIAL): port `calc_valuation` (mod_calc.f90:3503-3774,
   the n_bond valuation columns: nfa/Enfa/val/Eval/... → jx 138-151) + run the
   DISASTER ensemble (simulate_ensemble disaster=True already exists). table_7 is the
   `dis` branch of calc_moments (innovations = temp(2:end,jx_val_*) - temp(1:end-1,
   jx_val_E*); covariances with nfa innovations). Needs var_indices jx_val_* map.
4. **Table 8** (SUBSTANTIAL): swap-line MIT experiment (`collect_swap_moments.m` +
   the swap IRF in mod_results) — a counterfactual policy experiment, NOT the moment
   ensemble. Lower priority.
5. **Figures 2-20** (LARGEST): generalized IRF port (MIT shocks z/ω/p/zf/bg/fix from
   mod_results create_results IRF section) → extract_irfs.m → matplotlib (Python first
   pass per §6; Stata only if a fig becomes a big_cy deliverable). Figs 4/6/7/10 need
   the empirical SAMPLE-PATH sim (sample{1,2,3}_mat.txt + NAG e01zmf/e01znf →
   scipy RegularGridInterpolator/griddata). Headline = fig 10 (ω-shock vs swapped-Tbill).

**FLAGS:** (a) pulling .pt/figures to local needs SYNC_ALLOWED_big_cy=(data/output/
klreplication) in sherlock-agent-com (Angus). (b) The benchmark NFA-level residual
(−15.6 vs −23) is benign/localized — do NOT re-chase without KL's policy dump.

### 13.6 RESUME NOTE — simulation/moments for KL Table-2 (session 7, in progress)

**Goal (Angus): port simulation/moments → validate KL Table-2 on the solved
spec-1 → THEN Phase-2.** Do NOT launch the 9-spec array until Table-2 validates.

**DONE + committed:**
- `reference/table2_var_spec.md` — exact definition of every Table-2 series
  (column map jx=13+results_col, calc_bond_prices formula, extract_series
  de-trending). KEY: Table-2 uses REALIZED rates (consecutive periods); only
  `nfa` + `uip_pvt` use next-period info; only **moment 8** (div-price/carry corr)
  needs the recursive bond ladder; de-trend levels by `× exp(cumsum(z_shock))`.
- `src/klrep/simulate/simulate.py` (committed d955e81) — the SIMULATION ENGINE:
  `build_sim_coeffs` (fits Smolyak coeffs to the solved per-node transition,
  9 economic states, and 14 policies [s,q,l_h,l_f,c_h,c_f,nom_ih,spread,infl_h,
  infl_f,share_h,share_f,bF_h,bF_f]); `stochastic_ss`; `burn_in` (disaster ON);
  `simulate_ensemble` (no-disaster: product nodes only + clip [-1,1]; disaster:
  big_weight_vec + no clip). Validated on a coarse local solve (stoch-SS
  converges, states in bounds, policies sane). RNG ≠ NAG → moments statistical.

**⚠ REFINEMENTS discovered session 7b (read before coding):**
- **Moments 8, 9, AND 10 need the bond ladder** (not just 8). `rA` (levered equity
  return) = `(1+d2e)·exp(rk) − d2e·exp(rfh_lev)` where `rfh_lev` blends the
  **20-period** bond realized returns (`rq20_h/rq20_fh` from `q19/q20`). So
  `exc_retA` (→ moments 8,9,10) needs the bond ladder. **Deliverable WITHOUT bonds:
  the 12 moments {1,2,3,4,5,6,7,11,12,13,14,15}.** Defer 8/9/10 until calc_bond_prices
  (recursive bond pricing over maturities [1,2,3,4,19,20]) is ported.
- **Architecture decision (do it this way):** compute `results_vec` ON THE GRID via a
  new `compute_results_vec(const, g)` (per grid point: y_h,y_f,l,c_h,c_f, **inv =
  k_next_new−(1−δ)k**, s, q_real=q/P1, pi=pi_cur0/P1, P1,P2, qx=P1/P2, kappa_h,
  savings_h, h_sav=savings_h/P1, h_kap=kappa_h, nom_ih, nom_if, infl_h, infl_f) —
  reuse period_block + the converged g (k_next_new = Σ savings·(1−share)/q is clean
  on the grid). Then fit Smolyak coeffs and **interpolate along the simulated
  STANDARDIZED state path** (`sim['std_series']`, now stored). This matches KL
  (interpolate results_mat) and AVOIDS the de-trending bugs of computing inv/nfa
  from consecutive economic states. `nfa_t = h_sav_t − h_kap_{t+1}·exp(−dis_{t+1})·
  q_t` uses the next period's interpolated h_kap. De-trend levels by
  `× exp(cumsum(z_shock))`; rfh = log(nom_ih_{t-1}/infl_h_t) realized.

**REMAINING (next session) — pieces:**
0. `src/klrep/model/compute_results_vec.py` — grid-level Table-2 vars (see refinement
   above). Then in post: fit coeffs (interp_factor/solve on smol_polynom) +
   interpolate at `sim['std_series']` (basis@coeffs) → de-trended series.
1. (superseded by §0 + interpolation) `src/klrep/post/results_path.py` — from `sim` output (state_series 9, pol_series
   14, z_shock_series) compute the Table-2 economic vars per `table2_var_spec.md`:
   production (yh,yf,c_h,c_f,inv,w,pi) via the period_block formulas; qx =
   P_div_P_h(1)/P_div_P_h(2); portfolio (savings = wealth+w·l−c, h_kap =
   sav·(1−share)/q, h_bh_sav, h_bf_sav, h_sav); **nfa = h_sav_t −
   h_kap_{t+1}·exp(−dis_{t+1})·q_t**; realized returns from consecutive periods
   (rfh = home real bond return = (nom_ih_{t-1})/infl_h_t·…; rff_h; rA = capital
   return (pi_t+(1−δ)q_t)/q_{t-1}·exp corrections; exc_retA = rA−rfh; exc_rf);
   yields (1-period ≈ nominal); y_growth; infl_h/f. De-trend levels by
   exp(cumsum(z_shock)). (Reuse model/period_block.py + portfolio logic.)
2. `src/klrep/post/moments.py` — the 15 Table-2 moments + 3 regressions (see
   calc_moments.m, transcribed in table2_var_spec.md), averaged over n_sims
   (collect_moments). **Moment 8 needs div_price_smoothed (20-period bond ladder
   = calc_bond_prices port) — defer/flag; deliver the other 14 first.**
3. Cluster entry `run_klrep_moments.py` (or extend run_klrep): load
   `solution_spec1.pt` from oak, rebuild `const` from load_param_file(1) +
   calc_steady + build_* (deterministic), assemble SolverState, build_sim_coeffs,
   simulate (n_burn=10000, n_sims=100, n_sim_periods=400, both disaster regimes),
   compute Table-2, PRINT the table to stdout (so it returns via the synced .out
   log — no SYNC_ALLOWED needed). Add a `run_klrep_moments` block to batch_submit
   (CPU or GPU). cc-merge → sync_git → submit → read .out → compare to KL Table-2
   "Model" column (table_2.tex) → write Table-2 closeness note.
   **GATE: if Table-2 aligns → Phase-2 (specs=1-9, KLREP_CONV=1e-8, 7h). If a
   MATERIAL discrepancy → STOP + flag Angus before spending GPU on all 9.**

KL Table-2 "Model" targets (from table_2.tex, the bar to hit): ζ* 1.6; σ(Δlogc)
0.5%; σ(Δlogy*) 0.8%; ρ(y*/y) 0.5; σ(Δlogx) 1.6%; 4·Er 2.0%; nfa/4y −23%;
ρ_{-1}(r^e, r*+Δq−r) 0.5; 4·E[r^e−r] 5.2%; β(Δnfa/y, r^e−r) 0.6; b*_Hs/4y 3.8%;
ℓ 1.0; ℓ* 1.0; logP/P_{-1} 0.0%; logP*/P*_{-1} 0.0%.

### 13.5 PHASE-1 GATE PASSED — baseline solved on Sherlock GPU (2026-06-13)

`run_klrep` job 29449338_1 (spec 1, benchmark) on an **A100-80GB**: full grid
(589 states) converged to **diff 9.94e-9 in 1733 iters, 3h09m**, on `cuda:0`
(float64). Central node vs det. SS: q/s/ℓ at SS; v_h 0.826, 1+i 1.00697 (both
lower = disaster-risk pricing); **share_h −0.383 (equity bias)**. Baseline ALIGNS
qualitatively + at SS level; the GPU pipeline works end-to-end. **Full Table-2
moment comparison still needs the simulation port.** Details + table:
`klreplication/PHASE1_BASELINE.md`. STOPPED at the gate per dispatch (awaiting
Angus before the 9-spec array). Infra fixes en route: 16G→96G OOM; CPU→GPU node;
conda-solver grind→pip torch; conda-not-found→bootstrap+py3-assert; 4h→7h wall;
configurable KLREP_CONV (1e-6 gate / 1e-8 final).

### 13.4 RESUME NOTE — Sherlock migration (04b), state as of session 5

**DONE + pushed (origin/main 0c4a16b):**
- `klreplication/run_klrep.py` — Sherlock entry: array index N → param_file_N,
  full-grid solve (auto-CUDA) + SS validation, writes solution_spec<N>.pt +
  ss_check_spec<N>.txt to $big_cy/data/output/klreplication. Tested locally
  (coarse, 8 iters) — works end-to-end.
- `klreplication/inputs/` — all KL inputs committed (88K: params/sample/data) so
  the cluster gets them via sync_git (no local→oak data push needed). run_klrep
  reads params from here.
- `batch_controller.sh` — `.py` branch env fixed to `big_cy_klrep_env`.
- `batch_submit.sh` — `specs="${specs:-1}"` + `setup_klrep_env_sherlock` (env
  build, CPU) + `run_klrep` (GPU array `--array=${specs}`, --gpus=1) blocks.
- `setup_klrep_env_sherlock.sh` — now bootstraps conda from the group/home
  miniconda profile (compute nodes lack conda on PATH).
- Sherlock connectivity confirmed (squeue works).

**⚠ BLOCKER (flagged to Angus) — `sherlock-agent-com` doesn't know `big_cy`:**
`VALID_PROJECTS=("lx1" "rier" "ecb1")` only (line 13 of
~/code/claude_scripts/sherlock-agent-com). To submit, Angus must register big_cy:
add "big_cy" to VALID_PROJECTS; add a project-root case (~line 124-160) →
`PROJECT_REMOTE_ROOT="/oak/stanford/groups/maggiori/Lab/lewis"` (repo-name
appended ⇒ /oak/.../lewis/big_cy, matching project_strings.py); ensure a cluster
git checkout of big_cy exists at the path submit/sync_git expect AND at
batch_controller's `base_folder=$HOME/big_cy`; optionally
`SYNC_ALLOWED_big_cy=(data/output/klreplication)` to sync results back. This is
shared tooling + cluster setup = Angus's domain (convention: ask before using a
non-listed project).

**NEXT once unblocked:** `sherlock-agent-com big_cy sync_git` →
`sherlock-agent-com big_cy submit "setup_klrep_env_sherlock run_klrep"` (specs
defaults to 1 = baseline) → squeue/sacct → read ss_check_spec1.txt, compare to
KL SS/Table-2 → Phase-1 gate report to Angus. THEN port simulation (mod_results,
3 clipping regimes) + var_indices/moments for the full Table-2 moment comparison
(SS-level validation is in run_klrep already; Table-2 moments need simulation).
calc_bond_prices/calc_valuation needed for the NFA/valuation table rows.

### 13.3 RESUME NOTE — solver assembled + running (state as of session 4)

**DONE this session:** `model_const.py` (ModelConst bundle + big_weight),
`model/equilibrium_step.py` (the full 10-step batched equilibrium update incl.
nested batched Brent bond clearing), `solve/solve_model.py` (outer fit/eval/damp/
rebuild loop with the §13.2 damping + staged activation). **The solver RUNS and
CONVERGES stably on a coarse grid** (mus_dims=[1]*6+[0,0,1], 15 states): diff
decreases monotonically, bond clearing activates at iter 11 (shares 0→0.7, finite),
foreign trading at 100; q≈1.004/s≈0.99 sit at SS. 49 tests pass + smoke test.

**Timing gate (coarse-grid):** coarse 0.7s/iter (15 states); full grid (589
states) 3.1s/iter pre-bond-clearing (more once Brent active). Convergence to 1e-8
needs ~thousands of dampened iters → benchmark solve ~hours on CPU, all 9 calibs
~a day. Tensor-native ⇒ GPU later cuts this to minutes. Flagged to Angus.

**COARSE-GRID CONVERGENCE VALIDATION (800 iters, ~15 min):** diff fell
1.2e-2(it100) → 8.4e-6(it800) monotonically — the dampened fixed point converges,
no stalling/blowup. Central node vs deterministic SS: q 1.0007/1.0051, s
0.9986/0.9918, l_h 1.002/1.000 (real quantities AT SS); v_h 0.81 vs 1.0, nom_i
1.003 vs 1.011 (LOWER — correct: the global solution prices disaster risk that the
deterministic SS ignores; with γ=21 the safe rate + welfare fall), share_h = -0.51
(home leveraged into capital = equity bias, the model's central mechanism), bF_h
-0.11. These risk-adjusted deviations are economically RIGHT, not bugs. Definitive
gate = simulated moments vs KL tables (needs post-proc + full grid). Late
convergence rate ⇒ ~1300 coarse iters to 1e-8; full-grid benchmark ≈ few hours.

**NEXT (remaining):** (a) run benchmark to convergence (background/next session;
consider a 1e-5/1e-6 tolerance for a first pass) + finish the SS/Table-2 gate;
(b) port calc_bond_prices (3101-3501) + calc_valuation (3503-3774) — the post-solve
n_bond ladder + valuation terms (NFA decomposition etc.); (c) simulate (3 clipping
regimes, mod_results.f90) + generalized IRFs; (d) post-processing
(var_indices/read_series + de-trending, moments, figures, tables) + closeness note.
A `run/solve_calibration.py` entry point + saving the solved SolverState to disk
would let the long solve run once and post-processing iterate cheaply.

Possible refinements to check if convergence stalls: (i) the foreign-bond leverage
bound sign (calc_excess_foreign_bond_nom temp_r branch); (ii) whether the
bracket_expand step sizes (1e-2 home, 1e-3 spread) need widening for some states;
(iii) the w_choice log-damp only affects diagnostics here (w recomputed each step).

### 13.1 RESUME NOTE — equilibrium-step assembly (state as of session 3)

**Built + tested (47 tests green)** under `src/klrep/`: config; grid (Smolyak
elem/grid/basis/interp); quad (GH); params; setup (shock_grid, smolyak_setup,
state_grid); solve/steady_state; model/{calc_bundle, util_fun, portfolio
(portfolio_return + portfolio_foc), period_block (compute_current_period,
compute_q_new)}; solve/brent (batched_brent + batched_bracket_expand). The whole
SETUP pipeline + the equilibrium step's static algebra + the root-find machinery
are validated. **All the LEGO pieces for the equilibrium step now exist and are
tested.**

**Remaining = assemble `model/equilibrium_step.py` + `solve/solve_model.py`** from
those pieces, following the exact Fortran (now fully read; key facts below).

**Equilibrium step data shapes** (batch S=n_states, Q=n_quad, n_I=2):
- `big_weight_vec` (S,Q): `[:, :Q-1] = quad_weight_vec*(1-p_dis)`, `[:, Q-1] = p_dis`,
  `p_dis = exp(state_grid[:,IDX_DIS])`. NOTE the no-shock node already has
  quad_weight 0, so it gets weight 0 in expectations (it's only for the stochastic
  SS in simulation). Sum over Q = 1.
- STEP 2 next-period (mod_calc.f90:2316-2396), all (S,Q,·): unpack nxt_mat columns
  `[v_h,v_f,mc_h,mc_f,s,q,l1,l2,infl1,infl2,csp1,csp2]`; `ch/cf_spending_nxt =
  nxt_mat[...,10/11]*exp(dz_vec)`; `k_nxt=nxt_mat_2[...,0]`, `next_k=k_nxt[:,0]`;
  `zf_nxt2=zeta*exp(nxt_mat_2[...,1])`; `omg_nxt=exp(nxt_mat_2[...,2])+omg_shift-
  b_lmbd*shock_nxt`; `homegood_infl_nxt` (2349-2350); `kappa_nxt`; `y_next`,
  `w_next_choice`, `pi_nxt` (2361-2372); **`rk_vec=((1-δ)q_nxt+pi_nxt)/q_current`
  then `rk_vec[:,Q-1]*=exp(dz_vec[Q-1])`** (disaster node); `P_div_P_h_nxt`;
  `rf_vec[...,0]=nom_i1/homegood_infl_nxt[...,0]`, `rf_vec[...,1]=(nom_i1+nom_i2)/
  homegood_infl_nxt[...,1]`; seigniorage (2390-2396): if bg_yss>0,
  `seig=cf_spending_nxt*(bg_yss+shock_nxt)*omg_nxt` (the (ch+cf) cancels), transfer
  `[+seig,-seig]`.
- STEP 4 home bond clearing (outer_iter>10): batched Brent over S on `nom_i[:,0]`
  to zero `excess_b = sum_agent savings*share`. Residual `calc_excess_bond_nom`:
  per agent `savings=wealth+w_choice*l-c_spend` (c_spend clamped to
  [min_cons_sav, wealth+w*l-min_cons_sav], min_cons_sav=1e-8); natural-leverage
  bounds (3826-3842, note numerator uses RAW rf_vec cols, denom uses
  rf_home=rf_vec0/(1-omg)); `share` via batched Brent on `portfolio_foc` with
  r1=rf_home, r2=rk, `next_period_share=(savings*r_alpha+exp(dz_vec)*q_l_ss+seig)/
  tot_wealth_ss`, `r_alpha=portfolio_return(share,bF_share,rf_home,rf_foreign,rk)`;
  constraint-binding -> clamp to bound. `b_temp=savings*share`.
- STEP 5 foreign clearing (foreign_trading & outer_iter>100): batched Brent on the
  SPREAD `nom_i[:,1]` to zero `-sum bh_temp`; inner `calc_bond_portfolio_share`
  solves `bF_share` (Brent on portfolio_foc, r1=rf_foreign, r2=rf_home);
  `bf_temp=savings*bF_share`, `bh_temp=(share-bF_share)*savings`.
- STEP 6 value/θ (2818-2867): per agent `r_alpha` (with /(1-omg) on home),
  `next_period_share`, `EV=(Σ w·(v_temp·nps)^(1-γ))^(1/(1-γ))`, `labor_part`,
  `objf=(v_norm·labor_part^(1/ies)·c^((ies-1)/ies)+β·EV^((ies-1)/ies))^(ies/(ies-1))`,
  `v_new=objf·((wealth+q_l_ss)/tot_wealth_ss)^(-1)`; `θ_nxt=nxt_wealth1/Σ`;
  `k_next_new=Σ savings*(1-share)/q_current`. β includes `bbeta_adj`.
- STEP 7 terms of trade (2873-2903): damped (0.1) fixed point on s_new.
- STEP 8 consumption (2909-3000): damped (0.5) Euler per agent; updates
  `M_vec_mat[:,iii]=β·mc_temp/util_c_deriv·(v_temp/EV)^(1/ies-γ)·nps^(-γ)`,
  `cons_update=c/temp^ies`, `temp=Σ M·r_alpha_omg·c_cost/c_cost_nxt·w`; if phi_w==0
  also update labor (2974-2978).
- STEP 9 labor (3006-3064): if phi_w>0, damped (0.005) wage Phillips curve
  (3024-3037) then `l=kappa·((zf^(1-α)(1-α)/w/s_vec)^(1/α))`.
- STEP 10 inflation (3072-3074): Taylor-rule inversion using `nom_i_vec_in`
  (the INCOMING nom_i, saved at entry) + `ih_last/if_last`, rho_i, phi_h/f, phi_yh/f,
  tayl_ic.
- outer_iter<=10 fallback (3079-3097): share=0, nom_i from M_vec_mat averages.
- Outer loop `solve_model` (PLAN §13): fit coeffs (interp_factor/solve on
  smol_polynom), per-iteration build nxt via Smolyak basis @ coeffs, call
  equilibrium_step (batched over ALL states at once), damp each object, rebuild
  next_state_mat (clamp [-1,1]), recompute diff, until <1e-8 or max_iter=5000.
  Per-object damping weights: read from the update section of calc_sol (lines
  ~1100-1937 — not yet transcribed; CHECK there for the exact relaxation weights).
- Then `calc_bond_prices` (3101-3501), `calc_valuation` (3503-3774), then simulate
  (3 clipping regimes) + IRFs + post-proc.

**Next action:** transcribe the per-object damping weights from calc_sol's update
block (mod_calc.f90 ~1100-1937), then write equilibrium_step.py + solve_model.py,
run a COARSE-grid solve (reduce mus_dims), gate on SS/Table-2, then scale.

### 13.2 Damping weights + next_state rebuild + convergence (mod_calc.f90:524-709)

**Per-object damping** `x = x + w*(x_new - x)` (line refs 676-704):
v 1.0, mc 1.0, next_state_mat 0.1, q 0.2, s 0.1, l_aggr 0.2, c_spending 0.2,
nom_i(1) 0.1, nom_i(2) 0.1, infl 0.05, k_next 0.2, share 0.1; bF_share 0.0 for
outer_iter<101 else 0.05. **w_choice is LOG-damped**: `w = w*(1 + 0.2*log(w_new/w))`.

**next_state_mat_new** (per active dim, then clamp [-1,1]); exogenous dims
zf/dis/omg are NOT recomputed (kept from init, so damping leaves them fixed):
- k:  `(k_next_new/exp(dz_vec_adj) - k_grid_mean)/k_grid_dev`  (per quad node)
- θ:  `(theta_nxt - tht_h_grid_mean)/tht_h_grid_dev`
- wh: `(w_choice_new[0]/exp(dz_vec_adj) / (varsigma1+(1-varsigma1)s^(σ-1))^(1/(σ-1)) - wh_grid_mean)/wh_grid_dev`
- wf: `(w_choice_new[1]/exp(dz_vec_adj)*s / (varsigma2 s^(1-σ)+(1-varsigma2))^(1/(σ-1)) - wf_grid_mean)/wf_grid_dev`
- ih: `((nom_i[0]-1) - ih_grid_mean)/ih_grid_dev`
- if: `((nom_i[0]+nom_i[1]-1) - if_grid_mean)/if_grid_dev`
Also `k_next_mat_new = k_next_new/exp(dz_vec_adj)*exp(dz_vec)` (the next-iter k_nxt:
= k_next_new for non-disaster nodes, *exp(-disast_shock) at the disaster node).

**Convergence diff** (line 604): max over [abs log-ratio of v, mc, c_spending,
k_next, l_aggr, w_choice, q, infl ; abs diff of share, nom_i ; abs diff
next_state]. Stop when diff < 1e-8 (conv_crit) or outer_iter == max_iter (5000).

**Share solve bounds (calc_excess_bond_nom):** `temp = (-rk -
bF_share*(rf_vec1 - rf_vec0)) / (rf_home - rk + sqrt_eps)` with rf_home =
rf_vec0/(1-omg); share_low = max(-10, max_q{temp : temp<0}+eps); share_high =
min(5, min_q{temp : temp>0}-eps). If foc(share_low)<=0 -> share_low; elif
foc(share_high)>=0 -> share_high; else Brent. b_temp = savings*share;
excess_b = sum_agent b_temp. (Foreign analog uses safe bounds [-10,10] and
temp_r = rf_vec1 - rf_vec0/(1-omg); bh_temp = (share-bF_share)*savings;
excess = -sum bh_temp.) savings = wealth + w_choice*l - c_spend (indep of nom_i).
