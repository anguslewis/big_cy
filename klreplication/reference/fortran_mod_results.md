# Fortran `mod_results.f90` — simulation/output module digest

Source (read-only):
`…/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/Code/src/fortran/mod_results.f90` (1186 lines)

This module's single public entry point, `create_results`, takes an already-solved
model (policy/value matrices on a 9-dim Smolyak grid, written to disk by the
solver) and produces (a) a stochastic-steady-state report, (b) simulated panels
for business-cycle moments, (c) generalized impulse responses, (d) optional
long-lived "ll" IRFs, and (e) optional interpolated empirical-sample paths. All
output is written as plain-text `.txt` files into `results_folder`
(`../output/tmp/res_<i>/`, where `<i>` is the run index passed as command-line
arg 1). The MATLAB post-processing reads these files and computes the reported
moments.

---

## 0. Key parameters (from `mod_param.f90`)

| Symbol | Value | Meaning |
|---|---|---|
| `smolyak_d` | 9 | state-space dimension (`mod_param.f90:42`) |
| `n_I` | 2 | number of agents/countries (`:42`) |
| `n_interp` | 52+6*12 = 124 | number of interpolated "other variables" (`:42`) |
| `n_shocks` | 4 | normal shocks: z, zf, p, ω (`:55`) |
| `n_sample` | 100 | # IRF / sample starting points (`:43`) |
| `n_sims` | 100 | # simulation paths in the moment ensemble (`:43`) |
| `n_irf_periods` | 200 | IRF horizon (`:43`) |
| `n_sim_periods` | 400 | length of each simulation path (`:43`) |
| `n_burn` | 10000 | burn-in length (`:43`) |
| `n_sample_periods` | 94 | length of empirical-sample path (`:43`) |
| `n_valuation` | 500 | (used only in solver; not used here) (`:43`) |
| `sidx_z, sidx_zf, sidx_p, sidx_omg, sidx_dis` | 1,2,3,4,5 | shock indices (`:49`) |
| `sidx_bg` | 1 | (bond-supply ll index alias) (`:50`) |
| `n_irf` | 6 | # of MIT-shock IRFs (`:49`) |
| `irf_indices` | `[0,1,2,3,4,5]` | which shock each IRF perturbs (0 = baseline/no shock) (`:51`) |
| `irf_indices_ll` | `[1,1,6]` | long-lived IRF shock indices (`:51`) |
| `n_irf_ll` | 3 | # of long-lived IRFs (`:50`) |
| `irf_ll_length` | 40 | # of held-shock periods in ll IRFs (`:52`) |
| `irf_ll_hold_fixed` | `[0,1,1]` | per-ll-IRF "hold fixed" flag (`:52`) |
| `n_quad` | `n_GH_points + 2` | total quad nodes (GH + no-shock + disaster) (`:265`) |
| `no_shock_idx` | `n_quad-1` (scalar in `(1)`) | index of the zero-shock node (`:425`) |
| disaster node | `n_quad` | `shock_grid(n_quad,1) = -disast_shock`, rest 0 (`:422`) |
| `run_bg`, `run_samp` | 0/1 from param file 64/65 | toggle ll-IRFs and sample paths (`:241-242`) |

`irf_shock_sizes = [0, −5·σ_z, −2·σ_zf, +2·σ_dis, +2·σ_omg, 0]` (`mod_param.f90:254`).
`irf_ll_shock_sizes = [−σ_omg, −σ_omg, 0]` (`:257`).

### shock_grid layout (`mod_param.f90:389-425`)
- Rows `1 … n_quad-2`: the `n_GH_points` Gauss–Hermite quadrature nodes. Each row
  `i` holds the 4 correlated shock realizations
  `shock_grid(i,:) = cov_mat · quad_vec_temp` where `quad_vec_temp =
  [z, zf, p, omg]` standardized nodes (Cholesky `cov_mat`, `:404`). Built by a
  nested loop over `z × zf × p × omg` (`:391-414`), counter-major.
- Row `n_quad-1` (= `no_shock_idx`): all zeros — the deterministic/no-shock node,
  weight 0 (`:418-419`).
- Row `n_quad`: disaster node, only the z-component = `−disast_shock` (`:422`).
- `quad_weight_vec` has length `n_GH_points+1`; entries `1..n_GH_points` are the
  GH product weights, entry `n_GH_points+1`(= the no-shock slot) is 0 (`:407,418`).

---

## 1. Subroutine / function inventory

| Routine | Signature | Purpose | Lines |
|---|---|---|---|
| `create_results` | `subroutine create_results()` (no args; run index via `get_command_argument(1,…)`) | sole public routine; does everything below | `:20–1182` |

No other module-level routines. External library calls used:
- `Smolyak_Polynomial` (from `mod_smolyak`) — evaluate Smolyak basis at a state (`:22`, many call sites).
- `F07ABF` (LAPACK `dgesvx`) — solve linear systems to get polynomial coefficients (`:151,167,…`).
- `DGEMM` (BLAS) — basis × coeffs → interpolated values (`:198,…`).
- NAG RNG: `G05KFF` (seed), `g05saf` (uniform draws), `g05ndf` (sample integers w/o replacement) (`:237,246,304,315,481`).
- NAG interpolation: `e01zmf` / `e01znf` — multidim scattered-data interpolation, used only in the sample-path block (`:1139,1141`).

`create_results` is one long monolithic routine. Logical phases:
| Phase | Lines | Description |
|---|---|---|
| Setup / write grid metadata | `:95–119` | parse run index, write `grid.csv`, `grid_locs.csv` |
| Read solution files | `:124–135` | `next_state_mat.dat`, `results_mat.dat`, `valuation_mat.dat` |
| Build interpolation coeffs | `:140–187` | transition, other-vars, state, valuation coeffs |
| Stochastic steady state | `:189–227` | fixed-point iteration at no-shock node; prints SS |
| Burn-in | `:233–289` | 2·n_burn periods, collect last n_burn states |
| Sample starting points | `:297–306` | draw `n_sims` start states from collection |
| Simulation (no disaster) | `:308–390` | n_sims × n_sim_periods; write 4 `sim_*` files |
| Simulation (with disaster) | `:392–464` | same but disaster node enabled; write 4 `sim_dis_*` |
| IRF starting points | `:467–487` | draw n_sample starts; force start #1 = mean state |
| MIT-shock IRFs (generalized) | `:489–792` | loop over `irf_indices`; write `<tag>_irf_*` files |
| Long-lived IRFs | `:795–978` | only if `run_bg==1`; write `bg1/fix/fx0_irf_*` |
| Sample paths (interpolated) | `:982–1180` | only if `run_samp==1`; read `sample{1,2,3}_mat.txt`; write `samp{1,2,3}_*` |

---

## 2. Simulation

### 2.1 Inputs read from disk (solver output), `:124–135`
- `next_state_mat.dat` → `next_state_mat(n_quad, n_active_dims, n_states)`:
  next-state policy at each grid point for each quadrature node. Unformatted STREAM.
- `results_mat.dat` → `interp_input_mat(n_states, n_interp)`: the 124 "other
  variables" at each grid point.
- `valuation_mat.dat` → `valuation_input_mat(14, n_states)`: 14 valuation series.

These are turned into Smolyak polynomial coefficients via `F07ABF`:
- `smol_coeffs(qqq,:,:)` — per-quad next-state coeffs, `:147–155`.
- `interp_coeffs` — other-vars coeffs (`nrhs2=n_interp`), `:165–169`.
- `state_coeffs` — maps basis → the 9 raw state variables (`nrhs3=smolyak_d`), built from `state_grid`, `:172–179`.
- `value_coeffs` — 14 valuation series (`nrhs4=14`), `:181–187`.

### 2.2 Stochastic steady state, `:189–206`
Start at `smol_grid(1,:)`; iterate `current_state ← Smolyak(current_state) ·
smol_coeffs(no_shock_idx,:,:)` (the zero-shock transition) until
`max|Δ| < sqrt_eps`. Result `stochsteady_state`. Then interpolate other-vars and
state vars there and print (Capital, Wealth_h, wh, wf, disaster p, ζ/zf, ih, if),
`:208–227`.

### 2.3 Random numbers, `:233–246`
NAG generator seeded `SEED=712`, `LSTATE=633`, `GENID=3` (`:234–237`). Uniform
draws via `g05saf`. The uniforms are mapped to a quadrature node by walking the
cumulative weight vector — i.e. each period draws a discrete next-shock node from
the quadrature distribution (this is a **stochastic** simulation that integrates
out shocks via the quadrature nodes, NOT continuous shock draws).

### 2.4 Burn-in, `:239–289`
- `n_periods = 2*n_burn` (= 20000). Draw `simshock_series(n_periods,1)` uniforms.
- Start at `stochsteady_state`.
- Each period (`:250–289`):
  1. `polyn_points = Smolyak(current_state)`.
  2. Interpolate raw state to read `p_dis = exp(state_tmp(idx_dis))`.
  3. Build `big_weight_vec`: nodes `1..n_quad-1` get `quad_weight_vec·(1−p_dis)`,
     node `n_quad` (disaster) gets `p_dis` (`:264–265`). Disaster IS active in burn-in.
  4. Choose `q_nxt` = smallest node whose cumulative weight ≥ the uniform draw (`:267–272`).
  5. For `ttt > n_burn`, save `current_state` into `state_collection(:,ttt-n_burn)`
     (so `state_collection` is `n_active_dims × n_burn`) (`:274–276`).
  6. Advance: `next_state = Smolyak(current_state)·smol_coeffs(q_nxt,:,:)`,
     clipped to `[−1,1]` (`:278–287`).

### 2.5 Simulation starting points, `:297–306`
`g05ndf` samples `n_sims` integers (1..n_burn) without replacement →
`starting_locs`; `starting_states = state_collection(:,starting_locs)`.

### 2.6 Moment-ensemble simulation, `:308–390` (no disaster) and `:392–464` (with disaster)
`n_periods = n_sim_periods (=400)`. Draw `simshock_series(n_periods,n_sims)`.
Loop `sss = 1..n_sims`:
- `current_state = starting_states(:,sss)`; `shock_series(:,1,sss) =
  shock_grid(no_shock_idx,:)` (= zero) (`:328–330`).
- Each `ttt = 1..n_periods` (`:332–373`):
  1. `polyn_points = Smolyak(current_state)`.
  2. Interpolate & store: `state_series(:,ttt,sss)` (9 raw states),
     `other_vars_series(:,ttt,sss)` (124 vars), `value_series(:,ttt,sss)` (14).
  3. Pick `q_nxt` walking cumulative `quad_weight_vec` — but loop bound
     `q_nxt < n_quad-1`, so **disaster node excluded** (no-disaster run) (`:350–355`).
  4. Advance state via `smol_coeffs(q_nxt,…)`, clip to `[−1,1]`.
  5. Store realized shock for next period: `shock_series(:,ttt+1,sss) =
     shock_grid(q_nxt,:)` (`:367–369`).

The **with-disaster** version (`:397–447`) is identical except: weight vector is
`big_weight_vec` with disaster prob `p_dis = exp(state_series(idx_dis,ttt,sss))`
(`:421–425`), loop bound `q_nxt < n_quad` (disaster reachable), and next-state is
**not** clipped to `[−1,1]` (`:434–436`).

### 2.7 Distinguishing the three "samples" (sample1/2/3)
These are NOT the simulation ensemble. They are read only in the `run_samp==1`
block (`:982–1180`) from `../src/params/sample{1,2,3}_mat.txt`. They are
pre-drawn standardized **empirical** shock paths of length `n_sample_periods=94`,
used to produce a model-implied "sample path" matching the data sample. The three
files are three alternative empirical shock sequences (e.g. different
historical/counterfactual draws); each is run separately and written to
`samp1_*`, `samp2_*`, `samp3_*`. (⚠ UNCLEAR exactly what distinguishes 1 vs 2 vs
3 economically — the file only labels them by index; likely different
identified-shock series. Check MATLAB driver / data construction.)

### 2.8 `sampleK_mat.txt` format, `:1038–1071`
- Read with format `(N f100.0)`, `N = n_shocks*n_sample_periods = 4*94 = 376`,
  into `sample_mat_tmp(n_shocks, n_sample_periods)` (`:1041–1042`). Fortran
  column-major ⇒ the flat stream is read column-by-column: period 1's 4 shocks
  (z, zf, p, omg), then period 2's 4 shocks, … So each "row" of the text file
  (16 comma-separated values seen = 4 periods × 4 shocks) just continues filling
  columns. Confirmed: first values `0, 0, -8.294, -5.306, 0, 0, -0.662, 0.608,…`
  → z=0, zf=0 (inactive in this calibration), p and omg carry the action.
- **Role of the leading zeros:** columns `sidx_z` and `sidx_zf` are 0 because in
  the baseline calibration `n_uni_quad_z`/`n_uni_quad_zf` may be 1 (those shocks
  not active), so the empirical sample only varies p (sidx_p) and omg (sidx_omg).
  The code only copies *active* shocks (those with `n_uni_quad_* > 1`) into
  `sample_mat` (`:1045–1071`); p and omg are additionally clipped to the
  `shock_grid` support for periods 2..94 (`:1056–1069`).
- **Shock → state mapping in the sample path** (`:1083–1153`):
  - For `iii = 1..n_sample`, start at `starting_states(:,iii)`, but overwrite the
    p- and omg-state coordinates with the standardized period-1 sample values:
    `(sample_mat_tmp(sidx_p,1)−dis_grid_mean)/dis_grid_dev` and
    `(…omg…−omg_grid_mean)/omg_grid_dev` (`:1087–1098`); clip to `[−1+eps,1−eps]`.
  - Each period: interpolate other-vars & state; store realized shocks
    `shock_series(:,ttt,iii) = sample_mat_tmp(:,ttt)` (period 1 = zero) (`:1120–1124`).
  - Advance: compute next-state at all `n_quad-2` GH nodes, then **interpolate**
    (NAG `e01zmf`/`e01znf`) across the active-shock node grid `shock_grid_temp`
    to the desired continuous shock `sample_mat(:,ttt+1)`, giving the next state
    per active dimension (`:1126–1143`). Clip to `[−1+eps,1−eps]`.
  - Average over the `n_sample` starts → mean sample path (`:1159–1161`).

---

## 3. Impulse responses

### 3.1 IRF starting points, `:467–487`
`g05ndf` draws `n_sample` start indices from the burn-in collection →
`starting_states`. **Start #1 is forced to the mean state**
`sum(state_collection,2)/n_burn` (`:486–487`). So IRFs are computed from a sample
of stochastic-steady states (plus the mean), then averaged ⇒ **generalized IRFs**.

### 3.2 MIT-shock IRFs (generalized), `:489–792`
`n_periods = n_irf_periods (=200)`. Loop `fff = 1..size(irf_indices)` (6 IRFs).
For each:
- Re-read baseline solution (`next_state_mat.dat`, `results_mat.dat`,
  `valuation_mat.dat`) and the shock-transition matrix
  `next_state_shocktrans_mat_<shock_char>.dat` where `shock_char = irf_indices(fff)`
  as 1 digit (`:508–527`). Build coeffs incl. `smol_irftransition_coeffs` (1-node
  transition into the shocked state) and `smol_irf_coeffs` (post-shock dynamics).
- Parallel loop `iii = 1..n_sample` (OpenMP, `:562–664`):
  - Period 1 & 2: sit at the start state (stochastic steady) — both periods equal
    (`:578–594`).
  - Period 3: jump via `smol_irftransition_coeffs` into the **shocked** state;
    interpolate vars there using the *irf* coeffs; record the shock in
    `shock_series(:,3,iii)`. For a z-shock the realized shock value is overwritten
    with `irf_shock_sizes(fff)` (`:627–628`); for a disaster shock
    (`sidx_dis`) period 3 uses `shock_grid(n_quad,:)` (`:629–631`).
  - Periods 4..200: propagate with the no-shock transition `smol_coeffs(no_shock_idx)`,
    storing zero shocks (`:636–661`).
- **Generalized-IRF differencing** (`:668–691`):
  - When `irf_indices(fff)==0` (baseline/no-shock run): store the whole panel as
    `*_series_N` (the "no shock" reference), and build `*_series_N_helper` = the
    period-1 value broadcast across all periods. The reported IRF for the baseline
    is just the average level `sum(series_N,3)/n_sample`.
  - Otherwise: `irf = sum( shocked − series_N + series_N_helper , 3)/n_sample`,
    i.e. (shocked path − matched no-shock path + initial level), averaged over the
    n_sample starting states. This is the **average treatment-effect / generalized
    IRF**, re-centered at the initial level. NOT a single deterministic IRF.
- Written under a per-shock tag (see §4).

### 3.3 Long-lived ("ll") IRFs, `:795–978` (only if `run_bg==1`)
Loop `fff = 1..size(irf_indices_ll)` (3). These hold a shock fixed for
`irf_ll_length=40` periods. They read time-indexed solution files
`next_state_shocktrans_mat_<10+fff>_<1001>.dat`, then per period
`next_state_shock_mat_<10+fff>_<1000+ttt>.dat` and
`results_shock_mat_<10+fff>_<1000+ttt>.dat` (`:807–861`). Only `iii=1` (the mean
start) is used — **no sampling** (`:817`). After the 40 held periods, periods
`43..200` propagate with the no-shock transition (`:908–929`).
- Differencing reuses the previously stored `*_series_N` / `*_series_N_helper`
  from §3.2 baseline: `irf = path(:,:,1) − series_N(:,:,1) + series_N_helper(:,:,1)`
  (`:933–936`). (⚠ No value-series for ll IRFs — only state, vars, shock written.)
- File tags by `fff`: 1→`bg1`, 2→`fix`, 3→`fx0` (`:939–975`). Note the code
  switches on `fff` not `irf_indices_ll(fff)` (the latter is commented out, `:938`).

---

## 4. Output files (EXHAUSTIVE)

All paths are `trim(results_folder)//<name>` = `../output/tmp/res_<i>/<name>`.
Write format for all series files: `write(10,*) array` — list-directed
(whitespace-separated) ASCII, **column-major flatten** of the Fortran array, with
`ACCESS="STREAM"`, `STATUS="replace"`. The two CSVs use explicit formats.

### 4.1 Metadata (CSV)
| File | Format | Contents |
|---|---|---|
| `grid.csv` | `(7i6)` | 7 ints: `n_I, n_states, n_active_dims, n_interp, n_shocks, n_spread, smolyak_d` (`:108–111`) |
| `grid_locs.csv` | `(14f10.4)` | 14 reals: grid means/devs — `k,tht_h,wh,wf,ih,if,omg` (mean,dev pairs) (`:114–119`) |

### 4.2 Simulation ensemble — no disaster (`:376–390`)
Each array is `(dim, n_sim_periods=400, n_sims=100)`.
| File | Array shape | Content tag |
|---|---|---|
| `sim_state_series.txt` | `(smolyak_d=9, 400, 100)` | raw state series |
| `sim_vars_series.txt` | `(n_interp=124, 400, 100)` | other variables series |
| `sim_value_series.txt` | `(14, 400, 100)` | valuation series |
| `sim_shock_series.txt` | `(n_shocks=4, 400, 100)` | realized shock series |

### 4.3 Simulation ensemble — with disaster (`:450–464`)
| File | Array shape | Content tag |
|---|---|---|
| `sim_dis_state_series.txt` | `(9, 400, 100)` | state w/ disaster |
| `sim_dis_vars_series.txt` | `(124, 400, 100)` | vars w/ disaster |
| `sim_dis_shock_series.txt` | `(4, 400, 100)` | shocks w/ disaster |
| `sim_dis_value_series.txt` | `(14, 400, 100)` | valuation w/ disaster |

### 4.4 MIT-shock generalized IRFs (`:694–790`)
Each array is `(dim, n_irf_periods=200)` (already averaged over n_sample).
Tag prefix selected by `irf_indices(fff)`: `0→none`, `1(sidx_z)→g`,
`2(sidx_zf)→zf`, `3(sidx_p)→p`, `4(sidx_omg)→omg`, `5(sidx_dis)→dis1`.
For each tag `T` in {`g`,`p`,`omg`,`zf`,`none`,`dis1`}:
| File | Array shape | Content tag |
|---|---|---|
| `<T>_irf_state_series.txt` | `(9, 200)` | IRF raw states |
| `<T>_irf_vars_series.txt` | `(124, 200)` | IRF other vars |
| `<T>_irf_shock_series.txt` | `(4, 200)` | IRF shocks |
| `<T>_irf_value_series.txt` | `(14, 200)` | IRF valuation |

(So 6 tags × 4 files = 24 IRF files.) Note: `g` = the z/TFP shock IRF (named
"generalized"/g because `sidx_z`), `dis1` = disaster IRF.

### 4.5 Long-lived IRFs (`:940–975`, only if `run_bg==1`)
Each array `(dim, 200)`. Tags by `fff`: 1→`bg1`, 2→`fix`, 3→`fx0`. **No value
file.**
| File | Array shape | Content tag |
|---|---|---|
| `bg1_irf_state_series.txt`, `fix_…`, `fx0_…` | `(9, 200)` | ll-IRF states |
| `bg1_irf_vars_series.txt`, `fix_…`, `fx0_…` | `(124, 200)` | ll-IRF vars |
| `bg1_irf_shock_series.txt`, `fix_…`, `fx0_…` | `(4, 200)` | ll-IRF shocks |

(3 tags × 3 files = 9 files.)

### 4.6 Sample paths (`:1163–1176`, only if `run_samp==1`)
Each array `(dim, n_sample_periods=94)` (averaged over n_sample). Tag = `samp<fff>`,
`fff∈{1,2,3}`. **No value file.**
| File | Array shape | Content tag |
|---|---|---|
| `samp1_state_series.txt`, `samp2_…`, `samp3_…` | `(9, 94)` | sample-path states |
| `samp1_vars_series.txt`, `samp2_…`, `samp3_…` | `(124, 94)` | sample-path vars |
| `samp1_shock_series.txt`, `samp2_…`, `samp3_…` | `(4, 94)` | sample-path shocks |

(3 tags × 3 files = 9 files.)

### 4.7 Files READ (not written) by this module
- `next_state_mat.dat`, `results_mat.dat`, `valuation_mat.dat` (solver output).
- `next_state_shocktrans_mat_<d>.dat` (MIT IRFs).
- `next_state_shocktrans_mat_<dd>_<iter>.dat`, `next_state_shock_mat_<dd>_<iter>.dat`,
  `results_shock_mat_<dd>_<iter>.dat` (ll IRFs; `dd=10+fff`, `iter=1000+ttt`).
- `../src/params/sample{1,2,3}_mat.txt` (empirical shock paths).

---

## 5. Moment computation

**Almost none in Fortran.** This module only:
- Computes and *prints to stdout* the stochastic steady state (`:217–227`) — not
  saved to a moment file.
- Averages IRF/sample panels across starting points (the `sum(...)/n_sample`
  generalized-IRF means, `:681–690`, `:1159–1161`).

It does **not** compute business-cycle moments (means, std devs, correlations,
autocorrelations) — despite the "Business cycle moments" banner (`:320`,`:392`),
those banners just label the simulation-generation phase. All actual moment
computation (and IRF post-processing/plots) happens in **MATLAB**, which reads the
`.txt` series above. The Python port must reproduce these series files (or
equivalent in-memory arrays with identical shapes/orderings), then port the MATLAB
moment code separately.

---

## 6. Flags / uncertainties

- ⚠ **`g` prefix:** the z/TFP-shock IRF is written with prefix `g` (`g_irf_*`),
  selected by `irf_indices==sidx_z` (`:694`). Confirm the MATLAB side reads `g_*`
  as the productivity-shock IRF.
- ⚠ **sample1/2/3 economics:** what distinguishes the three empirical shock files
  is not documented in this module (see §2.7). Resolve from the MATLAB/data scripts.
- ⚠ **ll-IRF differencing** reuses `*_series_N(:,:,1)` and `*_series_N_helper`
  from the §3.2 baseline run; these must still be in scope/allocated (they are,
  not deallocated before this block). The Python port must keep that baseline
  reference path around.
- ⚠ **Disaster prob in burn-in vs no-disaster sim:** burn-in (`:264`) and the
  with-disaster sim (`:424`) enable the disaster node; the first ("no disaster")
  sim explicitly excludes it via the `q_nxt < n_quad-1` bound (`:352`). Preserve
  this distinction.
- ⚠ **State clipping asymmetry:** no-disaster sim clips next_state to `[−1,1]`
  (`:360`); with-disaster sim does NOT (`:434`); sample path clips to `[−1+eps,1−eps]`.
- The `value_series`/`valuation` arrays are 14-wide everywhere; the meaning of the
  14 columns is set by the solver (`valuation_mat.dat`), not defined here.
- `n_valuation=500` and several declared arrays (`kappa_mat`, `l_aggr_mat`,
  `q_mat`, `v_mat`, `state_spread_*`, `next_zf_mat`, etc.) are declared but unused
  in the active code path — likely vestigial. Do not port unless MATLAB needs them.
