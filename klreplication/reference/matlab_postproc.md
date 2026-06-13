# Kekre & Lenel (2024, AER) — MATLAB post-processing digest

Reference for porting the MATLAB post-processing of the Fortran solver output to
Python (no MATLAB runtime). Source (read-only):

`…/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/Code/src/matlab/`

All `file:line` citations are to that folder. The Fortran solver writes per-calibration
output into `../output/tmp/res_<i>/`; data for comparison lives in `../src/data/`.

---

## 0. Top-level numbers and the calibration set

From `main.m:14-29`:

| name | value | meaning |
|---|---|---|
| `n_params` | 65 | columns in each `param_file_<i>.csv` (header at `param_file_1.csv` row 1) |
| `n_sims` | 100 | number of simulation paths |
| `n_sim_periods` | 400 | periods per simulation path |
| `n_irf_periods` | 200 | periods per IRF |
| `n_sample_periods` | 94 | periods in `samp1/2/3` historical samples |
| `n_irf` | 8 | number of IRF *types* collected per calibration |

`n_comp` (number of calibrations actually run) is read at runtime from
`../output/tmp/n_comp.txt` (`main.m:38-39`). In the shipped package
`create_param_files.m:139` sets `n_comp = 9`.

**Calibration indices** (`main.m:22-30`) — these are the `res_<i>` folders and the
columns of every `table_*_moms` matrix:

| idx | name | description |
|---|---|---|
| 1 | `ix_bm` | benchmark / "Model" |
| 2 | `ix_no_omg` | no safety shock (σ^ω = 0) |
| 3 | `ix_symm` | symmetric preferences γ = γ* |
| 4 | `ix_symm_flex` | symmetric + flexible wages (χ^W = 0) |
| 5 | `ix_tayl_y0` | Taylor rule with φ^y = 0.5/4 |
| 6 | `ix_tayl_rho` | Taylor rule φ^y = 0.5/4, ρ^i = 0.5 |
| 7 | `ix_no_omg_symm` | no ω and symmetric |
| 8 | `ix_nocorr_nobg` | no corr(p,ω), no gov debt b^g = 0 |
| 9 | `ix_nocorr` | no corr(p,ω) |

`param_input_vec` (5 × n_comp; rows = [mean p, std p, mean ω, std ω, skew ω]) is
saved in `../src/params/additional_params.mat` (`create_param_files.m:266-267`) and
loaded via `read_results.m:15`. Used by `create_recession_fig` and `create_safety_fig`
to scale the data series, and by `print_params`.

---

## 1. `main.m` — driver flow

Order of operations (`main.m`):

1. **Preallocate** table matrices (`main.m:42-52`) and raw series arrays (`:55-67`),
   plus `collected_irfs = zeros(n_irf_periods-2, 37, n_irf, n_comp)` = `(198,37,8,n_comp)`
   (`:70`). The 37 is the number of derived IRF columns from `extract_irfs`; 198 = periods 3:200.

2. **Loop `ccc = 1:n_comp`** (`:72`):
   - Set `data_path = ../output/tmp/res_<ccc>/`, `param_file = ../src/params/param_file_<ccc>.csv` (`:76-77`).
   - `read_results` (`:81`) → reads params, grid, all series. `save data.mat` (refreshable via `new_results`).
   - Open num-check text file `../output/tables/num_checks/results_<ccc>.txt` (`:88`).
   - `print_params` (`:92`) → writes parametrization to that file.
   - **Without disaster** (`dis=false`, `:98`):
     - Pooled series across all sims → `extract_series` + `numerical_check` (`:101-103`).
     - Loop `sss=1:n_sims`: per-path `extract_series` + `calc_moments` (`:106-110`).
   - **With disaster** (`dis=true`, `:115`): same, on `sim_dis_series` (`:116-124`).
   - `collect_moments` (`:127`) → averages per-sim moments into `table_*_moms(:,ccc)`.
   - **IRFs** (`:133-142`): `extract_irfs` called on 8 IRF types, stored in
     `collected_irfs(:,:,k,ccc)` for k=1..8:
     | k | source series | shock |
     |---|---|---|
     | 1 | `p_irf_series` | disaster-probability p |
     | 2 | `omg_irf_series` | safety ω |
     | 3 | `g_irf_series` | productivity z (home) — note: `g_irf` files |
     | 4 | `dis1_irf_series` | disaster realization |
     | 5 | `zf_irf_series` | foreign productivity z_f |
     | 6 | `bg1_irf_series` | swap-line response |
     | 7 | `fix_irf_series` | swap line, fixed nominal rate |
     | 8 | `fx0_irf_series` | no-shock baseline, fixed nominal rate |
     - **`collected_irfs(:,:,7) = fix - fx0`** (`:142`): k=7 is differenced against k=8.
   - If `ccc == ix_bm` (`:144`): `collect_swap_moments`, `create_recession_fig`, `create_safety_fig`.
   - `extract_irfs` also returns `irf_idxs` (struct mapping names→columns) and `irf_titles`
     (LaTeX labels), set once from the first call (`:133`).

3. **After loop**: `create_tables` (`:153`), then `create_figures` (`:154`).

**Figure → script map:**
- `create_figures.m`: fig 2, 3, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 (+ many `addl_figs/`).
- `create_recession_fig.m`: fig 4, 6, 7, and `table_A1.tex`.
- `create_safety_fig.m`: fig 10.

**Table → script map:**
- `create_tables.m`: table 1–10 (`.tex`) + table A2.
- `create_recession_fig.m`: table A1.
- Moments fed by `calc_moments` / `collect_moments` / `collect_swap_moments`.

---

## 2. Reading layer — THE KEYSTONE COLUMN MAP

### 2.1 File set per `res_<i>/` and how `read_series.m` parses it

`read_series.m` is invoked once per `series_id` (set in `read_results.m:72-135`). For
each id it reads three or four whitespace-delimited text files and **vertically stacks**
them into one matrix `<series_id>_series` with rows = variables, cols = periods
(× n_series for sims).

Files and their reshape dims (`read_series.m:12-44`):

| file suffix | reshape `[rows, n_periods, n_series]` | rows = |
|---|---|---|
| `_state_series.txt` | `[smolyak_d, n_periods, n_series]` | `smolyak_d` state cols |
| `_vars_series.txt` | `[n_interp, n_periods, n_series]` | `n_interp` policy/aux vars |
| `_shock_series.txt` | `[n_shocks, n_periods, n_series]` | `n_shocks` shocks |
| `_value_series.txt` | `[14, n_periods, n_series]` | 14 valuation/accounting vars |

`n_periods`/`n_series` chosen by id prefix (`read_series.m:1-10`):
- `sim*` → `n_sim_periods=400`, `n_sims=100`
- `sam*` (samp1/2/3) → `n_sample_periods=94`, `n_series=1`
- else (IRFs) → `n_irf_periods=200`, `n_series=1`

`smolyak_d`, `n_interp`, `n_shocks` come from `grid.csv` (`read_results.m:46-52`,
`gridsize_vec`):
`n_I, n_states, n_active_dims, n_interp, n_shocks, n_spread, smolyak_d`.

**Stacking order** (`read_series.m:36-49`):
- For `samp*`, `bg*`, `fix`, `fx0` (no value file):
  `series = [shock_series; state_series; other_vars_series]` → **137 rows**.
- For everything else (incl. `sim*`, `*_irf` except the swap ones):
  `series = [shock_series; state_series; other_vars_series; value_series]` → **151 rows**.

This matches the preallocation in `main.m`: sim/irf arrays have 151 rows
(`main.m:55-61`), the swap/samp arrays 137 (`:62-67`).

### 2.2 `get_var_indices.m` — absolute row indices into the stacked matrix

CRITICAL: the `jx_*` indices in `get_var_indices.m` are **absolute row positions** in the
stacked `[shocks(4); states(?); vars; value(14)]` array, NOT positions within any single
file. The numbering implies the layout is:

- rows **1–4**: shocks (`n_shocks` effective = 4 used here)
- rows **5–13**: states (so the model uses **9 active state rows** here; `jx_k=5` … `jx_omg=13`,
  and `jx_yh=14` is the first vars-block row)  ⚠ VERIFY `smolyak_d == 9` for the benchmark grid
  (depends on `grid.csv`; the `extract_series` code assumes these exact positions).
- rows **14–137**: `other_vars_series` (the n_interp block) ⚠ VERIFY `n_interp == 124`
- rows **138–151**: 14-element value/valuation block.

**Full mapping** (transcribed verbatim from `get_var_indices.m:3-167`):

#### Shocks (rows 1–4)
| jx | row | meaning |
|---|---|---|
| `jx_z_shock` | 1 | home productivity growth shock (z; disaster lives here) |
| `jx_zf_shock` | 2 | foreign productivity shock |
| `jx_p_shock` | 3 | disaster-probability shock |
| `jx_omg_shock` | 4 | safety/convenience shock |

#### States (rows 5–13)
| jx | row | meaning |
|---|---|---|
| `jx_k` | 5 | capital stock (de-trended) |
| `jx_tht_h` | 6 | home wealth share θ_H |
| `jx_zf` | 7 | foreign productivity level z_f |
| `jx_whs` | 8 | home real wage chosen last period (state) |
| `jx_wfs` | 9 | foreign real wage chosen last period (state) |
| `jx_p` | 10 | log disaster prob (level recovered via `exp`) |
| `jx_rho_ih` | 11 | lagged home nominal rate term (Taylor inertia) |
| `jx_rho_if` | 12 | lagged foreign nominal rate term |
| `jx_omg` | 13 | log safety wedge ω (level via `exp(.)+omg_shift`) |

#### Vars block (rows 14–137) — `other_vars_series`
| jx | row | meaning |
|---|---|---|
| `jx_yh` | 14 | home output (goods units, de-trended) |
| `jx_yf` | 15 | foreign output |
| `jx_lh` | 16 | home labor |
| `jx_lf` | 17 | foreign labor |
| `jx_ch` | 18 | home aggregate consumption |
| `jx_cf` | 19 | foreign aggregate consumption |
| `jx_inv` | 20 | investment (capital units) |
| `jx_pi` | 21 | capital rental income / profit per unit (π) |
| `jx_pif` | 22 | foreign-side π (`pif`) ⚠ VERIFY (unused in extract_series) |
| `jx_P_P_hh` | 23 | P/P_h (home CPI over home-good price) |
| `jx_P_P_hf` | 24 | P*/P_h* (foreign CPI over home-good price) |
| `jx_qx` | 25 | real exchange rate q (RER) |
| `jx_s` | 26 | terms of trade s |
| `jx_Erfh` | 27 | E[r] home risk-free (home-currency, home-investor) |
| `jx_Erff` | 28 | E[r*] foreign risk-free |
| `jx_Erk` | 29 | E[r^k] return on capital |
| `jx_q` | 30 | price of capital q^k |
| `jx_sh_h` | 31 | home portfolio share (risky/total) |
| `jx_sh_f` | 32 | foreign portfolio share |
| `jx_ssh_h` | 33 | home safe-bond sub-share (dollar bond) |
| `jx_ssh_f` | 34 | foreign safe-bond sub-share |
| `jx_aggrw` | 35 | aggregate wealth |
| `jx_wh` | 36 | home wage (in home basket) |
| `jx_wf` | 37 | foreign wage |
| `jx_infl_h` | 38 | home gross inflation P/P_{-1} |
| `jx_infl_f` | 39 | foreign gross inflation |
| `jx_nom_i_h` | 40 | home gross nominal rate i |
| `jx_nom_i_f` | 41 | foreign gross nominal rate i* |
| `jx_h_ksav` | 42 | home savings in capital |
| `jx_h_kap` | 43 | home deployed capital κ |
| `jx_f_kap` | 44 | foreign deployed capital |
| `jx_Erfh_f` | 45 | E[r] home rf, foreign-investor units |
| `jx_Erff_f` | 46 | E[r*] foreign rf, foreign-investor units |
| `jx_Erk_f` | 47 | E[r^k] foreign-investor units |
| `jx_Ekpa` | 48 | E[κ] expected deployed capital |
| `jx_v1` | 49 | home value function V1 |
| `jx_v2` | 50 | foreign value function V2 |
| `jx_invh` | 51 | home investment spending (home basket) |
| `jx_invf` | 52 | foreign investment spending |
| `jx_h_bhsav` | 53 | home savings in home-currency bonds b_H |
| `jx_bg` | 54 | extra government safe debt (de-trended) |
| `jx_ce1` | 55 | home certainty-equiv CE1 |
| `jx_ce2` | 56 | foreign certainty-equiv CE2 |
| `jx_h_sav` | 57 | home total savings |
| `jx_f_sav` | 58 | foreign total savings |
| `jx_Erf_spread` | 59 | E[rf spread] |
| `jx_corr_rf_spread_mh` | 60 | corr(rf-spread, m_H) |
| `jx_corr_rf_spread_mf` | 61 | corr(rf-spread, m_F) |
| `jx_corr_rf_spread_mhB` | 62 | corr(rf-spread, m_H) variant B → Table 10 |
| `jx_corr_rf_spread_mfB` | 63 | corr(rf-spread, m_F) variant B → Table 10 |
| `jx_corr_p_mhB` | 64 | corr(p, m_H) variant B |
| `jx_corr_p_mfB` | 65 | corr(p, m_F) variant B |

**Bond price ladder** (zero-coupon prices q at maturities 1,2,3,4,19,20 quarters;
suffixes: `h` = home/Treasury-with-ω, `f` = foreign, `hw` = home "world"/no-ω). Rows 66–83:
| jx | row | | jx | row | | jx | row |
|---|---|---|---|---|---|---|---|
| `jx_q1h` | 66 | | `jx_q2hw` | 71 | | `jx_q4h` | 75 |
| `jx_q1f` | 67 | | `jx_q3h` | 72 | | `jx_q4f` | 76 |
| `jx_q1hw` | 68 | | `jx_q3f` | 73 | | `jx_q4hw` | 77 |
| `jx_q2h` | 69 | | `jx_q3hw` | 74 | | `jx_q19h` | 78 |
| `jx_q2f` | 70 | |  |  | | `jx_q19f` | 79 |

| jx | row | | jx | row |
|---|---|---|---|---|
| `jx_q19hw` | 80 | | `jx_q20f` | 82 |
| `jx_q20h` | 81 | | `jx_q20hw` | 83 |

**Expected bond returns E[q]** (rows 84–101): `jx_Eq{1,2,3,4,19,20}{h,f,hw}`:
| jx | row | jx | row | jx | row |
|---|---|---|---|---|---|
| `jx_Eq1h` 84 | `jx_Eq1f` 85 | `jx_Eq1hw` 86 | `jx_Eq2h` 87 | `jx_Eq2f` 88 | `jx_Eq2hw` 89 |
| `jx_Eq3h` 90 | `jx_Eq3f` 91 | `jx_Eq3hw` 92 | `jx_Eq4h` 93 | `jx_Eq4f` 94 | `jx_Eq4hw` 95 |
| `jx_Eq19h` 96 | `jx_Eq19f` 97 | `jx_Eq19hw` 98 | `jx_Eq20h` 99 | `jx_Eq20f` 100 | `jx_Eq20hw` 101 |

**Foreign-investor expected bond returns fE[q]** (rows 102–119): `jx_fEq{...}{h,f,hw}`:
| jx | row | jx | row | jx | row |
|---|---|---|---|---|---|
| `jx_fEq1h` 102 | `jx_fEq1f` 103 | `jx_fEq1hw` 104 | `jx_fEq2h` 105 | `jx_fEq2f` 106 | `jx_fEq2hw` 107 |
| `jx_fEq3h` 108 | `jx_fEq3f` 109 | `jx_fEq3hw` 110 | `jx_fEq4h` 111 | `jx_fEq4f` 112 | `jx_fEq4hw` 113 |
| `jx_fEq19h` 114 | `jx_fEq19f` 115 | `jx_fEq19hw` 116 | `jx_fEq20h` 117 | `jx_fEq20f` 118 | `jx_fEq20hw` 119 |

**"Who prices the bond" indicators Pq** (rows 120–137): `jx_Pq{...}{h,f,hw}`:
| jx | row | jx | row | jx | row |
|---|---|---|---|---|---|
| `jx_Pq1h` 120 | `jx_Pq1f` 121 | `jx_Pq1hw` 122 | `jx_Pq2h` 123 | `jx_Pq2f` 124 | `jx_Pq2hw` 125 |
| `jx_Pq3h` 126 | `jx_Pq3f` 127 | `jx_Pq3hw` 128 | `jx_Pq4h` 129 | `jx_Pq4f` 130 | `jx_Pq4hw` 131 |
| `jx_Pq19h` 132 | `jx_Pq19f` 133 | `jx_Pq19hw` 134 | `jx_Pq20h` 135 | `jx_Pq20f` 136 | `jx_Pq20hw` 137 |

#### Value block (rows 138–151) — `value_series` (only in non-swap, non-samp series)
Documented in `get_var_indices.m:148-167`. Two blocks of 7: realized vs. lagged-expectation.
| jx | row | meaning |
|---|---|---|
| `jx_val_nfa` | 138 | current NFA position using expected κ |
| `jx_val_EnfaEnd` | 139 | discounted final NFA position |
| `jx_val_Eval` | 140 | discounted valuation term |
| `jx_val_Enx` | 141 | discounted net-exports term (current nx, shifted one period) |
| `jx_val_Evalrf` | 142 | discounted valuation r component |
| `jx_val_Evalrff` | 143 | discounted valuation r* component |
| `jx_val_Evalseign` | 144 | discounted seigniorage valuation component |
| `jx_val_Enfa` | 145 | E_{t-1} of jx_val_nfa (lagged expectation) |
| `jx_val_EEnfaEnd` | 146 | E_{t-1} of jx_val_EnfaEnd |
| `jx_val_EEval` | 147 | E_{t-1} of jx_val_Eval |
| `jx_val_EEnx` | 148 | E_{t-1} of jx_val_Enx |
| `jx_val_EEvalrf` | 149 | E_{t-1} of jx_val_Evalrf |
| `jx_val_EEvalrff` | 150 | E_{t-1} of jx_val_Evalrff |
| `jx_val_EEvalseign` | 151 | E_{t-1} of jx_val_Evalseign |

These 14 value rows feed **Table 7** (external-adjustment covariance decomposition);
innovations = realized minus lagged-expectation (`calc_moments.m:114-120`).

### 2.3 Params, grid, extra_data (read_results.m)
- `param_file_<i>.csv` (`read_results.m:3-13`): 65 named params → struct `trg_prms`.
  Header order is the column list in `param_file_1.csv` row 1 (zeta, bbeta_h, bbeta_f,
  gma_h, gma_f, ies_h, ies_f, chi, sigma, varsigma_h, varsigma_f, l_target, delta, aalpha,
  chiX, inv_share_h, varphi_w, vareps_w, phi_pi_h, phi_pi_f, phi_y_h, phi_y_f, tayl_ic_h,
  tayl_ic_f, rho_i, b_lmbd, bg_yss, sig_z, std_zf, rho_zf, disast_p, disast_shock, p_rho,
  p_std, omg_mean, omg_std, omg_shift, omg_rho, corr_omg_dis, tht_trgt_h, theta_h_dev,
  k_dev_param, k_grid_adj, w_dev_param, w_grid_adj, ih_grid_mean, if_grid_mean, ih_grid_dev,
  if_grid_dev, mu_k, mu_tht, mu_zf, mu_wh, mu_wf, mu_p, mu_ih, mu_if, mu_omg, n_quad_z,
  n_quad_zf, n_quad_p, n_quad_omg, foreign_trading, run_bg, run_samp). Plus
  `trg_prms.debt_to_equity = 1.0` set manually (`:14`).
- `grid.csv` (7 numbers): `n_I, n_states, n_active_dims, n_interp, n_shocks, n_spread, smolyak_d` (`:46-52`).
- `grid_locs.csv` (14 numbers): grid means/devs for k, θ_H, w_h, w_f, i_h, i_f, ω (`:54-67`).
- `extra_data.csv` (10 numbers, `:32-42`): `chi0_vec(1:2)` (labor disutilities),
  `sig_dis`, `sig_omg`, `dis_grid_mean`, `omg_grid_mean`, `yh_ss`, `aggr_css`.

---

## 3. `extract_series.m` — derived series (one big script, operates on `temp_series`)

`temp_series` is a `[periods × 151]` matrix (one sim path, or pooled, or one IRF).
`extract_series.m:3` first runs `get_var_indices`. Then it builds named series. Key
transformations (cite `extract_series.m`):

**Trend handling.** `z_shocks = temp(:,1)`; cumulative growth `z_series = exp(cumsum(z_shocks))`
(`:13`). Most quantity series are stored de-trended and re-trended by multiplying by `z_series`
(e.g. `k_series = temp(:,jx_k).*z_series` `:22`; same for yh, yf, inv, c_h, c_f, wages, savings,
capital, aggregate wealth, bg).

**Disaster extraction** (`:15-18`): `dis_series` isolates z-shocks equal to `-disast_shock`
(within sqrt(eps)); `z_shocks_no_dis = z_shocks - dis_series`. `k_chosen_series = k/exp(dis)`
shifted one period (`:24-25`).

**Levels from logs.** `p_series = exp(temp(:,jx_p))` (`:38`);
`omg_series = exp(temp(:,jx_omg)) + omg_shift` (`:43`), then **adjusted for gov debt**:
`omg_series -= b_lmbd*bg_xtra_rel_series` (`:139`). `q_series = temp(:,jx_q)` (`:85`);
`q_f_series = q*P_Phh/P_Phf` (`:86`).

**Exchange rate / inflation** (`:88-100`): `P_h_series=cumprod(infl_h)`,
`P_f_series=cumprod(infl_f)`; `E_change = -log(qx_{t}/qx_{t+1} * infl_h_{t+1}/infl_f_{t+1})`,
nominal level `E_series = [1; exp(cumsum(E_change))]`; `qx_change = -log(qx_t/qx_{t+1})`.

**Consumption split** (`:113-122`): home/foreign consumption split into home-good and
foreign-good components via CES with `varsigma_*`, `sigma`, terms of trade `s`.

**Bond yields** (`:394-405`): `yieldN = -log(qN)/N` for N=1,2,3,4 each of h, hw, f.

**Realized bond returns** `rqN_*` (`:179-240`): holding-period log returns of N-maturity
bonds, e.g. `rq1_h = log(1/infl_h_{t+1}/q1_h_t)`; cross-currency variants `_fh`, `_hf`, `_hwf`
add the relative-price-change term `log(P_Phf_{t+1}/P_Phf_t / P_Phh_{t+1}*P_Phh_t)`.
All shifted/padded so length = `periods` (`:211-240`).

**Risk-free / equity returns** (`:263-301`):
- `rfh = log(nom_ih_{t}/infl_h_{t+1})` (home rf, log) (`:263`).
- `rfh_omg = log(nom_ih/infl_h/(1-omg))` (convenience-adjusted) (`:267`).
- `rff = log(nom_if/infl_f)`; `rff_h = rff + relative-price change` (`:271-274`).
- `rk = log(exp(dis)*(π + (1-δ)q)/q_{-1})` (`:275`); `rk_f` adds relative-price term.
- Levered Treasury/foreign-bond returns `rfh_lev`, `rff_lev` using ζ weights (`:280-281`).
- `rA = log((1+D/E)exp(rk) - (D/E)exp(rfh_lev))` levered equity return (`:285`);
  `rA_f` foreign analog. (`debt_to_equity=1`.)
- Excess returns: `exc_retA = rA - rfh`, `exc_retA_omg = rA - rfh_omg`,
  `exc_rf = rff_h - rfh`, `exc_ret = rk - rfh`, foreign variants (`:289-301`).

**Portfolios / NFA** (`:124-336`):
- Portfolio shares `sh_h, sh_f, ssh_h, ssh_f` (`:125-128`).
- Savings split: `h_b_sav = h_sav - h_ksav`; `h_bf_sav = h_b_sav - h_bh_sav` (`:151-153`).
- `nfa = h_sav_{t} - h_kap_{t+1}*exp(-dis_{t+1})*q_t` (`:314`); `nfa_k` capital-only analog.
- `nfa_rel = nfa/yh`; `nfa_rel_growth = Δnfa/yh` (`:319-321`); composition ratios
  `nfa_k_rel, nfa_bh_rel, nfa_bf_rel` etc. relative to savings or yh (`:323-331`).
- `nfa_return` portfolio return on NFA using mean shares × component returns (`:334`).

**Net exports** (`:409-414`): `imports = chf/s/P_Phh`; `exports = cfh + inv_h + capital
disinvestment terms`; `nx = exports - imports`; `nx_rely = nx/yh`.

**Valuation terms** (`:423-448`): `val_series`, `valb_series` and r/r* sub-components used
by IRF column 25 and the Table-7 logic (note Table 7 instead uses the jx_val_* columns).

**Dividends / price-dividend** (`:346-392`): two dividend definitions (`div_series_1` per-unit,
`div_series_2` aggregate), price-dividend ratios, and a smoothed PD ratio
`div_price_smoothed_series_1` (4-quarter `movsum`) used in the return-predictability regression.

**Regression inputs** (`:467-470`): `y_growth_series = 100*4Q log output growth`;
`uip_pvt_series = 100*(yield1_f - yield1_hw - E_change_{t+1})` (UIP deviation).

**SDFs** (`:472-477`): marginal-utility-of-consumption `mc_h/mc_f` with Epstein-Zin
adjustment via `v1/ce1`, `v2/ce2`, IES/RRA params → `sdf_h, sdf_f`.

**`collected_data`** (`:417-420`): a 17-column matrix snapshot (used elsewhere / sanity).

> Python port note: many series are built with a leading-element pad (`[x(1); x]`) to keep
> length = `periods` after differencing/shifting. Reproduce these paddings exactly or
> downstream means/stds will differ at the boundaries.

---

## 4. `extract_irfs.m` — IRF derived columns (function)

`function [irf_data, irf_idxs, irf_titles] = extract_irfs(temp_series, trg_prms)`
(`extract_irfs.m:7`). Internally calls `extract_series` (`:12`), then builds
`irf_data` of size `[198, 37]` (periods 3:200 of each series minus its period-1 value,
× 10000 to express in basis points). `E_series` is re-based to the no-shock path (`:16`).

The 37 columns (each line sets `irf_data(:,c)`, `irf_idxs.<name>=c`, `irf_titles{c}`):

| col | irf_idxs name | series (Δ vs t=1, ×10000) | LaTeX title |
|---|---|---|---|
| 1 | `p` | `p_series` | $p$ |
| 2 | `excA` | `exc_retA_series` | $r^e - r$ |
| 3 | `excrf` | `exc_rf_series` | $r^* - \Delta\log q - r$ |
| 4 | `excrk` | `exc_ret_series` | $r^k - r$ |
| 5 | `Erh` | `E_rfh_series` | $E[r]$ |
| 6 | `Erf` | `E_rff_f_series` | $E[r^*]$ |
| 7 | `log_qx` | `log(qx_series)` | $\log q$ |
| 8 | `log_E` | `log(E_series)` | $\log E$ |
| 9 | `thth` | `tht_h_series` | $\theta$ |
| 10 | `log_yh` | `log(yh_series)` | $\log y$ |
| 11 | `log_yf` | `log(yf_series)` | $\log y^*$ |
| 12 | `log_ch` | `log(c_h_series)` | $\log c$ |
| 13 | `log_cf` | `log(c_f_series)` | $\log c^*$ |
| 14 | `log_lh` | `log(lh_series)` | $\log \ell$ |
| 15 | `log_lf` | `log(lf_series)` | $\log \ell^*$ |
| 16 | `log_wh` | `log(wh_series)` | $\log w$ |
| 17 | `log_wf` | `log(wf_series)` | $\log w^*$ |
| 18 | `log_infl_h` | `log(infl_h_series)` | $\log P/P_{-1}$ |
| 19 | `log_infl_f` | `log(infl_f_series)` | $\log P^*/P^*_{-1}$ |
| 20 | `log_x` | `log(inv_series)` | $\log x$ |
| 21 | `omg` | `omg_series` | $\omega$ |
| 22 | `log_qk` | `log(q_series)` | $\log q^k$ |
| 23 | `log_kh` | `log(h_kap_series)` | $\log \kappa$ |
| 24 | `nfa_ybar` | `nfa_series / yh(1)` | $nfa/y_0$ |
| 25 | `val_ybar` | `valb_series / yh(1)` | $val/y_0$ |
| 26 | `nx_ybar` | `nx_series / yh(1)` | $nx/y_0$ |
| 27 | `log_z` | `log(z_series)` | $\log z$ |
| 28 | `log_zf` | `zf_series` | $\log z_f$ |
| 29 | `ih` | `nom_ih_series` | $i$ |
| 30 | `if` | `nom_if_series` | $i^*$ |
| 31 | `bg` | `bg_rel_series` | $-b^g_{H,s}/(c+\zeta^* q^{-1}c^*)$ |
| 32 | `matteo` | `0.1*omg_rho.^(0:.)` | $\varsigma$ (placeholder process) |
| 33 | `pog` | `0.1*omg_rho.^(0:.)` | $\log\gamma^*$ (placeholder) |
| 34 | `EexcA` | `E_excA_series` | $E[r^e-r]$ |
| 35 | `retA` | `rA_series` | $r^e$ |
| 36 | `nfa_y` | `nfa_rel_series` | $nfa/y$ |
| 37 | `nx_y` | `nx_rely_series` | $nx/y$ |

Note cols 24/25/26 divide by **period-1 yh** (constant), while 36/37 use the time-varying ratio.

---

## 5. Moment scripts (populate "Model" columns)

### 5.1 `calc_moments.m` — per-simulation (`sss`), per-disaster-flag

Stores into `table_*_moms_tmp(:,sss)`. Means across sims taken later in `collect_moments`.
Most Table 2/3/4/5/6/10 moments computed only when `~dis`; Table 7 only when `dis`.

**Table 2 (targeted moments)** `calc_moments.m:33-47` (each ×100 / ×4 as noted):
1. `mean(yf/(s*yh))` — relative output (`:33`).
2. `100*std(diff(log(c_h)))` (`:34`).
3. `100*std(diff(log(yf)))` (`:35`).
4. `corr(yf/yh_t, yf/yh_{t+1})` AR(1)-ish (`:36`).
5. `100*std(diff(log(inv)))` (`:37`).
6. `4*100*mean(rfh)` — annualized mean rf (`:38`).
7. `100*mean(nfa_rel/4)` (`:39`).
8. `corr(resid_dp, resid_carry)` — corr of residuals from (a) regression of
   `4*100*exc_retA(5:end)` on smoothed PD ratio (`:16-18`) and (b) the carry regression
   `uip_pvt` on past yield diff and output growth (`:21-25`) (`:40`).
9. `4*100*mean(rA - mean(rfh))` — equity premium (`:41`).
10. `reg_coeffs1(2)` — slope of `nfa_rel_growth` on `exc_retA` and `exc_rf` (3-var regression `:28-31`).
11. `100*bg_yss*mean((c_f/qx)/(4*yh))` (`:43`).
12. `mean(lh)` (`:44`); 13. `mean(lf)` (`:45`).
14. `100*mean(log(infl_h))` (`:46`); 15. `100*mean(log(infl_f))` (`:47`).

**Table 3 (comovements)** `calc_moments.m:52-72`:
1. slope of `uip_pvt` on lagged 4Q output growth (`reg_coeffs2(2)`, `:54-58`).
2. slope of `uip_pvt` on `100*exc_retA` (`reg_coeffs3(2)`, `:60-64`).
3. `reg_coeffs1(3)` — the `exc_rf` coefficient from the Table-2 NFA regression (`:67`).
4. `100*mean((h_ksav - h_kap)/(4*yh))` — net equity position (`:69`).
5. `100*mean(h_bh_sav/(4*yh))` — home-currency bonds (`:70`).
6. `100*mean(h_bf_sav/(4*yh))` — foreign-currency bonds (`:71`).

**Table 4 (second moments)** `calc_moments.m:79-84`:
1. `4*100*std(rfh)`; 2. `4*100*std(exc_retA)`; 3. `4*100*std(rff_h - rfh)`;
4. `100*std(qx_change)`; 5. `100*std(E_change)`;
6. `corr(Δlog qx, Δlog c_f - Δlog c_h)` (`:84`).

**Table 5 (output vol)** `calc_moments.m:92-93`: `100*std(log(yh_t/yh_{t-1}))`, same for yf.

**Table 6 (NFA)** `calc_moments.m:101-106`:
1. `100*std(nfa_rel_growth)`; 2. `100*std(nx_rely)`; 3. `100*std(nfa_rel_growth - nx_rely)`;
4. `100*mean(nfa_rel_growth)`; 5. `100*mean(nx_rely)`; 6. `100*mean(nfa_rel_growth - nx_rely)`.

**Table 7 (external adjustment, `dis` only)** `calc_moments.m:114-136`:
Innovations = `temp(:,jx_val_*)_{t} - temp(:,jx_val_EE*)_{t-1}` (realized minus lagged
expectation) for nfa, val, val_rf, val_rff, val_seign, nx, nfa_resid (`:114-120`). Then:
`varnfa = var(nfa_innov)`; `cov1 = cov(-nx_innov, nfa_innov)`; `cov2 = cov(-val_innov,…)`;
`cov2A/B/C` for rf/rff/seign sub-components; `cov3 = cov(nfa_resid_innov, nfa_innov)`
(`:122-128`). Stored as rows 1–6 of `table_7_moms_tmp` (`:131-136`).

**Table 10 (portfolios / risk premia)** `calc_moments.m:150-154`:
1. `100*mean(h_ksav/h_sav)`; 2. `100*mean(h_bh_sav/h_sav)`; 3. `100*mean(h_bf_sav/h_sav)`;
4. `mean(corr_rf_spread_m1B)` (jx 62); 5. `mean(corr_rf_spread_m2B)` (jx 63).

### 5.2 `collect_moments.m` — average across sims + parameter rows

- `table_K_moms(:,ccc) = mean(table_K_moms_tmp, 2)` for K ∈ {2,3,4,5,6,7,10} (`:10-18`).
- **Table 7 normalization** (`collect_moments.m:17`): after averaging, rescale so the
  decomposition sums to 100%: `table_7_moms(:,ccc) = 100*table_7/sum(table_7([1,2,6],ccc))`.
- `table_1_params(:,ccc)` — 16 externally-set params (`:22-39`): ies_h, sigma,
  varsigma_h-1/(1+ζ), chi, aalpha, delta, vareps_w, varphi_w, phi_pi_h, -disast_shock,
  100*exp(disast_p), p_rho, sig_dis, exp(omg_mean), omg_rho, corr_omg_dis.
- `table_2_params(:,ccc)` — 15 calibrated params (`:41-57`): zeta, sig_z,
  std_zf*sqrt(1-rho_zf²), rho_zf, chiX, bbeta_f, bbeta_h, sig_omg, gma_f, gma_h, bg_yss,
  chi0_vec(1), chi0_vec(2), 100*tayl_ic_h, 100*tayl_ic_f.

### 5.3 `collect_swap_moments.m` — Table 8 (benchmark only)

Operates on `collected_irfs(:,:,7,ccc)` (fix − fx0, with-rate-pegged) and
`(:,:,6,ccc)` (bg1, floating). Each is **scaled so the impact ω response = −14 bp**:
`scale = -14/tmp(1, irf_idxs.omg)` (`collect_swap_moments.m:4,18`). Then
`table_8_moms(:, j, ccc)` (j=1 fixed, j=2 floating):
1. `tmp(1, log_E)`; 2. `tmp(1, retA)+tmp(1, log_infl_h)`; 3. `max(tmp(:, log_yh))`;
4. `max(tmp(:, log_yf))`; 5. `tmp(1, log_qx)`; 6. `tmp(1, excA)`; 7. `tmp(1, retA)`;
8. `tmp(1, nfa_y)`; 9. `tmp(1, nfa_y) - tmp(1, nx_y)` (`:6-28`).

---

## 6. Figures

### 6.1 `make_figure.m` — helper

`fig_handle = make_figure(fig_size, series_set, series_set_2, x_vec_array, y_vec_array,
title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, datetick_dmy, thickness, colorscheme)`
(`make_figure.m:7`).
- `fig_size = [nrows, ncols]` subplot grid; iterates `nnn=1:nrows*ncols`.
- `x_vec_array{nnn}`, `y_vec_array{nnn}` are matrices `[T × n_lines]`. `series_set` = which
  columns to plot; `series_set_2` = which columns to include when computing y-limits.
- Colors keyed on number of lines (`:32-56`): 1→darkblue; 2→{lightblue,darkblue} (or swapped
  for `colorscheme='modelvsdata'`); 4→{lightblue,darkblue,red,green}; else 3-color ramp.
- Draws a zero line, then each series; title/ylabel in LaTeX; legend only on subplot
  `legend_loc1` at location `legend_loc2`; optional `datetick`; auto y-limits with 5% pad
  (NaN-safe) (`:88-102`).

### 6.2 `create_figures.m` — IRF figures (after the ccc loop)

Window: periods `1:29` (`start_t=1, len=29`, `create_figures.m:9-10`). Each block selects an
`irf_idx` (1=p, 2=ω, 3=z, 4=disaster, 5=zf — i.e. the k-index into `collected_irfs`) and a
`series_to_plot` list of `irf_idxs.*`, builds per-subplot y/x/title arrays from
`collected_irfs(start_t:len, col, irf_idx, :)` (the `:` = all calibrations).

| fig | irf_idx (shock) | grid | series_to_plot |
|---|---|---|---|
| **2** | 1 (p) | 2×3 | p, excA, excrf, log_qx, thth, log_yh (`:17`); printed for ccc=ix_bm (`:31-33`) |
| **11** | 1 (p) | 3×4 | p, Erh, Erf, excrk, log_qk, excrf, log_qx, log_E, thth, log_kh, nfa_ybar, val_ybar (`:41`) |
| **12** | 1 (p) | 3×4 | nx_ybar, log_ch, log_cf, log_x, log_infl_h, log_wh, log_lh, log_yh, log_infl_f, log_wf, log_lf, log_yf (`:65`) |
| **3** | 2 (ω) | 2×3 | same 6 as fig 2 but with `irf_idxs.omg` first; compares calibs {ix_symm_flex, ix_symm, ix_bm} with legend (`:94,115-121`) |
| **5** | 2 (ω) | 2×3 | ih, excA, excrf, log_qx, thth, log_yh; compares {ix_tayl_rho, ix_tayl_y0, ix_bm} (`:137,146-151`) |
| **19** | 2 (ω) | 3×4 | omg, Erh, Erf, excrk, log_qk, excrf, log_qx, log_E, thth, log_kh, nfa_ybar, val_ybar (`:157`) |
| **20** | 2 (ω) | 3×4 | nx_ybar block (same 12 as fig 12) (`:182`) |
| **15** | 3 (z) | 3×4 | log_z, Erh, Erf, excrk, log_qk, excrf, log_qx, log_E, thth, log_kh, nfa_ybar, val_ybar (`:213`) |
| **16** | 3 (z) | 3×4 | nx_ybar block (`:239`) |
| **13** | 4 (disaster) | 3×4 | log_z, … (same 12 as fig 15) (`:269`) |
| **14** | 4 (disaster) | 3×4 | nx_ybar block (`:295`) |
| **17** | 5 (zf) | 3×4 | log_zf, … (same 12) (`:325`) |
| **18** | 5 (zf) | 3×4 | nx_ybar block (`:351`) |

(For figs 2/11/12/19/20/13/14/15/16/17/18 a single-calibration version is also written to
`addl_figs/` for every ccc; the numbered figure is the `ccc==ix_bm` one.) Figs 3 and 5 are
the multi-calibration comparison plots.

### 6.3 `create_recession_fig.m` — figs 4, 6, 7 + table_A1 (benchmark only)

Reads `../src/data/modelvsdata.csv` (skip 1 header row). Columns used (`:24-34`): col4=p,
col5=ω(/100/100), col6=Δlog q, col7=Δlog c_h, col8=Δlog c_f, col9=Δlog y_h, col10=Δlog y_f,
col11=nfa/y, col12=Δnfa/y, col13=nx/y, col14=i_h, col15=i_f-i_h(→i_f). Date vector
`datenum(1995,1:3:309,1)` (quarterly). Data p rescaled by `param_input_vec(1,1)`; data ω
demeaned and rescaled to `param_input_vec(4,1)` std (`omg_scaling`, `:17-23`).

Model series: `temp_series = samp1_series` (and samp2, samp3) → `extract_series` (`:42,124,141`).
Window `start_idx..end_idx` = 2006Q2..2011Q3 (`:37-38`).

- **fig 4** (`:177-179`): 2×3, "Model" vs "Data" using samp1 (lines 2,1) for
  p, ω, Δlog y, Δlog y*, Δnfa/y, (Δnfa−nx)/y.
- **fig 6** (`:181-184`): 2×3, "only p / only ω / Model" using samp2(only-p), samp3(only-ω),
  samp1(model) (lines 3,4,2).
- **fig 7** (`:186-191`): 1×3, i, i*, θ — model vs data.
- **table_A1.tex** (`:200-220`): Q3'07–Q3'09 cumulative changes (Data vs Model) of Δlog y,
  Δlog y*, Δnfa/y, using `samp1_series` and cumulated data growth.

### 6.4 `create_safety_fig.m` — fig 10 (benchmark only)

Reads `../src/data/fig_10.csv` (skip header). 18 cols = 6 series × {point, lower CI, upper CI}
(`:14-36`): ω, log q, Δi (i*−i), excA, y_h, y*−y_h. Plots the ω-shock IRF
`collected_irfs(:,:,2,1)` overlaid on data with CI bands (2×3). The model IRF is rescaled by
`irf_scaling = omg_data(1)/irf_vec(1)` (and ω undone by `/omg_scaling`) so the ω panels line up
(`:56-58`); Δi panel uses `4*(if - ih)` (`:91`); last panel uses `log_yf - log_yh` (`:144`).
Output: `fig_10.pdf` (`:155`).

---

## 7. `create_tables.m` — table contents (1–10 + A2)

All write `../output/tables/table_<n>.tex`. `comp_idxs` selects which calibration columns
appear.

- **Table 1** (`:8-37`): external params, benchmark only. 16 rows from `table_1_params`
  (IES, trade elast σ, home bias, Frisch, α, δ, ε, χ^W, φ^π, disaster shock, disaster risk,
  ρ^p, σ^p, ω skewness, ρ^ω, ρ^{pω}). Cols: Description, Value, Notes.
- **Table 2** (`:40-142`): targeted moments, benchmark only. 15 rows; cols: param name,
  description, **Value** (`table_2_params`), moment name, **Target** (hardcoded
  `targeted_values`, `:99-115`), **Model** (`table_2_moms`).
- **Table 3** (`:145-174`): comovements. `comp_idxs=[ix_bm, ix_no_omg, ix_symm]`. Cols:
  Data (hardcoded `[-0.17,0.23,1.38]` with SEs), Model, No ω, γ=γ*. Rows = the 6
  `table_3_moms` (3 regression betas + 3 memo portfolio positions).
- **Table 9** (`:176-188`): SAME 6 moments as Table 3 but `comp_idxs=[ix_bm, ix_tayl_y0,
  ix_tayl_rho]` (alternative Taylor rules); no data column.
- **Table 4** (`:190-218`): `table_4_moms`, cols Data (hardcoded `[2.9,33.6,15.5,3.8,3.8,0.1]`)
  / Model / No ω / γ=γ*. 6 rows (5 std-dev moments + 1 correlation).
- **Table 5** (`:220-237`): `table_5_moms`, data `[0.59,0.81]`. 2 rows (σ Δlog y, σ Δlog y*).
- **Table 6** (`:240-264`): `table_6_moms`, data `[11.0,1.0,10.9,-2.8,-3.2,0.4]`. 6 rows
  (3 std + 3 mean NFA/nx moments).
- **Table 7** (`:267-283`): `table_7_moms` rows 1,2,6 (nx, val, residual; normalized to %),
  `comp_idxs=[ix_bm, ix_no_omg, ix_symm]`. 3 rows, shares of Var(Δnfa).
- **Table A2** (`:285-307`): full 6-row external-adjustment decomposition (`table_7_moms`
  rows 1–6: nx, val, val_rf, val_rff, val_seign, residual).
- **Table 8** (`:310-322`): swap lines, benchmark only, from `table_8_moms(:,:,ix_bm)`.
  Impact rows (log E, log P·rᵉ, Δnfa/y, Δ(nfa−nx)/y) with hardcoded data −72bp / +151bp on
  first two; peak rows (log y, log y*). Two model columns = fixed vs floating nominal rate.
- **Table 10** (`:325-338`): `table_10_moms`, `comp_idxs=[ix_no_omg_symm, ix_no_omg,
  ix_nocorr_nobg, ix_nocorr, ix_bm]` (5 columns, building up the model). 5 rows: k/a, b_H/a,
  b_F/a, ρ(rf-spread, m), ρ(rf-spread, m*+Δlog q).

---

## 8. `numerical_check.m` and `print_params.m`

- **`numerical_check.m`**: writes a "NUMERICAL VERIFICATION" block to the per-run text file
  comparing simulated state means/std vs. grid means/devs for k, θ_H, w_h, w_f, ρ·i_h, ρ·i_f,
  ω (`:4-33`). Diagnostic only — checks the simulation stays inside the Smolyak grid box. No
  figure/table output.
- **`print_params.m`**: dumps to the per-run text file (`:1-39`) the full 65-param vector
  (names + values), the steady-state-calibrated labor disutilities `chi0_vec`, steady-state
  values `yh_ss`/`aggr_css`, and the log-normal process parameters (ω/p innovation std, mean)
  plus their target inputs from `param_input_vec(:,ccc)`. Also assigns `trg_prms.chi0_h/_f/yh_ss`
  used later in `extract_series` SDF formulas. No figure/table output.

---

## 9. Reproduction gotchas for the Python port

1. **Stacked-row indexing.** `jx_*` are absolute rows of the concatenated
   `[shocks; states; vars; value]` matrix, in the order written by `read_series.m:38/46/48`.
   The Python reader must concatenate the four text files in exactly that order and the value
   block must be appended only for the 151-row series (not for samp/bg/fix/fx0, which are 137).
2. **De-trending convention.** Almost every quantity is `stored_value * exp(cumsum(z_shock))`.
   Returns instead use `exp(dis)` corrections. Replicate the `z_series` / `dis_series` logic
   (`extract_series.m:13-26`) before anything else.
3. **Boundary padding.** Differenced/shifted series are padded with their own first element
   (`[x(1); x]`). Means/stds in `calc_moments` slice off boundaries (e.g. `(2:end-1)`,
   `(3:end-1)`); index windows must match exactly.
4. **Table 7 expectations.** Uses the dedicated jx_val_* (138–151) realized-vs-lagged-
   expectation columns and is computed **only on the disaster sample** (`dis=true`); then
   renormalized to sum to 100% across rows {1,2,6} in `collect_moments.m:17`.
5. **IRF scalings.** Table 8 scales IRFs to a −14bp ω impact; fig 10 rescales the model IRF to
   match the data's first ω point and undoes `omg_scaling`; `collected_irfs(:,:,7)` is the
   `fix − fx0` difference. All IRFs are in basis points (×10000) and start at period 3.
6. **Regression residual-correlation moment** (Table 2 row 8) chains three regressions
   (`regress` with intercept) and correlates two residual vectors — reproduce the exact sample
   windows (`5:end`, `4:end-1`) and the dividend-smoothing offset.
