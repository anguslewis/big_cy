# Structural map of `mod_calc.f90` (Kekre & Lenel 2024 AER)

Source (read-only): `…/Safety/Quantitative Model/Code/src/fortran/mod_calc.f90` (5328 lines).
Authors: Rohan Kekre & Moritz Lenel. Public entry points: `calc_steady`, `calc_sol`.

This file is the reference for a tensor-native (PyTorch/JAX) Python re-implementation. All
line citations are `mod_calc.f90:line`.

---

## 0. Conventions, state space, dimensions

- **2 countries** `n_I = 2` (1 = Home/US, 2 = Foreign/RoW). `iii` indexes country/agent.
- **9-dim Smolyak state grid**, columns of `state_grid(idx_*, sss)`:
  1. `idx_k`   — aggregate capital `k`
  2. `idx_thh` — Home wealth share θ_h
  3. `idx_zf`  — relative TFP (foreign), log
  4. `idx_wh`  — Home wage state (in bundle units; see below)
  5. `idx_wf`  — Foreign wage state
  6. `idx_dis` — disaster probability `p`, stored as log → `p_dis = exp(state_grid(idx_dis,·))`
  7. `idx_ih`  — lagged Home nominal rate `i_h` (level − 1 stored)
  8. `idx_if`  — lagged Foreign nominal rate `i_f`
  9. `idx_omg` — convenience yield ω, log → `omg = exp(state_grid(idx_omg,·)) + omg_shift − b_lmbd·bond_supply_shock`
- `sss` = grid-point index (1…`n_states`). `n_active_dims` = number of stochastic dims actually
  carried in the transition (the others are deterministic/degenerate). `vector_mus_dimensions(d) > 0`
  flags an active dim — the `counter` bookkeeping that recurs everywhere maps active dims to columns
  of `next_state_mat`.
- **Quadrature**: `n_quad` nodes. Node `1…n_quad-1` are the non-disaster Gauss-Hermite nodes with
  weights `quad_weight_vec`; node `n_quad` is the **disaster node**. The combined weight vector is
  ```
  big_weight_vec(1:n_quad-1) = quad_weight_vec*(1 - p_dis)
  big_weight_vec(n_quad)     = p_dis
  ```
  (e.g. `:2256`). Expectations over t+1 are always `sum(big_weight_vec * x_nxt)`.
- **4 shocks** in `shock_grid`/`trans_shock_vec`: indices `sidx_z` (world TFP growth), `sidx_zf`
  (relative TFP), `sidx_p` (disaster prob), `sidx_omg` (convenience). `dz_vec` = growth-shock node
  vector; `dz_vec_adj` = its detrending adjustment. The disaster node carries an extra depreciation
  hit applied via `exp(dz_vec(n_quad))`.
- **Scaling**: value `v` and consumption are stored **scaled by representative-agent wealth μ̃**.
  Concretely `v_new = objf · ((wealth + q_l_ss)/tot_wealth_ss)^{-1}` (`:2860`) and the marginal
  utility `mc_new = util_c_deriv · ((q_l_ss + wealth)/tot_wealth_ss)^{1/IES}` (`:2957`). The constants
  `q_l_ss_vec` (PV of labor endowment) and `tot_wealth_ss_vec` (total SS wealth) come from steady state.
- `dp` = double precision (`base_lib`). `conv_crit = 1E-8` (module parameter, `:24`).

**Module-level saved state** (set by `calc_steady`, read by `calc_sol` and helpers), `:22-23`:
`v_ss_vec(n_I)`, `mc_ss_vec(n_I)`, `k_ss`, `s_ss`, `rf_ss`, `l_ss_vec(2)`, `q_ss`,
`tot_wealth_ss_vec(n_I)`, `q_l_ss_vec(n_I)`, `c_cost_ss_vec(n_I)`, `w_ss_vec(2)`.

---

## 1. Subroutine / function inventory (in file order)

### `calc_steady()` — `:28-200`
Args: none (operates on module globals + `mod_param`). Computes the deterministic steady state to
center the Smolyak grid. Sets module SS vars and the grid means/devs (`k_grid_mean`, `wh_grid_mean`,
`wf_grid_mean`, `k_grid_dev`, …) in `mod_param`, writes `extra_data.csv`. See §2.

### `calc_sol()` — `:204-1937`
Args: none. The master solver. Backward value-function iteration to a stationary recursive
equilibrium, then bond pricing, valuation accounting, and IRF construction. See §3. Key local arrays
declared `:226-292`. `n_nxt = 8 + 2*n_I = 12` interpolated next-period variables (`:224`).

### `calc_unexpected_transition(...)` — `:1940-2174`
Signature (`:1940-1958`):
```
in:    nxt_mat(1,n_nxt), n_nxt, nxt_mat_2(1,3), sss, c_spending_vec(n_I), s,
       l_aggr_vec(2), q_current, infl_vec(2), trans_shock_vec(n_shocks),
       bond_supply_shock, bond_supply_shock_nxt
inout: nom_i_vec(2), share_vec(n_I), bF_share_vec(n_I)
out:   theta_nxt(1)
```
A **single-node** ("1" instead of `n_quad`) version of the equilibrium step used to compute the
*impact* of an unexpected (MIT) shock `trans_shock_vec`: takes today's policies as given, applies the
shock to the realized state, and returns the implied next-period Home wealth share `theta_nxt`. Used
to build the IRF transition matrix `next_state_irftransition_mat`.

### `calc_equilibrium_and_update(...)` — `:2176-3099`
Signature (`:2176-2197`):
```
in:    nxt_mat(n_quad,n_nxt), n_nxt, nxt_mat_2(n_quad,3), sss, outer_iter,
       c_spending_vec(n_I), s, l_aggr_vec(2), q_current, infl_vec(2),
       bond_supply_shock, bond_supply_shock_nxt
inout: nom_i_vec(2), share_vec(n_I), bF_share_vec(n_I), M_vec_mat(n_quad,n_I)
out:   inv_h, inv_f, q_new, s_new, k_next_new, c_spending_vec_new(n_I), infl_new(2),
       v_new(n_I), mc_new(n_I), l_aggr_new(2), w_choice_new(2), theta_nxt(n_quad)
```
**The per-grid-point equilibrium step** — the economic core. Solves the within-period equilibrium and
produces updated policies/prices/values for one state `sss`, given interpolated continuation values
`nxt_mat`. The 10 documented steps live here. See §4.

### `calc_bond_prices(...)` — `:3101-3501`
Signature (`:3101-3120`): like `calc_equilibrium_and_update` but `share_vec`/`bF_share_vec`/`nom_i_vec`
are `intent(in)`, plus `in: nxt_bond_prices(n_quad,3)`, `out: who_prices_bond(3)`,
`bond_prices_out(3)`, `E_rb(6)`, `results_vec(n_interp-6*12)`. Given the **converged** policy
functions, recursively prices a ladder of nominal bonds (recursion `q^{(n+1)} = E[M · q'^{(n)}/π]`),
computes expected bond returns `E_rb`, identifies which agent is the marginal pricer
(`who_prices_bond`), and packs all reporting statistics (`results_vec`) for output. Pure
post-processing — does not feed the iteration.

### `calc_valuation(...)` — `:3503-3774`
Signature (`:3503-3522`): `out: valuation_vec(7)`, `in: nxt_valuation(n_quad,7)`. Accumulates external
balance-sheet / valuation accounting (NFA, net exports `nx`, return-on-assets vs return-on-liabilities
decompositions) by the same recursive-expectation trick. The 7 slots are valuation terms iterated to a
fixed point in `calc_sol` (`:807-884`); slot 4 (`nx`) is shifted to its conditional expectation in a
dedicated pass (`:913`). Pure post-processing.

### `calc_excess_bond_nom(...)` — `:3776-3868`
Signature (`:3776-3794`):
```
in:    aggr_wealth, omg, v_vecs(n_quad,n_I), rf_vec(n_quad,2), rk(n_quad),
       seignorage_transfer(n_quad,n_I), bbeta_adj(n_I), c_cost_vec(2),
       c_cost_nxt_vec(n_quad,2), q_current, w_choice_vec(2), l_vec(2), sss,
       big_weight_vec(n_quad), c_temp(n_I), c_spend_temp(n_I), mc_nxt_mat(n_quad,n_I)
inout: share_guess_vec(n_I), bF_share_guess_vec(n_I)
out:   k_temp_vec(n_I), b_temp_vec(n_I), excess_b, constraint_binding_vec(n_I)
```
**Aggregate Home-bond excess demand** at a candidate nominal rate. For each agent: holds
consumption/savings fixed, computes natural leverage bounds (no-negative-payoff), then calls
`calc_portfolio_share` to get the optimal capital-vs-bond share. Returns total bond demand
`excess_b = sum(b_temp_vec)` (zero at equilibrium i_h). This is the residual the Brent loop in
Step 4 zeros.

### `calc_portfolio_share(...)` — `:3871-4082`
Signature (`:3871-3885`): `inout: share_guess`; `out: constraint_binding`. Solves one agent's optimal
risky/safe split `share` (= fraction in capital vs. Home bond) by zeroing `portfolio_foc` with a
hand-rolled **Brent** root-finder (inverse-quadratic + secant + bisection fallback, `:4008-4073`).
First clips bounds and detects corner solutions (`constraint_binding = ±1`). `r1=rf_home/(1-omg)`,
`r2=rk`.

### `calc_excess_foreign_bond_nom(...)` — `:4086-4175`
Signature (`:4086-4103`): mirror of `calc_excess_bond_nom` for the **Foreign-currency bond**. Returns
`excess_b_home = -sum(bh_temp_vec)` — the net market-clearing residual for the spread `i*−i`. Calls
`calc_bond_portfolio_share` per agent. `share_guess_vec` (capital share) is `intent(in)` here; the
unknown is the foreign-bond share `bF_share`.

### `calc_bond_portfolio_share(...)` — `:4178-4360`
Signature (`:4178-…`): same Brent skeleton as `calc_portfolio_share` but zeros the FOC for the
**foreign-vs-home bond** margin: `portfolio_foc(..., r1=rf_foreign, r2=rf_home, ...)` (`:4212` etc.).
Unknown is `bond_share` (= `bF_share`).

### `util_fun(consumption, labor, iii, util, util_c_deriv, labor_part)` — `:4364-4386`
EZ-with-GHH-labor period utility. Handles `IES = 1` (log) as a special case. Returns `util`, its
consumption derivative `util_c_deriv`, and the labor-adjustment factor `labor_part`.

### `portfolio_foc(next_period_wealth, v, mc, big_weight_vec, r1, r2, c_cost_growth, gma, ies)` — `:4388-4400`
Function → scalar `foc`. The portfolio first-order condition residual:
```
sdf = v^(1/ies − gma) · mc
foc = Σ big_weight · sdf · next_period_wealth^(−gma) · (r1 − r2)/c_cost_growth   (scaled by |Σ big_weight·sdf|)
```
This is THE residual all four portfolio routines zero. **Critical for the Python port to replicate exactly.**

### `portfolio_return(bond_share, foreign_share, rh, rf, rk)` — `:4402-4410`
Function → `r_alpha(n_quad)`:
`(bond_share − foreign_share)·rh + foreign_share·rf + (1 − bond_share)·rk`.
Total portfolio return given capital share, foreign-bond share, and the three gross returns.

### `calc_bundle(c_h, c_f, varsigma)` — `:4414-4422`
CES home/foreign consumption aggregator with elasticity `sigma`, weight `varsigma`.

### `calc_equilibrium_ifixed(...)` — `:4424-5175`
Signature like `calc_equilibrium_and_update` plus `out: excess_b_out` (`:4424-4444`). Variant of the
core step for the **fixed-nominal-rate** scenario (ZLB / forward-guidance IRFs): the Home nominal rate
is *held fixed* and the Home-bond market does **not** clear via `i_h`; instead `excess_b_out` is
returned so an outer Brent loop in `calc_sol` can solve for inflation that restores clearing. Used only
on the long-lived-IRF "held fixed" path.

### `update_ifixed(...)` — `:5178-5324`
Signature `:5178-5193`. Companion to `calc_equilibrium_ifixed`: given a candidate Home inflation
`infl_h` and fixed `ih_fixed`, iterates the within-period block to convergence and returns the implied
bond-market residual `excess_b` (plus refreshed `v_new`, `mc_new`, `w_choice_new`). This `excess_b` is
the residual zeroed by the Brent loop at `calc_sol:1539-1618`.

---

## 2. `calc_steady` — deterministic steady state (`:28-200`)

Goal: solve the non-stochastic SS to (a) set module SS scalars/vectors and (b) set Smolyak grid
centers/half-widths in `mod_param`.

1. **Primitives** (`:45-62`): `zf_ss = [1, zeta]`; wealth-weighted average discount/IES
   `bbeta_avrg = Σ tht_trgt·bbeta`, `IES_avrg = Σ tht_trgt·IES`. SS returns
   `rk_ss = rf_ss = 1/bbeta_avrg − 1`. Labor targeted: `l_ss_vec = l_target`. Init `s_ss = 1`,
   `theta_ss_vec = tht_trgt_vec`.
2. **Outer fixed point on terms of trade `s_ss`** (`do while diff > sqrt_eps`, `:66-114`). Each pass,
   given `s_ss`, compute in closed form: capital price `q_ss` (`:70`), profit `pi_ss = (rk_ss+ddelta)·q_ss`,
   per-country deployed capital `kappa_h_ss`,`kappa_f_ss` (`:76-77`) and `k_ss = kappa_h+kappa_f`,
   investments `inv_h`,`inv_f`, outputs `y_h_ss`,`y_f_ss`, aggregate consumptions, `aggr_wealth`,
   wages `w_ss_vec` (`:94-97`), per-agent home/foreign consumption `c_h_ss_vec`,`c_f_ss_vec`
   (`:100-105`). Update `s_ss_new` from the foreign-good market-clearing ratio (`:107`); **damped**
   update `s_ss += 0.01·(s_ss_new − s_ss)` (`:111`). Hard stop after 10000 iters.
   ⚠ The damping is very slow (0.01) — Python can use a proper Newton/Brent on the same residual.
3. **Bundle cost** `c_cost_ss_vec` (CES price index, `:117-118`).
4. **Per-agent labor-disutility calibration `chi0`** (`:122-157`): inner fixed point per agent
   (`do while diff > 1E-12`, damped 0.5, `:128-143`) solving the labor FOC so SS labor hits
   `l_target`. Then set `v_normalization_vec`, `v_ss_vec`, `mc_ss_vec`, and store `c_bundle_ss_vec`.
   `chi0_vec` written into `mod_param`.
5. **Wealth aggregates** (`:159-161`): `q_l_ss_vec = w_ss_vec/(1−bbeta_vec)` (PV of labor endowment);
   `tot_wealth_ss_vec = q_l_ss_vec + theta_ss_vec·k_ss·(1+rk_ss)`. These two drive the μ̃-scaling.
6. **Grid centering** (`:181-187`): 
   `k_grid_mean  = k_grid_adj·(kappa_h_ss + kappa_f_ss)`,
   `wh_grid_mean = w_grid_adj·w_ss_vec(1)`, `wf_grid_mean = w_grid_adj·w_ss_vec(2)`,
   `k_grid_dev   = k_dev_param·k_grid_mean`, `wh_grid_dev = w_dev_param·wh_grid_mean`,
   `wf_grid_dev  = w_dev_param·wf_grid_mean`. (θ/zf/dis/ih/if/omg means+devs are set elsewhere in
   `mod_param` setup; `calc_steady` only sets k & wages.)
7. Print SS table; write `extra_data.csv` with `chi0(1..2), sig_dis, sig_omg, dis_grid_mean,
   omg_grid_mean, y_h_ss, aggr_css` (`:193-198`).

Root-finders used: two **damped fixed-point** loops (s_ss; chi0). No NAG.

---

## 3. `calc_sol` — outer backward iteration + post-processing (`:204-1937`)

### 3.1 Setup / initialization (`:294-429`)
Initialize policy/price/value matrices at SS (`l_aggr_mat=1`, `q_mat=q_ss`, `s_mat=s_ss`,
`infl_mat=1`, `share_mat=0`, `bF_share_mat=0` — "start from only dollar holdings", `v_mat=v_ss`,
`mc_mat=mc_ss`, etc.). Build the initial next-state transition `next_state_mat(n_quad, n_active_dims,
n_states)` dim-by-dim with the `counter`/`vector_mus_dimensions` bookkeeping (`:356-419`), clip to
`[-1,1]` (`:422-426`).

### 3.2 Main backward-iteration loop (`do while diff > conv_crit .and. outer_iter < max_iter`, `:437-709`)
Each iteration:
1. Pack the 12 interpolated variables into `interp_input_mat(n_states, n_nxt)` (`:449-460`):
   columns `[v_1, v_2, mc_1, mc_2, s, q, l_1, l_2, infl_1, infl_2, c_spend_1, c_spend_2]`.
2. **Solve Smolyak coefficients** via NAG `F07ABF` (linear solve `smol_polynom · coeffs = interp_input`,
   `:465`). `smol_coeffs(n_states, n_nxt)`.
3. **OpenMP parallel DO over grid points** `sss = 1…n_states` (`:468-593`, `SCHEDULE(static)`):
   - Evaluate Smolyak basis at the `n_quad` next-states: `polyn_points = Smolyak_Polynomial2(...)`
     (`:493`); interpolate continuation values `nxt_mat = polyn_points · smol_coeffs` via BLAS `DGEMM`
     (`:498`).
   - Append non-interpolated next-period quantities `nxt_mat_2 = [k_next, next_zf, next_omg]` (`:502-504`).
   - Load current guesses for `sss`, call **`calc_equilibrium_and_update`** (`:518`).
   - Store all `_new` outputs and rebuild this point's `next_state_mat_new` columns (`:524-589`).
4. Clip `next_state_mat_new` to `[-1,1]` (`:596-601`).
5. **Convergence** `diff` = max over relative-log changes of v, mc, c_spend, k_next, l, w, q, infl and
   absolute changes of share, nom_i, next_state (`:604-614`).
6. **Dampened update** of every object (`:676-704`):
   | object | weight |
   |---|---|
   | `v_mat`, `mc_mat` | 1.0 (full) |
   | `next_state_mat` | 0.1 |
   | `q_mat` | 0.2 |
   | `s_mat` | 0.1 |
   | `l_aggr_mat` | 0.2 |
   | `w_choice_mat` | 0.2 (in log space: `w·(1+0.2·log(new/old))`) |
   | `c_spending_mat` | 0.2 |
   | `nom_i_mat(1,·)`, `nom_i_mat(2,·)` | 0.1 |
   | `infl_mat` | 0.05 |
   | `k_next_mat` | 0.2 |
   | `share_mat` | 0.1 |
   | `bF_share_mat` | 0.0 for `outer_iter<101`, then 0.05 |

   ⚠ **Phase-in schedule (important):** the model is solved in stages controlled by `outer_iter` inside
   `calc_equilibrium_and_update`: bond-market clearing for `i_h` only runs when `outer_iter > 10`
   (else `nom_i` set from a no-trade Euler eq, `share=0`); foreign-bond/spread clearing only when
   `foreign_trading==1 .and. outer_iter > 100`; foreign-bond *damping* starts at iter 101. So a faithful
   Python port must replicate this staged activation, not just the damping weights.

### 3.3 Bond pricing pass (`:715-802`)
For each bond maturity `bbb = 1…n_bond`: solve coeffs for current bond-price slab via `F07ABF`
(`nrhs2=3`), then OpenMP loop calling **`calc_bond_prices`** to get `q_bond_mat(:,:,bbb+1)`,
`E_rb`, `who_prices_bond`, `results_vec`, packed into `results_mat`. Maturities `≤4` and the top two
are stored (`bbb_counter` logic).

### 3.4 Valuation pass (`:804-963`)
Iterate the 7 valuation terms to a fixed point over `bbb = 1…n_valuation-1` (`F07ABF` `nrhs3=7`,
OpenMP loop calling **`calc_valuation`**), then a pass that converts `nx` (slot 4) to its conditional
expectation (`:898-917`), and a final pass building `valuation_mat_store(14, n_states)` = the 7 current
terms plus their 7 one-period-ahead expectations (`:931-963`).

### 3.5 Output (`:965-982`)
Stream-write `next_state_mat.dat`, `results_mat.dat`, `valuation_mat.dat`.

### 3.6 IRFs (`:985-1204`) — standard one-period shocks
For each shock in `irf_indices`: copy the solution into `*_irf_*` arrays, build the **unexpected
transition** matrix `next_state_irftransition_mat`. For the disaster shock it just selects the disaster
node (`:1012`); otherwise it sets `trans_shock_vec`, re-solves coeffs, and iterates a damped loop
(weight 0.2) calling **`calc_unexpected_transition`** to fixed point (`:1049-1191`). Writes
`next_state_shocktrans_mat_<i>.dat`.

### 3.7 Long-lived IRFs (`:1207-1936`) — only if `run_bg == 1`
For each shock in `irf_indices_ll`, count **time backwards** `ttt = 1000+irf_ll_length … 1001`
(`:1239`). Two regimes:
- **Normal (`ttt > 1000 + irf_ll_hold_fixed`)**: damped inner loop
  (`do while diff>conv_crit .and. outer_iter<max_iter_irf`, `:1267-1503`) calling
  `calc_equilibrium_and_update` with `outer_iter=10000` (so all blocks active) and passing
  `bond_supply_shock`. Damping uses **clipped step sizes** (`where abs(Δ)>tol` → tiny signed step,
  else normal weight) for q/s/l/w (`:1471-1493`); nom_i 0.5, infl 0.2, k 0.2, share 0.5,
  bF_share 0.05.
- **Held-fixed (else branch, `:1504-…`)**: the **fixed-i_h / ZLB path**. Per grid point, an outer
  hand-rolled **Brent** loop (`initBrent`/`updateBrent` from `base_lib`, `:1539-1618`) brackets and
  zeros the bond-market residual `excess_b` returned by `update_ifixed`, solving for inflation `infl`
  given fixed `ih_fixed = 1.01`. (`calc_equilibrium_ifixed`/`update_ifixed` carry the economics.)

### OpenMP
Every grid-point loop is wrapped in `!$OMP PARALLEL DEFAULT(NONE) … !$OMP DO SCHEDULE(static)` with
explicit SHARED/PRIVATE lists (`:468`, `:728`, `:813`, `:889`, `:922`, `:1052`, `:1270`, `:1506`).
In Python these become **vectorized batch ops over the `n_states` axis** (the loop body is independent
across `sss`). `omp_set_num_threads` is called twice (commented context).

### Iteration counts
`max_iter` (main), `max_iter_irf` (long-lived IRF) come from `mod_param`. Convergence tolerance
`conv_crit = 1E-8` throughout. First 10 iterations of the IRF loops force `diff=1` to avoid premature exit.

---

## 4. `calc_equilibrium_and_update` mapped to the 10 steps (`:2176-3099`)

Preamble (`:2239-2310`): store incoming `nom_i_vec_in` (needed for Taylor rule); read current state
`k_aggr`, `zf_vec=[1, zeta·exp(zf)]`, `p_dis`, `big_weight_vec`, wage states `w_current`, lagged rates
`ih_last`,`if_last`; set `omg`, `bbeta_adj` (stationarity nudge on θ, `:2309-2310`).

| # | Step (code section) | Solves for | Held fixed | Mechanics / key lines |
|---|---|---|---|---|
| 1 | **Current production** (`:2273-2302`) | `kappa_vec` (capital allocation across countries), `pi_current`, `y_current`, `w_choice_vec`, `P_div_P_h`, `aggr_wealth`, `wealth_vec` | labor `l_aggr_vec`, `q_current`, `s` | capital split equalizes returns (`:2278`); profit `pi=(1/s)·α·(zf·l/κ)^{1-α}`; wages in home-good units; `aggr_wealth = Σ κ·((1−δ)q + π)`. |
| 2 | **Next production → capital returns** (`:2316-2396`) | continuation `v_temp,mc_temp`, `s_nxt,q_nxt,l_aggr_nxt,infl_nxt,c*_spending_nxt` (from `nxt_mat`), `kappa_nxt`, `y_next`, `pi_nxt`, `rk_vec`, `homegood_infl_nxt`, `rf_vec`, `seignorage_transfer` | the interpolated futures themselves | `rk = ((1−δ)q' + π')/q`, disaster node ×`exp(dz)` (`:2375-2376`). `homegood_infl_nxt` converts bundle inflation to local-good inflation (`:2346-2350`). `rf_vec(:,1)=i_h/π_h'`, `rf_vec(:,2)=(i_h+spread)/π_f'`. Seigniorage rebate when `bg_yss>0` (`:2390-2396`). |
| 3 | **Price of capital `q`** (`:2398-2405`) | `q_new` | everything else | `q_new = (CES inv price)·(next_k/k)^{chiX}` — adjustment-cost capital price. Only updated for next iter. |
| 4 | **Home bond market → nominal rate `i_h`** (`:2407-2616`) | `nom_i_vec(1)` | inflation, consumption, **within-country composition** (`bF_share`,`share`) | Only if `outer_iter>10`. Bracket `excess_b` (from `calc_excess_bond_nom`) by stepping `nom_i` ±1E-2 (`:2447-2504`), then **hand-rolled Brent** (IQI/secant/bisection, `:2530-2598`) zeroing aggregate Home-bond demand. Error if `|excess_b|>1E-2`. |
| 5 | **Foreign bond market → spread `i*−i`** (`:2618-2812`) | `nom_i_vec(2)` (the spread) | overall bond holdings, `share` | Only if `foreign_trading==1 .and. outer_iter>100`. Same bracket-then-Brent structure on `calc_excess_foreign_bond_nom`'s residual (step ±1E-3, `:2650-2804`). |
| 6 | **Value fn + θ transition + aggregate capital savings** (`:2814-2867`) | `v_new`, `theta_nxt`, `savings_vec`, `nxt_wealth_vec`, `k_next_new` | consumption guess | Per agent: clip `c_spend`, compute portfolio return `r_alpha_vec` (`:2841`), `savings`, next-period wealth share, `EV = (Σ big_weight·(v'·share')^{1−γ})^{1/(1−γ)}`, `objf` (EZ aggregator, `:2858`), μ̃-scaled `v_new` (`:2860`). `theta_nxt = nxt_wealth_h/Σnxt_wealth` (`:2864`). `k_next_new = Σ savings·(1−share)/q` (`:2867`). |
| 7 | **Terms of trade `s`** (`:2869-2903`) | `s_new`, `inv_h`, `inv_f` | consumption spending | Inner damped fixed point (weight 0.1, `:2879-2903`): given investment `next_k−(1−δ)k`, split into home/foreign investment via CES; `s_update` from goods-market clearing ratio (`:2889`). |
| 8 | **Consumption-savings FOC** (`:2905-3000`) | `c_spending_vec_new`, `mc_new`, `M_vec_mat` (per-agent SDF over quad nodes), `savings_vec`; labor if `phi_w==0` | shares, prices | Per agent inner damped loop (weight 0.5, `:2943-2990`): Euler residual `temp = Σ M·r_alpha_omg·(c_cost/c_cost')·big_weight`, update `c = c/temp^{IES}` (`:2963-2965`). `M_vec_mat(:,iii) = bbeta·mc'/util_c_deriv·(v'/EV)^{1/IES−γ}·share'^{−γ}` (`:2961`) — the **pricing kernel reused in steps 4,5,9**. If flexible wages (`phi_w==0`) labor solved jointly (`:2973-2979`). |
| 9 | **Labor supply / wage Phillips curve** (`:3002-3066`) | `l_aggr_new`, `w_choice_new` | — | If `phi_w>0`: per-agent damped loop (weight 0.005, `:3016-3043`) on the **wage Phillips curve** (`:3024-3037`) — wage markup vs. Rotemberg adjustment with `M_vec_mat` discounting; then `l = κ·((zf^{1−α}(1−α)/w/s)^{1/α})`. If `phi_w==0`, labor from step 8. |
| 10 | **Inflation via Taylor rules** (`:3068-3074`) | `infl_new(1)`, `infl_new(2)` | — | Invert the Taylor rules: `infl_h = ((i_h/(1+ih_last)^{ρ_i})^{1/(1−ρ_i)}·exp(−tayl_ic_h)/y_h^{φ_yh})^{1/φ_h}`; analogous for foreign with `i_h+spread`. Uses `nom_i_vec_in` (incoming guess), `ih_last`/`if_last`. |

Tail block (`:3076-3097`): if `outer_iter ≤ 10`, override `share_vec=0` and set `nom_i_vec(1)` from a
no-trade Euler equation (max over agents, `:3083`); else keep solved `share_temp`. Foreign analogue for
`nom_i_vec(2)`/`bF_share` (`:3090-3097`).

---

## 5. Root-finding / optimization inventory

**NAG / external linear algebra (not root-finding, but must be matched):**
- **`F07ABF`** (NAG = LAPACK `DGESVX`, expert general linear solve with equilibration & condition
  estimate). Solves `smol_polynom · coeffs = RHS` for Smolyak interpolation coefficients.
  Calls: `:465` (12 policy vars), `:725` (3 bond cols), `:810/887/919` (7 valuation cols),
  `:1044/1261` (IRF passes). **Python: `np.linalg.solve` / `scipy.linalg.lu_solve` on `smol_polynom`.**
- **`DGEMM`** (BLAS matrix multiply) — `nxt = polyn_points · coeffs`. 13 calls. **Python: `@`/`einsum`.**

**Root-finders that zero economic residuals (Python must match residual definitions exactly):**

| Routine / location | Unknown | Residual / equation zeroed | Bracket / guess | Method |
|---|---|---|---|---|
| Step 4 (`:2447-2598`) | Home nominal rate `nom_i_vec(1)` | `excess_b = Σ b_temp_vec` from `calc_excess_bond_nom` (aggregate Home-bond demand) | start = last `nom_i`, step ±1E-2 until sign change | hand-rolled Brent (IQI/secant/bisection, `brent_delta=5E-16`) |
| Step 5 (`:2650-2804`) | spread `nom_i_vec(2)` | `excess_b_home` from `calc_excess_foreign_bond_nom` (net Foreign-bond clearing) | start = last spread, step ±1E-3 | hand-rolled Brent |
| `calc_portfolio_share` (`:4008-4073`) | capital share `share` | `portfolio_foc(r1=rf_home/(1−omg), r2=rk)` = 0 | `[share_low, share_high]` natural-leverage bounds; corner check first | hand-rolled Brent |
| `calc_bond_portfolio_share` (`:4178-4360`) | foreign-bond share `bF_share` | `portfolio_foc(r1=rf_foreign, r2=rf_home)` = 0 | `[safe_low/high_guess_fixed]` clipped bounds | hand-rolled Brent |
| `calc_sol` ZLB else-branch (`:1539-1618`) | inflation `infl` (given fixed `ih_fixed=1.01`) | `excess_b` from `update_ifixed` (bond clearing under fixed rate) | bracket via `initBrent`, expand until bracketed | **library** Brent (`initBrent`/`updateBrent`, tol 1E-14, 1000 it) |

**Damped fixed-point (Newton-free) inner loops** — each is a simple successive-substitution with a
relaxation weight; a Python port can keep them or replace with a proper solver, but the *residual* must
match:
- `calc_steady` s_ss (`:66-114`, weight 0.01) and chi0 (`:128-143`, weight 0.5).
- Step 7 terms-of-trade `s_new` (`:2879-2903`, weight 0.1).
- Step 8 consumption Euler `c` (`:2943-2990`, weight 0.5) + labor (weight 0.2 when `phi_w==0`).
- Step 9 wage Phillips curve `w_temp` (`:3016-3043`, weight 0.005).
- `update_ifixed` block loop (`:5200-…`, ≥20 iters then conv).

⚠ The four hand-rolled "Brent" implementations are byte-for-byte identical control flow (IQI vs.
secant choice, mflag bookkeeping, A/B swap). A single shared `brentq`-style helper can replace all four
in Python; just feed the right residual closure. They are NOT NAG calls despite the NAG-style comments.

---

## 6. Interpolation (off-grid t+1 values)

Continuation values at the `n_quad` future states of grid point `sss` are evaluated by **anisotropic
Smolyak polynomial interpolation**:
1. Build the basis matrix `polyn_points = Smolyak_Polynomial2(next_state_mat(:,:,sss), n_active_dims,
   n_quad, n_states, max_smol_level, smol_elem_ani)` (from `mod_smolyak`). Inputs are the *normalized*
   next-state coordinates in `[-1,1]^{n_active_dims}` and the anisotropic element list `smol_elem_ani`.
2. Coefficients `smol_coeffs` are precomputed once per outer iteration by solving
   `smol_polynom · smol_coeffs = interp_input_mat` with `F07ABF` (`smol_polynom` = basis evaluated at
   the grid nodes themselves → effectively a change of basis from nodal values to spectral coeffs).
3. `nxt_mat = polyn_points · smol_coeffs` (`DGEMM`) gives the `n_quad × n_nxt` interpolated futures.

Non-interpolated futures (`k_next`, `next_zf`, `next_omg`) are carried directly in `nxt_mat_2` because
they are deterministic functions of the current state + shock (`:502-504` etc.).
Single-point version (`polyn_points_trans`, 1 row) is used in the unexpected-transition pass (`:1079`).

**Python:** reimplement `Smolyak_Polynomial2` (anisotropic Smolyak tensor of Chebyshev basis) and keep
the coefficient-solve / multiply structure. Coordinates are clipped to `[-1,1]` before evaluation
(`:422-426`, `:596-601`).

---

## 7. Glossary of key variables

**Prices / quantities (current period):**
- `s`, `s_mat`, `s_new` — terms of trade (price of foreign good in home good); `s_vec=[1,s]`.
- `q_current`, `q_mat`, `q_new` — price of capital. `q_ss` SS value.
- `kappa_vec(2)` — capital deployed in each country (allocation of `k_aggr`).
- `pi_current(2)` — profit per unit capital by country.
- `y_current(2)` — output by country.
- `w_choice_vec(2)`, `w_choice_mat`, `w_choice_new` — wages in home-consumption-good units.
- `w_current(2)` — lagged wage *state* (`state_grid(idx_wh/idx_wf)`), stored in bundle units `w·P_h/P`.
- `P_div_P_h(2)` / `P_div_P_h_nxt` — bundle price ÷ home-good price (CES price index) per country;
  converts spending↔consumption: `consumption = c_spend / c_cost`, `c_cost = P_div_P_h(iii)`.
- `aggr_wealth`, `wealth_vec(n_I)`, `wealth` — total and per-agent beginning-of-period wealth.
- `inv_h`, `inv_f` — home/foreign physical investment.

**Next period (length `n_quad`):**
- `v_temp(:,iii)`, `mc_temp(:,iii)` — interpolated continuation value & marginal utility.
- `s_nxt`, `q_nxt`, `l_aggr_nxt`, `infl_nxt`, `kappa_nxt`, `y_next`, `pi_nxt`, `w_next_choice`.
- `rk_vec` — gross capital return (disaster node scaled).
- `rf_vec(:,1)` — gross **Home** real bond return `= i_h / π_h'`; `rf_vec(:,2)` — **Foreign** bond
  `= (i_h + spread)/π_f'`.
- `homegood_infl_nxt(:,k)` — local-good gross inflation π_h', π_f' (derived from bundle inflation + Δs).
- `seignorage_transfer(:,n_I)` — convenience-yield seigniorage rebate (Home +, Foreign −).
- `omg`, `omg_nxt` — convenience yield ω (today / next), `= exp(state) + omg_shift − b_lmbd·bond_shock`.

**Policies / choices:**
- `share_vec(n_I)`, `share_mat` — **capital share** of savings (1−share in Home bond, roughly).
- `bF_share_vec(n_I)`, `bF_share_mat` — **foreign-currency-bond share** of savings.
- `c_spending_vec(n_I)`, `c_spending_mat` — aggregate consumption *spending* in home-good units
  (consumption = spending / `c_cost`).
- `c_vec(n_I)`, `ch_vec` — consumption (bundle / home-good component).
- `nom_i_vec(2)` — `[i_h, spread]` where spread = i* − i_h. `nom_i_mat(2,:)` is the spread, not i*.
- `savings_vec(n_I)` — savings `= wealth + w·l − c·c_cost`.
- `l_aggr_vec(2)`, `l_aggr_mat`, `l_aggr_new` — aggregate labor by country.
- `infl_vec(2)`, `infl_mat`, `infl_new` — gross bundle inflation by country.
- `k_next`, `k_next_new`, `k_next_mat` — chosen aggregate capital next period (`next_k` = no-disaster).
- `theta_nxt(n_quad)` — next-period Home wealth share (state-2 transition).

**Pricing kernel / value recursion:**
- `M_vec_mat(n_quad, n_I)` — per-agent stochastic discount factor over quadrature nodes
  (`= bbeta·mc'/util_c_deriv·(v'/EV)^{1/IES−γ}·share'^{−γ}`, `:2961`). **Reused across steps 4,5,8,9
  and bond pricing** — the central object.
- `EV` — certainty-equivalent continuation `(Σ big_weight·(v'·share')^{1−γ})^{1/(1−γ)}`.
- `next_period_share` — `(savings·r_alpha + dz·q_l_ss + seigniorage)/tot_wealth_ss` (μ̃-scaled wealth ratio).
- `r_alpha_vec` / `r_alpha_omg_vec` — total portfolio gross return (the `omg` version always includes
  the convenience wedge `/(1−omg)` on the Home bond; the plain version is for resource accumulation).
- `objf` — EZ utility aggregator value (pre-scaling); `v_new = objf·(scaling)^{-1}`.
- `mc_new`, `util_c_deriv`, `labor_part` — from the consumption/labor FOC (and `util_fun`).
- `bbeta_adj(n_I)` — small θ-dependent discount-factor nudge for stationarity (`bbeta_coeff`).

**Residuals / solver scratch:**
- `excess_b`, `excess_b_*` — bond-market excess demand (Home or net foreign), zeroed by Brent.
- `share_FOC*`, `portfolio_foc` output — portfolio FOC residual.
- `nom_i_A/B/C/D/S`, `share_guess_A/B/C/D/S`, `mflag` — Brent bracket bookkeeping.
- `constraint_binding(_vec)` — corner-solution flag (±1 = upper/lower leverage bound binds, 0 interior).
- `big_weight_vec` — disaster-augmented quadrature weights (see §0).

**Grid scaling constants (from `mod_param`, set by `calc_steady`):**
`k_grid_mean/dev`, `wh_grid_mean/dev`, `wf_grid_mean/dev`, `tht_h_grid_mean/dev`, `zf_grid_mean/dev`,
`dis_grid_mean/dev`, `omg_grid_mean/dev`, `ih_grid_mean/dev`, `if_grid_mean/dev` — used to map raw
state values into `[-1,1]` Smolyak coordinates (and back).

---

## 8. Things the Python port must get exactly right (flagged)

1. **`portfolio_foc` residual + the four Brent solvers** (§5) — the entire portfolio/asset-pricing
   block hinges on this scaled FOC and the natural-leverage bounds in `calc_excess_bond_nom`/
   `calc_excess_foreign_bond_nom`. Replicate `r1/r2` ordering and the `/(1−omg)` convenience wedge.
2. **μ̃-scaling of `v` and `mc`** (`q_l_ss_vec`, `tot_wealth_ss_vec`, `:2850-2860`, `:2957`) and the
   `M_vec_mat` definition (`:2961`) — get the EZ exponents (`1/IES`, `γ`, `1−γ`) and the `next_period_share`
   normalization right or nothing converges.
3. **Staged activation by `outer_iter`** (bond clearing >10, foreign trading >100, bF_share damping
   from 101) plus the per-object damping weights (§3.2) — convergence depends on this schedule.
4. **Disaster node handling** — last quadrature node carries `p_dis` weight + extra depreciation
   (`rk_vec(n_quad)*=exp(dz)`), and `big_weight_vec` construction must be exact.
5. **Inflation = inverted Taylor rule** (step 10) using the *incoming* `nom_i_vec_in`, and the
   bundle→local-good inflation conversion `homegood_infl_nxt` (`:2346-2350`).
