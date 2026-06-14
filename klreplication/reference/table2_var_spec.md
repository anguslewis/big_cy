# Kekre & Lenel (2024 AER) — Table 2 variable specification

Precise spec of every series the **Table-2 block** of `calc_moments.m` touches, for a Python
port of the moment pipeline. All file paths relative to
`…/big_cy/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/Code/src/`.

Source files (read-only):
- `matlab/calc_moments.m` — the Table-2 block (lines 12–49).
- `matlab/extract_series.m` — builds named `*_series` from `temp_series`.
- `matlab/get_var_indices.m` — column index map `jx_*` into `temp_series`.
- `matlab/read_series.m` — assembles `temp_series` from solver `.txt` output.
- `fortran/mod_calc.f90` — `calc_bond_prices` (3101–3501), `calc_valuation` (3503–3774), assembly (≈441–977).

---

## 0. The stacked series matrix (`temp_series`)

`read_series.m:36–49` stacks, **row-wise** (one column per variable after transpose in the
caller; `extract_series.m` indexes `temp_series(:, jx_*)` so rows = time, cols = vars):

```
temp_series = [ shock_series ; state_series ; other_vars_series ; value_series ]
                (4 rows)       (smolyak_d=9)   (n_interp=124)       (14)
```

So the `jx_*` integer indices (`get_var_indices.m`) decompose as:

| Block            | jx range | size | source array            |
|------------------|----------|------|-------------------------|
| shocks           | 1–4      | 4    | `sim_shock_series`      |
| states           | 5–13     | 9    | `sim_state_series`      |
| **other vars**   | 14–137   | 124  | `sim_vars_series` (`results_mat`) |
| valuation        | 138–151  | 14   | `sim_value_series`      |

⚠ The "124 other vars" = `n_interp`. The `jx_*` index for an other-var = **`13 + (results_mat column)`**.
e.g. `jx_yh = 14` ↔ `results_mat` col 1 = `y_current(1)`.

### Layout of the 124 `results_mat` columns (`mod_calc.f90:789–794`)

`results_vec` has length `n_interp - 6*12 = 124 - 72 = 52`. The remaining 72 are the bond
term-structure block, written in a loop over maturities `bbb_counter = 1..6` corresponding to
maturities **[1, 2, 3, 4, 19, 20]** (loop guard `bbb<=4 .or. bbb>n_bond-2`, `n_bond=20`,
`mod_calc.f90:718–719, 791–794`):

| results_mat cols | jx range | content                                           |
|------------------|----------|---------------------------------------------------|
| 1–52             | 14–65    | `results_vec` (scalars, see §1)                   |
| 53–70 (18=6×3)   | 66–83    | `q_bond_mat` prices: per mat × {h, f, hw}         |
| 71–106 (36=6×6)  | 84–119   | `E_rb(1:6)` per mat: {Eqh,Eqf,Eqhw, fEqh,fEqf,fEqhw} |
| 107–124 (18=6×3) | 120–137  | `who_prices_bond` per mat × {h, f, hw}            |

Within each maturity the 3-currency order is **(h, f, hw)** = `bond_prices_out(1:3)`
(`calc_bond_prices:3367–3386`): 1 = home private bond, 2 = foreign bond, 3 = home **omega-adjusted**
("hw") bond priced with `1/(1-omg)`.

### De-trending machinery (common pattern)

`z_shocks = temp_series(:,1)`; cumulative TFP trend `z_series = exp(cumsum(z_shocks))`
(`extract_series.m:7,12`). Any "level" quantity stored by the solver in **de-trended units** is
re-trended by `.*z_series` (output, consumption, capital, wages, investment, savings, wealth, bg).
Returns/ratios are scale-free and are NOT multiplied by `z_series`.

Disaster correction: `dis_series` isolates the disaster component of `z_shocks`
(`extract_series.m:15–17`); capital "chosen at t" divides it out via `./exp(dis_series)` and the
realized capital return `rk` multiplies `exp(dis_series(2:end))` (`:275`). `dz_vec` (Fortran) is the
quadrature-node TFP growth incl. the disaster node — already baked into `results_vec` next-period
expectations.

---

## 1. Per-variable spec

For each: **(a)** jx name + index and results_mat col; **(b)** the Fortran formula producing the
column; **(c)** the `extract_series.m` transform into the named `*_series`.

### Output — `yh_series`, `yf_series`
- **jx_yh = 14** (col 1), **jx_yf = 15** (col 2).
- Fortran (`calc_bond_prices:3210, 3449–3450`):
  `y_current = (zf_vec**(1-aalpha)) * kappa_vec**aalpha * l_aggr_vec**(1-aalpha)`,
  `results_vec(1)=y_current(1)`, `results_vec(2)=y_current(2)`. (`zf_vec=[1, zeta*exp(zf)]`.)
- extract (`:57–58`): `yh_series = temp_series(:,jx_yh).*z_series`; same for `yf` ×`z_series`.

### Terms of trade — `s_series`
- **jx_s = 26** (col 13). Fortran `results_vec(13)=s` (`:3461`); `s` = P_f-good / P_h-good rel. price.
- extract (`:54`): `s_series = temp_series(:,jx_s)` (no detrend).

### Aggregate consumption — `c_h_series`, `c_f_series`
- **jx_ch = 18** (col 5), **jx_cf = 19** (col 6).
- Fortran (`:3326–3327, 3453–3454`): `consumption = c_spend / P_div_P_h(iii)`, `c_vec(iii)`;
  i.e. real consumption in each country's **own basket** units. (`c_spend` is nominal spend in
  home-good units; `P_div_P_h` = basket cost in home-good units.)
- extract (`:113–114`): `c_h_series = temp_series(:,jx_ch).*z_series`; `c_f` ×`z_series`.

### Investment — `inv_series`
- **jx_inv = 20** (col 7). Fortran `results_vec(7) = k_nxt(1) - (1-ddelta)*k_aggr` (`:3455`) =
  aggregate investment in capital units.
- extract (`:107`): `inv_series = temp_series(:,jx_inv).*z_series`.

### Labor — `lh_series`, `lf_series`
- **jx_lh = 16** (col 3), **jx_lf = 17** (col 4). Fortran `results_vec(3)=l_aggr_vec(1)`,
  `(4)=l_aggr_vec(2)` (`:3451–3452`).
- extract (`:50–51`): `lh_series = temp_series(:,jx_lh)` (no detrend; labor is stationary).

### Inflation — `infl_h_series`, `infl_f_series`
- **jx_infl_h = 38** (col 25), **jx_infl_f = 39** (col 26).
- Fortran (`:3473–3474`): `results_vec(25)=infl_vec(1)`, `(26)=infl_vec(2)` = realized gross CPI
  inflation (country basket) this period. ⚠ jx 21/22 (`jx_pi`,`jx_pif`) are a *different* object —
  `jx_pi` = `pi_series` (per-unit capital profit, `results_vec(8)=pi_current(1)/P_div_P_h(1)`); do
  NOT confuse `jx_pi` with inflation. Inflation is jx_infl_h/f.
- extract (`:89–92`): `infl_h_series = temp_series(:,jx_infl_h)` (gross, no detrend);
  `P_h_series = cumprod(infl_h_series)` etc.

### Nominal policy rates — `nom_ih_series`, `nom_if_series`
- **jx_nom_i_h = 40** (col 27), **jx_nom_i_f = 41** (col 28).
- Fortran (`:3475–3476`): `results_vec(27)=nom_i_vec(1)` (gross home nominal rate),
  `(28)=nom_i_vec(1)+nom_i_vec(2)` (gross foreign nominal rate; `nom_i_vec(2)` is the spread).
- extract (`:143–144`): `nom_ih_series = temp_series(:,jx_nom_i_h)`.

### Price of capital — `q_series`
- **jx_q = 30** (col 17). Fortran `results_vec(17)=q_current/P_div_P_h(1)` (`:3465`) = real price of
  capital in home basket. extract (`:85`): `q_series = temp_series(:,jx_q)`.

### Real exchange rate — `qx_series`  ⚠TRICKY
- **jx_qx = 25** (col 12). Fortran `results_vec(12) = P_div_P_h(1) / P_div_P_h(2)` (`:3460`) =
  (home basket cost in home-good units) / (foreign basket cost in home-good units) = **real exchange
  rate** (home/foreign price-level ratio in common units). With
  `P_div_P_h(1)=(varsigma_1 + (1-varsigma_1)*s^(sigma-1))^(1/(1-sigma))`,
  `P_div_P_h(2)=s^{-1}(varsigma_2 s^{1-sigma}+(1-varsigma_2))^{1/(1-sigma)}` (`:3215–3217`).
- extract (`:93`): `qx_series = temp_series(:,jx_qx)` (no detrend).
- Derived: `qx_change_series = -log(qx_{t-1}/qx_t)` (`:96`, with first entry duplicated);
  nominal `E_change_series = -log(qx_{t-1}/qx_t · infl_h_t/infl_f_t)` (`:94–95`).

### Profit per unit capital — `pi_series` (feeds rk, div)
- **jx_pi = 21** (col 8). Fortran `results_vec(8)=pi_current(1)/P_div_P_h(1)` (`:3456`),
  `pi_current = s_vec^{-1}·aalpha·(zf_vec·l/kappa)^{1-aalpha}` (`:3205`).
- extract (`:145`): `pi_series = temp_series(:,jx_pi)` (no detrend; per-unit, ratio-like).

### Goods price ratios — `P_Phh_series`, `P_Phf_series`
- **jx_P_P_hh = 23** (col 10) = `P_div_P_h(1)`; **jx_P_P_hf = 24** (col 11) = `P_div_P_h(2)`
  (`:3458–3459`). extract (`:46–47`). Used to convert returns between baskets.

### Realized risk-free returns

**`rfh_series`** (home real safe rate, realized log):
- Built from `nom_ih_series` and `infl_h_series`, NOT a single column. extract (`:263–264`):
  `rfh_series = log( nom_ih_{t-1} / infl_h_t )`, then prepend first value. **Realized-from-consecutive-
  periods** (uses t and t-1). This is the *return on the home private/policy bond*, no omega.

**`rff_series`/`rff_h_series`** (foreign safe rate, in foreign then home basket):
- extract (`:271–274`): `rff_series = log( nom_if_{t-1}/infl_f_t )`;
  `rff_h_series = rff_{t} + log( P_Phf_t/P_Phf_{t-1} / P_Phh_t · P_Phh_{t-1} )` — converts the
  foreign-bond return into home-basket units via the realized real-exchange-rate change. Prepend first.

**`rfh_omg_series`** (omega-adjusted home rate, used as the wedge-free benchmark in valuations):
- extract (`:267–268`): `log( nom_ih_{t-1}/infl_h_t / (1 - omg_{t-1}) )`. Note `omg_series`
  (`:43,139`): `omg = exp(temp(:,jx_omg)) + omg_shift - b_lmbd*bg_xtra_rel`. **jx_omg = 13** (col... it
  is a *state*, col index 13, not an other-var).

### Realized capital / equity return — `rk_series`, `rA_series`  ⚠
- **`rk_series`** (`:275–276`): `rk = log( exp(dis_t)·(pi_t + (1-delta)·q_t)/q_{t-1} )`. Uses
  `pi_series`(jx_pi), `q_series`(jx_q), and the disaster correction. Realized, consecutive periods.
- **`rA_series`** = levered equity return (`:285`):
  `rA = log( (1+debt_to_equity)·exp(rk) - debt_to_equity·exp(rfh_lev) )`, where
  `rfh_lev_series` (`:280`) blends the 20-period (long) bond realized returns:
  `rfh_lev = log( 1/(1+zeta)·exp(rq20_h) + zeta/(1+zeta)·exp(rq20_fh) )`. `debt_to_equity`, `zeta` ∈ params.
  - `rq20_h_series` (`:183`): `log( q19_h_{t}/infl_h_{t}/q20_h_{t-1} )` (realized holding return on the
    20-per bond = buy 20-per at t-1, it becomes a 19-per at t). `rq20_fh` = foreign-20 in home basket
    (`:198`). These need `q19_h, q20_h, q19_f, q20_f` columns (jx 78/81, 79/82 etc.) — see bond block.

### Excess returns
- **`exc_retA_series`** (`:292`) `= rA_series - rfh_series` (levered equity excess return).
- **`exc_rf_series`** (`:300`) `= rff_h_series - rfh_series` (foreign-minus-home safe excess, home basket).
- (also `exc_ret_series = rk - rfh`, etc., not used by Table 2.)

### Bond yields — `yield1_h_series`, `yield1_f_series`, `yield1_hw_series`
- From the bond-price columns. **jx_q1h = 66** (col 53), **jx_q1f = 67** (col 54),
  **jx_q1hw = 68** (col 55). Fortran `bond_prices_out(1:3)` = (h, f, hw) 1-period prices
  (`calc_bond_prices:3367–3377`), written into the maturity-1 slot.
- extract (`:394–405`): `yield1_h_series = -log(q1_h_series)`; `yield1_f = -log(q1_f)`;
  `yield1_hw = -log(q1_hw)`. (n-period: `-log(qn)/n`.)
- ⚠ `q1_h` is the **private** home bond (no omega); `q1_hw` is the omega-adjusted ("convenience")
  home bond. Table-2 carry regression uses `yield1_hw` for UIP but `yield1_h` for the Fama spread.

### UIP residual — `uip_pvt_series`  ⚠TRICKY
- extract (`:469–470`): `uip_pvt = 100·( yield1_f_{t} - yield1_hw_{t} - E_change_{t+1} )`, then prepend
  0. So: **foreign 1-per yield − home 1-per (omega-adjusted "hw") yield − expected log nominal
  depreciation next period**. `E_change_series` (`:94`) = `-log(qx_{t-1}/qx_t·infl_h_t/infl_f_t)`,
  here shifted forward by one (`E_change_series(2:end)` ↔ time t's *next-period* realized change used
  as the UIP forward term). Uses **next-period** info (it is the realized depreciation used in a
  forward-looking UIP deviation). Annualized ×100 (already in pct, gross-of-frequency; the regressions
  multiply other terms by 4).

### Output growth (4-quarter) — `y_growth_series`
- extract (`:467–468`): `y_growth = 100·( log(yh_{t}) - log(yh_{t-4}) )`, zero-pad first 4.
  4-quarter log output growth, in pct. Uses `yh_series` (already z-trended).

### Dividend-price (smoothed) — `div_price_smoothed_series_1`  ⚠TRICKY
- Levered dividend on the equity claim. extract:
  - `div_series_1` (`:359–361`): per-unit-capital levered cash flow,
    `div_1_t = pi_t + (1-delta)q_t - debt_to_equity/(1+debt_to_equity)·q_{t-1}·exp(rfh_lev_t)
               - 1/(1+debt_to_equity)·q_t`. (Uses t and t-1 ⇒ realized, consecutive periods;
    prepend first value `:366`.)
  - `qE_series_1` (`:370`): `= q_t/(1+debt_to_equity)` (the equity slice of capital price).
  - `div_price_series_1 = div_1/qE_1` (`:385`).
  - **smoothed** (`:389–390`): `div_price_smoothed_series_1 = movsum(div_series_1,[3,0]) ./ qE_series_1`
    — trailing 4-quarter sum of dividends over current equity price (annualized D/P). First 3 entries
    set to `div_price_smoothed_series_1(1)*4` (`:390`). This is why the regression starts at t=4/5.

### NFA — `nfa_series`, `nfa_rel_series`, `nfa_rel_growth_series`  ⚠TRICKY
- **`nfa_series`** (`:314–315`): net foreign assets =
  `nfa_t = h_sav_{t} - h_kap_{t+1}·exp(-dis_{t+1})·q_{t}`,
  i.e. home savings at t minus the value of domestically *deployed* capital next period (disaster-
  adjusted, valued at this period's q). **Uses next-period** `h_kap` and `dis`. Then duplicate last.
  - `h_sav_series` = **jx_h_sav = 57** (col 44), `results_vec(44)=savings_vec(1)/P_div_P_h(1)`
    (`:3492`), ×`z_series` (`:148`).
  - `h_kap_series` = **jx_h_kap = 43** (col 30), `results_vec(30)=kappa_vec(1)` (`:3478`), ×`z_series`
    (`:156`). (`h_ksav` = jx_h_ksav = 42, col 29 = `savings·(1-share)/P` `:3477`.)
- **`nfa_rel_series`** (`:319): `nfa_series ./ yh_series`.
- **`nfa_rel_growth_series`** (`:320–321`): `(nfa_t - nfa_{t-1})/yh_t`, prepend first.
- ⚠ The Fortran `calc_valuation` also produces an *expectational* nfa (`valuation_vec`, jx 138–151)
  used only in **Table 7** (the `dis` branch); the Table-2 `nfa_*` use the **realized** construction
  above, not the valuation columns.

### Net exports (for completeness; Table 2 mom 11 uses c_f/qx not nx)
- `nx_series` (`:409–414`) = exports − imports from `chf/cfh/inv_h/h_kap`; `nx_rely = nx/yh`.

### Expected returns (NOT used by Table 2; listed to avoid confusion)
- `E_rfh_series=jx_Erfh=27` (col 14), `E_rff=jx_Erff=28` (col 15), `E_rk=jx_Erk=29` (col 16),
  `results_vec(14..16)=E_rfh,E_rff,E_rk` (`:3462–3464`). Table 2 uses **realized** rfh/rff/rk/rA, not these.

---

## 2. The 15 Table-2 moment formulas (`calc_moments.m:33–47`)

Let `series(a:b)` denote the MATLAB slice. `trg_prms` are calibrated params.

1. `mean( yf ./ (s · yh) )` — relative output level (foreign/home), tot-adjusted. (`:33`)
2. `100·std( diff(log(c_h)) )` — home consumption growth volatility, pct. (`:34`)
3. `100·std( diff(log(yf)) )` — foreign output growth volatility, pct. (`:35`)
4. `corr( (yf/yh)_{t}, (yf/yh)_{t+1} )` — autocorr of relative output. (`:36`)
5. `100·std( diff(log(inv)) )` — investment growth volatility, pct. (`:37`)
6. `4·100·mean( rfh )` — annualized mean home real safe rate, pct. (`:38`)
7. `100·mean( nfa_rel / 4 )` — mean NFA / (annual) output, pct (÷4 to annualize yh). (`:39`)
8. `corr( resid_dp, resid_carry )` — corr of the two regression residuals (see §3 R1,R2). (`:40`)
9. `4·100·mean( rA_{2:end-1} ) - 4·100·mean( rfh_{2:end-1} )` ⚠ — written as
   `4*100*mean( rA(2:end-1) - mean(rfh(2:end-1)) )`; algebraically the annualized levered equity
   premium (mean rA − mean rfh), ×4×100. (`:41`)
10. `reg_coeffs1(2)` — slope on excess equity return in the NFA-growth regression (see §3 R3). (`:42`)
11. `100·trg_prms.bg_yss · mean( (c_f/qx) ./ (4·yh) )` — government-debt / annual-output style
    scaling moment; `c_f/qx` converts foreign consumption to home-basket units. (`:43`)
12. `mean( lh )` — mean home labor. (`:44`)
13. `mean( lf )` — mean foreign labor. (`:45`)
14. `100·mean( log(infl_h) )` — mean home (net, log) inflation, pct. (`:46`)
15. `100·mean( log(infl_f) )` — mean foreign log inflation, pct. (`:47`)

(`sss` indexes the simulation/parameter draw; results averaged across `n_sims` in `collect_moments.m`.)

## 3. The three regressions (`calc_moments.m:15–31`)

**R1 — excess equity returns on (lagged, smoothed) dividend-price** (`:15–18`):
- `y = 4·100·exc_retA_series(5:end)`  (annualized levered equity excess return, starts t=5 b/c of
  dividend smoothing)
- `x = div_price_smoothed_series_1(4:end-1)`  (lagged smoothed D/P)
- `regress(y, [1, x])` → keep residual `resid_dp`.

**R2 — Fama / carry: excess foreign bond return on yield diff and past output growth** (`:20–25`):
- `fama_yield_pvt_series = 100·4·( yield1_f_series - yield1_h_series )`  (annualized foreign−home
  1-period **private** yield spread)
- `x = uip_pvt_series(5:end)`  (the UIP residual, dep. var here)
- `y = fama_yield_pvt_series(4:end-1)`  (lagged yield diff)
- `y2 = y_growth_series(4:end-1)`  (lagged 4-q output growth)
- `regress(x, [1, y, y2])` → keep residual `resid_carry`.
  (Moment 8 = `corr(resid_dp, resid_carry)`.)

**R3 — NFA growth on excess equity and excess foreign-bond returns** (`:27–31`):
- `x  = exc_retA_series(3:end-1)`     (levered equity excess return)
- `x2 = exc_rf_series(3:end-1)`       (foreign−home safe excess return, home basket)
- `y  = nfa_rel_growth_series(3:end-1)`
- `reg_coeffs1 = regress(y, [1, x, x2])`.
  - Moment 10 = `reg_coeffs1(2)` (loading on excess equity).
  - (`reg_coeffs1(3)`, loading on excess fbond, is Table 3 mom 3.)

---

## Appendix — which series need next-period / expectational vs realized-consecutive

- **Realized from consecutive periods (t, t-1 or t+1):** `rfh, rff/rff_h, rfh_omg, rk, rA, rq20_*`,
  all `exc_*`, `div_series_1`, `div_price_smoothed`, `qx_change`, `E_change`, `nfa` (uses t+1 `h_kap`),
  `nfa_rel_growth`, `y_growth` (t, t-4), `c/y diff` moments.
- **Forward-looking / next-period info:** `uip_pvt` (uses `E_change_{t+1}`), `nfa_series` (uses
  `h_kap_{t+1}`, `dis_{t+1}`).
- **Solver-internal expectations (already inside a single column, NOT recomputed in MATLAB):**
  `E_rfh/E_rff/E_rk` (jx 27–29) and the `E_rb`/valuation columns — these integrate over `big_weight_vec`
  quadrature nodes inside Fortran. Table 2 does not use these (uses realized analogues).
- **De-trending location:** entirely in `extract_series.m` via `.*z_series` on level vars and the
  disaster `exp(dis_series)` corrections; returns/ratios are left untouched.
