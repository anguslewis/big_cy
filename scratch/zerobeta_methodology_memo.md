# Documentation memo — How DiTella, Hébert, Kurlat & Wang (2025, JPE) estimate the zero-beta interest rate

**Author:** `big_cy-zerobeta` (Assignment 05 / E0) · **Date:** 2026-06-13
**Status:** Understanding & documentation only. No replication, no plan, no estimation code. (Per the scope boundary in the dispatch.)

**Sources read first-hand**
- Paper: `…/Dropbox/Angus Lewis/papers/DiTella_Hebert_Kurlat_Wang_2025_JPE_The_zero_beta_interest_rate.pdf` (main text, Secs. 1–7).
- Appendix: `…_Online_appendix.pdf` (Secs. A–M, 90 pp.).
- Public replication code, cloned to `scratch/TheZeroBetaRate/` from `github.com/bhebert/TheZeroBetaRate` (commit fetched 2026-06-13). Pipeline = Python `DataCode/` (data assembly) → MATLAB `MatlabCode/` (GMM + tests) → Stata `StataCode/MakeTables.do` (tables).

**Reference shorthand.** "MT §" = main-text section; "App §" = appendix section; code files are cited as `File.m:function` / `file.py`.

> **One-paragraph summary.** The zero-beta rate `R_{0,t}` is the conditional expected return of a stock portfolio whose return is orthogonal to the SDF — i.e., a unit-investment portfolio with zero loading on every priced factor. DHKW estimate a *time series* of it by GMM over **~130 pre-formed, value-weighted equity test portfolios** (not individual stocks): the modelled rate `R_{0,t}=γ'Z_t` is identified as the best linear predictor (in instruments `Z_t`) of the realised return on the **minimum-variance** zero-beta portfolio, where the covariance matrix used to form that portfolio is a **Ledoit-Wolf (2020) nonlinear-shrinkage** estimate. Estimation is **exactly identified** via a reduced-rank (rank-1, PSD) weighting matrix that collapses the N asset-pricing moments onto the single zero-beta portfolio combination; the **consumption Euler equation** supplies the over-identifying moments used only for a **Stock-Wright (2000) weak-identification-robust S-test** of the curvature parameter σ. The headline: the real zero-beta rate averages ~8.5%/yr, is volatile (SD ~9%/yr), runs ~7.6%/yr above the T-bill, and *falls* in bad times (co-moving with expected consumption growth) rather than spiking.

---

## 1. The estimand and the GMM procedure

### 1.1 The estimand
The **zero-beta interest rate** `R_{0,t}` is "the conditional expected return of a stock portfolio orthogonal to the stochastic discount factor" (MT abstract; MT §1). With nominal SDF `Λ_t = δ^t c_t^{-σ} P_t^{-1} ξ_t`, the object tested is implication (iii): `R_{0,t}^{-1} = E_t[ δ (c_{t+1}/c_t)^{-σ} (P_t/P_{t+1}) ]` (MT §1, eqs. 1–3; App §A). The contribution over Black-Jensen-Scholes (1972)/Shanken (1986) is producing a **time series**, not just an average level.

### 1.2 The three model layers (MT §2.1)
- **Factor model for test assets** (MT eq. 4): `R_{i,t+1} − R_{0,t} = α_i + Σ_j β_{ij} F_{j,t+1} + ε_{i,t+1}`, betas constant, excess returns defined **relative to `R_{0,t}`** (not the T-bill). The first factor is the market *excess of the zero-beta rate*, `F_{1,t+1}=R_{m,t+1}−R_{0,t}`; the rest are zero-investment.
- **Linear SDF** (MT eqs. 5–6) with time-varying `R_{0,t}` and prices of risk.
- **Predictive specification** (MT eq. 7): `R_{0,t}(γ) = γ'Z_t`, linear in `L+1` instruments with `Z_{0,t}=1` and the T-bill yield `R_{b,t}` included so the traditional view `R_{0,t}=R_{b,t}` is nested.

### 1.3 The infeasible three-step logic (MT §2.2), then jointly by GMM
(1) time-series regressions → (α, β); (2) pick a zero-beta portfolio `w` with `w'β=0`, `w'ι=1`; (3) OLS-predict the portfolio return with `Z_t` → γ. This is infeasible because step 1 needs `R_{0,t}` to form excess returns, so all parameters are estimated **jointly**.

### 1.4 The stacked moment system (MT §2.3; App §M p.85)
Parameter vector **(θ, δ, σ)** with **θ = (α, β, γ)**. Three moment blocks:
1. **Time-series / beta moments**: `E[ε̂_{t+1}(θ) ⊗ F_{t+1}] = 0` — N×(K+1) moments identifying (α, β).
2. **Instrumented zero-beta pricing moments**: with the symmetric orthogonal projector `H(β) = I_N − β(β'β)^+ β'`, `E[ H(β)·(R_{t+1} − ι·γ'Z_t) ⊗ Z_t ] = 0` — this is the predictive (γ) regression; `H` guarantees the priced combination is genuinely zero-beta.
3. **Consumption Euler moments**: `E[ (δ (c_{t+1}/c_t)^{-σ} (P_t/P_{t+1}) R_{0,t} − 1) ⊗ Z_t ] = 0`.

Blocks 1–2 estimate (α, β, γ, δ); **block 3 is held back as over-identifying** and used only to build the S-test (App §M; MT §4.2).

**In the code**, the concentrated objective (`InstMomentsConc.m`, `InstrumentGMMConcOpt.m`) targets exactly `K+2` moments — the `K+1` block-2 moments (`amoments = (w·H·(R−zbrate)).*[1;Z]`) plus the **single unconditional** block-3 Euler moment (`consmoments` row for the constant). The remaining `K` conditional Euler moments are the over-identifying set fed to the S-test (`InstrumentGMMWrapperConc.m:166–172`). Note the code's parameters are `[Rf; γ(K); ρ; σ]` where `ρ` ≈ the paper's `δ` (level/time-preference) and code-`γ` predicts the *spread* of the zero-beta rate over the bill yield (the paper's `γ` predicts the nominal level; `RunTest.m:295–305` patches the saved output to the paper convention).

### 1.5 Exact identification via a reduced-rank, PSD weighting matrix (MT §2.3; App §M p.85)
The GMM weight matrix is **block-diagonal and only positive semi-definite**:
```
W_T(θ) = blkdiag(  I_{N(K+1)} ,  w_T(θ)w_T(θ)' ⊗ I_{L+1} ,  e_0 e_0'  )
```
The middle block is **rank-1** (`w_T w_T'`): it projects the N asset-pricing moments onto the *single* zero-beta-portfolio combination before interacting with the `L+1` instruments. This "selects moments and achieves exact identification, as opposed to using a full-rank weighting matrix with over-identifying restrictions" (MT §2.3, citing Cochrane 2009); it is "akin to a GLS weight matrix" and is motivated by a Gaussian-MLE argument (App §L, relationship to Shanken 1986). The PSD (vs. PD) relaxation is the formal generalization of Stock-Wright they must make (App §M pp.84, 87–88). `W` depends on `T` because `w_T` embeds the (Ledoit-Wolf) covariance estimate.

### 1.6 The minimum-variance, unit-investment, zero-beta portfolio (MT §2.3 eq. 10; `PortfolioWeight.m`)
Among all `w` with `w'ι=1` and `w'β=0`, pick the one minimizing `w'Σ_R w`. Closed form:
```
w(θ) = Σ_R^{-1} A (A' Σ_R^{-1} A)^{-1} e_1,   A = [ι, β],   e_1 = [1; 0_M]
```
`PortfolioWeight.m` implements exactly this (`Sinv=inv(Sigma); mult=inv(Amat'*Sinv*Amat)*e1; weight=(Sinv*Amat*mult)'`). The weights are over the **N test portfolios**, and `Σ_R` is the Ledoit-Wolf estimate (§1.7). Minimizing variance maximizes the precision of the γ predictive regression. (`InstMomentsConc.m` additionally pre-multiplies by `H=I−β·pinv(β)`, but `w·H=w` already, so this is a numerical safeguard.)

### 1.7 The Ledoit-Wolf covariance estimator (App §K; `BetaSigma.m`, `ExternalCode/analytical_shrinkage.m`)
`Σ_R(θ)` is built (matching App §K and `BetaSigma.m:50–75`):
1. Factor-model preconditioner `Σ_F = β Σ_K β' + diag(Σ̂ − β Σ_K β')` (K-factor model implied covariance; `Σ_K` = factor covariance). This is their adaptation of LW (2017) §4.2, using a **K-factor** (not single-factor) preconditioner.
2. Symmetric square root `Σ_F^{1/2}` (`sqrtm`), whiten the de-meaned excess returns `e_P = Σ_F^{-1/2}(R − ι R_{0,t} − μ)`.
3. Apply the **Ledoit-Wolf (2020) analytical nonlinear shrinkage** estimator to `Cov(e_P)` (`analytical_shrinkage.m`, downloaded from Wolf's site; LW substitute the analytical 2020 estimator for the 2017 nonlinear one — "similar OOS performance, substantially faster").
4. Re-inflate: `Σ_R = Σ_F^{1/2} Σ̂_c Σ_F^{1/2}`.

The shrinkage *target* is implicit in the preconditioner (the factor-model structure); the appendix defers the shrinkage-intensity formula to the LW code rather than restating it. Alternatives (robustness, `RunAll.m`/App §J): `Sample` (raw sample cov — "not recommended"), `Diag`/`PCA`, and `I` (identity → the zero-beta portfolio closest to equal-weighted).

---

## 2. THE CORE QUESTION — the portfolio choice set

**Confirmed from code: the algorithm chooses portfolio shares over PRE-FORMED PORTFOLIOS, never individual securities.** The GMM never sees a single stock's return. Individual stocks enter only upstream, to *build* the test portfolios and their pre-sorting market beta.

### 2.1 What the weights range over
`w(θ) ∈ R^N` is a vector of shares over the **N test portfolios** (`RunTest.m:58`: `R = table2array(Portfolios(:,2:end))'`). In the **main spec, N = 130**. The constraints are `w'ι=1` (unit investment) and `w'β=0` (zero loading on all 7 factors); among those, min-variance under `Σ_R`. So the choice set is the affine subspace {unit-investment, factor-neutral combinations of the 130 portfolios}.

### 2.2 The exact test-asset universe (MT §3 "Equity Portfolios"; App §C.1 pp.14–17; `portfolio_construction.py`)
Built à la Fama-French (1993), monthly, **1973–2020**, CRSP×COMPUSTAT matched (closely follows Drechsler's WRDS Python replication, App fn. 54):

- **81 beta-sorted portfolios** = three independent 3×3×3 sorts (27 each), value-weighted:
  - market β × size (ME) × book-to-market,
  - market β × size × operating profitability,
  - market β × size × investment.
- **+ 49 industry portfolios** (4-digit SIC, value-weighted) from Ken French.
- **Total N = 130.** The FF3 robustness spec (`27_plus_Industry`, earlier-version main) = 27 (β×size×B/M) + 49 = **76** (App §J p.51).

**Construction details (verified in `portfolio_construction.py:generate_groups`):**
- **Breakpoints at 30th/70th percentiles**, computed on **NYSE stocks only** (`exchcd==1`) restricted to **US common shares** (`shrcd ∈ {10,11}`) with `bm>0, me>0`. (Code segment `[0, 0.3, 0.7, 1]`.)
- **Beta breakpoints are computed *within* each size×characteristic group**, not unconditionally. The code comment is explicit about *why*: *"the small stocks may have inaccurate betas due to liquidity reasons"* (`portfolio_construction.py:107–109`). This is a deliberate small-stock safeguard (see §3).
- **Value-weighting**: returns are `wavg(retadj, wt)` where `wt` is end-of-May market value compounded by the ex-dividend return ("to eliminate the effect of equity issuance on weights", App fn. 60).
- Annual June rebalancing (`ccm_june`), Fama-French timing.

### 2.3 The bottom-20% exclusion (MT §3; App §C.1 p.16; `portfolio_construction.py:59–63`)
*"At the end of each June, the bottom 20% of stocks, in terms of market value, are dropped."* In code: per-date 20th-percentile ME cutoff, keep `me >= cutoff`. Applies to the **beta-sorted** portfolios. The 49 industry portfolios are *not* filtered (they include the smallest stocks) but, being value-weighted and not beta-sorted, small stocks carry negligible weight there (MT fn. 13). The **NoDrop20** robustness spec re-runs everything with this filter off (App §J p.50, Table 7, Fig. 15).

### 2.4 The pre-sorting market beta (`A2_beta_estimation.py`)
Single-factor CAPM beta of each stock's monthly excess return on the contemporaneous CRSP market excess return, **5-year rolling window, minimum 24 monthly obs**. This β is used **only to assign stocks to beta terciles**, not as the GMM β. (Implementation notes to flag: the docstring says "one to three years (two minimum)" but the code uses a 60-month window with min 24; the return-censoring to (−50%,100%) mentioned in the docstring is **commented out**; a `beta_segment` mismatch exists between `00_create_data.py` `[0,0.3,0.7,1]` and `A3_summary_statistics.py` `[0,0.33,0.66,1]` — cosmetic for summary stats only.)

### 2.5 Data defining the universe
CRSP monthly stock file (`crsp.msf`), CRSP events (`crsp.mseall`: shrcd/exchcd/SIC), CRSP delisting (`crsp.msedelist`), COMPUSTAT fundamentals (`comp.funda`: BE, OP, INV inputs), CCM link (`crsp.ccmxpf_lnkhist`), CRSP market index (`crsp.msi`), Ken French 49-industry returns. Full provenance in §7.

---

## 3. Safeguards against small/illiquid idiosyncratic stocks — and an assessment

**The concern (restated).** A bundle of small, illiquid stocks whose *idiosyncratic* variance dominates their factor loadings can look like it has ~zero β not because loadings are genuinely zero but because illiquidity/stale pricing attenuates measured betas. Bundling many such names diversifies the idiosyncratic noise and produces a portfolio that *appears* zero-beta and low-variance — exactly what a min-variance, factor-neutral optimizer is drawn to.

**Safeguards actually present (paper + code):**
1. **Pre-formed portfolios, not individual stocks** (§2.1). The optimizer cannot cherry-pick individual micro-caps; it can only weight 130 diversified portfolios. Strong mitigant against the single-stock version of the concern.
2. **Bottom-20%-by-ME exclusion** from the beta-sorted set (§2.3). Motivated (MT fn. 12) by **Novy-Marx & Velikov (2022)**: the smallest deciles are illiquid and have betas *attenuated* relative to other stocks (this drives Frazzini-Pedersen 2014 results). Directly targets the attenuated-beta channel.
3. **Beta sorts formed within size groups** (§2.2) — explicitly because "small stocks may have inaccurate betas due to liquidity reasons" (code comment). Prevents liquidity-driven beta error from contaminating the cross-size beta ranking.
4. **Monthly (not daily) betas** (§2.4). The paper attributes the robustness of NoDrop20 to using monthly betas, which are far less liquidity-contaminated than daily betas.
5. **Value-weighting** within every portfolio (§2.2) — small names get small weight even inside their bucket.
6. **Minimum-variance + Ledoit-Wolf shrinkage** (§1.6–1.7). The paper's stated motivation (MT §3 p.14) is that forcing a zero-beta portfolio out of only FF25 size×value portfolios "would load heavily on poorly estimated residuals"; using many portfolios + LW shrinkage avoids over-fitting the residual covariance. MT fn. 14 cautions that if test assets are highly correlated conditional on factors, LW performs poorly, so the main spec includes enough factors (7) to explain co-movement.
7. **NoDrop20 robustness** (App §J Table 7, Fig. 15): re-running with the bottom two deciles *included* leaves the level (~8.3–9.2%), SD, S-test, and consumption tracking "very similar." This is the most direct empirical rebuttal of the small-stock-artifact story.

**Assessment of adequacy.**
- Safeguards 1–5 are well-aimed at the **individual-stock** and **measured-beta-attenuation** versions of the concern, and 7 provides genuine empirical reassurance that *excluding vs including* small stocks does not move the result.
- The **residual exposure** is at the **portfolio level**: even after value-weighting and the 20% cut, the beta-sorted portfolios in the *low-beta × small* corner are still tilted toward less-liquid names. If those portfolios' betas are attenuated by stale pricing, the min-variance optimizer — whose objective is precisely to find the lowest-variance factor-neutral combination — will happily overweight them, because diversifying across portfolios reduces the idiosyncratic variance that LW shrinkage further smooths. The min-var + LW step thus **mitigates over-fitting of the covariance** but does **not** correct *systematically biased (attenuated) portfolio betas*; nothing in the pipeline re-estimates portfolio betas with a Dimson/Scholes-Williams lead-lag correction. So the portfolio-level illiquidity channel is the part the safeguards address least directly.
- The covariance-weighting *does* matter: the **CovI** robustness (identity Σ → portfolio closest to equal-weighted, i.e. *least* concentrated in low-variance/illiquid corners) is the **one spec where the S-test behaves materially differently** (App §J Fig. 34), while CovSample and CovPCA preserve the headline. This is a useful internal diagnostic: it shows the result is partly a property of the min-variance/Ledoit-Wolf weighting, not just of the test-asset set.

---

## 4. Evaluating the hypothesis (Angus's prior)

**Prior:** the zero-beta rate's "very large movements in bad times" suggest the zero-beta portfolio loads on small, weird, idiosyncratic stocks.

**What the paper/code actually show — and a key factual correction.** The estimated real zero-beta rate **falls in bad times, it does not spike** (MT §3.1 p.19; MT §4.3). The predictive signs (MT Table 1): it rises >1-for-1 with the T-bill yield, falls with lagged inflation; high unemployment, inverted yield curve (low/negative term spread), and high excess bond premium predict a *low* zero-beta portfolio return — "the zero-beta rate falls when a recession is likely," co-moving with expected consumption growth. So the large bad-times *movements* are **downward** and are **explained by recession-predictor instruments**, not by a flight-to-safety idiosyncratic blow-up. A monetary-tightening shock pushes the zero-beta rate and consumption growth down together (MT §4.3, Fig. 4), while the real T-bill return rises.

**Does the methodology *permit* a small/weird-stock artifact?** Partly yes, as a *level*-and-*volatility* matter, not as the source of the cyclical pattern:
- The min-variance + factor-neutral construction **permits** concentration in the lowest-variance corner of the (post-filter) portfolio space, which can include less-liquid value-weighted portfolios (§3). The variation in `R_{0,t}` over time is *predictable* variation driven by `Z_t`, not raw realized-return noise — so it is not mechanically idiosyncratic. But the *amplitude* of that predictable variation depends on the portfolio's exposure, and the CovI diagnostic (§3) shows the weighting scheme affects it.
- The authors themselves flag the **genuinely open** piece: the zero-beta rate is occasionally estimated **below the T-bill yield** and its excess volatility is large (SD ~9%/yr; std of `E_t[pd]` ≈ 29% of the consumption-claim p/d ratio). They acknowledge this could reflect **linear-extrapolation / in-sample over-fitting** of the predictive regression or genuinely low expected returns on riskier/less-liquid assets (MT §7 p.47). Consistent with over-fitting: the **Ridge** spec (cross-validated penalty shrinking `R_{0,t}` toward the T-bill) cuts the excess-volatility share from 29% to 19% (MT §6.1) — i.e., a non-trivial slice of the in-sample variation is fragile.

**Bottom line for the prior.** The bad-times *direction* (down, recession-predicted) is not an idiosyncratic artifact and is well-identified. The residual concerns that *do* survive are: (i) the **portfolio-level illiquidity-attenuated-beta** channel that the safeguards address least (§3), and (ii) **in-sample over-fitting / extrapolation** of the predictive regression, which the authors concede and which Ridge partially deflates. Neither is *resolved* here — both are documented as legitimate open questions a euro-area replication (or a re-estimation with lead-lag-corrected portfolio betas and the CovI/Ridge variants) could probe. (Per scope: documented, not resolved empirically.)

---

## 5. Full implementation inventory (for a future replication)

### 5.1 Optimization / convergence (`InstrumentGMMWrapperConc.m`, App §D)
- **Solver:** MATLAB `fmincon` (interior-point default), `Display off`, `MaxIterations 1e5`, `MaxFunctionEvaluations 1e10`, `OptimalityTolerance 1e-10`, `StepTolerance 1e-14`, `FunctionTolerance 1e-10`, with `TypicalX` scaling.
- **Concentration:** the optimizer searches **only over γ** (code: `[Rf; γ; ρ]`); (α, β), `Σ_R`, and the weight `w(θ)` are computed analytically at each γ (`BetaSigma.m` → `PortfolioWeight.m`). σ is **fixed** in any one solve and swept on a grid for the S-test. App §D: "it is possible to analytically compute the (α,β) parameters and weighting matrix w(θ) given an estimate of γ, which greatly speeds up the computation."
- **Starting values:** `R_{0,t}=R_{b,t}` (γ_1=1, others 0); refined by two infeasible-OLS iterations (regress the initial zero-beta portfolio return on `Z_t`).
- **Convergence handling:** objective should hit ~0 (exact identification); if `obj > 1e-6`, fall back to `GlobalSearch` with **fixed seed `rng(12302018,'twister')`**, re-seeding from infeasible-OLS up to a retry cap; hard error if still not converged or `fmincon` flag < 1. Wide box bounds avoid singular matrices.
- **Runtime:** ~2 h for `RunAll.m` on an M1 (per README); the bootstrap is the expensive piece (§5.2).

### 5.2 Weighting matrix, standard errors, weak-ID test, bootstrap
- **Moment covariance `Ω`** (`Wopt.m`): main spec `har='COV'` → **plain sample covariance** of the moments (residuals assumed serially uncorrelated under the conditional moment restriction; Hansen-Singleton 1982 / Cochrane 2009 ch. 11.7). Alternative: **Newey-West** with `nma = ceil(0.75·T^{1/3})` lags (Lazarus-Lewis-Stock-Watson 2018 textbook tuning) — available but off by default.
- **Standard errors** (`InstGMMSE.m`; App §M pp.88–90): standard GMM sandwich with the **exact-identification weight matrix** `Wmat` (identity on the `Nm−K` estimation moments, zero on the `K` over-identifying Euler moments). `pvar = (1/T)·GWGiGW·Ω·GWGiGW'`, `GWGiGW = (J'WJ)^{-1}J'W`. **Estimation error in betas (and alphas) is propagated**: the Jacobian `J` is taken w.r.t. the *full* parameter vector `[θ; β(:); α]`, and the moment-side variance uses the residual-maker `R̃ = I − J·GWGiGW` so the test moments are purged of estimation error in (α, β, γ). The **Jacobian is mixed**: the projection (β) block analytic (`InstMomentsProjectionJacobian.m`), the asset-pricing/Euler block numeric (`jacobianest.m`, ExternalCode) — see `InstMomentsFullJacobian.m`. SEs that "account for estimation error in the other parameters" are reported throughout (MT Table 1; App Tables 7–10). OLS-comparison SEs (which ignore beta error) are only slightly smaller — the dominant uncertainty is the predictive regression itself.
- **Weak-identification-robust test** (`RunTest.m:508–586`; App §M; MT §4.2): **Stock-Wright (2000) S-statistic.** For each σ on the grid `sigs` (≈ 0.1…10), re-estimate (θ, δ) exactly, then `S(σ) = m_Euler' (V̄_Test)^{-1} m_Euler` with the `K` held-back conditional Euler moments; `S(σ) ~ χ²_L` (`thresh = chi2inv(0.95, K)`). Invert over the grid → 95% confidence **set** for σ (equivalently the IES `1/σ`). Conclusion: fails to reject for σ ≳ 5 (excl. 2020) / σ ≳ 2 (incl. 2020), i.e. IES ∈ ~[0.1, 0.2]. (Kleibergen/AR proper are *not* used; the implementation is specifically Stock-Wright with the PSD-weight-matrix generalization.) A separate **Wald test** of `R_{0,t}=R_{b,t}+const` (RF coeff = 1, others 0) is reported with asymptotic and bootstrap p-values; the **Olea-Pflueger (2013) effective F** is reported for weak-instrument pre-testing.
- **Bootstrap** (`Bootstrap.m`; App §H): a **VAR(1) residual (iid) bootstrap**, *not* block or wild. Fit `[F_t; Z_t; π_t] = μ + Φ[F_{t-1};Z_{t-1};π_{t-1}] + η_t` on the excl-2020 sample; **resample residuals iid with replacement** (`datasample`), preserving the contemporaneous outcome-regressor innovation correlation (the Stambaugh 1999 / Bauer-Hamilton 2018 finite-sample channel) while `Φ` preserves regressor persistence; reconstruct from the **Feb-1973** initial values; re-run the full estimation. **1000 reps × 2 samples** (574 obs incl. 2020 / 562 excl.). Used for (a) finite-sample p-values of the Wald and consumption-predictability (effective-F) tests, and (b) a **power analysis** against the alternative "zero-beta rate no better than the T-bill" calibrated to ρ(real ZB, E[Δc]) = 0.4442. **Fixed base seed 12302018**, `RandStream('mlfg6331_64')`. Runtime: ~3–6 min/rep/core (README) — days on a single core.
- **Ridge / K-fold** (`K_Fold.m`, App §E): used **only** when `RunRidge=true` (off in main). 10-fold contiguous-block CV selects the ridge penalty ψ (shrinking γ_1→1, γ_2..L→0; level γ_0 unpenalized) by out-of-sample squared surprise return.

### 5.3 Risk factors (7) — MT §3; App §C.2
1–5. **Fama-French (2015) five factors** — Mkt, SMB, HML, **RMW** (profitability), **CMA** (investment) — Ken French library.
6. **Treasury term factor** — Fama 6–10y Treasury portfolio (CRSP, maturity >60 & ≤120 mo) minus 1-mo T-bill (App fn. 62).
7. **Default factor** — ICE BofA 15+ yr US Corporate total return (FRED `BAMLCC8A015PYTRIV`) minus the 10+ Fama Treasury portfolio (App fn. 63).
(Consumption is **not** a factor in the main spec — added only in robustness — so consumption data does not enter `R_{0,t}` construction, sharpening the Euler test.)

### 5.4 Instruments / predictors (5) — MT §3; App §C.3. `Z_t` (lagged one extra month vs. returns, `RunTest.m:50–53`):
1. **T-bill yield** `RF` (1-mo, Fama-French).
2. **Rolling 12-month inflation** `CPI_rolling` (12-mo MA of monthly log-change in core CPI, FRED `CPILFESL`).
3. **Term spread** `TSP` (10y `GS10` − 1-mo bill).
4. **Excess bond premium** `EBP` (Gilchrist-Zakrajšek 2012, updated Favara et al. 2016).
5. **Unemployment** `UMP` (U3, FRED `UNRATE`).
Instruments are standardized to unit variance before estimation and de-scaled afterward (`RunTest.m:173–186`, `ZSDiag`). Robustness instruments: CAPE, shadow spread (Lenel et al. 2019 via GSW), BAA−AAA spread, lagged consumption, realtime Sahm index.

### 5.5 Universe, sample, frequency
- **Universe:** US common stocks (shrcd 10/11) in CRSP matched to COMPUSTAT, bottom-20%-by-ME excluded, formed into the 130 portfolios; **never individual stocks in the GMM**.
- **Frequency:** monthly. **Sample:** Jan 1973 start (first month all instruments available); returns predicted from Mar 1973. **Incl. 2020 → 574 months; excl. 2020 → 562 months.** Everything reported both ways because of the >6-SD (one ~17-SD) COVID consumption months; excl-2020 treated as the better guide.

### 5.6 Consumption data — App §C.4
Real per-capita **nondurables + services** growth (Jagannathan-Wang style, `jw_m.csv`), `Δc_t = 100·log(c_t/c_{t-1})`. NIPA Table 2.8.5 (nominal PCE), deflator 2.8.4, population 2.6. Robustness: **nondurables only** (Hall, `hall_m.csv`).

### 5.7 Other tuning constants / vintage choices
WRDS pull pinned to `wrds==3.2.0`; COMPUSTAT BE = `SEQ + TXDITC − PS` (PS = coalesce PSTKRV/PSTKL/PSTK; TXDITC = coalesce TXDITC / TXDB+ITCB / 0; BE→NaN if ≤0); OP = `(REVT−COGS−XINT−XSGA)/BE`; INV = `AT/lag_AT − 1`; delisting `retadj=(1+ret)(1+dlret)−1`; June/Dec Fama-French timing; σ grid capped at 10 (power declines for larger σ); MATLAB R2023b+ required.

---

## 6. Information inventory for another sample (e.g., euro area) — *checklist of required inputs, NOT a plan*

To port the methodology to another market, the following analogues are required. (Listed as inputs only; hardest gaps flagged ⚠.)

- **Equity universe + accounting match**: a CRSP-equivalent monthly return panel (with delisting returns and share/exchange codes) matched to a COMPUSTAT-equivalent fundamentals source, to build size, book-to-market, operating profitability, and investment characteristics. *(EU: Compustat Global / Worldscope / Refinitiv; ⚠ cross-country share-class, currency, and listing-venue handling; ⚠ a clean delisting-return analogue.)*
- **Pre-sorting market beta**: a market-index return for 5-yr monthly CAPM betas. *(EU: a euro-area market portfolio — definition choice: pan-EA vs. country-level.)*
- **Test portfolios**: either self-built 3×3×3 (β×size×{B/M, OP, INV}) value-weighted sorts on the universe above **plus** industry portfolios, or a pre-built library. ⚠ **No off-the-shelf Ken-French-style euro-area 3×3×3 beta-sorted set exists** — these must be constructed, and breakpoint conventions (which exchange defines breakpoints? the NYSE analogue?) must be chosen.
- **Risk factors (7 analogues)**: FF5 for the region (e.g., AQR/French international or self-built), a **term factor** (long- vs. short-Treasury — EA sovereign choice ⚠: which sovereign curve in a multi-country bloc?), and a **default factor** (long corporate − long sovereign total returns; ⚠ EA IG corporate index).
- **Instruments (5 analogues)**: short rate, 12-mo trailing core inflation (HICP), term spread, **an excess-bond-premium analogue** ⚠ (the hardest gap — GZ/EBP requires firm-level corporate bond + issuer data; an EA EBP would need either an existing study's series or de-novo construction à la GZ), and unemployment.
- **Covariance estimator**: Ledoit-Wolf (2020) `analytical_shrinkage.m` — language-portable, no data dependency.
- **Consumption + deflator + population**: NIPA-equivalent monthly real per-capita nondurables+services (⚠ EA: monthly consumption is largely **quarterly** in Eurostat national accounts — a monthly indicator/interpolation would be needed, a real gap), with an HICP-type deflator and population.
- **Monetary-shock series** (only for the MP exercise): an EA high-frequency policy shock (e.g., Altavilla et al. EA-MPD) and/or a narrative series.
- **Sample/frequency**: monthly; a 1999+ euro sample is much shorter than the US 1973–2020, which bears on the weak-ID S-test power.

**Hardest gaps:** (1) a euro-area **EBP analogue**; (2) **monthly consumption** (EA data are quarterly); (3) **self-built 3×3×3 beta-sorted portfolios** with defensible breakpoint conventions in a multi-country bloc; (4) the **term/default factors'** sovereign-reference choice.

---

## 7. Data provenance & programmatic accessibility

Two lines configure everything (`DataCode/path_variables.md`: `main_path`, `WRDS_USERNAME`). The repo ships all non-WRDS inputs in `Raw Data/Archive.zip`; **only WRDS is pulled at runtime**. There is **no FRED API, no `pandas_datareader`, no download URL, and no `requests`/`urllib` anywhere in the pipeline** — every non-WRDS series is a hand-staged file (filenames carry provenance suffixes like `_FRED`, `_Shiller`, `_CapitalIQ`).

### 7.1 Programmatically pulled — WRDS via the `wrds` package (subscription-gated)
`A1_CCM_download.py` opens `wrds.Connection(...)` and issues six `conn.raw_sql(...)` queries (the **only** programmatic pulls in the whole pipeline; `CCM.py` merely reads the resulting CSVs):

| WRDS table | Role | Output |
|---|---|---|
| `comp.funda` | Compustat annual fundamentals (BE/OP/INV inputs); `indfmt='INDL', datafmt='STD', popsrc='D', consol='C', datadate≥1959` | `FF5_comp.csv` |
| `crsp.msf` | Monthly stock returns/prices/shares; `date≤2020-12-31` | `crsp_m.csv` |
| `crsp.mseall` | Monthly events: shrcd, exchcd, SIC/NAICS, divamt; `exchcd∈[1,3]` | `mseall.csv` |
| `crsp.msedelist` | Delisting returns | `dlret.csv` |
| `crsp.ccmxpf_lnkhist` | CCM link history (drop NR/NP/NU links) | `ccm_link.csv` |
| `crsp.msi` | Monthly market index (vwretd/vwretx) | `msi.csv` |

**Access barrier:** a paid **WRDS institutional subscription** entitled to CRSP, Compustat, and CCM (plus possible Duo 2FA). Not a stable public URL.

### 7.2 Manual / special-access (bundled in `Archive.zip`, read from `Raw Data/`)
- **WRDS "Bond Returns by WRDS"** → `Raw Data/Factors/CRSP bond returns.csv`. README lists it as a WRDS requirement, but it is **not** in the `A1` pull — it is a **manual WRDS web-query export** (subscription-gated). Feeds the term and default factors.
- **FRED-website CSVs** (free, easily re-automatable via FRED API / `fredgraph.csv`, but currently hand-downloaded): `BAA`, `AAA`, `GS10`, `TB3MS`, `CPILFESL`, `UNRATE`, `SAHMREALTIME`, `ExcessBondPremium` (GZ/EBP), `BAMLCC8A015PYTRIV` (ICE corporate), `DTB3`.
- **Ken French Data Library** (free, stable Dartmouth URLs, hand-downloaded): `ff3.xlsx`, `F-F_Research_Data_5_Factors_2x3.csv`, `49_Industry_Portfolios.CSV`, (`6_Portfolios_2x3.xlsx` — present but its caller is commented out).
- **Shiller CAPE** `CAPE_Shiller.xls` (free, Yale/Shiller site).
- **Federal Reserve GSW yield curve** `feds200628.csv` (free, Fed FEDS; used for the shadow-spread instrument).
- **BEA NIPA Excel** tables 2.8.5 / 2.8.4 / 2.6 (free, but **layout-dependent** hard-coded `skiprows/usecols/nrows` make re-download fragile → effectively manual).
- **Nakamura-Steinsson** `NS.xlsx` and **Romer-Romer** `RR_monetary_shock_monthly.dta` (+ Greenbook inputs and a Stata `.do`) — authors' replication archives, no API (MP exercise only).
- **S&P Capital IQ LIBOR** `LIBOR_3M_CapitalIQ.xls` — **licensed** (hardest barrier) — but **currently unread** (dormant input). Same for **Lettau CAY** `CAY_Lettau.csv` (free, unread).

### 7.3 Automatability verdict
- **Central hard barrier:** WRDS (CRSP/Compustat/CCM + the manual Bond-Returns export). Everything genuinely needed beyond WRDS is **free public data** whose only "barrier" is that the shipped pipeline reads it from hand-staged files rather than fetching it — all are straightforwardly re-automatable via FRED / Ken French / Shiller / Fed / BEA endpoints. A future automated rebuild would replace the `Raw Data/*` staging with FRED-API + French/Shiller/Fed/BEA downloaders and keep only WRDS as a credentialed pull.
- **Dormant/orphaned items** (note for a clean rebuild): `CAY_Lettau.csv`, `LIBOR_3M_CapitalIQ.xls` (unread); `generate_value_spread` + `6_Portfolios_2x3.xlsx` (commented out); and standalone `shadow_spread.py` / `consumption.py` / `monetary_shocks.py` / `price_dividend_ratio.py` are superseded duplicates of the live functions in `auxiliary_functions.py` (the standalone shadow-spread variant even uses a different maturity and RF source).

---

## Appendix — code-file → role map (quick reference)

| File | Role |
|---|---|
| `DataCode/A1_CCM_download.py` | WRDS SQL pulls → `Raw Data/CCM/` |
| `DataCode/A2_beta_estimation.py` | 5-yr rolling CAPM pre-sorting betas |
| `DataCode/CCM.py`, `00_create_data.py`, `auxiliary_functions.py` | CCM merge, characteristics, instruments/factors/consumption assembly → `Input/` |
| `DataCode/portfolio_construction.py` | 3×3×3 beta-sorted + 49-industry test portfolios; bottom-20% cut |
| `MatlabCode/RunAll.m` | Driver: main spec + all robustness specs |
| `MatlabCode/RunTest.m` | Per-spec orchestration: load, lag, standardize, estimate, S-test sweep, figures, save |
| `MatlabCode/BetaSigma.m` | OLS (α,β); `Σ_R` via Ledoit-Wolf preconditioning + analytical shrinkage |
| `MatlabCode/PortfolioWeight.m` | Min-variance, unit-investment, zero-beta weights (eq. 10) |
| `MatlabCode/InstrumentGMMWrapperConc.m` | GMM solve (concentrated over γ), SEs, S-statistic |
| `MatlabCode/InstrumentGMMConcOpt.m` | Concentrated GMM objective (`K+2` exactly-identified moments) |
| `MatlabCode/InstMomentsConc.m` | The three moment blocks |
| `MatlabCode/InstGMMSE.m`, `InstMomentsFullJacobian.m`, `InstMomentsProjectionJacobian.m`, `Wopt.m` | Sandwich SEs, mixed Jacobian, moment covariance (`COV`/NW) |
| `MatlabCode/K_Fold.m` | 10-fold CV for ridge penalty (ridge spec only) |
| `MatlabCode/Bootstrap.m`, `BootstrapGraphs.m` | VAR(1) residual bootstrap (1000×2), Wald/F p-values, power |
| `MatlabCode/MonetaryShock.m`, `EquityAnalysis.m`, `Model/monetary_shock_model.m` | MP-shock and equity-premium/NPV exercises (downstream of the core estimate) |
| `ExternalCode/analytical_shrinkage.m` | Ledoit-Wolf (2020) nonlinear shrinkage |
| `ExternalCode/jacobianest.m` | Numerical Jacobian (Euler block of SE) |
| `StataCode/MakeTables.do` | Assembles output `.csv` → paper tables |

---

# Addendum C — How the betas are estimated (Assignment 05c)

**All confirmed from code, not inferred.** There are **two distinct beta objects** in the pipeline, computed in different files, with different windows, factor sets, and units. They are **not the same object**.

## C-1. The two betas at a glance

| | **Sorting beta** | **GMM / estimation beta `β`** |
|---|---|---|
| Code | `DataCode/A2_beta_estimation.py` | `MatlabCode/BetaSigma.m` |
| Unit | **Individual stock** | **Test portfolio** (the ~130) |
| Window | **Rolling 60-month** (5-yr), min 24 mo | **Full sample**, all `T` months (static) |
| Frequency | Monthly | Monthly |
| Factors | **Market only** (single-factor CAPM) | **Full factor set** (7 in main spec) |
| Benchmark for "excess" | Stock return − risk-free (`eret = ret − rf`) | Returns in excess of the **zero-beta rate** (`R − ι·R_{0,t}`) |
| Constant? | Yes (`eret ~ mktrf`) | Yes (`alphas`) |
| Thin-trading correction | **None** | **None** |
| Role | Sorts stocks into beta terciles | Factor-loading matrix inside the GMM + min-variance portfolio |

## C-2. Static or dynamic? (Q1)
- **Sorting beta — rolling, but used annually.** `A2_beta_estimation.py:129` sets the window `min_dt = current_dt + MonthEnd(-60)` (a **60-month / 5-year trailing window**), requires **≥ 24 monthly observations** (`:134`, else NaN), and recomputes the beta for *every* month of each stock's history. In portfolio construction the beta is then read off only at the **June** formation date (`portfolio_construction.py` merges `betas` onto `ccm_june`), Fama-French style. So it is a rolling-window estimate, refreshed at each annual rebalancing.
- **GMM beta `β` — static (full-sample, constant).** `BetaSigma.m:28` computes `BetaFull = Xv \ X'`, a single OLS over **all `T` months** (`T = size(R,2)`, `:6`). The betas are **constant over the sample** — there is no rolling/expanding/conditional window. `β` *does* get recomputed at every GMM iteration, but only because the dependent/independent variables are taken in excess of the current parameter-implied zero-beta rate (`zbrate`, `:9`), not because the betas are time-varying. The **only** time-varying-beta option is the **`VaryingBetas`** robustness spec (`RunTest.m:189–205`), which appends factor×instrument interactions as extra factors — making loadings linear in the instruments `Z_t` — and is **off in the main spec** (`opts.VaryingBeta=false`). No DCC / GARCH / state-dependent specification anywhere.

## C-3. The regression itself (Q2)
- **Sorting beta:** `result = smf.ols(formula='eret ~ mktrf', data=est_df)` (`A2:145`) — monthly stock **excess** return on the **contemporaneous market excess return only**, with an intercept. Single-factor CAPM; β = `params['mktrf']` (`:146`). No SMB/HML/etc. (The docstring's mention of censoring returns to (−50%, 100%) is **commented out** in shipped code; the docstring also says "one to three years (two minimum)" but the code uses a 60-month window with min 24.)
- **GMM beta:** `BetaSigma.m:23–30` — portfolio returns in excess of the zero-beta rate (`X = R − ι·zbrate`) regressed on the **full factor set** in excess of the zero-beta rate (`Xm = Rm − ι_M·zbrate`) plus a constant (`Xv = [Xm; ones]`); `Beta = BetaFull(1:end-1,:)`, `alphas = BetaFull(end,:)`. In the main spec the factor set is the 7 factors (Mkt, SMB, HML, RMW, CMA, term, DEF); only the market is treated as a unit-investment factor (`iotaM`). **No lead-lag / Dimson (1979) / Scholes-Williams (1977) correction** — a repo-wide grep for `dimson|scholes|shift(|lead|lag…beta|nonsync` returns nothing for either beta.

## C-4. Where in the pipeline (Q3)
- **(a) Sorting:** the `A2` stock-level CAPM betas define the **beta terciles** of the 3×3×3 test portfolios. Critically, the beta breakpoints are formed **within each size×characteristic cell**, not unconditionally (`portfolio_construction.py:107–163`), with the in-code rationale: *"the small stocks may have inaccurate betas due to liquidity reasons"* (`:107–109`).
- **(b) Estimation:** the `BetaSigma.m` portfolio-level betas are the factor-loading matrix `β` used (i) in the orthogonal projector `H(β)` and the zero-beta moment (`InstMomentsConc.m`), and (ii) in the minimum-variance, factor-neutral portfolio weights `w(θ) = Σ_R^{-1}[ι,β](…)^{-1}e_1` (`PortfolioWeight.m`).
- **Same object? No.** The sorting beta is a *single-factor, rolling, stock-level* CAPM beta used only to *bucket* stocks; the estimation beta is a *full-factor, full-sample, portfolio-level* loading measured in excess of the *zero-beta rate*. They share neither factor set, window, unit, nor benchmark.

## C-5. Implication for the illiquidity-attenuated-beta concern (Q4)
Ties directly to §3–§4 of the main memo. The estimation choices push **against** the worst form of attenuation but **do not eliminate** the residual portfolio-level channel:
- **Monthly (not daily) returns — the main mitigant.** Non-synchronous-trading attenuation is severe at daily frequency and modest at monthly. Both betas use monthly returns; the paper explicitly attributes the **NoDrop20 ≈ Main** robustness to monthly betas reducing the liquidity impact (MT fn. 12). This is the single biggest reason attenuation is *less* likely to drive the result.
- **5-year sorting window** gives enough observations for a stable stock beta, at the cost of spanning regimes.
- **Within-size-group beta tercile sort *is* the explicit compensation.** Forming beta breakpoints inside each size cell prevents liquidity-attenuated small-stock betas from contaminating the *cross-stock* beta ranking — so a "low-beta" bucket reflects genuine low beta within a liquidity class, not merely small/illiquid names. The code comment names exactly this motivation.
- **What remains uncorrected.** Neither beta applies a Dimson/Scholes-Williams lead-lag adjustment. So *residual* attenuation in the monthly betas is not corrected, and — because attenuation from systematic stale pricing does **not** fully wash out under aggregation — the **portfolio-level `β` can still be biased toward the cross-sectional average** for the less-liquid (small × low-characteristic) portfolios. The minimum-variance optimizer, whose objective is to find the lowest-variance factor-neutral combination, can then overweight portfolios that *appear* factor-neutral partly because their loadings are attenuated. This is precisely the residual concern flagged in main-memo §3–§4; the within-size sort and the 20% cut blunt it but the GMM `β` itself carries no thin-trading correction.

**Decision to flag for Angus.** The natural robustness probe — should a future replication estimate the **portfolio-level `β` with a Dimson/Scholes-Williams (lead-lag) correction** (and/or re-run **CovI** and **Ridge**) to bound how much of the zero-beta rate's variation could be attenuation/over-fitting? The current code does not, and this is the cleanest test of whether the illiquidity channel matters for the headline. (Documentation only — not attempted here.)
