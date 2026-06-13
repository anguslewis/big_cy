# Companion memo — Spec-by-spec data inputs & Ken-French-only feasibility

**Author:** `big_cy-zerobeta` (Assignment 05b, follow-up to E0) · **Date:** 2026-06-13
**Status:** Documentation / analysis only. No replication, no estimation code, no implementation plan.
**Reads with:** `scratch/zerobeta_methodology_memo.md` (E0) and the cloned repo `scratch/TheZeroBetaRate/`. Spec list verified against `MatlabCode/RunAll.m`; Ken French availability verified against the live KF data library (June 2026).

> **Bottom line up front.**
> 1. **Every** GMM specification in the paper builds its test assets from the **same custom security-level construction** — 81 beta-sorted 3×3×3 portfolios (market-beta terciles formed *within* size groups, crossed with size and one of B/M, OP, INV; value-weighted; bottom-20%-by-market-cap excluded) **plus the 49 Ken French industry portfolios**. Only the 49-industry block is KF-published; the 81 beta-sorted portfolios are not. So **no spec is *exactly* reproducible from Ken French portfolios alone.**
> 2. **But a *rough* KF-only US version is more feasible than expected**, because KF *does* publish **"25 Portfolios Formed on Size and Market Beta"** (5×5) and univariate beta-sorted portfolios. A KF-only US test-asset set — {25 Size×Beta} ∪ {25 Size×B/M} ∪ {25 Size×OP} ∪ {25 Size×INV} ∪ {49 Industry} — captures the *beta dispersion* DiTella argue is essential, at the cost of the joint beta×characteristic sorting and the 20% small-cap exclusion. **Rough verdict (b) for all US specs.**
> 3. **The equity-factor side is KF-only for the FF3Only / MktOnly specs** (Mkt, SMB, HML all from KF); the main 7-factor spec additionally needs the **two bond factors** (term, default), which are *not* KF (proxiable from FRED/CRSP).
> 4. **Instruments `Z` are never KF** (except the risk-free rate, which lives in the KF factor files) — they are the universal binding non-KF input (FRED/Shiller/Fed). **Consumption is never KF** — but, crucially, **consumption is *not* needed to recover the zero-beta rate path or the excess-volatility / safe-rate-puzzle result** (see §B-2.1); it is required only for the σ/IES S-test and the consumption-factor specs.
> 5. **Other regions:** KF publishes FF5 + momentum + size×{B/M, OP, INV} portfolios for Europe / Japan / Asia-Pacific / North America / Developed / Emerging — but **no international beta-sorted portfolios and no international industry portfolios**, from **July 1990** only. So a ported estimate loses exactly the two test-asset blocks that matter most (beta dispersion + industries) and runs on a much shorter sample. This is the crux to flag for Angus.

---

## B-1. Q1 — The GMM specifications and a spec × data-input matrix

All specifications are `RunTest` calls in `MatlabCode/RunAll.m` with a varied `opts` struct. Every spec is automatically estimated on **both** the incl-2020 (574 mo) and excl-2020 (562 mo) samples inside `RunTest.m`; `Pre2020Only` is a separate flag. The estimation machinery (concentrated GMM over γ, min-variance zero-beta portfolio, Stock-Wright S-test, sandwich SEs) is identical across specs — they differ only along six axes:

- **(A) test-asset set**, **(B) factor set**, **(C) instrument set `Z`**, **(D) consumption factor**, **(E) covariance estimator**, **(F) ridge on/off**.

### B-1.1 The distinct test-asset sets (only three are actually used)
| Asset file | N | Composition | Used by |
|---|---|---|---|
| `FF5_plus_Industry_Portfolios_Nominal.csv` | **130** | 81 beta-sorted (β×size×B/M, β×size×OP, β×size×INV; 27 each) + 49 industry; bottom-20%-ME **excluded** | **Main** + all factor/instrument/cov/cons/ridge variants (default) |
| `NoDrop20/FF5_plus_Industry…csv` | 130 | same sorts, bottom-20%-ME **included** | **NoDrop20** only |
| `27_plus_Industry_Portfolios_Nominal.csv` | **76** | 27 (β×size×B/M) + 49 industry | **FF3Industry** only |
| *(AltData BLLP / Paul / InformativeFactors)* | — | — | commented out in shipped code |

The headline point for Q3: **all three sets contain custom beta-sorted portfolios**; even the 76-portfolio FF3 set uses the custom β×size×B/M sort.

### B-1.2 The full specification list and the data-input matrix
Legend for data sources: **KF** = Ken French library; **FF/KF-RF** = risk-free in the KF factor file; **WRDS** = CRSP/Compustat security-level + WRDS bond returns; **FRED**, **Shiller**, **GSW** (Fed feds200628), **BEA** (NIPA). "Test assets = WRDS" because the 81 beta-sorted portfolios require security-level CRSP/Compustat to build, *plus* the 49 industry from KF.

| Spec (`opts.Name`) | Test assets | Factors | Instruments `Z` | Cov | Ridge | Cons factor | Distinct data inputs vs Main |
|---|---|---|---|---|---|---|---|
| **Main** | 130 (drop20) [WRDS+KF] | FF5 + term + DEF (7) [KF + WRDS/FRED] | RF, UMP, EBP, TSP, CPI_roll [KF-RF + FRED] | LW | no | none | — (baseline) |
| **NoDrop20** | 130 (no drop) [WRDS+KF] | 7 | base 5 | LW | no | none | small-cap filter off |
| **FF3Industry** | 76 [WRDS+KF] | 7 | base 5 | LW | no | none | drops OP/INV beta-sorts |
| **MktOnly** | 130 | **Mkt only (1)** [KF] | base 5 | LW | no | none | factor set |
| **FF3Only** | 130 | **Mkt,SMB,HML (3)** [KF] | base 5 | LW | no | none | factor set |
| **VaryingBetas** | 130 | 7 × (1+5 instr.) = **42** | base 5 | LW | no | none | instrument×factor interactions |
| **WithConsSigma1/5/10** | 130 | 7 + **nonlinear cons factor** [BEA] | base 5 | LW | no | NL (σ=1/5/10) | adds consumption factor |
| **LinearCons** | 130 | 7 + **linear cons factor** [BEA] | base 5 | LW | no | linear | adds consumption factor |
| **NDOnly** | 130 | 7 | base 5 | LW | no | none | consumption = nondurables only [BEA] |
| **AltBAAS** | 130 | 7 | RF, UMP, **BAAS**, TSP, CPI_roll [FRED] | LW | no | none | BAA−AAA spread replaces EBP |
| **AltCAPE** | 130 | 7 | base 5 + **CAPE** (6) [Shiller] | LW | no | none | adds CAPE |
| **AltCAPERidge** | 130 | 7 | base 5 + CAPE (6) | LW | **yes** | none | CAPE + ridge |
| **AltDP** | 130 | 7 | base 5 + **DP_ratio + CAPE** (7) [CRSP idx + Shiller] | LW | no | none | adds dividend-price + CAPE |
| **AltSahm** | 130 | 7 | RF, **SAHM**, EBP, TSP, CPI_roll [FRED] | LW | no | none | Sahm index replaces unemployment |
| **LagCons** | 130 | 7 | base 5 + **lagged Δc** (6) [BEA] | LW | no | none | adds lagged consumption growth |
| **Shadow** | 130 | 7 | base 5 + **shadow_spread** (6) [GSW + FRED] | LW | no | none | adds shadow-rate spread |
| **AltCPI** | 130 | 7 | RF, UMP, EBP, TSP, **CPI** (lagged, not rolling) [FRED] | LW | no | none | inflation instrument variant |
| **Ridge** | 130 | 7 | base 5 | LW | **yes** | none | ridge penalty (10-fold CV) |
| **CovSample** | 130 | 7 | base 5 | **Sample** | no | none | raw sample covariance |
| **CovI** | 130 | 7 | base 5 | **I** (≈ equal-wt) | no | none | identity covariance |
| **CovDiag** | 130 | 7 | base 5 | **Diag** | no | none | factor + diagonal residual |
| **CovPCA** | 130 | 7 | base 5 | **PCA (3)** | no | none | 3-PC residual structure |
| *(InfFactors)* | 130 (pre-2020) | Amman et al. 5 + term + DEF | base 5 | LW | no | none | optional; commented out |

**Reading the matrix:** the *test-asset block* and the *macro-instrument/consumption blocks* are the data-source bottlenecks. The factor set is the one axis where some specs (MktOnly, FF3Only) become **fully KF on factors**. Cov/Ridge/instrument variants change *no* data-source class relative to Main (they re-weight or re-instrument the same inputs), except instrument swaps that pull a different FRED/Shiller/GSW series.

---

## B-2. Q2 — What differs and what is robust across specifications

From the paper's robustness battery (App §J, Tables 7–10; MT §6) and the E0 reading:

### B-2.1 Spec-invariant conclusions (the headline survives)
- **High average real zero-beta rate (~8–9%/yr)** and **~7.6%/yr spread over the T-bill** — stable across NoDrop20, FF3Industry, NDOnly, the factor variants, and the instrument variants (App Tables 7–9 levels cluster near the Main ~8.55%).
- **Large excess volatility / the "safe-rate puzzle"** — the zero-beta rate is far more volatile and far better at explaining the consumption-claim price/dividend ratio than the T-bill, across specs.
- **Cyclicality:** the zero-beta rate **falls in bad times** and co-moves with expected consumption growth — direction is robust.
- **S-test:** fails to reject the consumption Euler equation for high σ (low IES) across the LW-covariance specs.
- **Not a small-stock artifact:** **NoDrop20** (include the bottom two deciles) ≈ Main — the single most direct robustness check for Angus's prior.

### B-2.2 Spec-sensitive results (what hinges on a choice)
- **Covariance estimator** is the most consequential axis for *inference*: **CovI** (identity → portfolio closest to equal-weighted, the least min-variance-concentrated) is the one spec where the **S-test behaves materially differently** (flatter, lower S-statistic; App Fig. 34), while **CovSample** and **CovPCA** preserve the headline. So the min-variance/Ledoit-Wolf weighting matters for the *shape* of the σ confidence set — though the *level* and *excess-volatility* results are broadly cov-robust.
- **Ridge** deflates the in-sample variation: the excess-volatility share of the consumption-claim p/d ratio falls from **29% → 19%** (MT §6.1) — still substantial, but a signal that part of the predictive variation is fragile / over-fit.
- **Including vs excluding 2020** shifts the σ confidence set (fail-to-reject for σ ≳ 2 incl. 2020 vs σ ≳ 5 excl.), because of the extreme COVID consumption months; excl-2020 is the authors' preferred guide.
- **Factor parsimony (MktOnly / FF3Only):** broadly preserve the level and consumption-tracking, but with fewer factors the zero-beta portfolio is neutralized on fewer dimensions (more residual exposure) — these are robustness illustrations, not the preferred spec.

### B-2.3 A pivotal fact for KF-only feasibility
**The zero-beta rate path and the excess-volatility/safe-rate result do *not* require consumption data.** In the concentrated GMM (`InstMomentsConc.m`), the rate `R_{0,t} = γ'Z_t + R_f + R_{b,t}` is pinned by the **K+1 predictive moments** (zero-beta-portfolio surprise × instruments), which involve only the portfolios, factors, and instruments. Consumption enters only the **one unconditional Euler moment** (identifying the SDF level parameter δ/ρ, which does *not* appear in the rate formula) and the **K conditional Euler moments** used for the σ/IES S-test. **Implication:** the level/path of the zero-beta rate and the excess-volatility magnitude are recoverable from **test assets + factors + instruments only** — BEA consumption (and hence the IES test) can be deferred. This materially widens what is feasible without NIPA.

---

## B-3. Q3 — Ken-French-only feasibility, by input class and by spec

### B-3.1 Test-asset portfolios
**What KF publishes (US):** "Portfolios Formed on Market Beta" (univariate deciles), **"25 Portfolios Formed on Size and Market Beta"** (5×5), "25 Size × B/M", "25 Size × OP", "25 Size × INV" (all 5×5; also 2×3 6-portfolio versions), and the **5/10/12/17/30/38/48/49 industry** portfolios — all value-weighted, NYSE breakpoints.

**What DiTella build (custom, security-level):** 81 portfolios from three 3×3×3 sorts that put **market-beta terciles *inside* each size×characteristic cell**, with a **bottom-20%-by-market-cap exclusion**, value-weighted — plus the 49 KF industry portfolios.

**Mapping:**
| DiTella block | KF-published analog? | Verdict |
|---|---|---|
| 49 industry portfolios | **Yes — identical** (DiTella *use* the KF 49-industry file) | **(a) exact** |
| 81 beta-sorted 3×3×3 (β within size × {B/M, OP, INV}) | **No.** KF has Size×Beta (5×5) and Size×{char} (5×5) *separately*, but not beta crossed with a second characteristic, not beta-tercile-within-size-cell, and never with a 20%-cap exclusion | **(c) not exact**; **(b) approximable** via separate KF Size×Beta + Size×char sets |
| bottom-20%-ME exclusion | **No** — not a feature of any KF portfolio (KF uses NYSE breakpoints, so micro-caps sit in the small-extreme buckets but are not dropped) | **(c)** |

**Rough KF-only US test-asset set (the closest buildable analog, prose not code):**
> {25 Size×Market-Beta} ∪ {25 Size×B/M} ∪ {25 Size×OP} ∪ {25 Size×INV} ∪ {49 Industry} ≈ 130–149 value-weighted KF portfolios.
- **What it preserves:** wide **market-beta dispersion** (the property DiTella stress as essential — a portfolio with β≈1 is insensitive to the zero-beta *level*), characteristic dispersion across all FF5 dimensions, and the full industry block. The Ledoit-Wolf min-variance machinery runs unchanged on any test-asset matrix.
- **What it loses / its compromises:** (i) **no joint beta×characteristic sorting** — KF's Size×Beta set has no B/M-OP-INV control and vice versa, so the precise factor-loading spread of DiTella's 3×3×3 cells is not reproduced; (ii) **no bottom-20% exclusion** — micro-cap influence is only mitigated (value-weighting + NYSE breakpoints), not removed (note this makes the analog conceptually closest to the **NoDrop20** spec, which DiTella show ≈ Main); (iii) 5×5 = 25 cells per sort vs DiTella's 3×3×3 = 27; (iv) breakpoint conventions differ slightly.

### B-3.2 Factors
| Factor | KF? | Note |
|---|---|---|
| Mkt-RF, SMB, HML, RMW, CMA | **Yes (KF)** | US (1963/1963/…) and international |
| Momentum (not used by DHKW, available) | Yes (KF) | — |
| Treasury **term** factor (Fama 6–10y − 1-mo bill) | **No** | WRDS/CRSP bond returns; *proxiable* from FRED Treasury yields/returns (rough) |
| **Default** factor (ICE 15y+ corp − long Treasury) | **No** | FRED ICE index + CRSP bond returns |
**Consequence:** **MktOnly** and **FF3Only** are **fully KF on the factor side**. The Main / 7-factor specs need the two bond factors (non-KF, FRED-proxiable). The bond factors also matter because DHKW include them partly to justify the term-spread and EBP *instruments* — dropping them is a modelling change, not just a data swap.

### B-3.3 Instruments `Z`
| Instrument | KF? | Source |
|---|---|---|
| RF (T-bill yield) | **Yes** — in the KF factor file | KF/FF |
| Rolling/lagged CPI inflation | No | FRED `CPILFESL` |
| Term spread (10y − 3m/1m) | No | FRED `GS10`, `TB3MS` |
| Excess bond premium (GZ) | No | FRED EBP |
| Unemployment / Sahm | No | FRED `UNRATE` / `SAHMREALTIME` |
| CAPE | No | Shiller |
| Shadow spread | No | Fed GSW `feds200628` + FRED `DTB3` |
| BAA−AAA | No | FRED |
| Dividend-price ratio | No (≈) | CRSP index; loosely proxiable from KF market return w/ and w/o dividends |
**Consequence:** **the instrument block is the universal binding non-KF input** (RF excepted). No zero-beta specification can be estimated from KF portfolios *alone* — every spec needs a macro/financial instrument set from FRED/Shiller/Fed. This is unavoidable and spec-independent.

### B-3.4 Consumption / macro
BEA NIPA (nondurables+services, deflator, population) — **never KF**. Needed **only** for the consumption-factor specs and the σ/IES **S-test**, **not** for the zero-beta rate path or the excess-volatility result (§B-2.3).

### B-3.5 Per-specification KF-only verdict
No spec is category **(a) exactly reproducible from KF portfolios alone** (all use custom beta-sorted test assets + non-KF instruments). Layered verdict:

| Spec(s) | Test assets | Factors | Instruments | Consumption | Overall KF-only verdict |
|---|---|---|---|---|---|
| **MktOnly**, **FF3Only** | (b) rough-KF (Size×Beta+char+industry) | **KF-only** | non-KF (FRED) | not needed for rate/excess-vol | **(b)** — the *most KF-leaning* specs: equity factors + test assets all KF (rough); only the macro instruments [FRED] are non-KF if you skip the IES test |
| **Main, NoDrop20, FF3Industry, VaryingBetas, Cov*, Ridge, Alt-instrument specs** | (b) rough-KF | KF equity factors **+ 2 bond factors (non-KF)** | non-KF | not needed for rate/excess-vol | **(b)** — rough-KF test assets + factors, but bond factors and instruments require FRED/CRSP |
| **WithCons*, LinearCons, NDOnly, LagCons** | (b) rough-KF | KF + cons factor (BEA) | non-KF | **required (BEA)** | **(b/c)** — additionally require BEA consumption as a *factor*; the S-test/IES content is intrinsically non-KF |
| *(any spec, the IES/σ S-test component)* | — | — | — | **required (BEA)** | **(c)** for the Euler/IES test specifically — needs consumption, never KF |

**Practical synthesis for Angus:** the maximally-KF path is **FF3Only (or MktOnly) on a KF Size×Beta + Size×char + industry test-asset set**, producing the **zero-beta rate level/path and the excess-volatility / safe-rate-puzzle magnitude** from **KF (test assets + factors) + FRED/Shiller/Fed (instruments)**, with **BEA consumption deferred** (only needed to add the IES S-test). The irreducible non-KF inputs are (i) the macro instrument set and (ii) consumption-for-the-Euler-test.

---

## B-4. Q4 — Other-region extension

### B-4.1 What KF publishes internationally (verified)
For **Europe, Japan, Asia-Pacific ex Japan, North America, Developed, Developed-ex-US** (23 developed countries) and a separate **Emerging** set, from **July 1990**:
- **FF3 + FF5 factor sets** and a **momentum factor** (Emerging: 5-factor + momentum).
- **Size × B/M, Size × OP, Size × INV** portfolios (2×3 and 5×5; Emerging 2×3 only) and **32-portfolio three-way sorts**.
- **No international beta-sorted portfolios** (the Size×Beta set is **US-only**).
- **No international industry portfolios** (5/…/49-industry are **US-only**).
- Returns are **USD-denominated** and the factor files use the **US one-month T-bill as RF**.

### B-4.2 Portability of the KF-only forms from Q3
A ported, KF-leaning estimate (FF5 international factors + KF international size×char portfolios as test assets) is feasible in *form*, but loses exactly the two test-asset blocks that carry the most identifying weight in DHKW's design:

| Required block | US (KF) | Europe/Japan/AsiaPac (KF) | Consequence |
|---|---|---|---|
| Beta-sorted test assets | Size×Beta available | **Absent** | **Lose the market-beta dispersion** DHKW call essential for identifying the zero-beta *level* via the market channel; cannot rebuild without security-level data |
| Industry portfolios | 49 available | **Absent** | Lose ~38% of the 130-portfolio diversification block |
| Size×{B/M,OP,INV} | 5×5 | 5×5 (2×3 Emerging) | Available — characteristic dispersion preserved |
| 3-way 32-portfolio sorts | — | available | Partial substitute for joint sorting, but not beta-targeted |
| FF5 / momentum factors | yes | yes (from 1990) | Available |
| RF | KF (US bill) | KF files use **US** bill; a region-specific short rate is non-KF | Currency/RF convention decision needed |

### B-4.3 Binding non-KF inputs for a non-US estimate
- **Region instruments `Z`:** regional short rate, trailing inflation (e.g. HICP), term spread, **an EBP analog** (hardest — no off-the-shelf Gilchrist-Zakrajšek for Europe/Japan), unemployment — all from regional central banks / FRED / OECD, none from KF.
- **Region consumption (for the IES test only):** monthly real per-capita nondurables+services — **the binding macro gap**: euro-area / Japan consumption in national accounts is largely **quarterly**, so a monthly indicator or interpolation is required. (Not needed for the rate path itself — §B-2.3.)
- **Currency / RF convention:** KF international returns are USD with US RF; a credible region-specific zero-beta rate would want own-currency returns and the region's own short rate — a modelling decision, and own-currency single-stock returns are exactly the messy security-level data the KF route is meant to avoid.
- **Sample length:** international KF starts **July 1990** (~35 yrs) vs DHKW's 1973–2020 (~47 yrs); the Stock-Wright S-test power (already modest) degrades with the shorter sample.

### B-4.4 Most portable specification
A **FF5, LW-covariance** spec on **KF international size×{B/M,OP,INV} (5×5) + 32-portfolio sorts** as test assets, with **regional FRED/OECD instruments**, deferring consumption (hence the IES test) — would **port** to Europe / Japan / Asia-Pacific / Developed. It would deliver a regional zero-beta *rate path* and an *excess-volatility / safe-rate-puzzle* magnitude, **but**:
- without the **beta-sorted** test assets (US-only in KF) it departs most from DHKW's identification design — the single most important caveat;
- without **industry** portfolios it loses a large diversification block;
- it inherits the USD/US-RF convention unless adjusted; and
- the IES/Euler validation requires region-specific (likely interpolated) monthly consumption.

---

## B-5. Decisions to flag for Angus

1. **Is a KF Size×Beta + Size×char + industry substitute an acceptable "rough" US test-asset set?** It preserves beta dispersion but loses the joint beta×characteristic sorting and the 20% small-cap cut. (Conceptually closest to the **NoDrop20** spec, which DHKW show ≈ Main — mildly reassuring.)
2. **For other regions, is losing *both* beta-sorted *and* industry portfolios acceptable?** This is the crux: DHKW argue beta dispersion is essential for identifying the zero-beta *level*. If it is judged essential, a credible non-US estimate likely **cannot** avoid security-level construction after all — which would defeat the motivation for the KF route. If a FF5-characteristic-portfolio set is judged "good enough," the KF route is viable post-1990.
3. **Currency / risk-free convention** for international KF data (USD returns + US RF) — accept, or insist on own-currency returns (→ back to messy security-level data)?
4. **Scope of the target object:** a **rate-path + excess-volatility** estimate is achievable without consumption (KF + FRED/regional macro). Adding the **σ/IES Euler S-test** requires monthly consumption (US: BEA; international: a quarterly→monthly gap). Decide whether the IES test is in-scope for a first non-US pass.

**Sources (KF availability):** [Ken French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html); [Description of Fama/French International Factors](http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/details_global.html).
