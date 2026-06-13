# VBDG options measure vs. the zero-beta rate — deep dive & reconciliation

*Memo for big_cy (T1a follow-up, 01b). Close reading of two papers:*
- *van Binsbergen, Diamond & Grotteria (2022, JFE 143(1):1–29), "Risk-Free
  Interest Rates"* — `vanbinsbergendiamondgrotteria2022riskfree` (note: JFE, not
  JF; the working-paper version is "van Binsbergen et al. 2019").
- *Di Tella, Hébert, Kurlat & Wang (2025, JPE), "The Zero-Beta Interest Rate"* —
  `ditella2025zerobeta`.

*Neutral framing. The T1a headline tension: two ABSOLUTE convenience-yield
measures disagree ~10× — VBDG's option box-spread gives ~40 bp, DiTella's
equity-based zero-beta rate gives ~2–4%. This memo asks why, and whether either
paper reconciles them. **Bottom line up front:** the gap is largely NOT a
contradiction — the two measure convenience for **different marginal investors /
pricing kernels**, and DiTella (who post-date and explicitly cite VBDG)
reconcile the two via **leverage**: the box rate is a leveraged intermediary's
near-riskless funding rate, the zero-beta rate is the unlevered required return
on riskless equity, and the two are bridged by the bank leverage ratio (~8×).
This is a substantive calibration decision for big_cy — see §5.*

---

## 1. VBDG methodology deep dive

### 1.1 The box-spread construction (what it prices, and why it is convenience-free)

VBDG infer a risk-free rate from **European index options** via put–call parity,
*without* using any safe-asset price. For each date `t`, maturity `T`, and strike
`K_i`, put–call parity for European options states

```
p_{i,t,T} − c_{i,t,T} = (P_{t,T} − S_t) + exp(−r_{t,T} T) · K_i        (VBDG eq. 1)
```

where `p`, `c` are put and call prices, `S_t` the index level, and `P_{t,T}` the
present value of dividends paid to `T` (the "dividend strip"). Regressing the
put-minus-call price on the strike `K_i` across the cross-section of strikes:
- the **slope** is `β = exp(−r_{t,T} T)` ⇒ `r_{t,T} = −(1/T) ln β`;
- the **intercept** is `P_{t,T} − S_t` (so dividends need not be known a priori).

Two estimators:
- **Estimator 1** — cross-sectional OLS of `(p−c)` on `K` at each `(t,T)`.
- **Estimator 2** — the **"Box trade"**: for every strike pair `(K_i, K_j)`,
  `r = −(1/T) ln[((p_i−c_i)−(p_j−c_j))/(K_i−K_j)]`, then take the **Theil–Sen
  median** over all pairs (robust to outliers). The box trade — buy puts/write
  calls below the average strike, write puts/buy calls above — delivers a
  **riskless payoff `K_i − K_j`** for a known cost, i.e. a synthetic zero-coupon
  bond manufactured purely from (risky) options.

**Why it is free of Treasury-style convenience.** The box payoff is a
collateralized derivative position, not a money-like safe asset. It cannot be
posted as HQLA, used as repo collateral the way Treasuries are, or otherwise
perform the monetary/regulatory/safety roles that earn Treasuries their
convenience premium. So the box rate is the rate of return an investor *forgoes
by holding a convenience-bearing safe asset instead of the convenience-free box*.
The **convenience yield = box rate − Treasury yield**.

The cross-sectional fit is essentially perfect (R² = 0.9999998), so put–call
parity holds in the data and the implied rate is a single well-identified number
across all strikes — VBDG read this as evidence that option prices are *not*
distorted by violations unrelated to convenience.

### 1.2 Data, estimation, and the benchmark that makes the number ~40 bp

- **Data:** all CBOE option trades and quotes on **SPX** (and **DJX** for
  robustness), European-style, **2004–2018**, at **minute** frequency; bid, ask,
  strike, maturity. Compute the rate each minute from midpoints, then take the
  **daily median** across minutes. Maturities span **~1 month to ~3 years**.
- **Robustness instruments (Section 3):** precious-metals futures (a storage
  cost-of-carry arbitrage) and SPX **futures** (an option–futures arbitrage,
  eq. 14–16) — both yield similar ~40 bp average convenience, confirming the rate
  is not an SPX-options artifact.
- **The benchmark and the magnitude.** VBDG fit Nelson–Siegel–Svensson curves to
  Box rates and to Treasuries (Gürkaynak-Sack-Wright). The headline:
  - **Convenience yield (Box − government bond) ≈ 35–40 bp on average** across
    maturities beyond 3 months (Fig. 6: 38 bp at the short end, ~35 bp at the
    long end; the term structure of convenience is **roughly flat**).
  - **Larger at short maturities** — an extra ~25 bp for T-bills vs. notes/bonds
    (so the very short Treasury convenience is ≈ 40 bp + 25 bp).
  - **Quadruples in the 2008 crisis** — ~97–100 bp at 12–30 months and ~160 bp
    below 3 months (Fig. 7).
  - Box sits **above** government bonds, **OIS**, and **GC repo** (Table 3, 3-mo:
    Box 1.74%, gov 1.42%, OIS 1.39%, GC repo 1.41% ⇒ Box−gov ≈ 32 bp, Box−OIS ≈
    35 bp, Box−repo ≈ 33 bp) — i.e. OIS and repo *themselves* carry convenience.
  - Box ≈ **LIBOR** pre-2008, **below** LIBOR after (LIBOR then embeds bank
    credit risk, not convenience).

So VBDG's "convenience yield" is measured **relative to Treasuries**, with the
**Box rate as the convenience-free benchmark** that sits above all the
near-money rates (Treasury, OIS, repo).

### 1.3 Identifying assumptions and caveats

**Assumptions for the rate to be risk-free:**
1. **Put–call parity holds** (it does, R²≈1) — so a unique discount rate prices
   all strikes; deviations would otherwise contaminate the rate.
2. **Midpoints reflect true option values** — needed so the box cost is the true
   riskless cost. VBDG show bid–ask spreads have negligible explanatory power for
   the rate level (SPX vs. the much less liquid DJX give the same average rate;
   regressions of the DJX–SPX spread on bid–ask measures have R²≈4%).
3. **The box trade is genuinely riskless** — no counterparty risk: CBOE/OCC
   margin and variation-margin rules, OCC super-senior status and AAA rating,
   Fed liquidity backstop. The CBOE's margin rules *recognize* the box as
   riskless, so it can be cheaply levered (this matters in §3).

**Stated caveats:**
- **Short maturities only (≤ ~3 yr).** Convenience beyond 2.5 yr cannot be
  measured without extrapolation; whether 10–30 yr convenience behaves similarly
  is "an open question."
- **The rate may be specific to option-market intermediaries.** VBDG are explicit
  that their rate "might measure the risk-free cost of capital (or funding cost)
  for expert financial intermediaries alone, because our rates are inferred from a
  sophisticated asset class." The primary option market-makers are leveraged
  intermediaries (Table 2: Susquehanna, Citadel, hedge funds).
- **Option-market demand pressure.** The rate could in principle be distorted by
  leverage/demand for derivatives (Garleanu-Pedersen-Poteshman; Hazelkorn et al.
  find their futures basis correlates 0.8 with the VBDG convenience measure),
  though VBDG argue cross-market consistency (options, futures, metals) limits this.

### 1.4 What ABSOLUTE object VBDG actually measure (Appendix A)

VBDG's Appendix A models a **hedge fund** (the primary option market-maker)
choosing a portfolio subject to financial frictions (debt overhang, costly equity
issuance, special repo collateral). The convenience yield `s_i` of asset `i` is
the *extra amount the prime broker will lend against it* (special collateral —
e.g. government debt at low "special repo rates"). For any set of risky assets
`j` that replicate a riskless payoff (`Σ w_j δ_j = 1` — a box),

```
s_i = p_i − Σ_j w_j p_j           (VBDG eq. A.17)
```

i.e. the specialness of asset `i` = its price minus the price of the
synthetic-riskless replication. The box rate is the **risk-free cost of capital
for the hedge fund's pricing kernel `m_h`** — *"the correct risk-free rate for a
pricing kernel for risky asset markets in which our hedge [fund] is active, even
though that pricing kernel is distorted by financial frictions."*

**So VBDG measure: a tradable, no-arbitrage, near-riskless borrowing/lending rate
— the funding/cost-of-capital rate of a leveraged intermediary active in options —
free of Treasury money/collateral/regulatory convenience.** It is an *absolute*
level (priced against a convenience-free tradable claim), but it is the
*intermediary's debt-pricing* rate, short-horizon, and the convenience it
recovers is that of Treasuries **over a tradable riskless claim**, not over the
household's intertemporal consumption price.

VBDG flag this themselves in the **conclusion**: *"if the risk-free cost of
capital in options markets differs from the time value of money for households,
[this] would imply the convenience yield on a safe asset depends on which
investor buys it. Understanding how the costs of capital and convenience yields
provided by safe assets vary across investors is a promising direction for future
research."* — This is exactly the seam to the zero-beta rate.

---

## 2. What object does each measure capture? (side-by-side)

| | **VBDG Box rate** | **DiTella zero-beta rate** |
|---|---|---|
| **Object** | Tradable near-riskless borrowing/lending rate, from options via put-call parity | Intertemporal price of consumption = inverse mean growth of the SDF that prices equities, `R_{0,t} = E_t[Λ_{t+1}/Λ_t]^{-1}` |
| **What "the rate" is** | The cost of manufacturing a synthetic riskless payoff (box) | Expected return on a unit-investment, **minimum-variance** equity portfolio orthogonal to the SDF |
| **Marginal investor / kernel** | Option-market **intermediary** (hedge-fund) kernel `m_h` — prices **debt / arbitrage** | The SDF `Λ` that prices **equities** (≈ household / representative equity investor) — `η=0` for equities, `η>0` for bonds |
| **Identification** | **No-arbitrage** (put-call parity); model-free; needs no risk model | Requires a **linear factor SDF** (5 FF + bond + default), estimated betas, min-var zero-beta portfolio, GMM with macro instruments |
| **Tradable / leverable?** | **Yes** — box is legally/practically zero-risk; cheaply levered (CBOE margin ≈ a margin deposit) | **No** — a zero-beta equity portfolio is not recognized as riskless and cannot be cheaply levered |
| **Horizon** | Short (≤ ~3 yr); flat term structure | Long-horizon intertemporal price |
| **Convenience measured vs Treasuries** | Box − Treasury ≈ **35–40 bp** (≤160 bp in crisis) | Zero-beta − T-bill ≈ **2–4%** (real zero-beta ~8.3%; ~7.6%/yr over real T-bill) |
| **Level of the benchmark itself** | Box is **above** Treasury/OIS/repo but **far below** zero-beta | Zero-beta is the **highest** of the three |
| **Robustness of the *level*** | Very high (R²≈1; cross-market consistent) | Authors call the **level their least-robust result** (an omitted high-priced factor uncorrelated with the market could bias it); the **Euler-equation fit** is the robust result |

**The ordering (Diamond 2020, quoted by DiTella fn. 48):** there are *multiple*
risk-free rates, ordered —
```
Treasury / safe-asset rate   <   intermediary's risk-free rate   <   household's equity-kernel rate
(household demand for safety)     (≈ Box / OIS; prices debt)          (= zero-beta rate; prices equities)
```
The two "convenience yields" are `(Box − Treasury)` ≈ 40 bp and
`(ZeroBeta − Treasury)` ≈ 2–4%. The first is the *lower* rung of the second; they
are nested, not competing estimates of one number.

---

## 3. Does either paper reconcile the gap? (Yes — explicitly, DiTella → VBDG)

### 3.1 VBDG (2022, earlier) on the zero-beta rate
VBDG predate DiTella's JPE paper and engage the zero-beta idea via **Frazzini &
Pedersen (2014)** (intro, p. 3). Their position:
- Measuring a zero-beta rate **requires taking a stance on the risk factor** the
  beta is computed against. *"If the factor does not capture all risks relevant to
  investors, the zero-beta rate includes a risk-premium component."*
- By contrast, the box rate *"does not require specifying any particular risk
  model and implies a smaller spread than the spread estimates in their paper."*
- *"The spread we estimate therefore measures the tightness of leverage
  constraints in any multi-factor generalization of their model."*

So VBDG frame their ~40 bp as the **clean, model-free** number and imply a larger
zero-beta-type spread is **partly risk premium** from an incomplete factor model —
**and** is a measure of **leverage-constraint tightness**.

### 3.2 DiTella (2025, later) on VBDG — the explicit reconciliation (§7, pp. 48–49)
DiTella **cite van Binsbergen et al. [2019]** (the VBDG working paper) and Diamond
(2020), and reconcile the two numbers directly. Three moves:

1. **Different objects; no reason to be equal.** *"A spread of 75 basis points
   between a AAA-rated corporate bond and Treasury bills is perfectly compatible
   with a spread of 6% between the equity-based zero-beta rate and Treasury bills.
   Both spreads could be called 'convenience,' but there is no reason to expect
   them to be the same"* (citing Diamond 2020's multiple-risk-free-rates result
   and Angeletos et al. 2023).

2. **Leverage is the bridge.** A bank that can lever a riskless investment ~8×
   requires only a small spread on the *levered* position even if its required
   return on *unlevered* riskless equity is high. Formally (fn. 47): if the bank
   borrows at `X` bp over T-bills and needs 600 bp on equity, it requires
   `(7/8)·X + (1/8)·600` bp on a leverable riskless investment ≈ 75 bp. **So a
   ~40–75 bp leveraged-arbitrage/AAA spread and a ~4–6% unlevered zero-beta spread
   are two views of the same friction**, related by the leverage ratio.

3. **The two are *connected* — a lower bound.** *"One could use arbitrage
   opportunities together with knowledge of banks' leverage constraints to
   estimate a lower bound on the zero-beta rate as perceived by banks."*
   Boyarchenko et al. (2018) do exactly this with the JPY–USD CIP violation: a
   bank levering it earns *"an average safe return on equity of over 4% above safe
   rates"* — a lower bound on the bank-perceived zero-beta rate. DiTella note *"a
   similar calculation is in principle possible for ... the spread between Box
   rates and OIS rates [van Binsbergen et al. 2019]."*

**Footnote 49 is the crispest statement of *why the box number is small*:** *"The
key difference between a zero-beta equity portfolio and ... a Box rate trade is
that the latter is recognized in law and practice as zero-risk and hence can be
easily levered by a bank ... From a bank's perspective lending via the Box trade
is equivalent, in terms of risk and leverage, to a margin account deposit at the
CBOE."* The box rate is therefore the bank's **near-riskless funding rate** (≈
OIS), **not** its required return on (unlevered) riskless equity. The zero-beta
equity portfolio is **not** recognized as riskless and cannot be levered the same
way, so its rate ≈ the unlevered required equity return (high).

### 3.3 Net: the reconciliation
The ~10× gap is **largely not a contradiction**. The two measure convenience for
different marginal investors and at different leverage:
- **Box (~40 bp)** = Treasury convenience over a *tradable riskless intermediary
  funding rate* (the bank's leveraged near-riskless cost of capital).
- **Zero-beta (~2–4%)** = the wedge between safe bonds and the *household /
  equity-SDF intertemporal price* (the unlevered required return on riskless
  equity).
- They are bridged by **leverage** (~8×: 40–75 bp × 8 ≈ 3–6%) and ordered by
  **Diamond's (2020) ladder** of multiple risk-free rates. Both papers ultimately
  attribute the gap to **leverage constraints / safe-asset demand**; neither
  number is "wrong."

The residual disagreement is about **how much of the zero-beta wedge is
convenience vs. an omitted risk premium** (VBDG's worry) — which DiTella's
min-variance, multi-factor-orthogonal construction is *designed* to suppress, and
which they argue is second-order for the Euler-equation test even if it bites for
the *level*.

---

## 4. Candidate reconciliations, as open questions (neutral)

Even granting §3, several distinct (non-exclusive) channels remain genuinely open;
each implies a different "true" absolute level:

1. **Marginal investor / pricing kernel (Diamond 2020 ladder).** Box = the
   intermediary's debt-pricing kernel; zero-beta = the equity-pricing/household
   kernel. *Open:* which kernel is the relevant "intertemporal price" for a given
   modeling question — and how far apart the two kernels' riskless rates are.

2. **Leverage scaling (DiTella's explicit channel).** Box ≈ levered-arbitrage
   funding rate; zero-beta ≈ unlevered cost of riskless equity; bridged by the
   bank leverage ratio. *Open:* the leverage ratio (and hence the implied
   household-kernel wedge) is calibrated, not measured; Boyarchenko et al. give a
   *lower* bound (>4%), not a point estimate.

3. **Risk-premium contamination of the zero-beta level (VBDG's channel).** If the
   factor SDF omits a high-priced factor uncorrelated with the market, the
   zero-beta *level* is biased up. *Open:* DiTella concede the level is their
   least-robust result; the Euler-equation fit is robust, the level less so.

4. **Horizon / term structure.** Box ≤3 yr with a flat convenience term
   structure; zero-beta is a long-horizon intertemporal price. *Open:* whether
   convenience (and any term/duration premium bundled into the equity wedge) rises
   with horizon — VBDG cannot measure convenience beyond ~2.5 yr.

5. **Tradability / leverability (fn. 49).** The box lacks exactly the
   money/collateral/regulatory convenience of Treasuries (a narrow object); the
   zero-beta wedge additionally reflects the gap between the household's
   equity-kernel rate and *any* tradable riskless rate (a broad object). *Open:*
   how to decompose the broad wedge into "Treasury-specific convenience" vs.
   "stock–bond SDF segmentation / missing-risk compensation."

6. **What's bundled into each "convenience" label.** Box−Treasury is narrowly
   Treasury vs. a tradable riskless claim; ZeroBeta−Treasury may bundle
   convenience + a missing-risk-compensation/segmentation component (DiTella's
   §6.4 multiple-SDF point; their list of stock-vs-bond-SDF facts). *Open:* no one
   has cleanly decomposed the *absolute* equity-based wedge.

---

## 5. Flags for Angus (decision points)

**This is a calibration-target decision, not a "which paper is right" question.**
For big_cy's bonds-in-utility convenience wedge, the relevant absolute anchor
depends on *whose* intertemporal margin the model's safe-bond convenience is meant
to capture:

- **If the convenience wedge is the gap between safe bonds and the HOUSEHOLD's
  intertemporal (consumption-Euler) price** — i.e. the object a representative-agent
  bonds-in-utility model most naturally maps to — then **DiTella's ~2–4% (equity
  zero-beta) is the right anchor**, and the small CIP/Treasury-basis/Box numbers
  are *lower rungs* (intermediary funding rates), consistent with a large absolute
  household wedge.

- **If the convenience wedge is the gap between safe bonds and a tradable riskless
  intermediary funding rate** — relevant if big_cy's mechanism runs through
  *leveraged intermediaries* (as in the Gabaix-Maggiori / Kekre-Lenel-GFC
  lineage) — then **VBDG's ~40 bp is the right anchor**, and DiTella's number is
  the *unlevered* cost of equity, related by leverage.

- **DiTella's own reconciliation supports a large absolute wedge for leveraged
  agents:** Boyarchenko et al.'s levered-CIP calculation gives a bank-perceived
  zero-beta rate of **>4% above safe rates** — i.e. even the *intermediary's*
  unlevered riskless-equity return is several %, not 40 bp. So a DiTella-magnitude
  absolute wedge is defensible **for the household/equity kernel**, while the
  ~40 bp box number is specifically the *levered funding* rate.

**Recommended framing for the model (for T1c / Angus to rule on):** treat the
absolute convenience yield as a **household (equity-kernel) object** and anchor to
DiTella magnitudes, while explicitly acknowledging that (i) the *tradable*
near-riskless rate intermediaries face carries only ~40 bp of convenience, and
(ii) the ~10× gap is the leverage ratio between the two, per DiTella §7. The
neutral survey must continue to carry the ~40 bp box number as the principal
*intermediary-side* counterpart wherever the DiTella number is cited — they are
two rungs of one ladder, and big_cy's choice of rung should be stated explicitly
and defended on which marginal investor prices the model's safe bonds.

*(No reconciliation requires discarding either estimate; the open empirical
questions in §4 — especially the leverage ratio and the convenience-vs-missing-risk
decomposition — are the live ones.)*

---

## 6. References added to `references.bib`
All citation details taken from DiTella et al.'s reference list (no web inference):
- `boyarchenko2018bankintermediated` — Boyarchenko, Eisenbach, Gupta, Shachar &
  Van Tassel, "Bank-Intermediated Arbitrage," FRBNY Staff Report 858 (2018) — the
  levered-CIP lower bound (>4%) on the bank-perceived zero-beta rate.
- `bakerwurgler2015capital` — Baker & Wurgler, "Do Strict Capital Requirements
  Raise the Cost of Capital?...", AER 105(5):315–20 (2015) — the zero-beta-rate
  cost-of-capital interpretation.
- `bakerhoeyerwurgler2020leverage` — Baker, Hoeyer & Wurgler, "Leverage and the
  Beta Anomaly," JFQA 55(5):1491–1514 (2020).
- `herrenbrueck2019interest` — Herrenbrueck, "Interest Rates, Moneyness and the
  Fisher Equation," SFU WP (2019) — a liquidity-based safe-rate-puzzle theory.

**Correction (also applied):** `diamond2020safety` is **published in the Journal
of Finance** (75(6):2973–3012, 2020) — DiTella's reference list confirms it; it
was previously a "verify-before-use" WP entry. Updated to `@article` and removed
from the verify list.
