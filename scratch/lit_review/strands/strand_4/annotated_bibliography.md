# Strand 4 — Relative vs. Absolute Measurement of Convenience Yields

**Annotated bibliography.** Compiled for big_cy. Neutral, descriptive survey.

The organizing question of this strand is *what each empirical convenience-yield
measure actually identifies*. The recurring distinction:

- A **RELATIVE** (gap) measure prices one near-safe asset against another
  (dollar vs. euro; Treasury vs. AAA corporate; Treasury vs. foreign govt bond;
  nominal vs. TIPS). Whatever convenience is **common to both legs nets out**, so
  the measure recovers only the *difference* in convenience between the two
  legs, not the absolute level on either.
- An **ABSOLUTE** (level) measure prices a safe asset against a constructed
  **convenience-free** benchmark — a synthetic risk-free rate backed out of
  risky-asset prices (option box spreads; a zero-beta equity portfolio). The
  spread of the safe yield below that benchmark is the *level* of the
  convenience yield.

Each entry records: WHAT two legs the spread/basis is between; whether it is
RELATIVE or ABSOLUTE and WHY; the magnitude; what it does NOT identify; and a
proximity score (1–5) to big_cy's relative-vs-absolute question.

Citekeys match `references.bib`.

---

## Category A — Cross-currency convenience (CIP basis): RELATIVE measures across currencies

### A1. Du, Tepper & Verdelhan (2018, *Journal of Finance*) — `dutepperverdelhan2018cip`
**Proximity: 5.**
- **What it measures:** The **cross-currency basis** — the gap between the
  direct dollar interest rate and the synthetic dollar rate obtained by
  borrowing in foreign currency and swapping into dollars via FX forwards. Legs:
  dollar money-market/cash rate vs. FX-swap-implied dollar rate (Libor, repo, and
  the credit-risk-free KfW-bond variant).
- **Relative or absolute:** **RELATIVE** — a pure currency-pair gap. CIP nets
  the dollar leg against the swapped foreign leg; any convenience common to both
  safe legs cancels. It identifies the *cross-currency* convenience/funding
  wedge, not the absolute convenience of either currency's safe asset.
- **Magnitude:** Average annualized |basis| ≈ 24 bps at 3M, 27 bps at 5Y
  (2010–2016, G10). Five-year JPY basis ≈ −90 bps end-2015. Larger at quarter-ends
  (balance-sheet/regulatory cost). Persists in the credit-risk-free KfW/repo basis.
- **Does NOT identify:** the absolute level of convenience on the dollar (or on
  any single currency's safe asset); it nets that out. Also abstracts from why
  the gap exists in levels — attributes it to post-GFC intermediary balance-sheet
  costs + imbalances in funding supply/demand across currencies.
- **Bridge to big_cy:** the canonical statement of the "relative" measure. The
  basis is *the gap*, not the level. Central exhibit for the strand's thesis.

### A2. Du & Schreger (2020, *Handbook of International Economics* chapter) — `duschreger2020handbook`
**Proximity: 5.**
- **What it measures:** Surveys CIP-deviation measurement (cross-currency basis
  across benchmark rates, currencies, maturities) and its drivers/implications.
  Also reviews **government-bond CIP** as a sovereign-credit / market-segmentation
  object.
- **Relative or absolute:** **RELATIVE** throughout — explicitly frames CIP
  deviations as an intermediation *fee* (a gap) charged to supply dollar funding.
  Clients pay because they face asymmetric funding conditions across currencies.
- **Magnitude:** Reviews the DTV facts; first PC explains 51% of 5Y G10 basis
  variation, −61% corr with broad dollar, −33% with VIX, +56% with S&P.
- **Does NOT identify:** absolute convenience levels; frames everything as
  near-arbitrage gaps tied to intermediary leverage constraints (IOER–fed funds,
  triparty–GCF repo, swap spread, Treasury cash-futures, CDS-bond).
- **Bridge to big_cy:** authoritative menu of relative cross-currency measures;
  enumerates the co-moving family of fixed-income "near-arbitrage" gaps.

### A3. Du, Im & Schreger (2018, *Journal of International Economics*) — `duimschreger2018treasury`
**Proximity: 5.**
- **What it measures:** The **U.S. Treasury Premium** — the CIP deviation
  *between government bond yields*: the yield gap between US Treasuries and
  currency-hedged foreign (G10) government bonds. Legs: US Treasury vs. swapped
  foreign sovereign.
- **Relative or absolute:** **RELATIVE** — a *difference* in safe-asset
  convenience between the US and foreign sovereigns. Convenience common to all
  advanced-economy safe govt bonds nets out; only the US-minus-foreign wedge
  survives.
- **Magnitude:** 5Y premium ≈ 21 bps pre-GFC, up to ~90 bps in the crisis,
  ≈ −8 bps post-crisis (it disappeared/reversed after the GFC).
- **Does NOT identify:** the absolute level of US (or foreign) convenience —
  only the gap between them, which is small and even negative post-crisis. A key
  cautionary fact for big_cy: relative measures can be near zero even if the
  common absolute level is large (2–4%).
- **Bridge to big_cy:** the cleanest illustration that the "Treasury premium"
  measured relative to other safe govt bonds is *not* the absolute convenience
  level; it is the cross-country differential.

### A4. Jiang, Krishnamurthy & Lustig (2021, *Journal of Finance*) — `jiangkrishnamurthylustig2021foreign`
**Proximity: 5.**
- **What it measures:** The **Treasury basis** — yield gap between US govt bonds
  and currency-hedged foreign govt bonds — as a proxy for the convenience yield
  foreigners derive from US safe assets; links it to the dollar exchange rate.
- **Relative or absolute:** **RELATIVE** — same construction as the Treasury
  Premium (US vs. swapped-foreign sovereign). Identifies *foreign-minus-US*
  convenience differential, interpreted as the marginal foreign investor's
  convenience on Treasuries.
- **Magnitude:** Treasury-basis variation accounts for up to 41% of quarterly
  dollar variation; basis widening → immediate USD appreciation then depreciation.
- **Does NOT identify:** the absolute convenience level; it measures the wedge
  that the *marginal cross-currency* investor perceives. Explicitly a relative
  object that the authors then map to the exchange rate.
- **Bridge to big_cy:** shows the relative basis *is informative for the dollar*
  but is still a gap; big_cy's claim is that the absolute level (DiTella) can be
  far larger than this measured gap.

### A5. Avdjiev, Du, Koch & Shin (2019, *AER: Insights*) — `avdjievdukochshin2019dollar`
**Proximity: 3.**
- **What it measures:** The triangular co-movement of the **broad dollar**, the
  cross-currency basis, and cross-border dollar bank lending. Basis is the same
  relative cross-currency gap as DTV.
- **Relative or absolute:** **RELATIVE** — uses the CIP basis (gap) and relates
  it to dollar strength as a barometer of bank risk-bearing capacity / shadow
  price of leverage.
- **Magnitude:** Stronger broad dollar ⇒ wider (more negative) basis; robust at
  daily frequency in levels and changes.
- **Does NOT identify:** absolute convenience; it interprets the *gap's*
  time-variation via the leverage-constraint channel.
- **Bridge to big_cy:** supply-side mechanism for why the relative gap moves;
  tangential to the level question.

### A6. Anderson, Du & Schlusche (2021, FEDS / NBER w28658) — `andersonduschlusche2021arbitrage`
**Proximity: 2.**
- **What it measures:** How global banks finance near-risk-free arbitrage
  (CIP basis, IOER–fed funds) with unsecured wholesale funding; response to the
  2016 MMF reform funding shock.
- **Relative or absolute:** Studies the **RELATIVE** CIP/IOER gaps as arbitrage
  spreads; about the *supply of arbitrage capital* that keeps the gaps from
  closing, not the level.
- **Magnitude:** Banks cut arbitrage positions (not loans) in response to the
  wholesale-funding shock — evidence that limited arbitrage capital sustains the gaps.
- **Does NOT identify:** absolute convenience.
- **Bridge to big_cy:** background on why relative gaps persist (limits to
  arbitrage); not a level measure.

---

## Category B — Domestic relative measures: safe vs. near-safe within one currency

### B1. Krishnamurthy & Vissing-Jorgensen (2012, *Journal of Political Economy*) — `krishnamurthyvissingjorgensen2012aggregate`
**Proximity: 5.**
- **What it measures:** The **convenience yield on Treasuries** via spreads of
  Treasuries against assets that differ in *either* liquidity *or* safety but not
  both (AAA/Aaa corporates, commercial paper, FDIC-insured deposits). Recovers a
  downward-sloping demand curve for Treasury convenience as a function of
  Treasury supply/GDP.
- **Relative or absolute:** **RELATIVE** — the AAA–Treasury (and similar)
  spread is the gap between Treasuries and a *next-best near-safe* asset. It
  identifies the convenience *of Treasuries over high-grade corporates*, not the
  convenience common to both. The "absolute" demand-curve estimate is still
  anchored to a near-safe comparison asset.
- **Magnitude:** Treasury convenience averages ≈ 73 bps (1926–2008); decomposes
  into liquidity and safety components, both downward-sloping in supply.
- **Does NOT identify:** the convenience common to AAA corporates and Treasuries
  alike (i.e., the part of the *absolute* level shared by both legs). If AAA bonds
  themselves carry substantial convenience, the spread understates the absolute level.
- **Bridge to big_cy:** the canonical domestic relative measure and the explicit
  source of the "Treasury convenience yield" language. big_cy's point: KVJ's
  spread is Treasury-minus-AAA, not Treasury-minus-(convenience-free rate).

### B2. Krishnamurthy (2002, *Journal of Financial Economics*) — `krishnamurthy2002bond`
**Proximity: 3.**
- **What it measures:** The **on-the-run/off-the-run ("bond/old-bond") spread**
  on 30Y Treasuries — a pure *liquidity* premium between two near-identical
  Treasury securities.
- **Relative or absolute:** **RELATIVE** — gap between two Treasuries differing
  only in liquidity. Identifies a liquidity component of convenience; both legs
  are Treasuries, so safety nets out entirely.
- **Magnitude:** ≈ 12 bps new-vs-old at issuance, collapsing toward ≈ 3 bps over
  the auction cycle; on-the-run premium varies with supply and aggregate liquidity.
- **Does NOT identify:** the safety component, or any absolute level.
- **Bridge to big_cy:** demonstrates the finest-grained relative measure (within
  Treasuries); shows convenience has separable liquidity vs. safety pieces.

### B3. Longstaff (2004, *Journal of Business*) — `longstaff2004flight`
**Proximity: 4.**
- **What it measures:** The **flight-to-liquidity premium** in Treasuries via
  the **Treasury–Refcorp spread**. Refcorp bonds are Treasury-guaranteed (so
  credit-risk-matched) but less liquid. Legs: Treasury vs. Treasury-guaranteed
  Refcorp.
- **Relative or absolute:** **RELATIVE** — a credit-risk-free *liquidity* gap.
  Because credit/safety is held fixed by the Treasury guarantee, it isolates
  pure liquidity convenience, not absolute level.
- **Magnitude:** Liquidity premium can exceed 15% of some Treasury values
  (large at long maturities); co-moves with consumer confidence, Treasury supply,
  mutual-fund flows.
- **Does NOT identify:** the safety component or the absolute level; only the
  Treasury-over-guaranteed-agency liquidity wedge.
- **Bridge to big_cy:** clean credit-controlled relative measure of liquidity
  convenience; complements KVJ's safety/liquidity decomposition.

### B4. Fleckenstein, Longstaff & Lustig (2014, *Journal of Finance*) — `fleckensteinlongstafflustig2014tips`
**Proximity: 4.**
- **What it measures:** The **TIPS–Treasury mispricing** — nominal Treasury vs.
  an inflation-swapped TIPS replicating identical nominal cash flows. Legs: nominal
  Treasury vs. synthetic nominal Treasury built from TIPS + inflation swaps.
- **Relative or absolute:** **RELATIVE** — a gap between two Treasury-issued
  instruments. Identifies the *extra* convenience of nominal Treasuries over
  TIPS (KVJ note TIPS are perceived as less money-like). Convenience shared by
  both Treasury types nets out.
- **Magnitude:** Mispricing averages 54.5 bps in yield (can exceed 200 bps);
  total ≈ $56bn, ≈ 8% of TIPS outstanding; Treasuries always *over*priced vs. TIPS.
- **Does NOT identify:** absolute level; only nominal-minus-TIPS convenience.
  Authors note Treasury–TIPS differentials therefore can't cleanly back out
  inflation expectations.
- **Bridge to big_cy:** another within-issuer relative gap; shows even
  same-issuer safe assets differ in money-likeness.

### B5. Klingler & Sundaresan (2019, *Journal of Finance*) — `klinglersundaresan2019negative`
**Proximity: 3.**
- **What it measures:** **Negative swap spreads** — the gap between the fixed
  rate on an interest-rate swap and the Treasury yield of the same maturity
  (30Y swap rate < 30Y Treasury yield since 2008). Legs: swap rate vs. Treasury yield.
- **Relative or absolute:** **RELATIVE** — swap-minus-Treasury gap. Negative
  spreads are usually read as evidence of a *positive* Treasury convenience
  (Treasuries rich), but K-S show pension duration demand + dealer balance-sheet
  limits drive it, complicating the convenience reading.
- **Magnitude:** 30Y swap spread persistently negative post-2008; tracks
  aggregate defined-benefit pension underfunding.
- **Does NOT identify:** a clean absolute (or even relative) convenience level —
  the spread is contaminated by hedging-demand and intermediary frictions.
- **Bridge to big_cy:** cautionary case — a popular "convenience" spread is
  partly a demand/segmentation artifact; relative measures are noisy proxies.

---

## Category C — ABSOLUTE-level measures: safe asset vs. a convenience-free benchmark

### C1. van Binsbergen, Diamond & Grotteria (2022, *Journal of Financial Economics*) — `vanbinsbergendiamondgrotteria2022riskfree`
**Proximity: 5.**
- **What it measures:** A **convenience-free risk-free rate** inferred from
  **box spreads on S&P 500 index options** (put-call parity / option-replicated
  riskless payoffs), at minutely frequency up to ~2.5–3Y. The implied **Treasury
  convenience yield** = (option-implied risk-free rate) − (Treasury yield).
- **Relative or absolute:** **ABSOLUTE (level).** This is the strand's key
  level-measure: the benchmark leg is a synthetic riskless rate from *risky*
  (option) markets that carries *no* Treasury-style convenience, so the spread to
  Treasuries is the *level* of Treasury convenience, not a gap to another safe asset.
- **Magnitude:** Treasury convenience ≈ **40 bps** on average; larger below 3M
  maturity; roughly **quadruples** during the financial crisis.
- **Does NOT identify:** convenience at long maturities (limited to ≈3Y); and the
  ~40 bps is an order of magnitude *below* the DiTella 2–4% zero-beta wedge —
  highlighting that option-box "absolute" measures and equity-based "absolute"
  measures disagree on the level.
- **Bridge to big_cy:** the cleanest *option-based* absolute measure. Its modest
  ~40 bps vs. DiTella's 200–400 bps is the central tension big_cy must engage:
  which "absolute" benchmark is the right intertemporal price.

### C2. Diamond & Van Tassel (2024, *Journal of Finance*) — `diamondvantassel2024riskfree`
**Proximity: 5.**
- **What it measures:** **Option-implied risk-free rates and convenience yields
  across 10–11 (G11) currencies**, inferring each currency's convenience-free
  rate from index-option prices, then the safe-bond convenience as the level
  spread in each country.
- **Relative or absolute:** **ABSOLUTE (level), country by country** — then
  *compares* the levels across countries. Provides both an absolute level per
  currency and a re-derivation of the cross-currency (relative) wedge.
- **Magnitude:** Convenience yields rise ~linearly with the level of a country's
  interest rate; US convenience is 5th-largest among G11; in crises convenience
  rises but the **US-minus-foreign difference does not**; option-implied CIP
  deviations are similar in size US-vs-each-other-country.
- **Does NOT identify:** a structural/macro intertemporal price (it is an
  option-implied financial rate, like vBDG); long maturities.
- **Bridge to big_cy:** the international extension of the absolute box-spread
  method — directly tests "absolute level per currency" vs. "relative gap," and
  finds the relative US premium small while absolute levels track rate levels.
  Core reference for big_cy's relative-vs-absolute decomposition internationally.

### C3. Di Tella, Hébert, Kurlat & Wang (2025, *Journal of Political Economy*) — `ditellahebertkurlatwang2025zerobeta`
**Proximity: 5.** *(Strand-1 anchor; included here for the bridge.)*
- **What it measures:** The **zero-beta interest rate** — expected return on an
  equity portfolio orthogonal to the SDF — as the true intertemporal price; the
  **safe-rate puzzle** is the wedge between this and the T-bill yield.
- **Relative or absolute:** **ABSOLUTE (level).** The convenience-free benchmark
  is built entirely from *equity* returns, so it embeds no safe-asset convenience.
  The zero-beta-minus-safe-rate gap is the *absolute* convenience/wedge on safe bonds.
- **Magnitude:** Zero-beta rate sits ≈ **2–4%** above safe bond yields; high,
  volatile, persistent; fits the consumption Euler equation where the T-bill does not.
- **Does NOT identify:** *which* friction generates the wedge or how it splits
  across asset types (deliberately model-free about the wedge's composition); it
  is a single aggregate level, not currency- or instrument-specific.
- **Bridge to big_cy:** the large-absolute-level claim that motivates big_cy. The
  gap to vBDG/Diamond-Van Tassel (~40 bps option-based vs. 200–400 bps equity-based)
  is the open empirical question this survey foregrounds.

### C4. Acharya & Laarits (2023, NBER w31863) — `acharyalaarits2023when`
**Proximity: 4.**
- **What it measures:** When Treasuries *earn* the convenience yield, from a
  **hedging** perspective — decomposes the stock-bond covariance into convenience,
  frictionless-rate, and default-risk terms.
- **Relative or absolute:** Engages the **absolute** convenience object (uses a
  frictionless-rate decomposition à la vBDG) but its contribution is about the
  *conditional* (state-dependent, beta/hedging) behavior of convenience, not a
  new static level.
- **Magnitude:** Convenience is low when Treasury-equity covariance is high;
  eroded by inflation expectations (substitution toward gold) and before
  debt-ceiling episodes; falls with Treasury supply.
- **Does NOT identify:** a single absolute level; focuses on time-variation and
  the hedging/negative-beta channel.
- **Bridge to big_cy:** connects the absolute convenience level to the *negative
  beta of safe/dollar assets* — directly relevant to KL's "dollar's negative beta"
  result that big_cy wants to reproduce via large convenience.

---

## Category D — Supply/demand & term-structure determinants (level vs. spread)

### D1. Greenwood, Hanson, Stein & Sunderam (2020, WP / later AER) — `greenwoodhansonsteinsunderam2020quantity`
**Proximity: 3.**
- **What it measures:** A quantity-driven (preferred-habitat, Vayanos-Vila)
  two-currency term-structure model: specialized investors absorb supply shocks
  to long-term bonds in two currencies; links bond term premia, FX, and CIP.
- **Relative or absolute:** Model of **spreads/premia driven by relative supply**
  — explains CIP deviations (relative) and term premia as compensation for
  bearing supply, not an absolute convenience level.
- **Magnitude:** Qualitative; matches QE→exchange-rate co-movement and
  post-2008 CIP deviations.
- **Does NOT identify:** an absolute convenience level; it is a supply-based
  theory of relative premia.
- **Bridge to big_cy:** supply-side micro-foundation for why convenience *spreads*
  vary with quantities; relevant to big_cy's "not enough dollar debt" mechanism.

### D2. Nenova (2025, WP / BIS) — `nenova2025global`
**Proximity: 4.**
- **What it measures:** **Demand elasticities and substitution patterns** for
  global government & corporate bonds (Koijen-Yogo demand system) using ~5,000
  granular mutual-fund portfolios; asks whether safe assets are *global* or
  *regional*.
- **Relative or absolute:** About **substitutability/segmentation** rather than a
  convenience-yield level per se; low *own* demand elasticities of US Treasuries /
  German Bunds proxy the *existence* of convenience benefits, and cross-elasticities
  reveal which assets are substitutes (relative positioning).
- **Magnitude:** US Treasuries substituted by *global* (corporate/EM) bonds;
  German Bunds substituted only within a narrow euro-area safe set; substitutability
  deteriorates in stress.
- **Does NOT identify:** a convenience-yield level; identifies elasticities and
  substitution geography.
- **Bridge to big_cy:** evidence on imperfect substitutability of dollar safe
  assets — supports big_cy's "insufficient dollar debt to satisfy US
  diversification demand" reinterpretation, and the global-vs-regional safe-asset
  distinction.

### D3. Du & Schreger (2016, *Journal of Finance*) — `duschreger2016local`
**Proximity: 3.**
- **What it measures:** The **local-currency sovereign credit spread** in EMs —
  local-currency sovereign bond yield over the synthetic local-currency risk-free
  rate built from cross-currency swaps. Legs: EM local-currency sovereign vs.
  swap-implied local risk-free.
- **Relative or absolute:** **RELATIVE** — a gap (default/segmentation) measured
  via the swap curve; methodologically the same swap-based synthetic-rate logic
  used in Treasury-premium work.
- **Magnitude:** Local-currency credit spreads positive and sizable; lower mean,
  lower cross-country correlation, lower global-risk sensitivity than FX-debt spreads.
- **Does NOT identify:** an absolute convenience level; it isolates sovereign
  credit/segmentation in local currency.
- **Bridge to big_cy:** shows the swap-synthetic-rate machinery in another
  setting; relevant to constructing convenience-adjusted foreign-currency yields.

---

## Category E — Currency wedge / theory context

### E1. Engel (2016, *American Economic Review*) — `engel2016exchange`
**Proximity: 4.**
- **What it measures:** Documents the joint **UIP puzzle** (high-rate currencies
  earn excess returns short-run) and the **level puzzle** (high-real-rate
  currencies are *stronger* than UIP warrants), arguing they are contradictory
  unless the risk premium / currency wedge has a particular structure.
- **Relative or absolute:** Conceptual — about the **currency risk-premium /
  convenience wedge** in the UIP relation (a relative, cross-currency object), and
  what models must deliver to reconcile both puzzles.
- **Magnitude:** cov(E_t ρ_{t+1}, r*_t − r_t) > 0 (UIP) coexists with stronger
  high-rate currencies in levels.
- **Does NOT identify:** a convenience level; it is a diagnostic of the joint
  restriction on the currency wedge.
- **Bridge to big_cy:** motivates a currency-specific convenience/safety wedge in
  the UIP block — exactly the object big_cy puts in the utility function.

### E2. Jiang, Krishnamurthy & Lustig (2023, *Review of Economic Studies*) — `jiangkrishnamurthylustig2023dollar`
**Proximity: 4.**
- **What it measures:** **Dollar safety** and the global financial cycle —
  extends the Treasury-basis/convenience-yield framework to dollar-asset safety
  premia and global risk.
- **Relative or absolute:** Uses the **RELATIVE** dollar/Treasury basis as the
  convenience proxy but embeds it in a model of the dollar's special role.
- **Magnitude:** Dollar convenience co-moves with the global financial cycle;
  drives the dollar's safe-haven appreciation in risk-off states.
- **Does NOT identify:** an absolute level; relies on the relative basis.
- **Bridge to big_cy:** structural use of the relative convenience measure to
  generate the dollar's negative beta — the result big_cy wants to obtain from a
  *large absolute* convenience instead.

### E3. Coppola, Krishnamurthy & Xu (2023, NBER w30984) — `coppolakrishnamurthyxu2023liquidity`
**Proximity: 3.**
- **What it measures:** A **liquidity-based theory of currency dominance**:
  endogenous search frictions make firms denominate debt in the most-liquid
  (dollar) unit; convenience yield (KVJ) makes it profitable to issue safe dollar claims.
- **Relative or absolute:** Theoretical; treats convenience yield as a level
  benefit driving issuance, but does not measure it. Provides the
  *positive-feedback* microfoundation for dollar dominance.
- **Magnitude:** Qualitative; rationalizes dollar's outsized debt share.
- **Does NOT identify:** any empirical level.
- **Bridge to big_cy:** supply-side microfoundation of dollar convenience and
  home/dollar bias — complements big_cy's demand-side convenience mechanism.

---

## Coverage note
20 entries spanning: cross-currency relative measures (A1–A6), domestic relative
measures (B1–B5), absolute-level measures (C1–C4), supply/demand & term-structure
(D1–D3), and currency-wedge/theory context (E1–E3). Anchor cross-reference:
DiTella et al. (C3) is the Strand-1 absolute-level anchor reproduced here for the
relative-vs-absolute bridge.
