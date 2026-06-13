# big_cy — Literature Synthesis

*Deliverable for T1a. Ties together the four strands — (1) the zero-beta rate and
the absolute level of convenience; (2) international models with safe-asset /
convenience demand; (3) mechanisms and microfoundations; (4) relative-vs-absolute
measurement. Companion to `frontier_map.md` (taxonomy) and
`annotated_bibliography.md` (entries). This is a neutral map of the field and its
open questions, not an argument for any one mechanism.*

---

## 1. The organizing distinction: a gap is not a level

Almost every empirical convenience-yield number in circulation is a **relative**
measure — the *difference* between two near-substitute safe assets. The
cross-currency basis (`dutepperverdelhan2018cip`) differences the dollar money rate
against an FX-swapped foreign rate (~25 bp). The US Treasury Premium
(`duimschreger2018treasury`) and the Treasury basis (`jiangkrishnamurthylustig2021foreign`)
difference US Treasuries against currency-hedged foreign sovereigns (~21 bp
pre-GFC). The Aaa–Treasury spread (`kvj2012aggregate`) differences Treasuries
against high-grade corporates (~73 bp). On/off-the-run (`krishnamurthy2002bond`),
Treasury–Refcorp (`longstaff2004flight`), and TIPS–Treasury
(`fleckensteinlongstafflustig2014tips`) difference Treasuries against still-closer
substitutes.

A difference cancels whatever convenience is **common to both legs**. If two safe
assets each carry a large convenience benefit, the spread between them can be
small, zero, or negative even when the *level* on each is large. This is not a
hypothetical: the US Treasury Premium is roughly **+21 bp before the GFC and −8 bp
after** (`duimschreger2018treasury`) — not because Treasury convenience vanished,
but because foreign safe sovereigns came to carry comparable convenience, so the
*gap* shrank and flipped. Diamond–Van Tassel (`diamondvantassel2024riskfree`)
confirm the point directly with absolute, option-implied, per-currency measures:
convenience *levels* rise with each country's interest-rate level and are
substantial everywhere, while the US-minus-foreign *difference* does not grow in
crises even as the levels do. **Relative and absolute convenience move
differently**, and a small relative number is not evidence against a large
absolute one.

Only two families of measure recover an absolute level, and both price the safe
asset against a **convenience-free benchmark built from risky assets**:

- **Option box spreads** (`vanbinsbergendiamondgrotteria2022riskfree`;
  `diamondvantassel2024riskfree`) extract a synthetic riskless rate from index
  options via put-call parity. The Treasury-vs-box gap is an absolute level —
  and it is **modest: ~40 bp** (≤ ~3Y maturity, ~4× in crises).
- **The zero-beta equity rate** (`ditella2025zerobeta`) extracts an
  intertemporal price from an equity portfolio orthogonal to the SDF. The
  zero-beta-vs-safe-rate gap is an absolute wedge — and it is **large: ~2–4%**
  (their real zero-beta rate averages ~8.3%, a ~7.6%/yr spread over the real
  T-bill).

These two absolute estimates **disagree by roughly an order of magnitude**. That
disagreement — not the familiar relative-vs-absolute point — is the live empirical
tension at the center of this literature, and the one big_cy must engage head-on
(§5).

---

## 2. Strand 1: what the zero-beta rate is, and why it is contested

Di Tella, Hébert, Kurlat and Wang (`ditella2025zerobeta`) revive Black's zero-beta
CAPM (`black1972capital`; `blackjensenscholes1972capm`) as a *time-varying*
object. Their construction: posit a linear-factor SDF (5 Fama-French factors + a
Treasury-bond factor + a default factor), build a minimum-variance equity
portfolio orthogonal to it (Ledoit-Wolf shrinkage, `ledoitwolf2017nonlinear`),
and project its return on macro instruments (T-bill, inflation, term spread, the
Gilchrist-Zakrajšek excess bond premium `gilchrist2012credit`, unemployment), all
by GMM. The resulting zero-beta rate (i) **satisfies the aggregate consumption
Euler equation** where the T-bill does not; (ii) is roughly equal to the average
market return (so the equity premium *relative to the zero-beta rate* is ≈ 0); and
(iii) is volatile enough to generate price-dividend-ratio volatility without
time-varying risk premia. They reinterpret the equity-premium puzzle
(`mehra1985equity`) and the risk-free-rate puzzle (`weil1989equity`) as a single
**"safe-rate puzzle"**: it is the *bonds* that are anomalous — too low, too
stable, disconnected from consumption — because safe bonds yield a large
non-pecuniary (convenience) return.

What is **settled**: a wedge exists. The security-market line is flat with a high
intercept (`blackjensenscholes1972capm`), the high zero-beta return survives rich
factor controls (`lopezlira2020common`), and Treasuries are demonstrably special
relative to near-substitutes (`kvj2012aggregate`; `vanbinsbergendiamondgrotteria2022riskfree`).
The wedge is time-varying and countercyclical (`gilchrist2012credit` as
instrument; vBDG's ~4× crisis widening).

What is **contested** (the authors are scrupulous that their *facts* survive but
the *level interpretation* does not): (a) the ~15× gap between the option-based
~40 bp and the equity-based ~2–4% — different "risk-free" notions, or an omitted
high-priced factor inflating the zero-beta level (DiTella et al.'s own stated
caveat)?; (b) whether the zero-beta rate is the intertemporal price, the SDF of
leverage-constrained investors (`frazzini2014betting`), or an incomplete-markets
multiple-SDF artifact; (c) no clean safety-vs-liquidity decomposition of the
*absolute* wedge exists (only of the relative spread, `kvj2012aggregate`,
`mota2021corporate`); (d) the source of the safe-rate puzzle itself (liquidity /
incomplete markets — `ditella2024aggregateeuler`, `kocherlakota2025public` — vs.
bond-stock SDF segmentation). The companion theory (`ditella2024aggregateeuler`)
and the fiscal "so-what" papers (`reis2021constraint`, `kocherlakota2025public`,
`angeletoscollarddellas2023public`) take the large wedge as input and show it
implies seigniorage and fiscal capacity "much higher than traditionally thought."

---

## 3. Strands 1↔3: what such a wedge would be *made of*

Strand 1 measures a wedge but is deliberately agnostic about its composition.
Strand 3 supplies the menu of economic sources — and the methodological choice
big_cy faces.

**The reduced-form route (bonds-in-utility / money-in-utility).** The Sidrauski
MIU lineage (`sidrauski1967rational`) drops a convenience term `v(B/P)` into
utility; its marginal value *is* the convenience yield. This is the device KL2024
uses, and it recurs across macro (`fisher2015structural` shows the workhorse DSGE
"risk-premium" shock is observationally a safe-asset demand shock;
`mianstraubsufi2025goldilocks`, `michaillatsaez2022economical`). It buys
tractability, a single Euler-equation wedge, an easily *time-varying* and
*crankable* level, and an estimable demand curve. It costs interpretability: the
term cannot separate liquidity from safety, is silent on why premia collapse in
crises, and is not policy-invariant (Lucas critique).

**The deep-microfoundation route.** The premium can instead be derived from
primitives: money/search frictions (`lagoswright2005unified`, `lagos2010asset`,
`vayanosweill2008search` — which microfounds the on-the-run fact
`krishnamurthy2002bond`); pledgeability/collateral and public-debt-as-liquidity
(`woodford1990public`, `holmstromtirole1998private`, `aiyagarimcgrattan1998optimum`,
`angeletoscollarddellas2023public`); information-insensitivity/safety
(`gortonpennacchi1990financial`, `danggortonholmstrom2017banks`); and
regulatory/private-money-creation externalities (`stein2012monetary`,
`ghs2015comparative`, `kvj2015impact`). These deliver structural interpretation,
state-contingency, and policy comparative statics, but are model-specific and
typically price only one leg.

An early structural antecedent of the safe-rate
puzzle reading is Bansal-Coleman (`bansalcoleman1996monetary`), in which money
services lower safe-bond returns enough to address the equity-premium and
risk-free-rate puzzles; the institutional/**regulatory** source is represented by
**d'Avernas-Vandeweyer (`davernasvandeweyer2024treasury`)** (T-bill scarcity and
bank reserve/repo demand pricing short safe assets), and the *quantity* of safe
assets by **Gorton-Lewellen-Metrick (`gortonlewellenmetrick2012safe`)** (the
safe-asset share is stable at ~33% of total US assets).

**The empirical anchors that discipline either route** are
Krishnamurthy-Vissing-Jorgensen (`kvj2012aggregate`) — the convenience yield is
real, ~70 bp on average, splits into a **liquidity** and a **safety** leg, and
traces a **downward-sloping demand curve** in Treasury/GDP (falling supply raises
it) — and Nagel (`nagel2016liquidity`) — the liquidity leg moves **~one-for-one
with the short nominal rate** (its opportunity cost). KVJ pins the *quantity*
slope; Nagel pins the *price/opportunity-cost* slope; both imply strong
countercyclicality. The tension big_cy inherits: these are **relative** spreads,
so they bound the *near-substitute gap* (tens of bp), and may badly understate the
*absolute* level DiTella infers. The mechanism literature tells you what a large
absolute wedge would be made of and how it should move — even if its level is
identified elsewhere.

---

## 4. Strand 2: how the international models are built, and where they part ways

The open-economy literature splits into **three modeling families**, plus a
cross-cutting estimation approach (`maggiori2022handbook` is the field's own map).

1. **Bonds-in-utility / exogenous-convenience structural NK.** Kekre-Lenel
   (`kekrelenel2024flight`) is the template: a two-country New Keynesian model
   with a time-varying BiU **safety shock** on dollar bonds *and* cross-country
   **heterogeneity in risk-bearing capacity** (the US is more risk-tolerant). A
   flight to safety appreciates the dollar, makes it a negative-beta hedge, and
   leads the US to hold a levered long-capital/short-dollar portfolio — the
   "world insurer." Choi-Kirpalani-Perez (`choikirpalaniperez2022marketpower`)
   add monopoly safe-asset supply (markup ~2/3 of the convenience yield). This is
   big_cy's home family.

2. **Intermediary / financier portfolio-balance.** Gabaix-Maggiori
   (`gabaixmaggiori2015liquidity`), Maggiori (`maggiori2017reservecurrencies`),
   Gourinchas-Ray-Vayanos (`gourinchasrayvayanos2022habitat`),
   Greenwood-Hanson-Stein-Sunderam (`greenwoodhansonsteinsunderam2023quantity`),
   Kekre-Lenel's GFC model (`kekrelenel2025gfc`), Bianchi-Bigio-Engel
   (`bianchibigioengel2023scrambling`). A constrained intermediary's limited
   risk-bearing prices currency/term/liquidity risk; gross flows move the
   exchange rate; convenience is implicit in a risk premium (or, in BBE, an
   endogenous settlement-liquidity yield).

3. **Safe-asset supply / shortage / dominance.** Caballero-Farhi-Gourinchas
   (`caballerofarhigourinchas2008imbalances`, `caballerofarhigourinchas2021policywars`),
   Caballero-Farhi (`caballerofarhi2018safetytrap`), Farhi-Maggiori
   (`farhimaggiori2017ims`), He-Krishnamurthy-Milbradt (`hekrishnamurthymilbradt2019safeasset`),
   Coppola-Krishnamurthy-Xu (`coppolakrishnamurthyxu2023liquidity`), Gopinath-Stein
   (`gopinathstein2021banking`). Emphasis on who *supplies* safety, scarcity, ZLB
   safety traps, monopoly rents, and *why the dollar* (coordination, fiscal
   capacity, invoicing complementarities, self-reinforcing liquidity).

Cross-cutting: a **demand-system estimation** approach (`koijenyogo2019demandsystem`,
`koijenyogo2020globaldemand`, `jiangrichmondzhang2022portfolio`) *measures* the
absolute US convenience yield (~2.15% on long-term debt, 1.70% on equity) and the
drivers of NFA composition — the empirical bridge to big_cy's quantitative aims
and the closest cousin to the in-house `rier` backbone.

Two further pieces sit especially close to big_cy. **Valchev (`valchev2020bond`)**
is the closest existing *structural* model to the project's mechanism: a
time-varying **bond convenience yield** in the households' Euler equation drives
the exchange rate and produces the UIP/forward-premium reversal — convenience-in-FX
without a standard currency risk premium. **Engel-Wu (`englewu2023liquidity`)**
supplies the matching empirics (the dollar tracks the liquidity/convenience yield
on dollar assets). And the *facts* big_cy must explain are the
exorbitant-privilege accounting of **Gourinchas-Rey (`gourinchasrey2007banker`)** —
the US as a "world venture capitalist," long high-yield equity and short low-yield
debt — while **Bocola-Lorenzoni (`bocolalorenzoni2020dollarization`)** microfounds
currency-composition (home/dollar) bias from a self-fulfilling savers' insurance
motive, and **Caballero-Krishnamurthy (`caballerokrishnamurthy2009fragility`)** is
the direct antecedent of the safe-asset-shortage lineage.

**Where they agree.** The dollar appreciates in bad times; the US earns a rent
("exorbitant privilege"); gross portfolio composition, not just net positions,
pins the exchange rate and premia.

**Where they disagree — the tension big_cy targets.** Is the engine of privilege
**convenience seigniorage** or **risk-bearing capacity**? JKL
(`jiangkrishnamurthylustig2024dollar`) and the convenience family argue the
privilege and the resolution of Maggiori's **reserve-currency paradox**
(`maggiori2017reservecurrencies` — a pure insurance story counterfactually
predicts the dollar *depreciates* in crises) come from seigniorage on convenience
demand and do **not** require the US to be more risk-tolerant. KL2024 uses *both*
channels, and its assumption that the US has **greater risk-bearing capacity** is
exactly the counterfactual lever (US risk tolerance is not empirically larger)
that big_cy proposes to replace with a large absolute convenience yield.
Devereux-Engel-Wu (collateral), Bianchi-Bigio-Engel (settlement), and
He-Krishnamurthy-Milbradt (coordination) offer yet other resolutions of the same
paradox — big_cy must position against this whole cluster.

---

## 5. The synthesis: the four strands meet at one unoccupied cell

Putting the strands together (see the plane in `frontier_map.md`):

- **Strand 4** establishes that the convenience numbers used to *discipline* the
  international models (KL2024, JKL) are **relative** — the Treasury/CIP basis —
  and so are, by construction, near-substitute gaps of tens of basis points.
- **Strand 1** establishes that the **absolute** level may be far larger (~2–4%),
  with the strong caveat that the equity-based number is contested and an
  order of magnitude above the option-based ~40 bp.
- **Strand 3** establishes that a wedge of either size can be represented as a
  bonds-in-utility term, and tells us how to *discipline* that term's slope (KVJ)
  and cyclicality (Nagel) against micro evidence — while warning that BiU is not
  policy-invariant.
- **Strand 2** establishes that the leading structural model (KL2024) currently
  leans on **risk-bearing heterogeneity** — a counterfactual assumption — to
  pin down US portfolio behavior, and that a *relative*-calibrated convenience
  shock sits alongside it.

The unoccupied cell is therefore precise: **no structural open-economy model has
asked whether an *absolute*, DiTella-magnitude convenience yield — disciplined as
a bonds-in-utility level rather than a relative basis, and let to operate
symmetrically across currencies (US households value safe dollars; foreign
households value safe bonds in both their own currency and dollars) — can carry
the weight currently borne by risk-bearing heterogeneity, and thereby microfound
home-currency bias.** The "quantity of safe debt" channel (JKL2024, GHSS2023,
Coppola-Krishnamurthy-Xu, Jiang-Richmond-Zhang) and the "absolute level" channel
(DiTella, Koijen-Yogo) have also never been combined in one quantitative GE model
— and big_cy's "not enough dollar debt to satisfy US diversification demand"
reading lives exactly at their intersection.

---

## 6. Open questions and disagreements (for downstream tasks)

1. **The order-of-magnitude disagreement on the absolute level.** Option-box
   (~40 bp) vs. zero-beta-equity (~2–4%). Candidate reconciliations: horizon
   (box spreads ≤ 3Y financial rates vs. a long-horizon intertemporal price);
   the *object* priced (a tradable riskless payoff vs. the consumption-Euler
   margin); composition (the equity wedge may bundle convenience with broader
   missing-risk-compensation). Whichever big_cy adopts is a *choice of benchmark*,
   and the option-based number is the principal countervailing evidence a neutral
   framing must present. Diamond–Van Tassel (`diamondvantassel2024riskfree`) is
   the natural international disciplining anchor (absolute, per-currency).
   **See `vbdg_zerobeta_reconciliation.md`** for the full deep dive: DiTella (§7)
   explicitly reconcile the two via *leverage* — the box rate (~40 bp) is a
   leveraged intermediary's near-riskless funding rate, the zero-beta rate (~2–4%)
   the unlevered required return on riskless equity, bridged by the bank leverage
   ratio (~8×) and ordered on Diamond's (2020) ladder of multiple risk-free rates.
   They are different objects, not competing estimates of one number; big_cy's
   choice of rung is a calibration decision about *which marginal investor* prices
   the model's safe bonds.

2. **Convenience vs. risk-bearing as the engine of exorbitant privilege.**
   Unresolved across `jiangkrishnamurthylustig2024dollar`,
   `maggiori2017reservecurrencies`, `kekrelenel2024flight`, `devereuxenglewu2023collateral`,
   `bianchibigioengel2023scrambling`. big_cy's bet is on convenience; the
   literature has not adjudicated.

3. **Steady-state level vs. time-varying shock.** KL and Strand-1 work emphasize
   *time-varying* safety demand (the cyclical wedge), whereas the home-bias /
   carry / NFA-composition facts are about the *average* (steady-state) level.
   Whether the average level or its variation does the work is conceptually open
   and is the explicit subject of the in-house theory task T1c.

4. **Reduced-form vs. microfounded convenience.** Is a BiU term an adequate
   stand-in, or does the policy-invariance/Lucas-critique problem bite for the
   counterfactuals big_cy wants (e.g., changing the supply of dollar debt)? The
   bridge papers (`kvj2015impact`, `ghs2015comparative`, `stein2012monetary`,
   `fisher2015structural`) are the templates for a *disciplined* BiU.

5. **Safety vs. liquidity composition of the absolute wedge.** Cleanly decomposed
   only for the *relative* spread (`kvj2012aggregate`, `mota2021corporate`,
   `nagel2016liquidity`); the absolute wedge's composition is unknown — and it
   matters for whether the convenience is currency-specific (dollar) or
   safe-asset-generic.

6. **Imperfect substitutability and the "quantity" reinterpretation.** Nenova
   (`nenova2025global`) finds US Treasuries are substituted *globally* while Bunds
   are substituted only *regionally*, and substitutability deteriorates in stress
   — consistent with big_cy's "insufficient dollar debt" reading, but the mapping
   from substitution elasticities to a convenience *level* is not established.

### Citations to verify before manuscript use
A handful of entries were assembled from citing-paper reference lists rather than
read in full, or have ambiguous published-vs-WP status. Verify before use:
`andersonduschlusche2021arbitrage` (NBER WP 28658 vs. a JF version — confirm
venue/volume/pages); `diamondvantassel2024riskfree` (JF forthcoming/2025, cited as
2024); `kocherlakota2025public`, `mota2021corporate`, `mota2025financially`,
`lopezlira2020common` (working-paper venue/year). Flagged in `references.bib`
header. (`diamond2020safety` is now resolved — published JF 75(6), 2020 — per the
VBDG/zero-beta deep dive, `vbdg_zerobeta_reconciliation.md`.)

---

## 7. Coverage note

~88 unique references across the four strands (79 in the initial sweep + 9 added
after the `librarian-critic` review — Valchev; Gourinchas-Rey; Bocola-Lorenzoni;
Caballero-Krishnamurthy; Engel-Wu "Liquidity and Exchange Rates"; d'Avernas-Vandeweyer;
Gorton-Lewellen-Metrick; Du-Hébert-Li; Bansal-Coleman), spanning top-5 generals (AER,
ECMA, JPE, QJE, REStud), top finance journals (JF, JFE, RFS), field journals
(JIE, JME, JMCB, JEP, JBES, J. Business), NBER/BIS working papers, and two
Handbook chapters. Foundational theory (Black 1972; Sidrauski 1967; Diamond-Dybvig
1983; Woodford 1990; Holmström-Tirole 1998; Gorton-Pennacchi 1990; Lagos-Wright
2005) is included alongside the 2018–2025 frontier. The collection is deliberately
**broad and neutral** — supply-side and demand-side, exogenous and endogenous,
relative and absolute, reduced-form and microfounded — so that big_cy's eventual
positioning rests on a faithful map rather than a curated one. The most important
deliberate inclusion against interest is the **option-based ~40 bp absolute
estimate** (`vanbinsbergendiamondgrotteria2022riskfree`,
`diamondvantassel2024riskfree`), which is the strongest evidence *against* the
large-absolute-convenience premise and must travel with the DiTella number
wherever it is cited.
