# big_cy — Master Annotated Bibliography

*Deliverable for T1a. Self-contained master organized by the four strands (~88 entries);
canonical citekeys match `references.bib`. Each entry gives a compact annotation
(question / method / finding) plus a **Rel/Abs** tag (does it work with a relative
gap or an absolute level?) and a **Prox** score (1–5; 5 = directly competes /
keystone). The four `strands/strand_*/annotated_bibliography.md` files hold the
longest per-entry treatments; this master consolidates and de-duplicates them.
See `frontier_map.md` for the cross-strand taxonomy and `synthesis.md` for the
essay.*

**Reading guide.** Strand 1 = the absolute level (the zero-beta rate and its
rivals). Strand 2 = international GE models with safe-asset/convenience demand.
Strand 3 = the economic sources of convenience and the reduced-form-vs-microfounded
split. Strand 4 = what each empirical measure does and does not identify. Papers
that span strands (KVJ 2012; DiTella; JKL; vBDG; Krishnamurthy 2002; FLL;
Coppola-K-X; GHSS; Angeletos et al.) are placed in their primary strand with a
cross-pointer.

---

## STRAND 1 — The zero-beta rate and the ABSOLUTE level of convenience yields

### Keystone

**`ditella2025zerobeta`** — Di Tella, Hébert, Kurlat & Wang, "The Zero-Beta
Interest Rate," *JPE* (2025, forthcoming; NBER WP 32464). **Rel/Abs: Absolute.
Prox 5.** Constructs a time-varying zero-beta rate — the expected return on a
minimum-variance equity portfolio orthogonal to a linear-factor SDF (5 FF + bond +
default factors), estimated by GMM with macro instruments. The real zero-beta rate
averages ~8.3% (~7.6%/yr over the real T-bill), fits the aggregate consumption
Euler equation (where the T-bill fails), makes the equity premium relative to it
≈ 0, and generates PD-ratio volatility without time-varying risk premia.
Reinterprets the equity-premium and risk-free-rate puzzles as one **"safe-rate
puzzle"**: safe bonds yield a large non-pecuniary (convenience) return. The source
of big_cy's large-absolute-CY premise; authors flag the *level* estimate as less
robust than the Euler test. (Project anchor paper.)

### Other estimates of the absolute level / true risk-free rate

**`vanbinsbergendiamondgrotteria2022riskfree`** — van Binsbergen, Diamond &
Grotteria, "Risk-Free Interest Rates," *JFE* 143(1):1–29 (2022). **Rel/Abs:
Absolute. Prox 5.** Extracts a convenience-free risk-free rate from S&P 500 option
**box spreads** (put-call parity), model-free, to ~3Y. Implied **Treasury
convenience ≈ 40 bp** on average, larger at short maturities, ~4× in the crisis.
The "small-number" absolute counterpoint to DiTella's "big-number" wedge — the
~15× gap between them is the strand's central open question. (Also Strand 4 / C1.)

**`reis2021constraint`** — Reis, "The Constraint on Public Debt when r<g but g<m,"
BIS WP 939 (2021). **Rel/Abs: Absolute-ish. Prox 4.** GE model endogenizing a
safety/liquidity "bubble premium" on public debt; frames the absolute wedge as
r-vs-m (the marginal product of capital), several %, and shows it generates large
"debt revenue"/seigniorage and fiscal space. Connects the absolute CY to fiscal
capacity. (Cross-ref Strand 4.)

**`mota2021corporate`** — Mota, "The Corporate Supply of (Quasi-)Safe Assets," WP
(2021); and **`mota2025financially`** — Mota & Siani, "Financially Sophisticated
Firms," WP (2025). **Rel/Abs: Relative (safety leg). Prox 3.** Build a near-riskless
synthetic bond (corporate + CDS) so the residual vs. Treasuries isolates a
**safety** premium; it explains a large share (R²≈60%) of spread variation, and
firms manufacture quasi-safe assets when it is high. Speaks to the safety-vs-liquidity
decomposition; a lower bound on the absolute safety premium.

### The zero-beta CAPM lineage

**`black1972capital`** — Black, "Capital Market Equilibrium with Restricted
Borrowing," *J. Business* 45(3) (1972). **Abs (theory). Prox 4.** Derives the
two-factor zero-beta CAPM: with no riskless borrowing, expected returns are linear
in beta with intercept the (high) zero-beta return — the conceptual origin of a
zero-beta rate above safe yields.

**`blackjensenscholes1972capm`** — Black, Jensen & Scholes, "The CAPM: Some
Empirical Tests," in *Studies in the Theory of Capital Markets* (1972). **Abs (avg
level). Prox 3.** First documentation of a flat security-market line with a high
positive intercept — the empirical zero-beta return above the T-bill that DiTella
makes time-varying.

**`frazzini2014betting`** — Frazzini & Pedersen, "Betting Against Beta," *JFE*
111(1) (2014). **Abs. Prox 3.** Leverage-constrained investors bid up high-beta
assets, flattening the SML and raising the zero-beta rate; the BAB factor is large
globally. DiTella's leading *alternative* interpretation of a high zero-beta rate
(a constraint multiplier, not the consumption intertemporal price).

**`lopezlira2020common`** — Lopez-Lira & Roussanov, "Do Common Factors Really
Explain the Cross-Section of Stock Returns?", WP (2020). **Abs (indirect). Prox 2.**
A diversified portfolio orthogonal to common factors still earns a high return;
cited to show the high zero-beta return is hard to factor away.

### Theory of why safe bonds are special / the bond–stock SDF wedge

**`diamond2020safety`** — Diamond, "Safety Transformation and the Structure of the
Financial System," WP (2020). **Abs (theory). Prox 4.** GE in which household
safe-asset demand makes the risk-free rate pricing *household* claims lie strictly
below the rate implied by pricing *risky* assets — multiple coexisting "risk-free
rates." The theoretical hinge for why the equity-based wedge and the AAA–Treasury
spread are both "convenience" yet need not be equal.

**`angeletoscollarddellas2023public`** — Angeletos, Collard & Dellas, "Public Debt
as Private Liquidity: Optimal Policy," *JPE* 131(11) (2023). **Abs (theory). Prox 3.**
Convenience/liquidity premium on public debt justifies issuing more than the
tax-smoothing benchmark; the safe-rate-vs-social-rate wedge reshapes optimal debt.
(Also Strand 3 / E.)

**`kocherlakota2025public`** — Kocherlakota, "Public Debt Bubbles and the Welfare
Effects of Higher Interest Rates," WP (2025). **Abs (uses wedge). Prox 3.** Takes
the DiTella configuration as input; raising rates can improve welfare when the safe
rate is below the zero-beta rate; seigniorage is much larger than usually measured.

**`ditella2024aggregateeuler`** — Di Tella, Hébert, Kurlat & Wang, "Aggregate Euler
Equations with Liquid and Illiquid Assets," WP (2024). **Abs (microfoundation).
Prox 4.** Kaplan-Violante-style HA model with liquid/illiquid accounts giving
conditions under which the aggregate Euler equation holds for the illiquid (equity)
return ≈ zero-beta rate but fails for liquid safe bonds. The companion theory
behind the keystone empirics; a template for microfounding a large absolute CY.

### Foundational puzzles being reinterpreted

**`mehra1985equity`** — Mehra & Prescott, "The Equity Premium: A Puzzle," *JME*
15(2) (1985). **Prox 2.** The ~6% equity premium over safe bonds needs implausible
risk aversion; DiTella reinterprets it as a safe-rate puzzle (the wedge to absorb).

**`weil1989equity`** — Weil, "The Equity Premium Puzzle and the Risk-Free Rate
Puzzle," *JME* 24(3) (1989). **Prox 2.** Matching the low risk-free rate needs
implausible preferences — the precursor framing of DiTella's safe-rate puzzle.

**`campbellshiller1988dividend`** / **`campbell1991variance`** — Campbell & Shiller
(*RFS* 1988); Campbell (*EJ* 1991). **Prox 2.** Present-value/variance
decompositions showing valuation-ratio volatility must reflect discount-rate
variation; DiTella repurpose the tool to show a volatile zero-beta rate generates it
without time-varying risk premia.

### Methods

**`ledoitwolf2017nonlinear`** (*RFS* 2017, **Prox 2**) — nonlinear covariance
shrinkage for the minimum-variance zero-beta portfolio. **`olea2013robust`**
(*JBES* 2013, **Prox 1**) — weak-instrument-robust inference for the GMM
instruments. **`gilchrist2012credit`** (*AER* 2012, **Prox 2**) — the excess bond
premium, a key DiTella instrument linking the wedge to credit conditions /
cyclicality.

---

## STRAND 2 — International models with safe-asset / convenience demand

### Directly related: structural open-economy convenience models

**`kekrelenel2024flight`** — Kekre & Lenel, "The Flight to Safety and International
Risk Sharing," *AER* 114(6):1650–1691 (2024). **Demand: exogenous BiU safety shock.
Open 2-ctry NK. Rel (Treasury basis). Prox 5.** big_cy's structural template.
Time-varying nonpecuniary convenience on safe **dollar** bonds + cross-country
**risk-bearing heterogeneity** (US more risk-tolerant). A flight to safety
appreciates the dollar, makes it a negative-beta hedge, and leads the US to hold a
levered long-capital/short-dollar portfolio. Calibrated to the Du-Im-Schreger
basis. Resolves Maggiori's reserve-currency paradox. The risk-tolerance assumption
is the counterfactual lever big_cy aims to replace. (Project anchor paper.)

**`jiangkrishnamurthylustig2024dollar`** — Jiang, Krishnamurthy & Lustig, "Dollar
Safety and the Global Financial Cycle," *ReStud* 91(5):2878–2915 (2024). **Demand:
exogenous, quantity-dependent + endogenous supply. Open US/RoW. Rel→abs framing.
Prox 5.** The GFC as a dollar cycle: foreign convenience demand for safe dollar
bonds decreases in their quantity; US monetary tightening cuts dollar-bond supply,
raises the convenience yield, appreciates the dollar. Argues seigniorage on
convenience — **not** risk-bearing — drives privilege and resolves the paradox.
The closest competitor/cousin to big_cy's mechanism.

**`jiangkrishnamurthylustig2021foreign`** — JKL, "Foreign Safe Asset Demand and the
Dollar Exchange Rate," *JF* 76(3):1049–1089 (2021). **Demand: exogenous (basis).
Rel. Prox 5.** Convenience-yield theory of the exchange rate: the dollar equals the
PV of rate differentials minus currency risk premia minus the foreign convenience
yield, measured by the **Treasury basis**. Foreigners earn ~2% convenience, ~90%
attributable to dollar denomination; convenience explains 16–41% of quarterly
dollar variation. big_cy's relative-measure point of departure. (Also Strand 4 / A4.)

**`maggiori2017reservecurrencies`** — Maggiori, "Financial Intermediation,
International Risk Sharing, and Reserve Currencies," *AER* 107(10) (2017). **Demand:
endogenous (intermediary asymmetry). Open 2-ctry. Prox 4.** Countries differ in
financial development; the US bears risk and insures RoW, its currency appreciating
in crises. States the **reserve-currency paradox** that a pure insurance story
counterfactually predicts dollar depreciation in crises — the puzzle the
convenience family claims to resolve.

**`devereuxenglewu2023collateral`** — Devereux, Engel & Wu, "Collateral Advantage,"
NBER WP 31164 (2023). **Demand: endogenous (collateral). Open 2-ctry NK.
Abs-leaning. Prox 4.** US government debt is superior bank collateral; in stress its
collateral value rises, appreciating the dollar despite falling US absorption —
resolving the paradox via a collateral channel.

**`bianchibigioengel2023scrambling`** — Bianchi, Bigio & Engel, "Scrambling for
Dollars," NBER WP 29457 (2023). **Demand: endogenous (settlement liquidity). Open
banking. Abs (endogenous level). Prox 4.** Interbank settlement frictions generate
an endogenous dollar convenience yield; when dollar funding risk rises, banks
demand dollars and the dollar appreciates.

### Intermediary / financier portfolio-balance

**`gabaixmaggiori2015liquidity`** — Gabaix & Maggiori, "International Liquidity and
Exchange Rate Dynamics," *QJE* 130(3) (2015). **Demand: exogenous flow. Open 2-ctry.
Rel. Prox 3.** Canonical financier/portfolio-balance FX model: limited-risk-bearing
global financiers intermediate currency imbalances; gross flows move the exchange
rate. Nests UIP and autarky; generates disconnect and carry. The methodological
ancestor of the KL lineage.

**`gourinchasrayvayanos2022habitat`** — Gourinchas, Ray & Vayanos, "A
Preferred-Habitat Model of Term Premia, Exchange Rates, and Monetary Policy
Spillovers," NBER WP 29875 (2022). **Demand: exogenous habitat. Open 2-ctry cts.
Rel. Prox 3.** Currency traders + bond clienteles + limited arbitrageurs who are
marginal everywhere, tightly linking UIP and term-premium violations; QE spillovers.

**`greenwoodhansonsteinsunderam2023quantity`** — Greenwood, Hanson, Stein &
Sunderam, "A Quantity-Driven Theory of Term Premia and Exchange Rates," *QJE*
138(4):2327–2389 (2023). **Demand: exogenous bond supply. Open 2-ctry. Rel
(cross-ccy qty). Prox 3.** Relative supplies of long bonds in two currencies drive
term premia and the exchange rate; matches QE→FX and post-2008 CIP deviations.
Connects "quantity of safe assets" to FX. (Cross-ref Strand 3/4.)

**`kekrelenel2025gfc`** — Kekre & Lenel, "A Model of U.S. Monetary Policy and the
Global Financial Cycle," NBER WP (Nov 2025). **Demand: exo habitat + endo
intermediary wealth. Open multi NK. Prox 4.** Intermediaries are short the dollar,
so US tightening appreciates the dollar, erodes their net worth, and raises the
global price of risk. Anticipates that lower future foreign dollar demand would
dampen US spillovers — engaging the absolute-level/quantity theme.

**`kekrelenel2022redistribution`** — Kekre & Lenel, "Monetary Policy,
Redistribution, and Risk Premia," *ECMA* 90(5) (2022). **Closed HANK. Prox 2.**
Heterogeneity in the marginal propensity to take risk; expansionary MP redistributes
to risk-tolerant households and lowers the risk premium. The methodological engine
(Epstein-Zin, heterogeneous risk-bearing) behind KL2024/2025.

**`kekrelenel2023swaplines`** — Kekre & Lenel, "The High Frequency Effects of Dollar
Swap Lines," WP (Dec 2023). **Empirical. Prox 2.** Swap-line expansions reduce
liquidity premia, compress CIP deviations, depreciate the dollar — the identified
moments calibrating KL2024's safe-dollar-supply experiment.

### Reserve-currency / dominance models turning on safe-asset supply

**`farhimaggiori2017ims`** — Farhi & Maggiori, "A Model of the International
Monetary System," *QJE* 133(1) (2018). **Demand: exogenous reserve demand +
strategic supply. Open multi. Abs level *as an endogenous monopoly rent* — distinct
from Koijen-Yogo's estimated level. Prox 4.** Monopoly/oligopoly
issuers of reserve assets earn an endogenous safety premium; self-fulfilling
confidence crises; rationalizes the Triffin dilemma and multipolar instability.
The supply-side counterpart to big_cy's demand-side story.

**`choikirpalaniperez2022marketpower`** — Choi, Kirpalani & Perez, "The
Macroeconomic Implications of US Market Power in Safe Assets," WP (2022). **Demand:
exogenous BiU + monopoly supply. Open 2-ctry. Abs (level + markup). Prox 4.** US
public debt yields a nonpecuniary convenience; the US is a monopoly provider facing
downward-sloping demand. Rejects price-taking; monopoly → underprovision (~half)
and a **markup ~2/3** of the convenience yield. Directly about the *level* and
safe-asset scarcity — close to "not enough dollar debt."

**`coppolakrishnamurthyxu2023liquidity`** — Coppola, Krishnamurthy & Xu, "Liquidity,
Debt Denomination, and Currency Dominance," NBER WP 30984 (2023). **Demand:
endogenous (search/settlement). Open multi-issuer. Abs (liquidity premium). Prox 4.**
Firms denominate debt in the most-liquid unit (largest safe short-term float),
yielding a self-reinforcing single dominant currency. Microfounds dollar-debt
dominance. (Also Strand 4 / E3.)

**`gopinathstein2021banking`** — Gopinath & Stein, "Banking, Trade, and the Making
of a Dominant Currency," *QJE* 136(2) (2021). **Demand: endogenous
(invoicing↔deposit↔bank). Open multi. Abs-leaning. Prox 3.** Dollar invoicing makes
dollar deposits safest for consumption, so banks fund and firms borrow in dollars —
a complementarity producing dominance and exorbitant privilege.

### Safe-asset shortage / safety trap

**`caballerofarhigourinchas2008imbalances`** — CFG, "An Equilibrium Model of 'Global
Imbalances' and Low Interest Rates," *AER* 98(1) (2008). **Demand: exogenous
asset-supply capacity. Open multi. Quantity/shortage. Prox 2.** Regions differ in
capacity to *supply* (safe) assets; the US supplies them, running deficits with
secularly low rates. The lineage of the safe-asset-shortage view.

**`caballerofarhi2018safetytrap`** — Caballero & Farhi, "The Safety Trap," *ReStud*
85(1) (2018). **Demand: exogenous safety. Closed/global ZLB. Abs/shortage. Prox 2.**
Excess safe-asset demand at the ZLB → a persistent demand-driven "safety trap";
QE-style safe-for-risky swaps and safe public debt are expansionary.

**`caballerofarhigourinchas2021policywars`** — CFG, "Global Imbalances and Policy
Wars at the ZLB," *ReStud* 88(6) (2021). **Open multi ZLB. Abs/shortage. Prox 3.**
Safety traps spill across borders; the reserve currency bears a disproportionate
share. KL2024's closest antecedent for the effects of safety shocks.

**`hekrishnamurthymilbradt2019safeasset`** — He, Krishnamurthy & Milbradt, "A Model
of Safe Asset Determination," *AER* 109(4) (2019). **Demand: endogenous
(coordination + fiscal). Open/global. Abs. Prox 3.** *Which* asset is "safe" is a
coordination outcome plus fiscal capacity; the US can be the global safe asset even
with comparable fundamentals. Microfounds why-the-dollar — the step big_cy's
exogenous-convenience shortcut skips.

### Demand-system / quantitative estimation

**`koijenyogo2020globaldemand`** — Koijen & Yogo, "Exchange Rates and Asset Prices in
a Global Demand System," NBER WP 27342 (2020). **Estimated latent demand. Open
36-ctry. Abs (estimated level). Prox 4.** A holdings-based demand system over
short/long debt and equity across 36 countries; macro variables explain 55–69% of
FX/yield/equity variation; estimates an **absolute** US convenience yield of
**2.15% (LT debt), 1.70% (equity)**. Highly relevant to the magnitude question.

**`koijenyogo2019demandsystem`** — Koijen & Yogo, "A Demand System Approach to Asset
Pricing," *JPE* 127(4) (2019). **Estimated. Closed (US eq.). Prox 1.** The
characteristics-based, IV-estimated method underpinning the global system and
Jiang-Richmond-Zhang.

**`jiangrichmondzhang2022portfolio`** — Jiang, Richmond & Zhang, "A Portfolio
Approach to Global Imbalances," WP (2022). **Estimated demand. Open multi.
Abs/quantity. Prox 4.** Decomposes US NFA and the privilege reversal; a demand shift
toward US *equities* drove the post-2010 reversal; documents a novel **quantity
dimension** (the US can issue far more long-term debt per unit yield rise).
Empirically closest to big_cy's NFA-composition target facts (cousin to `rier`).

### Cross-cutting reference / strategic dimension

**`maggiori2022handbook`** — Maggiori, "International Macroeconomics with Imperfect
Financial Markets," *Handbook of Int'l Economics* vol. 6, ch. 5 (2022). **Survey.
Prox 2.** The field's own organizing map: financier structure, risk-bearing
capacity, UIP-vs-CIP, carry, multi-asset extensions, the US-special branch.
Essential orientation.

**`claytonmaggiorischreger2024coercion`** — Clayton, Maggiori & Schreger, "A Theory
of Economic Coercion and Fragmentation," NBER WP 33309 (2024). **Strategic. Open
multi. Prox 1.** Dollar/safe-asset centrality as geoeconomic leverage; the
durability dimension of dollar dominance.

### Additions after librarian-critic review

**`valchev2020bond`** — Valchev, "Bond Convenience Yields and Exchange Rate
Dynamics," *AEJ: Macroeconomics* 12(2):124–166 (2020). **Demand: exogenous
(convenience-in-Euler). Open. Rel→abs. Prox 5.** The closest existing *structural*
model to big_cy's mechanism: a time-varying **bond convenience yield** enters the
households' bond Euler equation and drives exchange-rate dynamics, generating the
UIP/forward-premium reversal (short-run UIP failure, long-run reversal) without a
standard currency risk premium. A direct neighbor of `jiangkrishnamurthylustig2021foreign`
and `engel2016exchange`, and the natural prior art for a convenience-driven FX block.

**`gourinchasrey2007banker`** — Gourinchas & Rey, "From World Banker to World
Venture Capitalist: U.S. External Adjustment and the Exorbitant Privilege," in
Clarida (ed.), *G7 Current Account Imbalances*, U. Chicago Press, pp. 11–66 (2007).
**Empirical (NFA accounting). Open. Prox 4.** The foundational measurement of the
US external position: the US is long high-yield equity/FDI and short low-yield
debt ("world venture capitalist"), earning a return differential — the **empirical
statement of big_cy's target facts** (equity-biased external assets; exorbitant
privilege) and the antecedent of the in-house `rier` backbone.

**`bocolalorenzoni2020dollarization`** — Bocola & Lorenzoni, "Financial Crises,
Dollarization, and Lending of Last Resort in Open Economies," *AER* 110(8):2524–2557
(2020). **Demand: endogenous (savers' insurance motive). Open. Prox 4.** Liability
dollarization arises *self-fulfillingly*: because crises bring depreciations,
domestic savers demand a premium on local-currency savings, making local-currency
debt expensive and pushing borrowers into foreign currency. A microfoundation of
**currency-composition bias** directly relevant to big_cy's home-currency-bias
question (from the EM/borrower side).

**`caballerokrishnamurthy2009fragility`** — Caballero & Krishnamurthy, "Global
Imbalances and Financial Fragility," *AER P&P* 99(2):584–588 (2009). **Demand:
exogenous safe-asset demand. Open. Quantity/shortage. Prox 2.** The direct
antecedent (with `caballerofarhigourinchas2008imbalances`) of the safe-asset-shortage
lineage: foreign demand for US safe assets funds risky US intermediation, a source
of fragility. Background for the shortage view big_cy inherits.

**`englewu2023liquidity`** — Engel & Wu, "Liquidity and Exchange Rates: An Empirical
Investigation," *ReStud* 90(5):2395–2438 (2023). **Rel (liquidity/convenience
wedge). Open. Prox 4.** Empirically links the dollar to the *liquidity/convenience*
yield on dollar assets in an integrated UIP framework — the empirical companion to
the convenience-FX models (`valchev2020bond`, `jiangkrishnamurthylustig2021foreign`).
Distinct from the Devereux-Engel-Wu "Collateral Advantage" model
(`devereuxenglewu2023collateral`). (Cross-ref Strand 4.)

---

## STRAND 3 — Mechanisms & microfoundations of convenience yields

### Money/liquidity services in utility (MIU/BiU reduced form)

**`sidrauski1967rational`** — Sidrauski, "Rational Choice and Patterns of Growth in
a Monetary Economy," *AER P&P* (1967). **Source: liquidity. RF origin. Prox 3.** The
founding money-in-utility paper; the device that becomes bonds-in-utility — a
convenience term `v(B/P)` whose marginal utility is the convenience yield,
decreasing in real holdings. The root of the reduced-form branch big_cy will use.

**`fisher2015structural`** — Fisher, "On the Structural Interpretation of the
Smets-Wouters 'Risk Premium' Shock," *JMCB* 47 (2015). **Source: safety+liquidity.
RF. Prox 3.** The workhorse DSGE "risk-premium" shock is observationally a shock to
demand for safe/liquid short-term assets — a countercyclical, time-varying
convenience wedge. Directly relevant to parameterizing a safety shock.

**`mianstraubsufi2025goldilocks`** — Mian, Straub & Sufi, "A Goldilocks Theory of
Fiscal Deficits," *AER* 115(12) (2025). **Source: liquidity/safe demand. RF (macro
GE). Prox 3.** Government debt supplies a convenience service; an interior optimal
debt level equates marginal convenience to marginal fiscal cost (r<g). Bridges
KVJ-style micro estimates into a macro fiscal model with a downward-sloping demand.

**`michaillatsaez2022economical`** — Michaillat & Saez, "An Economical Business-Cycle
Model," *Oxford Econ. Papers* 74(2) (2022). **Source: wealth-in-utility. RF. Prox 2.**
Safe wealth in utility yields a store-of-value premium supporting a permanent ZLB —
a macro BiU cousin.

### Aggregate Treasury demand (keystone) and Treasury supply

**`kvj2012aggregate`** — Krishnamurthy & Vissing-Jorgensen, "The Aggregate Demand for
Treasury Debt," *JPE* 120(2):233–267 (2012). **Source: liquidity + safety. RF
estimation (micro-motivated). Prox 5.** *Keystone of the strand.* Using ~80 years
and the Aaa-Treasury / CP-bill spreads, estimates a **downward-sloping demand
curve**: falling Treasury/GDP raises the convenience yield. Decomposes it into a
**liquidity** and a **safety** leg (~70 bp combined average). The empirical anchor
for "how large and how variable," and the template for treating convenience as a
demand-curve object. Note: a **relative** spread, so a lower bound on the absolute
level. (Cross-ref Strands 1, 4.)

**`kvj2015impact`** — Krishnamurthy & Vissing-Jorgensen, "The Impact of Treasury
Supply on Financial Sector Lending and Stability," *JFE* 118(3) (2015). **Source:
liquidity + safety. RF/micro bridge. Prox 4.** Government safe-asset supply
**crowds out** private short-term money-like issuance; the moneyness premium falls
in total Treasury supply. Bridges the demand estimate to a microfounded supply side.

**`ghs2015comparative`** — Greenwood, Hanson & Stein, "A Comparative-Advantage
Approach to Government Debt Maturity," *JF* 70(4) (2015). **Source: liquidity/money
premium. Micro. Prox 4.** A money premium on riskless short-term claims; the
government trades it against rollover risk and crowds out private money creation.
Microfounds the dependence of the safe convenience yield on the maturity structure.

**`greenwoodhansonsteinsunderam2023quantity`** — (see Strand 2). **Source:
liquidity/segmentation supply-demand. Micro. Prox 3** in this strand: the
international supply-of-safe-assets channel relevant to big_cy's two-country mapping.

### Near-money / on-the-run liquidity premia (price evidence)

**`nagel2016liquidity`** — Nagel, "The Liquidity Premium of Near-Money Assets," *QJE*
131(4) (2016). **Source: liquidity. RF/micro bridge. Prox 5.** The near-money
liquidity premium (T-bills, GC repo) moves **~one-for-one with the short nominal
rate** (its opportunity cost): near zero at the ZLB, tens of bp when rates are high.
With KVJ, the second key empirical anchor — KVJ gives the supply slope, Nagel the
price/opportunity-cost slope and the cyclicality driver.

**`krishnamurthy2002bond`** — Krishnamurthy, "The Bond/Old-Bond Spread," *JFE*
66(2–3) (2002). **Source: liquidity/search. RF empirical. Prox 4.** The on/off-the-run
30Y spread (≈6–12 bp) reflects a preference for liquid bonds plus repo specialness;
widens in liquidity crises. Foundational evidence of a pure-liquidity premium within
Treasuries; motivates the search models. (Cross-ref Strand 4 / B2.)

**`fleckensteinlongstafflustig2014tips`** — Fleckenstein, Longstaff & Lustig, "The
TIPS-Treasury Bond Puzzle," *JF* 69(5) (2014). **Source: safety/liquidity. RF
empirical. Prox 3.** Nominal Treasuries overpriced vs. inflation-swapped TIPS by
~55 bp (to 200+) — a large convenience premium slow-moving capital fails to
arbitrage; narrows as capital flows in. (Cross-ref Strand 4 / B4.)

### Safety, collateral, information-insensitivity (microfoundations)

**`gortonpennacchi1990financial`** — Gorton & Pennacchi, "Financial Intermediaries
and Liquidity Creation," *JF* 45(1) (1990). **Source: safety/liquidity. Micro.
Prox 4.** Intermediaries create **information-insensitive** debt that trades at par
as a transaction medium; microfounds safety-as-information-insensitivity. The
premium can collapse when debt turns information-sensitive (crises).

**`danggortonholmstrom2017banks`** — Dang, Gorton, Holmström & Ordoñez, "Banks as
Secret Keepers," *AER* 107(4) (2017). **Source: safety/info-insensitivity. Micro.
Prox 3.** Debt is the least information-sensitive security; banks stay opaque so
debt trades like money. Convenience = the value of not producing information; lost
exactly when information is triggered (runs). Microfoundation behind KVJ's safety leg.

### Public debt as private liquidity / collateral (GE microfoundations)

**`woodford1990public`** — Woodford, "Public Debt as Private Liquidity," *AER P&P*
(1990). **Source: liquidity/collateral. Micro. Prox 4.** Seminal: government debt
relaxes private liquidity constraints and earns a liquidity premium even absent
transaction-service assumptions.

**`holmstromtirole1998private`** — Holmström & Tirole, "Private and Public Supply of
Liquidity," *JPE* 106(1) (1998); book-length **`holmstromtirole2011inside`** (MIT
Press, 2011). **Source: liquidity/collateral. Micro. Prox 4.** Imperfect
pledgeability creates a private liquidity shortage; safe public claims earn a
premium because they back continuation of investment. The premium is an
aggregate-state object — high when aggregate liquidity is scarce. A core
microfoundation for a convenience yield on safe public debt.

**`aiyagarimcgrattan1998optimum`** — Aiyagari & McGrattan, "The Optimum Quantity of
Debt," *JME* 42(3) (1998). **Source: liquidity/self-insurance. Micro. Prox 3.** In an
incomplete-markets economy, government debt provides self-insurance to
constrained households; an interior optimal debt quantity, with convenience falling
in debt — ties the KVJ slope to household heterogeneity.

**`angeletoscollarddellas2023public`** — (see Strand 1). **Source: liquidity/
collateral. Micro. Prox 3** here: the modern optimal-policy version mapping debt
quantity to the liquidity premium.

**`diamonddybvig1983bank`** — Diamond & Dybvig, "Bank Runs, Deposit Insurance, and
Liquidity," *JPE* 91(3) (1983). **Source: liquidity-insurance. Micro. Prox 2.** Banks
insure idiosyncratic liquidity shocks via demandable, run-prone claims — background
for why money-like claims exist and why their convenience is fragile.

### Search / OTC liquidity (New Monetarist & microstructure)

**`lagoswright2005unified`** — Lagos & Wright, "A Unified Framework for Monetary
Theory and Policy Analysis," *JPE* 113(3) (2005). **Source: liquidity/search. Micro.
Prox 3.** The workhorse search model in which liquid assets have value because of
search/commitment frictions — the deep benchmark for what BiU/MIU stands in for.

**`lagos2010asset`** — Lagos, "Asset Prices and Liquidity in an Exchange Economy,"
*JME* 57(8) (2010). **Source: liquidity/search. Micro. Prox 3.** Extends the logic to
asset pricing: assets earn a liquidity premium to the extent they serve as media of
exchange.

**`vayanosweill2008search`** — Vayanos & Weill, "A Search-Based Theory of the
On-the-Run Phenomenon," *JF* 63(3) (2008). **Source: liquidity/search. Micro.
Prox 4.** Microfounds the on-the-run premium (the Krishnamurthy 2002 fact) via
search frictions and short-selling.

### Regulatory / private money-creation; macro safe-asset framing

**`stein2012monetary`** — Stein, "Monetary Policy as Financial-Stability
Regulation," *QJE* 127(1) (2012). **Source: regulatory/money-creation externality.
Micro. Prox 3.** Intermediaries over-issue short-term money-like debt because they
ignore the crisis externality; the moneyness premium is what they harvest. The
externality-laden microfoundation for private safe-asset supply. (Adjacent HQLA/LCR
and reserve-demand plumbing — d'Avernas-Vandeweyer, Diamond-Jiang-Schrimpf,
Acharya-Rajan — should be pulled in if an institutional layer is needed; flagged as
a coverage extension.)

**`caballerofarhigourinchas2017safe`** — Caballero, Farhi & Gourinchas, "The Safe
Assets Shortage Conundrum," *JEP* 31(3) (2017). **Source: safety (macro framing).
Survey. Prox 3.** Synthesizes the shortage view: safe-asset demand outstrips supply,
depressing r* and raising the safety premium; the macro stakes of a large absolute CY.

### Additions after librarian-critic review

**`davernasvandeweyer2024treasury`** — d'Avernas & Vandeweyer, "Treasury Bill
Shortages and the Pricing of Short-Term Assets," *JF* 79(6):4083–4141 (2024).
**Source: regulatory/reserve & money-market plumbing. Micro. Prox 3.** Populates the
*institutional* source of convenience the frontier map names (HQLA / short-safe
demand): a structural model of the money market in which T-bill scarcity and
bank reserve/repo demand price short-term safe assets, matching the bill premium
and its interaction with the Fed's balance sheet. The plumbing layer beneath the
reduced-form "liquidity" leg of KVJ/Nagel.

**`gortonlewellenmetrick2012safe`** — Gorton, Lewellen & Metrick, "The Safe-Asset
Share," *AER P&P* 102(3):101–106 (2012). **Source: safety (quantity). Empirical.
Prox 3.** Documents that the share of safe assets in total US assets is remarkably
**stable at ~33%** over 1952–2010, even as the supplier mix shifts (government vs.
financial-sector "shadow money"). The canonical *quantity* fact for the
safe-asset literature — anchors big_cy's "how much safe debt is there?" /
"not enough dollar debt" theme. (Cross-ref Strand 4.)

**`bansalcoleman1996monetary`** — Bansal & Coleman, "A Monetary Explanation of the
Equity Premium, Term Premium, and Risk-Free Rate Puzzles," *JPE* 104(6):1135–1171
(1996). **Source: liquidity/money-services. Micro. Prox 2.** Money/liquid assets
yield non-pecuniary transaction services, so safe/liquid bonds carry a lower
equilibrium return — a monetary (convenience) resolution of the equity-premium and
risk-free-rate puzzles. An early structural antecedent of the DiTella "safe-rate
puzzle" reading (cited in the Strand-1 frontier notes). (Cross-ref Strand 1.)

---

## STRAND 4 — Relative-vs-absolute measurement

*Cross-currency relative measures (A), domestic relative measures (B),
absolute-level measures (C), supply/term-structure determinants (D), currency-wedge
context (E). The absolute-level entries C1/C3 (vBDG, DiTella) are detailed in
Strand 1; reproduced here only by pointer.*

### A. Cross-currency convenience (CIP basis) — RELATIVE across currencies

**`dutepperverdelhan2018cip`** — Du, Tepper & Verdelhan, "Deviations from Covered
Interest Rate Parity," *JF* 73(3) (2018). **Rel (x-ccy). Prox 5.** The canonical
**cross-currency basis**: dollar cash rate vs. FX-swap-implied dollar rate.
Average |basis| ≈ 24 bp (3M), 27 bp (5Y); JPY 5Y ≈ −90 bp; larger at quarter-ends;
persists in credit-risk-free (KfW/repo) form. Differences out convenience common to
both safe legs — the canonical relative measure.

**`duschreger2020handbook`** — Du & Schreger, "CIP Deviations, the Dollar, and
Frictions in International Capital Markets," *Handbook of Int'l Econ.* vol. 6
(working chapter, 2020/2022). **Rel. Prox 5.** Authoritative menu of cross-currency
relative measures, framing CIP deviations as an intermediation *fee* (a gap) and
enumerating the co-moving family of fixed-income near-arbitrage spreads.

**`duimschreger2018treasury`** — Du, Im & Schreger, "The U.S. Treasury Premium,"
*JIE* 112 (2018). **Rel (US vs. foreign sovereign). Prox 5.** The CIP deviation
*between government bonds*: US Treasuries vs. hedged foreign sovereigns. ≈ 21 bp
pre-GFC, up to ~90 bp in crisis, ≈ **−8 bp post-GFC**. The cleanest illustration
that the "Treasury premium" measured relative to other safe sovereigns is a
cross-country *differential*, not the absolute level — and can vanish/reverse even
if the common absolute level is large.

**`jiangkrishnamurthylustig2021foreign`** — (Strand 2). **Rel (Treasury basis).
Prox 5** here: the relative basis is informative for the dollar yet is still a gap;
big_cy's claim is the absolute level can far exceed it.

**`avdjievdukochshin2019dollar`** — Avdjiev, Du, Koch & Shin, "The Dollar, Bank
Leverage, and Deviations from Covered Interest Parity," *AER: Insights* 1(2) (2019).
**Rel. Prox 3.** A stronger broad dollar ⇒ wider basis; interprets the gap's
time-variation via banks' leverage-constraint shadow price.

**`andersonduschlusche2021arbitrage`** — Anderson, Du & Schlusche, "Arbitrage
Capital of Global Banks," NBER WP 28658 (2021; JF version exists — verify). **Rel.
Prox 2.** Limited arbitrage capital sustains the CIP/IOER gaps; banks cut arbitrage
(not loans) after the 2016 MMF-reform funding shock. Background on why relative gaps
persist.

### B. Domestic relative measures — safe vs. near-safe within one currency

**`kvj2012aggregate`** — (Strand 3 keystone). **Rel (Treasury vs. AAA). Prox 5**
here: the AAA–Treasury spread is Treasury-over-AAA, not Treasury-over-(convenience-free
rate); if AAA bonds are money-like, the ~73 bp understates the absolute level.

**`krishnamurthy2002bond`** — (Strand 3). **Rel (liquidity only). Prox 3** here: the
finest-grained relative measure (within Treasuries); safety nets out entirely.

**`longstaff2004flight`** — Longstaff, "The Flight-to-Liquidity Premium in U.S.
Treasury Bond Prices," *J. Business* 77(3) (2004). **Rel (liquidity,
credit-matched). Prox 4.** The Treasury–Refcorp spread (Refcorp is
Treasury-guaranteed but less liquid) isolates a pure **liquidity** premium — large
at long maturities, co-moving with confidence and Treasury supply.

**`fleckensteinlongstafflustig2014tips`** — (Strand 3). **Rel (within issuer). Prox 4**
here: nominal-over-TIPS convenience (~54.5 bp); shared Treasury convenience nets out.

**`klinglersundaresan2019negative`** — Klingler & Sundaresan, "An Explanation of
Negative Swap Spreads," *JF* 74(2) (2019). **Rel (contaminated). Prox 3.** Negative
30Y swap spreads, often read as Treasury richness, are driven by pension duration
demand + dealer balance-sheet limits — a cautionary case that a popular "convenience"
spread is partly a demand/segmentation artifact.

### C. ABSOLUTE-level measures — safe asset vs. a convenience-free benchmark

**`vanbinsbergendiamondgrotteria2022riskfree`** — (Strand 1 / B1). **Abs. Prox 5.**
Box-spread risk-free rate ⇒ Treasury convenience ≈ **40 bp** (≤3Y); the cleanest
option-based absolute measure and the principal countervailing evidence to DiTella.

**`diamondvantassel2024riskfree`** — Diamond & Van Tassel, "Risk-Free Rates and
Convenience Yields Around the World," *JF* (forthcoming 2025; cited as 2024).
**Abs (per currency). Prox 5.** Option-implied riskless rates and convenience yields
across G11; convenience *levels* rise ~linearly with each country's rate (US 5th
largest), while in crises levels rise but the **US-minus-foreign difference does
not**. The international extension of the box-spread method and the natural anchor
for disciplining big_cy's US-vs-RoW convenience parameters.

**`ditella2025zerobeta`** — (Strand 1 keystone). **Abs. Prox 5.** The equity-based
absolute wedge (~2–4%); the large-level claim motivating big_cy.

**`acharyalaarits2023when`** — Acharya & Laarits, "When do Treasuries Earn the
Convenience Yield? — A Hedging Perspective," NBER WP 31863 (2023). **Abs
(conditional). Prox 4.** Decomposes stock-bond covariance into convenience,
frictionless-rate, and default terms; convenience is low when Treasury-equity
covariance is high, eroded by inflation expectations and debt-ceiling episodes.
Connects the absolute level to the **negative beta** of safe assets — directly
relevant to KL's "dollar's negative beta" result big_cy wants from large convenience.

### D. Supply/demand & term-structure determinants

**`greenwoodhansonsteinsunderam2023quantity`** — (Strand 2). **Rel (supply-driven
premia). Prox 3** here: a supply-based theory of *relative* premia and CIP
deviations; the "not enough dollar debt" micro-foundation.

**`nenova2025global`** — Nenova, "Global or Regional Safe Assets: Evidence from Bond
Substitution Patterns," BIS WP (2025). **Segmentation/elasticities. Prox 4.** A
Koijen-Yogo demand system over ~5,000 fund portfolios: US Treasuries are substituted
by *global* (corporate/EM) bonds, Bunds only within a narrow euro-area safe set, and
substitutability deteriorates in stress. Supports big_cy's "insufficient dollar
debt / imperfect substitutability" reading.

**`duschreger2016local`** — Du & Schreger, "Local Currency Sovereign Risk," *JF*
71(3) (2016). **Rel. Prox 3.** EM local-currency sovereign credit spreads via the
swap curve; shows the swap-synthetic-rate machinery in another setting (relevant to
building convenience-adjusted foreign-currency yields).

### E. Currency wedge / theory context

**`engel2016exchange`** — Engel, "Exchange Rates, Interest Rates, and the Risk
Premium," *AER* 106(2) (2016). **Rel (currency wedge). Prox 4.** The UIP puzzle
(high-rate currencies earn excess returns) and the level puzzle (high-real-rate
currencies are *stronger*) are jointly contradictory unless the currency
risk/convenience wedge has a particular structure — motivating a currency-specific
convenience term in the UIP block, exactly the object big_cy puts in utility.

**`jiangkrishnamurthylustig2024dollar`** — (Strand 2). **Rel→abs. Prox 5** here: a
structural use of the relative convenience measure to generate the dollar's negative
beta — the result big_cy seeks from a *large absolute* convenience instead.

**`coppolakrishnamurthyxu2023liquidity`** — (Strand 2). **Abs (liquidity premium).
Prox 3** here: the supply-side microfoundation of dollar convenience and home/dollar
bias, complementing big_cy's demand-side mechanism.

### Additions after librarian-critic review

**`duhebertli2023intermediary`** — Du, Hébert & Li, "Intermediary Balance Sheets and
the Treasury Yield Curve," *JFE* 150(3):103746 (2023). **Rel (intermediary-priced).
Prox 3.** Dealers shifted from net-short to net-long Treasuries after the GFC; the
Treasury yield curve (and its convenience component) is priced off the intermediary
balance-sheet constraint, which moved the curve from a "net-short" to a "net-long"
configuration. The intermediary-constraint determinant of the *relative* Treasury
convenience that complements the demand-curve (KVJ) and opportunity-cost (Nagel)
views. (Cross-ref Strand 3 plumbing: `davernasvandeweyer2024treasury`.)

**`englewu2023liquidity`** — (see Strand 2). **Rel (liquidity/convenience wedge).
Prox 4** here: the empirical measurement of the dollar's liquidity/convenience yield
in the UIP relation — a *relative* (cross-currency) convenience object.

**`gortonlewellenmetrick2012safe`** — (see Strand 3). **Quantity. Prox 3** here: the
stable ~33% safe-asset share is the quantity counterpart to the relative price
measures — how much safe asset exists, not how special it is.

---

*End of master annotated bibliography. ~88 unique entries; full per-entry detail in
`strands/strand_*/annotated_bibliography.md` plus the post-review additions above.*
