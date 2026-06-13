# Strand 1 — The Zero-Beta Rate and the ABSOLUTE Level of Convenience Yields

**Annotated bibliography.** Organized by category. Each entry: full cite; one-paragraph
summary (question / method / data / finding incl. sign & magnitude); a line
`RELEVANCE TO ABSOLUTE-CY QUESTION`; a line `IDENTIFICATION OF THE LEVEL` (does it pin
down an *absolute* level, or only a *relative* gap?); and a proximity score (1–5) to the
big_cy question (5 = directly competes / is the keystone; 1 = foundational background).

The organizing distinction for this strand: estimates of the **absolute** level of the
safe-asset convenience yield / safety premium (a wedge between the safe rate and *some
notion of the true intertemporal price or true risk-free rate*) vs. the much larger
literature on **relative** convenience yields (gaps between two near-substitute safe
assets: Treasuries vs. AAA corporates, dollar vs. foreign safe bonds via the CIP basis).

---

## A. The keystone: the zero-beta rate as the intertemporal price

### A1. Di Tella, Hébert, Kurlat & Wang (2025) — "The Zero-Beta Interest Rate" — *JPE* (forthcoming/2025) — **PROXIMITY 5**

Di Tella, Sebastian, Benjamin Hébert, Pablo Kurlat, and Qitong Wang. "The Zero-Beta
Interest Rate." *Journal of Political Economy* (2025).

The central paper of the strand. **Question:** is the safe (T-bill) rate the right measure
of the intertemporal price of consumption, and if not, what is? **Method:** they construct
a *time-varying* measure of the **zero-beta rate** — the conditional expected return on a
unit-investment, minimum-variance equity portfolio orthogonal to the SDF that prices
stocks. They postulate a linear-factor SDF (5 Fama-French factors + a Treasury-bond term
factor + a corporate default factor), estimate betas, form a minimum-variance zero-beta
portfolio (Ledoit-Wolf covariance shrinkage), and project its return on five macro
instruments (T-bill yield, trailing inflation, term spread, excess bond premium of
Gilchrist-Zakrajšek, unemployment). All parameters are estimated jointly by **GMM**
(exactly identified, reduced-rank GLS-type weight matrix) so the zero-beta rate need not
be known to construct excess returns. **Data:** CRSP/Compustat equities (1973–2019/2020),
130 test portfolios (81 triple-sorted beta/size/value/investment/profitability + 49
industry), NIPA nondurables+services consumption. **Findings (magnitudes):** the average
*real* zero-beta rate is **~8.3%** (8.55% excluding 2020), with SD **~9.3%**, vs. a low,
stable real T-bill rate — implying a **spread of ~7.6%/yr** of the zero-beta rate over the
T-bill (nominal zero-beta ~12.0%, ~3.7% inflation). The zero-beta rate (i) fits the
aggregate consumption Euler equation remarkably well, unconditionally (corr with expected
consumption growth 0.56 full / 0.93 ex-2020) and conditional on Nakamura-Steinsson
monetary shocks, with IES in [0.1, 0.2]; (ii) is roughly equal to the average market
return → small/≈zero average equity premium *relative to the zero-beta rate*; (iii) is
volatile and persistent enough to generate price-dividend volatility (std of implied
PD ~29%, vs. 8% using the T-bill) without time-varying risk premia. The new **"safe-rate
puzzle"**: why is the T-bill yield so low, stable, and disconnected from both consumption
and the zero-beta rate? They explicitly *do not* resolve it but speculate (Section 7) that
liquidity/convenience-based theories (Di Tella et al. 2024; Kaplan-Violante; Bansal-Coleman;
Lagos; Herrenbrueck; Kocherlakota 2025) are the natural candidates. Crucially they argue
the equity-based ~6%+ "convenience"/wedge and the ~75bp AAA-Treasury spread are *both*
"convenience" but **there is no reason they should be equal** (citing Diamond 2020,
Angeletos et al. 2023), and that the absolute wedge can be far larger than relative spreads.

RELEVANCE TO ABSOLUTE-CY QUESTION: This is THE source of the "large absolute convenience
yield" claim. The ~6–8% gap between the equity-implied intertemporal price and safe bond
yields is an *absolute*, model-light estimate of how special safe bonds are — orders of
magnitude bigger than the relative (CIP / AAA) spreads big_cy contrasts it against.

IDENTIFICATION OF THE LEVEL: Pins down an absolute level (the intertemporal price) using
*only equity returns* + a linear-factor SDF assumption, without taking a stand on the
theory of the bond wedge or assuming the wedge is the same across safe assets. The level
identification rests on: (i) the linear factor SDF spans the relevant SDF innovations for
equities; (ii) ι (vector of ones) not in the span of the betas (so a zero-variance
portfolio cannot be built); (iii) the minimum-variance zero-beta portfolio loads only
weakly on omitted factors. The authors flag the level estimate as *less robust* than the
Euler-equation test (an omitted high-priced factor uncorrelated with the market could bias
the level). Section 6.4 (multiple-SDF / incomplete markets) is the key caveat: many SDFs
price equities and attribute the spread differently — the procedure recovers the R0 that
does not require an additional factor to which all assets are equally exposed.

---

## B. Other estimates of the ABSOLUTE level / true risk-free rate

### B1. van Binsbergen, Diamond & Grotteria (2022) — "Risk-Free Interest Rates" — *JFE* — **PROXIMITY 5**

van Binsbergen, Jules H., William F. Diamond, and Marco Grotteria. "Risk-Free Interest
Rates." *Journal of Financial Economics* 143, no. 1 (2022): 1–29. (NBER WP 26138, 2019.)

**Question:** what is the true risk-free rate, unaffected by the convenience yield on safe
assets? **Method:** put-call parity / **box spreads** on S&P 500 index options — a fully
collateralized options position that pays a fixed amount with no default risk and (unlike
Treasuries) no monetary/safety convenience — extracts a risk-free rate purely from risky
asset prices, model-free. **Data:** S&P 500 option quotes, minutely, maturities to ~2.5–3
years, 2004–2018ish. **Findings:** the implied (box) risk-free rate sits *above* Treasury
yields; the **Treasury convenience yield ≈ 40 bp on average**, larger at short maturities
(below 3 months), and quadruples during the 2008 crisis. Conventional and unconventional
monetary stimulus reduce the convenience yield in high-frequency event studies.

RELEVANCE TO ABSOLUTE-CY QUESTION: A direct, market-based estimate of the *absolute*
convenience yield of Treasuries relative to a genuinely risk-free (but unspecial) asset.
It is the natural "small-number" counterpoint to Di Tella et al.'s "big-number" (~6%)
wedge — the gap between these two absolute estimates is itself a central open question in
the strand (different notions of "risk-free": no-default vs. zero-beta intertemporal price).

IDENTIFICATION OF THE LEVEL: Yes — identifies an absolute risk-free rate from options, so
the Treasury–box gap is a genuine absolute convenience yield, not a relative spread between
two safe bonds. Rests on no-arbitrage in options + the box position being free of the
convenience benefits Treasuries carry. Measures convenience as *liquidity/safety of
Treasuries vs. a collateralized derivative*, a narrower notion than the zero-beta wedge.

### B2. Krishnamurthy & Vissing-Jorgensen (2012) — "The Aggregate Demand for Treasury Debt" — *JPE* — **PROXIMITY 5**

Krishnamurthy, Arvind, and Annette Vissing-Jorgensen. "The Aggregate Demand for Treasury
Debt." *Journal of Political Economy* 120, no. 2 (2012): 233–267.

**Question:** how much do investors value the liquidity and safety of US Treasuries, and
how does the convenience yield vary with supply? **Method:** time-series regressions of
yield spreads (Treasury vs. AAA/BAA corporates, CD-Treasury, etc.) on the Treasury debt/GDP
ratio; decomposition into a "safety" and a "liquidity" component using assets that share
one attribute but not the other. **Data:** 1926–2008 US Treasury supply and corporate yield
spreads. **Findings:** the convenience yield (liquidity + safety) on Treasuries averages
**~73 bp** over 1926–2008 and is strongly decreasing in Treasury supply (more debt → lower
convenience yield); both safety and liquidity attributes are priced. This is the canonical
"levels" estimate of Treasury specialness — and the explicit *relative* spread the big_cy
thesis says understates the absolute wedge.

RELEVANCE TO ABSOLUTE-CY QUESTION: The reference point for the "relative" view of
convenience yields. Its ~73 bp is what Di Tella et al.'s ~6% absolute wedge is meant to
dwarf; it defines the gap-vs-level tension at the heart of the strand.

IDENTIFICATION OF THE LEVEL: Measures a *relative* gap (Treasury vs. near-substitute
corporates/CDs). It does NOT pin down an absolute level: AAA corporates themselves carry
convenience, so KVJ's number is a lower bound on the absolute Treasury convenience yield.
The supply-elasticity result is the cleanest piece of "level" evidence (a downward-sloping
demand curve for safety), but the intercept (absolute level) is not identified.

### B3. Reis (2021) — "The Constraint on Public Debt when r<g but g<m" — BIS WP 939 / LSE — **PROXIMITY 4**

Reis, Ricardo. "The Constraint on Public Debt when r < g but g < m." Working paper, BIS WP
939 / LSE, 2021.

**Question:** if the safe rate r is below growth g, but the marginal product of capital m
exceeds g, how much fiscal space / seigniorage does the safe-asset wedge create? **Method:**
a general-equilibrium model that *endogenizes* a "bubble premium" on public debt arising
from its safety and liquidity (convenience), then derives debt sustainability and "debt
revenue." **Data:** calibration to US r, g, m. **Findings:** the wedge between m (or the
true return to capital) and the safe rate r — i.e. the absolute convenience/safety premium
on government debt — is large and is the source of substantial "debt revenue"/seigniorage;
more spending requires a larger bubble premium, but demand for debt falls, capping spending.
Inflation reduces fiscal space; financial repression raises it.

RELEVANCE TO ABSOLUTE-CY QUESTION: Frames the safe-asset wedge as r vs. m (the marginal
product of capital / true return), which is conceptually close to Di Tella et al.'s
zero-beta wedge and far larger than relative spreads. Connects the absolute CY directly to
fiscal capacity and seigniorage — a key implication Di Tella et al. flag (seigniorage "much
higher than traditionally thought").

IDENTIFICATION OF THE LEVEL: Conceptual/structural rather than an empirical point estimate:
it defines the absolute wedge as r minus m and shows it has first-order fiscal consequences,
but the magnitude is taken from calibration (the m–r gap, several %), not separately
identified from a portfolio/options procedure.

### B4. Mota (2025, WP) / Mota & Siani — corporate (quasi-)safe asset supply & the safety premium — **PROXIMITY 3**

Mota, Lira. "The Corporate Supply of (Quasi-)Safe Assets." Working paper (JMP, 2021;
revised). (And Mota & Siani, "Financially Sophisticated Firms," WP 2025, in the library.)

**Question:** how large is the *safety* premium specifically (as opposed to liquidity), and
how do corporations respond by manufacturing quasi-safe assets? **Method:** constructs a
near-risk-free synthetic bond from a corporate bond + a maturity-matched CDS, so the
residual yield difference vs. Treasuries isolates a safety premium; demand-based asset
pricing for the corporate supply response. **Findings:** the US Treasury *safety* premium
explains a large fraction (R² ~43–90%, mean ~60%) of time-series variation in relevant
spreads; firms expand quasi-safe issuance when the safety premium is high.

RELEVANCE TO ABSOLUTE-CY QUESTION: Speaks to the *decomposition* of the convenience yield
into safety vs. liquidity components — relevant for interpreting *what* the large absolute
wedge is (Di Tella et al. lean on liquidity/segmentation). Establishes safety alone is a
big, time-varying piece.

IDENTIFICATION OF THE LEVEL: Identifies a *relative* safety spread (corporate-vs-Treasury,
CDS-cleaned), not an absolute level. The CDS-stripped bond still carries some convenience,
so again a lower bound on the absolute safety premium.

### B5. Jiang, Krishnamurthy & Lustig (2021) — "Foreign Safe Asset Demand and the Dollar Exchange Rate" — *JF* — **PROXIMITY 4**

Jiang, Zhengyang, Arvind Krishnamurthy, and Hanno Lustig. "Foreign Safe Asset Demand and
the Dollar Exchange Rate." *Journal of Finance* 76, no. 3 (2021): 1049–1089.

**Question:** how does the convenience yield foreigners derive from US safe assets drive the
dollar? **Method:** identify the (relative) Treasury convenience yield from the **Treasury
basis** — the yield gap between US Treasuries and currency-hedged foreign government bonds
— and relate it to the dollar exchange rate in a portfolio model. **Data:** G10 government
yields, FX forwards, 1988–2017ish. **Findings:** Treasury-basis variation accounts for up
to **41%** of quarterly dollar variation; a wider basis → immediate dollar appreciation,
later depreciation. The convenience yield earned on dollar safe assets is the mechanism.

RELEVANCE TO ABSOLUTE-CY QUESTION: This is the international/relative-convenience benchmark
that big_cy contrasts with the absolute level — the "dollar vs. foreign safe bond" gap. It
also operationalizes the idea that foreign demand for dollar safety has macro/FX
consequences, central to big_cy's international application.

IDENTIFICATION OF THE LEVEL: Relative — the Treasury basis is dollar-vs-foreign-safe-bond,
a near-substitute gap. Does not pin down the absolute dollar convenience yield (foreign govt
bonds also carry convenience).

### B6. Jiang, Krishnamurthy & Lustig (2023) — "Dollar Safety and the Global Financial Cycle" — *REStud* — **PROXIMITY 3**

Jiang, Zhengyang, Arvind Krishnamurthy, and Hanno Lustig. "Dollar Safety and the Global
Financial Cycle." *Review of Economic Studies* (2023).

**Question:** how does time-varying demand for dollar safe assets propagate to global asset
prices and the financial cycle? **Method:** a model where a global convenience yield on
dollar safe assets co-moves with risk premia; tested against dollar, Treasury basis, and
global asset returns. **Findings:** a rise in dollar-safety demand raises the dollar and the
convenience yield and is associated with the global financial cycle; the dollar earns a
safety/convenience premium that is large and time-varying.

RELEVANCE TO ABSOLUTE-CY QUESTION: Extends the relative dollar-convenience machinery to
global comovement — the international-finance payoff big_cy is after. Useful as
context/mechanism, not as an absolute-level estimate.

IDENTIFICATION OF THE LEVEL: Relative (basis-based), as in B5.

---

## C. The zero-beta CAPM / flat security-market-line lineage (level of zero-beta return)

### C1. Black (1972) — "Capital Market Equilibrium with Restricted Borrowing" — *Journal of Business* — **PROXIMITY 4**

Black, Fischer. "Capital Market Equilibrium with Restricted Borrowing." *Journal of
Business* 45, no. 3 (1972): 444–455.

**Question:** what does asset-pricing equilibrium look like without riskless
borrowing/lending? **Method:** theory — derives the two-factor "zero-beta CAPM," replacing
the risk-free rate with the expected return on a portfolio uncorrelated with the market
(the zero-beta portfolio). **Findings:** equilibrium expected returns are linear in beta
with an intercept equal to the (high) zero-beta return, not the T-bill — the conceptual
origin of a zero-beta rate well above safe yields.

RELEVANCE TO ABSOLUTE-CY QUESTION: The theoretical lineage of Di Tella et al.'s object. The
high zero-beta intercept *is* the equity-side manifestation of a large wedge between the
intertemporal price and safe bond yields.

IDENTIFICATION OF THE LEVEL: Theoretical — defines the zero-beta return as the relevant
intertemporal price; later empirical work estimates its level.

### C2. Black, Jensen & Scholes (1972) — "The Capital Asset Pricing Model: Some Empirical Tests" — *Studies in the Theory of Capital Markets* — **PROXIMITY 3**

Black, Fischer, Michael C. Jensen, and Myron Scholes. "The Capital Asset Pricing Model:
Some Empirical Tests." In *Studies in the Theory of Capital Markets*, ed. M. C. Jensen,
79–121. New York: Praeger, 1972.

**Question:** does the empirical security market line match the standard CAPM? **Method:**
beta-sorted portfolios, time-series and cross-section tests. **Findings:** the SML is
**flatter than predicted with a high intercept** — the empirical zero-beta return is
significantly positive and well above the T-bill, the classic anomaly Di Tella et al. revive
and reinterpret.

RELEVANCE TO ABSOLUTE-CY QUESTION: First empirical documentation that the zero-beta return
(≈ intertemporal price implied by equities) is far above safe rates — the original "average
level of the zero-beta rate" finding Di Tella et al. extend to a time series.

IDENTIFICATION OF THE LEVEL: Estimates the *average* level of the zero-beta return
(absolute, equity-based) but not its time variation, and pre-dates the convenience-yield
interpretation.

### C3. Frazzini & Pedersen (2014) — "Betting Against Beta" — *JFE* — **PROXIMITY 3**

Frazzini, Andrea, and Lasse Heje Pedersen. "Betting Against Beta." *Journal of Financial
Economics* 111, no. 1 (2014): 1–25.

**Question:** why is the SML flat? **Method:** model + empirics where leverage-constrained
investors bid up high-beta assets, flattening the SML; a "betting against beta" (BAB)
factor earns high risk-adjusted returns across asset classes. **Findings:** constrained
investors → a high zero-beta rate and a flat SML; BAB is large and significant globally.

RELEVANCE TO ABSOLUTE-CY QUESTION: Di Tella et al. (Sec 6.4) use BAB as the leading
*alternative* interpretation of a high zero-beta rate: leverage-constrained agents value
"convenience" of safe bonds (the constraint multiplier enters their bond Euler equation but
not their stock Euler equation), generating a high zero-beta rate without consumption. This
is a competing micro-foundation for the same absolute wedge.

IDENTIFICATION OF THE LEVEL: Identifies a high zero-beta intercept and attributes it to
leverage constraints (a convenience-like force), but does not model consumption, so it
cannot speak to whether the zero-beta rate is the intertemporal price.

### C4. Lopez-Lira & Roussanov (2020, WP) — high-return portfolios orthogonal to factors — **PROXIMITY 2**

Lopez-Lira, Alejandro, and Nikolai Roussanov. "Do Common Factors Really Explain the
Cross-Section of Stock Returns?" Working paper, 2020.

**Question:** can a diversified portfolio orthogonal to common factors still earn a high
return? **Method:** strip common factors from returns and form diversified portfolios.
**Findings:** even removing almost all common factors, one can build a diversified portfolio
with a high average return — i.e. a high zero-beta-type return survives factor controls.

RELEVANCE TO ABSOLUTE-CY QUESTION: Supports the robustness of a high zero-beta return to
rich factor controls (Di Tella et al. cite it for exactly this), bolstering the claim that
the large absolute wedge is not an artifact of omitted standard factors.

IDENTIFICATION OF THE LEVEL: Indirect — shows the high level is hard to factor away; not a
structural level estimate.

---

## D. Theory of why safe bonds are special / the bond–stock SDF wedge

### D1. Diamond (2020, WP) — convenient safe assets in general equilibrium — **PROXIMITY 4**

Diamond, William. "Safety Transformation and the Structure of the Financial System."
Working paper (cited by Di Tella et al. as Diamond 2020).

**Question:** how does demand for safe/convenient assets shape the financial system and the
relation among different pricing kernels? **Method:** GE model with households valuing safe
assets and intermediaries supplying them. **Findings:** the household's safe-asset demand
makes the risk-free rate that prices *household* claims lie strictly *below* the risk-free
rate implied by both the household's and the intermediary's pricing kernels for *risky*
assets — i.e. multiple "risk-free rates" coexist and the bond convenience wedge differs by
who/what is being priced. (Di Tella et al. quote this directly to justify why the equity
wedge and AAA-Treasury spread need not coincide.)

RELEVANCE TO ABSOLUTE-CY QUESTION: The core theoretical statement that *the convenience
wedge is not a single number* — relative spreads (AAA-Treasury, ~75bp) and the equity-based
wedge (~6%) are both "convenience" but legitimately different. This is the intellectual
hinge of big_cy's "absolute vs. relative" thesis.

IDENTIFICATION OF THE LEVEL: Theoretical — explains why an absolute, equity-based level can
greatly exceed relative spreads; does not itself estimate the level.

### D2. Angeletos, Collard & Dellas (2023) — public debt as convenient safe asset / optimal debt — **PROXIMITY 3**

Angeletos, George-Marios, Fabrice Collard, and Harris Dellas. "Public Debt as Private
Liquidity: Optimal Policy." *Journal of Political Economy* (2023). (Di Tella et al. cite as
Angeletos et al. 2023.)

**Question:** what is optimal public debt when government debt provides liquidity/convenience
services? **Method:** incomplete-markets GE with a convenience/liquidity premium on public
debt. **Findings:** a sizable convenience premium justifies issuing more public debt than
the standard tax-smoothing benchmark; the wedge between the safe rate and the "social" rate
matters for optimal policy.

RELEVANCE TO ABSOLUTE-CY QUESTION: Reinforces Diamond's point that the convenience wedge is
structurally distinct across assets and large enough to reshape optimal debt policy —
context for the fiscal/seigniorage implications of a large absolute CY.

IDENTIFICATION OF THE LEVEL: Theoretical/quantitative via calibration of a convenience
premium; not an independent empirical level estimate.

### D3. Kocherlakota (2025) — fiscal policy welfare when safe rates carry a liquidity premium — **PROXIMITY 3**

Kocherlakota, Narayana. "Public Debt Bubbles in Heterogeneous-Agent Models" / fiscal-policy
welfare paper (cited by Di Tella et al. as Kocherlakota 2025).

**Question:** what are the welfare implications of raising interest rates when the safe-rate
spread is part risk premium, part liquidity premium? **Method:** incomplete-markets model
(Kaplan-Violante spirit) with a liquidity premium on safe debt. **Findings:** raising rates
*improves* welfare (in a zero-growth economy) precisely when the safe rate is negative and
the zero-beta rate satisfies the aggregate Euler equation — i.e. exactly the Di Tella et al.
configuration. Seigniorage (gap between actual interest expense and what it would be at the
zero-beta rate) is much larger than usually measured.

RELEVANCE TO ABSOLUTE-CY QUESTION: Takes the Di Tella et al. facts as input and derives
fiscal-welfare and seigniorage consequences of a *large absolute* convenience wedge — a key
"so what" for big_cy.

IDENTIFICATION OF THE LEVEL: Uses the absolute wedge as a primitive; does not estimate it.

### D4. Di Tella, Hébert, Kurlat & Wang (2024, WP) — liquid/illiquid accounts micro-foundation — **PROXIMITY 4**

Di Tella, Sebastian, Benjamin Hébert, Pablo Kurlat, and Qitong Wang. "Aggregate Euler
Equations with Liquid and Illiquid Assets." Working paper, 2024. (Companion to the JPE
paper.)

**Question:** can a heterogeneous-agent model rationalize an aggregate Euler equation for
illiquid (equity) but not liquid (bond) assets? **Method:** Kaplan-Violante-style model with
uninsurable idiosyncratic risk, liquid (bond-backed) and illiquid (stock) accounts, and
rebalancing frictions; extends Werning (2015) aggregation. **Findings:** conditions under
which the aggregate consumption Euler equation holds for the *illiquid* asset return (≈
zero-beta rate) but fails for *liquid* safe bonds — providing a micro-foundation for the
keystone empirics.

RELEVANCE TO ABSOLUTE-CY QUESTION: The companion theory that turns the empirical "safe-rate
puzzle" into a liquidity wedge — directly relevant to how big_cy might microfound a large
absolute CY in a structural model.

IDENTIFICATION OF THE LEVEL: Structural; rationalizes the level rather than re-estimating it.

---

## E. Classic puzzles being reinterpreted (foundational)

### E1. Mehra & Prescott (1985) — "The Equity Premium: A Puzzle" — *JME* — **PROXIMITY 2**

Mehra, Rajnish, and Edward C. Prescott. "The Equity Premium: A Puzzle." *Journal of
Monetary Economics* 15, no. 2 (1985): 145–161.

**Question:** can consumption covariance explain the equity premium over T-bills? **Method:**
calibrated Lucas-tree Euler equation. **Findings:** the historical ~6% equity premium over
safe bonds requires implausibly high risk aversion — the equity premium puzzle. Di Tella et
al. reinterpret this as a *safe-rate* puzzle: had Mehra-Prescott used the spread of the
market over the *zero-beta* rate (≈0) rather than over T-bills, there would be no puzzle —
"it is the bonds that create the puzzle."

RELEVANCE TO ABSOLUTE-CY QUESTION: Defines the ~6% number that the absolute convenience
wedge is meant to *absorb*. The reinterpretation is the conceptual bridge from "equity
premium" to "absolute convenience yield on safe bonds."

IDENTIFICATION OF THE LEVEL: Foundational; not a CY estimate.

### E2. Weil (1989) — "The Equity Premium Puzzle and the Risk-Free Rate Puzzle" — *JME* — **PROXIMITY 2**

Weil, Philippe. "The Equity Premium Puzzle and the Risk-Free Rate Puzzle." *Journal of
Monetary Economics* 24, no. 3 (1989): 401–421.

**Question:** with Epstein-Zin preferences, can one jointly match the equity premium and the
risk-free rate? **Method:** recursive-utility calibration. **Findings:** matching the low
risk-free rate requires implausible time-preference/IES parameters — the **risk-free rate
puzzle** (the safe rate is "too low"). This is the precursor framing of Di Tella et al.'s
safe-rate puzzle.

RELEVANCE TO ABSOLUTE-CY QUESTION: The "low safe rate" puzzle is exactly what a large
absolute convenience yield would explain — the safe rate is low because safe bonds yield a
big non-pecuniary (convenience) return.

IDENTIFICATION OF THE LEVEL: Foundational; frames the level question, no estimate.

### E3. Campbell & Shiller (1988) / Campbell (1991) — excess volatility & PD-ratio decomposition — **PROXIMITY 2**

Campbell, John Y., and Robert J. Shiller. "The Dividend-Price Ratio and Expectations of
Future Dividends and Discount Factors." *Review of Financial Studies* 1, no. 3 (1988):
195–228. Campbell, John Y. "A Variance Decomposition for Stock Returns." *Economic Journal*
101 (1991): 157–179.

**Question:** what drives price-dividend-ratio volatility — expected dividends, expected
returns, or rates? **Method:** log-linear present-value (Campbell-Shiller) decomposition;
VAR. **Findings:** safe rates and dividends are not volatile enough → variation in valuation
ratios must reflect time-varying *excess* returns (the excess-volatility puzzle). Di Tella
et al. show that using the volatile zero-beta rate instead of the T-bill, valuation-ratio
volatility can be generated by discount-rate variation *without* time-varying risk premia.

RELEVANCE TO ABSOLUTE-CY QUESTION: The volatility side of the keystone's reinterpretation —
a high *and volatile* zero-beta rate (hence a volatile absolute convenience wedge) replaces
time-varying risk premia. Supports that the absolute wedge is large *and* time-varying.

IDENTIFICATION OF THE LEVEL: Foundational/methodological; the decomposition is the tool Di
Tella et al. repurpose.

---

## F. Methods papers (econometric tools the level estimates rely on)

### F1. Ledoit & Wolf (2017) — nonlinear shrinkage covariance estimation — **PROXIMITY 2**

Ledoit, Olivier, and Michael Wolf. "Nonlinear Shrinkage of the Covariance Matrix for
Portfolio Selection: Markowitz Meets Goldilocks." *Review of Financial Studies* 30, no. 12
(2017): 4349–4388.

**Question:** how to estimate large covariance matrices for minimum-variance portfolios
without overfitting? **Method:** nonlinear shrinkage of eigenvalues. **Findings:**
outperforms alternatives on out-of-sample portfolio variance. Di Tella et al. use it to
build the minimum-variance zero-beta portfolio — a load-bearing ingredient in their level
estimate.

RELEVANCE TO ABSOLUTE-CY QUESTION: Methodological underpinning of the keystone's
absolute-level construction (the min-variance zero-beta portfolio).

IDENTIFICATION OF THE LEVEL: Tool, not estimate.

### F2. Olea & Pflueger (2013) — weak-instrument-robust inference — **PROXIMITY 1**

Montiel Olea, José Luis, and Carolin Pflueger. "A Robust Test for Weak Instruments."
*Journal of Business & Economic Statistics* 31, no. 3 (2013): 358–369.

**Question:** how to test for weak instruments under heteroskedasticity/autocorrelation?
**Method:** an effective F-statistic. **Findings:** Di Tella et al. report it (effective F
~30 GMM) to show their macro instruments meaningfully predict the zero-beta portfolio return.

RELEVANCE TO ABSOLUTE-CY QUESTION: Supports the statistical validity of the time-series
zero-beta-rate estimate.

IDENTIFICATION OF THE LEVEL: Tool, not estimate.

### F3. Gilchrist & Zakrajšek (2012) — the excess bond premium — *AER* — **PROXIMITY 2**

Gilchrist, Simon, and Egon Zakrajšek. "Credit Spreads and Business Cycle Fluctuations."
*American Economic Review* 102, no. 4 (2012): 1692–1720.

**Question:** what part of credit spreads reflects a non-default "excess bond premium"
(EBP)? **Method:** decompose corporate spreads into expected default and a residual EBP.
**Findings:** the EBP predicts business cycles and risk appetite. Di Tella et al. use the
EBP as a key instrument; it negatively predicts the zero-beta portfolio return (zero-beta
rate falls when recession is likely).

RELEVANCE TO ABSOLUTE-CY QUESTION: A key instrument linking the absolute wedge to credit
conditions/risk appetite — informs the *cyclicality* of the convenience wedge.

IDENTIFICATION OF THE LEVEL: Instrument/tool.
