# Strand 3 — What Gives Rise to Convenience Yields: Mechanisms & Microfoundations

Annotated bibliography. For each entry: **source of convenience** (liquidity /
safety / collateral / regulatory / search), whether it is **reduced-form (BiU/MIU)**
or **microfounded**, implications for the **level** and **cyclicality** of the
convenience yield, and **proximity** to big_cy (1–5).

Tags per entry: `[SOURCE] [RF or MICRO] [proximity P]`.

Note: this is a CREATOR document — collect and organize, no quality scoring.

---

## A. Money / liquidity services in utility (MIU / Sidrauski lineage); BiU as reduced form

### Sidrauski (1967) — Rational Choice and Patterns of Growth in a Monetary Economy
`[liquidity] [RF/MICRO origin] [P3]`
The founding money-in-the-utility-function (MIU) paper. Real balances enter the
representative household's utility directly as a stand-in for the unmodeled
transaction/liquidity services of money. Establishes the modeling device that
later becomes "bonds-in-the-utility" (BiU): a convenience term `v(B/P)` whose
marginal utility is the convenience yield. **Level/cyclicality:** silent on
empirics; the device makes the convenience yield equal to the marginal liquidity
service, which is decreasing in the real quantity held. The intellectual root of
the reduced-form branch big_cy plans to use.

### Fisher (2015) — Structural Interpretation of the Smets-Wouters "Risk Premium" Shock
`[safety+liquidity] [RF] [P3]`
Shows the workhorse DSGE "risk premium" shock is observationally equivalent to a
structural shock to the demand for safe and liquid short-term assets (Treasuries).
Effectively relabels a wedge in the Euler equation as a time-varying convenience
yield — a clean macro example of BiU/preference-shock reduced form. **Level/
cyclicality:** the wedge is time-varying and recessionary (demand for safety
spikes in downturns), giving the convenience yield a strongly countercyclical
flavor. Directly relevant to how big_cy parameterizes a time-varying safety shock.

### Michaillat & Saez (2022) — An Economical Business-Cycle Model
`[liquidity/wealth-in-utility] [RF] [P2]`
Wealth (a safe store of value) enters utility directly, generating a positive
"liquidity"/store-of-value premium that lets the model sit at a permanent ZLB.
A macro BiU cousin: convenience accrues to safe wealth rather than transactions
balances. **Level/cyclicality:** premium rises when desired saving outruns
investment; relevant as a demand-side rationalization of low safe rates.

### Mian, Straub & Sufi (2025, AER; NBER 2022) — A Goldilocks Theory of Fiscal Deficits
`[liquidity/safe-asset demand] [RF] [P3]`
Government debt supplies a convenience/liquidity service demanded by households;
the model pins down a debt level where the marginal convenience value equals the
marginal fiscal cost (r < g region). Macro general-equilibrium use of a BiU-style
convenience term to study optimal deficits. **Level/cyclicality:** convenience
yield is decreasing in the supply of public debt — the central downward-sloping
demand curve big_cy relies on. Bridges KVJ-style micro estimates into a macro
fiscal model.

---

## B. Aggregate demand for Treasury debt (keystone) and Treasury supply

### Krishnamurthy & Vissing-Jorgensen (2012, JPE) — The Aggregate Demand for Treasury Debt  *(KEYSTONE)*
`[liquidity + safety] [RF estimation, MICRO-motivated] [P5]`
The central empirical paper of the strand. Uses ~80 years of US data and the
**Aaa-Treasury** and **commercial-paper-bill** spreads as price measures of a
Treasury "convenience yield," then estimates a **downward-sloping demand curve**:
when the Treasury-debt-to-GDP ratio falls, the convenience yield rises (and vice
versa). Decomposes convenience into a **liquidity** component (Treasuries are
money-like, traded cheaply) and a **safety** component (absolute safety of
principal), identifying the two by comparing assets that differ in one dimension
but not the other (e.g., Aaa corporates are safe-ish but not as liquid; T-bills
vs. very-safe near-substitutes). **Level:** the combined convenience yield
averages on the order of ~70 basis points historically, with liquidity and safety
each contributing a meaningful share and the split shifting with the
debt-to-GDP regime. **Cyclicality:** the convenience yield is **decreasing in
Treasury supply** and **rises in crises / when safe-asset demand spikes**.
This is the paper big_cy must engage most directly: it is the empirical anchor
for "how large and how variable is the convenience yield," and the methodological
template for treating it as a demand-curve object rather than a fixed wedge.

### Krishnamurthy & Vissing-Jorgensen (2015, JFE) — Impact of Treasury Supply on Financial Sector Lending and Stability
`[liquidity + safety] [RF/MICRO bridge] [P4]`
Companion to KVJ 2012. Nonfinancial-sector demand for safe/liquid assets drives a
moneyness premium; the financial sector exploits it by issuing short-term money-like
claims against risky assets. Central prediction (confirmed): government safe-asset
supply **crowds out** private short-term debt issuance. **Level/cyclicality:** the
moneyness premium is declining in total Treasury supply; links the convenience
yield to financial-stability externalities. Bridges the demand estimate to a
microfounded supply-side story.

### Greenwood, Hanson & Stein (2015, JF) — A Comparative-Advantage Approach to Government Debt Maturity
`[liquidity/money-premium] [MICRO] [P4]`
Investors derive monetary services from riskless short-term claims (a money
premium on T-bills). The government trades the money premium it can harvest by
issuing short-term debt against rollover/refinancing risk; private intermediaries
compete to supply money-like claims, and with negative externalities of private
money creation the government should tilt issuance short to crowd them out.
**Level/cyclicality:** short-rate money premium is a falling function of the
quantity of short safe paper. A microfoundation for why the safe convenience
yield depends on the maturity structure of supply.

### Greenwood, Hanson, Stein & Sunderam (2020 WP) — A Quantity-Driven Theory of Term Premia and Exchange Rates
`[liquidity/segmentation supply-demand] [MICRO] [P4]`
Specialized bond investors absorb supply/demand shocks for long-term bonds in two
currencies; FX and long bonds load on the same short-rate factor. Quantity shifts
(e.g., QE) move term premia and the exchange rate; an extension rationalizes
post-2008 CIP deviations. **Level/cyclicality:** premia driven by relative
bond supply across currencies — directly on-theme for big_cy's two-country
convenience-yield mapping (the international supply-of-safe-assets channel).

---

## C. Near-money / on-the-run liquidity premia (empirical price evidence)

### Nagel (2016, QJE) — The Liquidity Premium of Near-Money Assets
`[liquidity] [RF/MICRO bridge] [P5]`
Shows the convenience/liquidity premium on near-money assets (T-bills, repo
relative to general collateral) moves **one-for-one with the short-term nominal
interest rate**: when the policy rate is high, the opportunity cost of holding
money-like assets is high, so the liquidity premium they command is high. Provides
a tight, almost arbitrage-style measurement of the liquidity component of the
convenience yield and a parsimonious driver (the level of nominal rates) for its
**cyclicality**. **Level:** premium scales with nominal rates (near zero at the
ZLB, tens of bps when rates are high). Together with KVJ 2012 this is the most
important empirical anchor for the strand: KVJ gives the supply/quantity slope,
Nagel gives the price/opportunity-cost slope.

### Krishnamurthy (2002, JFE) — The Bond/Old-Bond Spread
`[liquidity/search] [RF empirical, search-motivated] [P4]`
Documents the on-the-run/off-the-run 30-year Treasury spread (≈6 bps average,
1995–1999) and shows it reflects investors' preference for liquid (on-the-run)
bonds plus repo "specialness" carry costs; the spread covaries with Treasury
supply and aggregate liquidity factors. Foundational price evidence that even
within Treasuries a pure-liquidity convenience premium exists. **Cyclicality:**
spread widens in liquidity crises (e.g., Fall 1998). Motivates the search models
below.

### Fleckenstein, Longstaff & Lustig (2014, JF) — The TIPS-Treasury Bond Puzzle
`[safety/liquidity] [RF empirical] [P3]`
Nominal Treasuries are persistently overpriced relative to synthetic
(inflation-swapped TIPS) Treasuries by an average ~55 bps (up to 200+ bps),
implying a large convenience premium on the most liquid/safe instrument that
slow-moving capital fails to arbitrage away. **Level:** a direct, large estimate
of the Treasury convenience yield from a near-arbitrage. **Cyclicality:** mispricing
narrows as capital flows in; widens in stressed, capital-constrained states.

---

## D. Safety, collateral value, and information-insensitivity (microfoundations)

### Gorton & Pennacchi (1990, JF) — Financial Intermediaries and Liquidity Creation
`[safety/liquidity] [MICRO] [P4]`
Foundational "why safe debt is special" paper. Intermediaries create
**information-insensitive** debt (deposits) that is safe enough that uninformed
traders are not exploited, so it circulates at par as a transaction/liquidity
medium. Microfounds the safety-as-information-insensitivity source of convenience.
**Level/cyclicality:** convenience accrues to assets immune to adverse selection;
the premium can collapse when debt becomes information-sensitive (crises).

### Dang, Gorton, Holmström & Ordoñez (2017, AER) — Banks as Secret Keepers
`[safety/information-insensitivity] [MICRO] [P3]`
Formalizes the Dang-Gorton-Holmström information-insensitivity program: debt is
the least information-sensitive security; banks stay deliberately opaque so their
debt trades at par like money. Convenience = the value of not needing to produce
information. **Cyclicality:** money-like debt loses its convenience precisely when
information production is triggered (runs/crises) — a sharp account of why safe
premia are state-contingent. Microfoundation behind the "safety" leg of KVJ.

---

## E. Public debt as private liquidity / collateral (general-equilibrium microfoundations)

### Woodford (1990, AER P&P) — Public Debt as Private Liquidity
`[liquidity/collateral] [MICRO] [P4]`
Seminal GE statement that government debt relaxes private borrowing/liquidity
constraints and thus commands a liquidity premium (carries a convenience yield)
even absent transaction-service assumptions. Establishes "public debt as collateral
for private trade." **Level/cyclicality:** premium positive whenever private
liquidity is scarce; rises when constraints tighten.

### Holmström & Tirole (1998, JPE) — Private and Public Supply of Liquidity
`[liquidity/collateral] [MICRO] [P4]`
Imperfect pledgeability of corporate income generates a private shortage of
liquid instruments; firms hoard claims to weather liquidity shocks, and government
securities can fill the gap when private collateral is insufficient. The LAPM
(liquidity asset pricing) logic: safe public claims earn a liquidity premium
because they back continuation of investment. **Level/cyclicality:** premium is an
**aggregate-state** object — high when aggregate liquidity is scarce. Core
microfoundation for a convenience yield on safe public debt. (Book-length
treatment: Holmström-Tirole 2011, *Inside and Outside Liquidity*.)

### Aiyagari & McGrattan (1998, JME) — The Optimum Quantity of Debt
`[liquidity/self-insurance] [MICRO] [P3]`
In a Bewley-Aiyagari incomplete-markets economy, government debt provides
self-insurance/liquidity to borrowing-constrained households, so there is an
interior optimal debt quantity. Microfounds a downward-sloping convenience-yield-
vs-debt relationship from precautionary demand. **Level/cyclicality:** convenience
falls as debt rises; ties the KVJ demand slope to household heterogeneity.

### Angeletos, Collard & Dellas (2023, JPE) — Public Debt as Private Liquidity: Optimal Policy
`[liquidity/collateral] [MICRO] [P3]`
Modern optimal-policy version of Woodford/Holmström-Tirole: public debt eases a
financial friction (serves as collateral/liquidity buffer), but more debt lowers
the liquidity premium and raises the government's borrowing cost — pinning down a
unique long-run optimal debt level. **Level/cyclicality:** explicit mapping from
debt quantity to the liquidity premium; the policy-relevant frontier of the
public-liquidity literature.

### Diamond & Dybvig (1983, JPE) — Bank Runs, Deposit Insurance, and Liquidity
`[liquidity-insurance] [MICRO] [P2]`
The canonical liquidity-provision/maturity-transformation model: banks insure
households against idiosyncratic liquidity shocks by issuing demandable claims,
which are run-prone. Background microfoundation for why money-like claims exist
and why their convenience is fragile. **Cyclicality:** run equilibria destroy the
liquidity service in bad states. (The "relevant Diamond" for this strand.)

---

## F. Search / OTC liquidity (New Monetarist & market-microstructure microfoundations)

### Lagos & Wright (2005, JPE) — A Unified Framework for Monetary Theory and Policy Analysis
`[liquidity/search] [MICRO] [P3]`
The workhorse New Monetarist (search-theoretic) model in which money/liquid assets
have value because bilateral exchange suffers search and commitment frictions. The
explicit microfoundation that BiU/MIU stands in for: it derives, rather than
assumes, a liquidity value for money-like assets. **Level/cyclicality:** liquidity
premium depends on trading frictions and the supply of acceptable media of
exchange; provides the deep "what BiU means" benchmark for the methodological split.

### Lagos (2010, JME) — Asset Prices and Liquidity in an Exchange Economy
`[liquidity/search] [MICRO] [P3]`
Extends the New Monetarist logic to asset pricing: assets command a liquidity
premium to the extent they serve as media of exchange / can be readily traded,
generating return differences across otherwise-similar assets. Microfounds a
liquidity-based convenience yield in an asset-pricing setting. **Cyclicality:**
premium varies with the value of liquidity services and asset acceptability.

### Vayanos & Weill (2008, JF) — A Search-Based Theory of the On-the-Run Phenomenon
`[liquidity/search] [MICRO] [P4]`
Microfounds the on-the-run premium (the Krishnamurthy 2002 fact) via search
frictions and short-selling: liquidity concentrates in one security, raising its
price and lowering its yield endogenously. The deep model behind the on-the-run
convenience premium. **Level/cyclicality:** premium rises with search frictions /
shorting demand and in scarce-collateral states.

---

## G. Regulatory / institutional safety-production (banks as safety producers)

### Stein (2012, QJE) — Monetary Policy as Financial-Stability Regulation
`[regulatory/money-creation externality] [MICRO] [P3]`
Models private money creation: intermediaries issue too much short-term money-like
debt because they do not internalize the crisis externality, and the convenience/
moneyness premium is what they harvest. Provides the externality-laden
microfoundation for why private safe-asset supply responds to the convenience
yield (and why it should be regulated). **Level/cyclicality:** the money premium
governs the incentive to over-issue; complements GHS 2015 and KVJ 2015 on the
public-vs-private safety-supply margin. Closest "regulatory-demand /
safety-production" entry available in this strand's reading; the HQLA/LCR and
reserve-demand literature (d'Avernas-Vandeweyer, Diamond-Jiang-Schrimpf,
Acharya-Rajan) sits adjacent and should be pulled in if the survey needs the
institutional-plumbing layer.

### Caballero, Farhi & Gourinchas (2017, JEP) — The Safe Assets Shortage Conundrum
`[safety] [survey/RF framing] [P3]`
Synthesizes the "safe asset shortage" view: global demand for safe stores of value
outstrips supply, depressing safe rates and raising the safe-asset convenience
yield, with macro consequences (low r*, the safety trap). **Level/cyclicality:**
the shortage — and hence the safety premium — intensifies in flights to safety.
Frames the macro stakes of a large absolute convenience yield central to big_cy.

---

## Cross-strand pointers
- Relative-vs-absolute measurement (CIP basis, Aaa-Treasury) is Strand 1/2 territory;
  this strand supplies the *mechanisms* those spreads are pricing.
- Time-varying safety shocks (Fisher 2015; KVJ 2012 cyclicality) connect to
  Kekre-Lenel's safety shock and DiTella et al.'s large absolute wedge — see the
  big_cy anchor papers.
