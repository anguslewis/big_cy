# Strand 3 — Frontier Notes: Mechanisms → Models, and the RF-vs-Microfounded Tradeoff

## 1. Mechanism families (the competing economic SOURCES of convenience)

| Family | Core idea | Source tag | Representative models |
|---|---|---|---|
| **Money/liquidity services** | Safe/short claims yield transaction & store-of-value services | liquidity | Sidrauski (1967) MIU; Fisher (2015); Michaillat-Saez (2022); Lagos-Wright (2005), Lagos (2010) |
| **Public debt as private liquidity / collateral** | Gov debt relaxes private borrowing & liquidity constraints | liquidity/collateral | Woodford (1990); Holmström-Tirole (1998, 2011); Aiyagari-McGrattan (1998); Angeletos-Collard-Dellas (2023) |
| **Safety / information-insensitivity** | Debt safe enough to trade at par without adverse selection | safety | Gorton-Pennacchi (1990); Dang-Gorton-Holmström-Ordoñez (2017) |
| **Search / OTC liquidity** | Premia arise from trading/search frictions | search | Lagos-Wright (2005); Lagos (2010); Vayanos-Weill (2008); Krishnamurthy (2002) (evidence) |
| **Regulatory / private money-creation** | Institutions demand/supply safe assets; externalities | regulatory | Stein (2012); GHS (2015); KVJ (2015) |
| **Aggregate safe-asset demand (macro framing)** | Global scarcity of safe stores of value | safety | Caballero-Farhi-Gourinchas (2017); Mian-Straub-Sufi (2025) |

These are **not mutually exclusive**: KVJ (2012) explicitly decompose the
empirical convenience yield into a **liquidity** leg (money-likeness) and a
**safety** leg (absolute safety of principal), and most microfoundations target
one leg.

## 2. The reduced-form (BiU/MIU) vs deep-microfoundation split

**Reduced-form (BiU/MIU).** Put a convenience term `v(B/P)` (or wealth/safe-asset
holdings) directly in utility. Examples: Sidrauski lineage, Fisher (2015),
Michaillat-Saez (2022), Mian-Straub-Sufi (2025), and the way KVJ (2012) /
Nagel (2016) estimate a demand curve without committing to one micro mechanism.
- **Buys you:** tractability; a single Euler-equation wedge; easy to crank the
  *level* (to DiTella magnitudes) and to make the wedge time-varying/cyclical;
  estimable slope (KVJ) and price-driver (Nagel).
- **Costs:** the term is a black box — it cannot tell liquidity from safety from
  collateral, is silent on *why* the premium collapses in crises, and is not
  policy-invariant (Lucas critique): the convenience function may shift with
  supply, regulation, or information regime.

**Deep microfoundations.** Derive the premium from primitives:
search/commitment frictions (Lagos-Wright, Lagos, Vayanos-Weill), pledgeability/
collateral constraints (Holmström-Tirole, Woodford, Angeletos et al.,
Aiyagari-McGrattan), information-insensitivity (Gorton-Pennacchi, Dang et al.),
or money-creation externalities (Stein, GHS).
- **Buys you:** structural interpretation, state-contingency (why convenience
  vanishes in runs / scarce-collateral states), policy-relevant comparative
  statics, and discipline on which spreads identify which mechanism.
- **Costs:** model-specific, harder to take to multi-asset/international data,
  and each typically captures only one leg.

**Bridges.** KVJ (2015), GHS (2015), and Stein (2012) connect the estimated
demand curve to a microfounded supply side (private vs public safety production).
Fisher (2015) connects the macro preference-shock wedge to a structural
safe-asset demand shock. These are the natural templates for a project that wants
a BiU term that is *disciplined* by micro evidence.

## 3. What KVJ (2012) and Nagel (2016) establish (level & composition)

- **KVJ (2012):** the Treasury convenience yield is a real, sizable, and
  **supply-sensitive** object. Measured via Aaa-Treasury and CP-bill spreads, it
  averages on the order of ~70 bps over ~80 years, splits into **liquidity** and
  **safety** components (both material; the mix shifts with the debt regime), and
  traces out a **downward-sloping demand curve** in Treasury/GDP — falling supply
  raises the convenience yield. **Cyclicality:** rises in crises / flights to
  safety. This is the empirical *quantity* anchor.
- **Nagel (2016):** the **liquidity** component of the near-money premium moves
  essentially **one-for-one with the short-term nominal rate** (opportunity-cost
  of money). Near zero at the ZLB, tens of bps when rates are high. This is the
  empirical *price/opportunity-cost* anchor and pins the convenience yield's
  cyclicality to the monetary-policy stance.

Together they bound the *level* (tens of bps on the most liquid US safe assets in
relative terms) and the *drivers* (supply quantity; nominal-rate opportunity
cost; crisis safety demand) of the convenience yield. Note the tension with the
big_cy thesis: these are **relative/near-substitute** spreads, so they may
understate the **absolute** convenience yield that DiTella et al. (2025) infer
(2–4%). The mechanism literature here tells you *what such a wedge would be made
of* and *how it should move*, even if the absolute level comes from elsewhere.

## 4. Open gaps relevant to big_cy
- Most microfoundations price *one* leg (liquidity OR safety) and are
  single-country/closed-economy; little maps the mechanism to an **absolute**,
  **international**, two-currency convenience yield.
- The BiU term big_cy plans to use is the reduced form of these mechanisms;
  the disciplined move is to choose the convenience function's slope (KVJ) and
  cyclicality (Nagel / Fisher / crisis safety demand) to match the micro evidence,
  while acknowledging the policy-invariance caveat.
