# Strand 2 — Frontier notes: how these models relate

## The big organizing axis: where does safe-asset demand come from?

Two questions split the field. (1) Is the special demand for the safe/dollar asset
**exogenous** (preference/utility shock, bonds-in-utility, reduced-form clientele or
noise demand) or **endogenous** (collateral usefulness, intermediary/financial
constraints, search/settlement liquidity, coordination)? (2) Is the convenience
notion **relative** (a *gap* between near-substitutes — the Treasury basis, the
CIP/LIBOR basis, cross-currency spreads) or **absolute** (a *level* on the safe
asset itself, e.g. an estimated US convenience yield or a monopoly markup)?

These two questions are the spine of big_cy's contribution: nearly every
*structural* model that takes convenience as a primitive (KL2024, JKL2021/2024)
**disciplines it against a relative basis**, even though the underlying object is
conceptually an absolute level. The papers that work with an *absolute* notion
tend to be either demand-system estimates (Koijen-Yogo 2020) or supply-side /
monopoly models (Choi-Kirpalani-Perez, Farhi-Maggiori).

## Three modeling families

1. **Bonds-in-utility / exogenous-convenience structural NK** (Kekre-Lenel 2024;
   Choi-Kirpalani-Perez 2022; the demand block of JKL). The safe asset enters
   utility (or a nonpecuniary value), demand is a primitive shock, and the model
   traces GE implications for output, portfolios, and risk premia. **This is
   big_cy's template.** Cleanest place to "crank the absolute level."

2. **Intermediary / financier portfolio-balance** (Gabaix-Maggiori 2015;
   Maggiori 2017; Gourinchas-Ray-Vayanos 2022; Greenwood-Hanson-Stein-Sunderam
   2023; Kekre-Lenel 2025; Bianchi-Bigio-Engel 2023). Limited risk-bearing
   capacity of a constrained intermediary prices currency/term/liquidity risk.
   Convenience is implicit (a risk premium or, in BBE, an endogenous liquidity
   yield). Exchange rate pinned by gross flows + balance-sheet capacity.

3. **Safe-asset supply / shortage / dominance** (Caballero-Farhi-Gourinchas 2008,
   2021; Caballero-Farhi 2018; Farhi-Maggiori 2017; He-Krishnamurthy-Milbradt
   2019; Coppola-Krishnamurthy-Xu 2023; Gopinath-Stein 2021). Emphasis on who
   *supplies* safety, scarcity, ZLB/safety traps, monopoly rents, and why the
   dollar (network/complementarity/coordination/fiscal capacity).

Plus a cross-cutting **demand-system estimation** approach (Koijen-Yogo 2019/2020;
Jiang-Richmond-Zhang 2022) that *measures* the absolute convenience yield and
NFA-composition drivers — the empirical bridge to big_cy's quantitative aims.

## Comparative table

| Paper | Demand: exo/endo | Open/closed | Source of dollar/safe convenience | Implies for US carry / home bias / dollar | Rel/Abs convenience | Prox |
|---|---|---|---|---|---|---|
| Kekre-Lenel 2024 (AER) | Exo (BIU safety shock) + supply | Open 2-ctry NK | Nonpecuniary value of safe dollar bonds | US long-capital/short-dollar; dollar = neg-beta hedge; +avg NFA return from risk-tolerance gap | Relative (Treasury basis) | 5 |
| JKL 2024 (ReStud) | Exo demand (qty-dep) + endo supply | Open US/RoW | Foreign convenience on safe dollar claims | US issues cheap safe dollar debt (privilege); GFC=dollar cycle; resolves paradox w/o risk tolerance | Relative→absolute framing | 5 |
| JKL 2021 (JF) | Exo (measured via basis) | Open (kernel) | Foreign safety/liquidity on dollar; 90% is dollar per se | Convenience drives dollar level & disconnect | Relative (basis) | 5 |
| Maggiori 2017 (AER) | Endo (intermediary asymmetry) | Open 2-ctry | US financial-development depth (risk-bearing) | US=world insurer; reserve-ccy paradox raised | Absolute-ish risk price; no basis | 4 |
| Devereux-Engel-Wu 2023 | Endo (collateral) | Open 2-ctry NK | Collateral value of US govt debt for banks | Dollar appreciates in stress; resolves paradox | Absolute-leaning | 4 |
| Bianchi-Bigio-Engel 2023 | Endo (settlement liquidity) | Open banking | Dollar interbank settlement/reserve liquidity | Dollar appreciates with funding risk | Absolute (endo level) | 4 |
| Gabaix-Maggiori 2015 | Exo (noise/flow) | Open 2-ctry | none (segmentation + Gamma) | Carry & UIP from financier risk-bearing | Relative (ccy premium) | 3 |
| Gourinchas-Ray-Vayanos 2022 | Exo (habitat) | Open 2-ctry cont. | none (habitat segmentation) | Joint term+ccy premia; carry maturity-dependent | Relative | 3 |
| Greenwood-Hanson-Stein-Sunderam 2023 | Exo (bond supply) | Open 2-ctry | none (relative qty of long bonds) | Quantity→FX & term premia; QE, CIP | Relative (cross-ccy qty) | 3 |
| Kekre-Lenel 2025 (GFC) | Exo (habitat) + endo interm. wealth | Open multi NK | Implicit via dollar funding demand | US MP→global risk price; spillovers fall if dollar demand falls | Relative risk premia | 4 |
| Farhi-Maggiori 2017 (QJE) | Exo reserve demand + strategic supply | Open multi | No-devaluation commitment + fiscal capacity | Hegemon issues safe debt at premium; Triffin fragility | Absolute (safety premium) | 4 |
| Choi-Kirpalani-Perez 2022 | Exo (BIU) + monopoly supply | Open 2-ctry dyn | Nonpecuniary value of US public debt | Underprovision; markup ~2/3 of conv. yield; privilege | Absolute (level+markup) | 4 |
| Coppola-Krishnamurthy-Xu 2023 | Endo (search/settlement) | Open multi-issuer | Liquidity of dollar safe short-term float | Dollar debt dominance; self-reinforcing | Absolute (liquidity premium) | 4 |
| Gopinath-Stein 2021 (QJE) | Endo (invoicing↔deposit↔bank) | Open multi | Dollar unit-of-account → safest deposit | Dollar banking; EM mismatch; privilege | Absolute-leaning | 3 |
| Caballero-Farhi-Gourinchas 2008 | Exo (asset-supply capacity) | Open multi | US capacity to manufacture (safe) assets | US deficits; low rates; US-tilted portfolios | Quantity/shortage | 2 |
| Caballero-Farhi 2018 (safety trap) | Exo safety demand | Closed/global ZLB | Safety per se | Safe-asset shortage → ZLB recession | Absolute/shortage | 2 |
| Caballero-Farhi-Gourinchas 2021 | Exo safety demand | Open multi ZLB | Safety; reserve ccy absorbs trap | Cross-border traps; reserve-ccy disadvantage | Absolute/shortage | 3 |
| He-Krishnamurthy-Milbradt 2019 | Endo (coordination+fiscal) | Open/global | Self-fulfilling safe-asset status + fiscal capacity | Determines identity of global safe asset | Absolute | 3 |
| Koijen-Yogo 2020 (global demand) | Estimated latent demand | Open 36-ctry | US-asset convenience as estimated characteristic | Conv. yield 2.15% US LT debt, 1.70% US equity | **Absolute (estimated level)** | 4 |
| Jiang-Richmond-Zhang 2022 | Estimated demand | Open multi | US-asset quality/convenience | Drivers of US NFA; privilege reversal; qty dimension | Absolute/quantity | 4 |
| Koijen-Yogo 2019 (JPE) | Estimated demand | Closed (US eq.) | n/a (method) | Inelastic demand drives prices | n/a | 1 |
| Clayton-Maggiori-Schreger 2024 | n/a (strategic) | Open multi | Dollar network centrality as leverage | Fragmentation; durability of dominance | n/a | 1 |
| Maggiori 2022 (Handbook) | survey | open | maps all sources | comprehensive on FX/carry/portfolios | rel & abs | 2 |

## Where the models AGREE

- **The dollar appreciates in bad times / flights to safety.** Every safe-asset
  model delivers a countercyclical dollar — whether via convenience demand
  (KL2024, JKL, DEW), intermediary wealth (Maggiori, KL2025), or settlement
  liquidity (BBE). This is the shared empirical anchor.
- **The US/dollar is special and earns a rent ("exorbitant privilege").** Either
  as a safety/convenience premium (JKL, Choi-Kirpalani-Perez, Farhi-Maggiori,
  Gopinath-Stein, Koijen-Yogo) or as a risk-bearing insurance premium (Maggiori,
  KL2024). All agree the US issues safe claims cheaply and/or holds risky assets.
- **Gross positions / portfolio composition matter, not just net.** Portfolio-
  balance and demand-system families both insist the exchange rate and premia are
  pinned by who-holds-what, decoupled from macro fundamentals (disconnect).

## Where they DISAGREE (the live tension big_cy targets)

- **Convenience vs. risk-bearing as the engine of privilege.** JKL2024 and the
  convenience family argue the privilege/paradox-resolution comes from *seigniorage
  on convenience demand* and does **not** require the US to be more risk-tolerant.
  Maggiori 2017 and KL2024's *risk-bearing-heterogeneity* channel argue US
  **greater risk tolerance** is what makes it the world insurer with a positive
  average NFA return. KL2024 actually uses *both* — and its heterogeneity-in-
  risk-bearing assumption (US more risk-tolerant) is the counterfactual lever
  big_cy wants to *replace* with a large absolute convenience yield.
- **The reserve-currency paradox.** Maggiori 2017 (and CFG 2021) frame a pure
  insurance story as paradoxical (predicts dollar depreciation in crises). KL2024,
  JKL2024, DEW, and BBE each claim to *resolve* it — but via different mechanisms
  (BIU safety + nominal rigidity; convenience seigniorage; collateral; settlement
  liquidity). big_cy must position against this cluster of competing resolutions.
- **Relative vs. absolute convenience.** The structural BIU/convenience models
  *calibrate to a relative basis* (Treasury or LIBOR basis), so the magnitude of
  the wedge is, by construction, the near-substitute gap. The demand-system and
  monopoly papers (Koijen-Yogo 2020: ~2% on US debt; Choi-Kirpalani-Perez:
  markup ~2/3) and DiTella et al. (Strand 1) put the *absolute* level much higher.
  **This gap is big_cy's central opening:** crank the KL2024 convenience to the
  absolute DiTella level and re-examine carry / home bias / NFA composition.
- **Endogenous vs. exogenous, and quantity.** Coppola-Krishnamurthy-Xu,
  He-Krishnamurthy-Milbradt, and Gopinath-Stein *microfound* dollar dominance
  (search, coordination, invoicing), implying the convenience level and the
  *quantity* of dollar safe assets are jointly endogenous. JKL2024 and
  Greenwood-Hanson-Stein-Sunderam make the quantity of (dollar) safe debt a key
  driver — directly relevant to big_cy's claim that there isn't enough dollar
  debt outstanding to satisfy US diversification demand.

## What is missing (frontier / gap for big_cy to occupy)

- No structural open-economy model has yet asked whether a **DiTella-magnitude
  absolute convenience yield** (2–4%), rather than a relative basis, can carry the
  weight currently borne by **risk-bearing heterogeneity** in KL2024. That is the
  unoccupied cell.
- The convenience models put the wedge on *dollar* bonds; few let *both* US and
  foreign households value safe bonds in *multiple* currencies (own + dollar),
  which is what big_cy needs to microfound home-currency bias symmetrically.
- The quantity-of-safe-debt channel (JKL2024, GHSS, CKX) and the absolute-level
  channel (Koijen-Yogo, DiTella) have not been combined in one quantitative
  GE model — a natural big_cy synthesis.
