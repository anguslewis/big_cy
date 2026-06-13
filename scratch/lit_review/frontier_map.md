# big_cy — Frontier Map & Taxonomy

*Deliverable for T1a. Companion to `synthesis.md` (the essay) and
`annotated_bibliography.md` (the entries). Neutral survey: this map classifies
the field; it does not argue for a mechanism.*

The literature on convenience yields and safe-asset demand can be organized along
**six axes**. The first three are about *what the convenience yield is and how
it is measured*; the last three are about *how it is modeled in (international)
general equilibrium*. The same paper often sits at a different point on each
axis, so the taxonomy is multi-dimensional, not a single tree.

---

## The six classifying axes

| Axis | Pole A | Pole B | Why it matters |
|---|---|---|---|
| **1. Measured object** | **Relative** — a *gap* between two near-substitute safe assets (CIP basis; Treasury–AAA; Treasury–foreign-sovereign; nominal–TIPS) | **Absolute** — a *level* vs. a convenience-free benchmark built from risky assets (option box spreads; zero-beta equity rate) | Relative measures net out convenience common to both legs; only absolute measures recover the level. This is big_cy's organizing distinction. |
| **2. Magnitude** | **Small** (~0.4–0.75%: KVJ ~73bp; vBDG ~40bp; CIP ~25bp) | **Large** (~2–4%: DiTella zero-beta; Koijen-Yogo ~2%) | The headline disagreement. Small estimates are all relative or option-based; large estimates are equity-based or demand-system. |
| **3. Source of convenience** | **Liquidity / money-services** (transactions, near-money, search) | **Safety / collateral** (information-insensitivity, pledgeability, default-remoteness) | KVJ decompose the yield into both; most microfoundations price only one leg. Regulatory / reserve-demand (HQLA, T-bill/repo) is a third, institutional, source (`davernasvandeweyer2024treasury`). |
| **4. Demand origin** | **Exogenous** — preference shock / bonds-in-utility / reduced-form clientele or noise | **Endogenous** — collateral constraint, intermediary balance sheet, search/settlement, coordination | Determines policy-invariance and whether the wedge can be "cranked." big_cy's KL template is exogenous-BiU. |
| **5. Economy** | **Closed** / single-asset | **Open** / multi-currency, two-country | The international target facts (home bias, carry, NFA composition) require open-economy, multi-currency structure. |
| **6. What pins the wedge in GE** | **Quantity** — supply/scarcity of safe assets drives the premium | **Price/preference** — a demand-curve intercept or risk price drives it | "Not enough dollar debt" (big_cy) is a quantity story; a large absolute BiU level is a price/preference story. Several models combine them. |

---

## Master classification table

Proximity (Prox) = closeness to big_cy's question (5 = directly competes/keystone;
1 = background). "Rel/Abs" = the convenience notion the paper *works with*
(axis 1). "Demand" = axis 4. "Src" = axis 3 (liq/safe/coll/search/reg/—).

### Strand 1 — Zero-beta & the absolute level

| Paper | Object | Rel/Abs | Magnitude | Identifies level? | Prox |
|---|---|---|---|---|---|
| `ditella2025zerobeta` | Zero-beta rate vs. T-bill | **Abs** | **~2–4% (ZB ~8.3% real; ~7.6%/yr over T-bill)** | Yes — equity-only intertemporal price | 5 |
| `vanbinsbergendiamondgrotteria2022riskfree` | Treasury vs. option box rate | **Abs** | **~40 bp** (≤3Y), ~4× in crisis | Yes — option-implied risk-free | 5 |
| `kvj2012aggregate` | Treasury vs. corporates/CDs | Rel | ~73 bp avg | No — Treasury-minus-near-safe gap | 5 |
| `reis2021constraint` | r vs. m (MPK) wedge | Abs-ish | several % | Conceptual (calibrated bubble premium) | 4 |
| `mota2021corporate` | Treasury *safety* premium (CDS-cleaned) | Rel | large; R²≈60% of spread var. | No (lower bound on safety leg) | 3 |
| `black1972capital` / `blackjensenscholes1972capm` | Zero-beta return (flat SML) | Abs (avg) | high, > T-bill | Avg level only | 4 / 3 |
| `frazzini2014betting` | Zero-beta rate via leverage constraints | Abs | high | Alt. interpretation (constraint, not consumption) | 3 |
| `mehra1985equity`, `weil1989equity` | Equity-premium / risk-free-rate puzzles | — | ~6% premium | Foundational (the wedge to absorb) | 2 |
| `ditella2024aggregateeuler` | Liquid/illiquid HA microfoundation | Abs | — | Rationalizes the level | 4 |
| `diamond2020safety`, `angeletoscollarddellas2023public`, `kocherlakota2025public` | Why the wedge is asset-specific / fiscal | Abs | — | Theory of multiple risk-free rates | 4/3/3 |

### Strand 2 — International models with safe-asset / convenience demand

| Paper | Demand | Economy | Src of $ convenience | Implies for carry / home bias / dollar | Rel/Abs | Prox |
|---|---|---|---|---|---|---|
| `kekrelenel2024flight` | **Exo (BiU safety shock)** + supply | Open 2-ctry NK | Nonpecuniary value of safe $ bonds | US long-capital/short-$; $ = neg-beta hedge; +avg NFA return from **risk-tolerance gap** | Rel (Treasury basis) | **5** |
| `jiangkrishnamurthylustig2024dollar` | Exo demand (qty-dep) + endo supply | Open US/RoW | Foreign convenience on safe $ claims | US issues cheap safe $ debt; GFC = $ cycle; resolves paradox **w/o risk tolerance** | Rel→abs framing | **5** |
| `jiangkrishnamurthylustig2021foreign` | Exo (measured via basis) | Open (kernel) | Foreign $ safety; ~90% is $ per se; ~2% conv. | Convenience drives $ level & disconnect | Rel (basis) | **5** |
| `maggiori2017reservecurrencies` | Endo (intermediary asymmetry) | Open 2-ctry | US financial-development depth | US = world insurer; raises the **reserve-currency paradox** | Risk price | 4 |
| `devereuxenglewu2023collateral` | Endo (collateral) | Open 2-ctry NK | Collateral value of US debt for banks | $ appreciates in stress; resolves paradox | Abs-leaning | 4 |
| `bianchibigioengel2023scrambling` | Endo (settlement liquidity) | Open banking | $ interbank settlement liquidity | $ appreciates with funding risk | Abs (endo) | 4 |
| `gabaixmaggiori2015liquidity` | Exo (noise/flow) | Open 2-ctry | none (segmentation + Γ) | Carry & UIP from financier risk-bearing | Rel | 3 |
| `gourinchasrayvayanos2022habitat` | Exo (habitat) | Open 2-ctry cts | none (habitat) | Joint term+ccy premia; carry maturity-dependent | Rel | 3 |
| `greenwoodhansonsteinsunderam2023quantity` | Exo (bond supply) | Open 2-ctry | none (relative qty of long bonds) | Quantity → FX & term premia; QE; CIP | Rel (cross-ccy qty) | 3 |
| `kekrelenel2025gfc` | Exo habitat + endo interm. wealth | Open multi NK | Implicit via $ funding demand | US MP → global risk price; spillovers fall if $ demand falls | Risk premia | 4 |
| `farhimaggiori2017ims` | Exo reserve demand + strategic supply | Open multi | No-devaluation commitment + fiscal cap. | Hegemon issues safe debt at a premium; Triffin fragility | **Abs (safety premium)** | 4 |
| `choikirpalaniperez2022marketpower` | Exo (BiU) + **monopoly supply** | Open 2-ctry dyn | Nonpecuniary value of US public debt | Underprovision; **markup ~2/3** of conv. yield; privilege | **Abs (level+markup)** | 4 |
| `coppolakrishnamurthyxu2023liquidity` | Endo (search/settlement) | Open multi-issuer | Liquidity of $ safe short-term float | $ debt dominance; self-reinforcing | Abs (liq premium) | 4 |
| `gopinathstein2021banking` | Endo (invoicing↔deposit↔bank) | Open multi | $ unit-of-account → safest deposit | $ banking; EM mismatch; privilege | Abs-leaning | 3 |
| `caballerofarhigourinchas2008imbalances` | Exo (asset-supply capacity) | Open multi | US capacity to manufacture safe assets | US deficits; low rates; US-tilted portfolios | Quantity/shortage | 2 |
| `caballerofarhi2018safetytrap` | Exo safety demand | Closed/global ZLB | Safety per se | Safe-asset shortage → ZLB recession | Abs/shortage | 2 |
| `caballerofarhigourinchas2021policywars` | Exo safety demand | Open multi ZLB | Safety; reserve ccy absorbs trap | Cross-border traps; reserve-ccy disadvantage | Abs/shortage | 3 |
| `hekrishnamurthymilbradt2019safeasset` | Endo (coordination + fiscal) | Open/global | Self-fulfilling safe status + fiscal cap. | Determines *which* asset is the global safe asset | Abs | 3 |
| `koijenyogo2020globaldemand` | Estimated latent demand | Open 36-ctry | US-asset convenience as characteristic | **Conv. ≈ 2.15% US LT debt, 1.70% US equity** | **Abs (estimated)** | 4 |
| `jiangrichmondzhang2022portfolio` | Estimated demand | Open multi | US-asset quality/convenience | Drivers of US NFA; privilege reversal; **quantity dimension** | Abs/quantity | 4 |
| `koijenyogo2019demandsystem` | Estimated demand | Closed (US eq.) | — (method) | Inelastic demand drives prices | — | 1 |
| `kekrelenel2022redistribution` | — (risk-bearing heterogeneity) | Closed HANK | — | The MPR/risk-bearing engine of the KL toolkit | — | 2 |
| `maggiori2022handbook` | survey | open | maps all sources | comprehensive on FX/carry/portfolios | rel & abs | 2 |
| `claytonmaggiorischreger2024coercion` | — (strategic) | Open multi | $ network centrality as leverage | Durability of $ dominance | — | 1 |
| `valchev2020bond` *(added)* | Exo (conv-in-Euler) | Open | Bond convenience yield | Convenience drives FX; UIP reversal — closest model to big_cy | Rel→abs | **5** |
| `gourinchasrey2007banker` *(added)* | — (NFA accounting) | Open | — | US long-equity/short-debt; exorbitant privilege — the target facts | — | 4 |
| `bocolalorenzoni2020dollarization` *(added)* | Endo (savers' insurance) | Open | — (FX-risk premium) | Self-fulfilling liability dollarization → currency-composition bias | — | 4 |
| `englewu2023liquidity` *(added)* | Exo (liquidity wedge) | Open | $ liquidity/convenience | Liquidity yield priced into UIP; empirical | Rel | 4 |
| `caballerokrishnamurthy2009fragility` *(added)* | Exo safe demand | Open | US safe-asset supply capacity | Funds risky US intermediation; fragility | Quantity/shortage | 2 |

### Strand 3 — Mechanisms & microfoundations

| Paper | Source | RF or Micro | Implication for level / cyclicality | Prox |
|---|---|---|---|---|
| `kvj2012aggregate` | liquidity + safety | RF estimation (micro-motivated) | ~70bp avg; **downward-sloping in supply**; rises in crises | 5 |
| `nagel2016liquidity` | liquidity | RF/micro bridge | premium moves **~1:1 with the short nominal rate** | 5 |
| `sidrauski1967rational` | liquidity (MIU) | RF origin | conv. = marginal liquidity service, decreasing in real qty held | 3 |
| `fisher2015structural` | safety+liquidity | RF | DSGE "risk-premium" shock = safe-asset demand shock; countercyclical | 3 |
| `mianstraubsufi2025goldilocks` | liquidity/safe demand | RF (macro GE) | conv. falls in public-debt supply; pins optimal deficit | 3 |
| `michaillatsaez2022economical` | wealth-in-utility | RF | store-of-value premium; permanent ZLB | 2 |
| `woodford1990public`, `holmstromtirole1998private`, `holmstromtirole2011inside`, `aiyagarimcgrattan1998optimum`, `angeletoscollarddellas2023public` | liquidity/collateral | Micro | public debt relaxes private constraints; premium high when liquidity scarce | 4/4/—/3/3 |
| `gortonpennacchi1990financial`, `danggortonholmstrom2017banks` | safety / info-insensitivity | Micro | conv. = value of no adverse selection; collapses when info is produced (runs) | 4/3 |
| `lagoswright2005unified`, `lagos2010asset`, `vayanosweill2008search` | search/OTC liquidity | Micro | premium from trading frictions; the deep "what BiU means" | 3/3/4 |
| `krishnamurthy2002bond` | liquidity/search (evidence) | RF empirical | on/off-the-run ≈ 6–12 bp; widens in liquidity crises | 4 |
| `fleckensteinlongstafflustig2014tips` | safety/liquidity (evidence) | RF empirical | TIPS–Treasury ≈ 55 bp (to 200+); slow-moving-capital limit | 3 |
| `stein2012monetary`, `ghs2015comparative`, `kvj2015impact` | regulatory / money-creation | Micro / bridge | private safety supply responds to the premium; crowding-out | 3/4/4 |
| `diamonddybvig1983bank` | liquidity-insurance | Micro | banks create demandable claims; convenience fragile (runs) | 2 |
| `caballerofarhigourinchas2017safe` | safety (macro framing) | survey | shortage → low r*, safety trap | 3 |
| `davernasvandeweyer2024treasury` *(added)* | **regulatory/reserve plumbing** | Micro | bill scarcity + reserve/repo demand price short safe assets | 3 |
| `gortonlewellenmetrick2012safe` *(added)* | safety (quantity) | Empirical | safe-asset share **stable ~33%** of total assets | 3 |
| `bansalcoleman1996monetary` *(added)* | liquidity/money-services | Micro | money services lower safe-bond returns; resolves EP/RF puzzles | 2 |

### Strand 4 — Relative-vs-absolute measurement

| Measure (paper) | Two legs | Rel/Abs | Magnitude | Prox |
|---|---|---|---|---|
| CIP basis (`dutepperverdelhan2018cip`, `duschreger2020handbook`) | USD cash vs. FX-swap-implied USD | **Rel** (x-ccy) | ~24bp 3M, ~27bp 5Y; JPY 5Y ~−90bp | 5/5 |
| US Treasury Premium (`duimschreger2018treasury`) | US Treasury vs. hedged foreign govt | **Rel** | ~21bp pre-GFC; **~−8bp post** | 5 |
| Treasury basis (`jiangkrishnamurthylustig2021foreign`) | US govt vs. hedged foreign govt | **Rel** | up to 41% of $ variation | 5 |
| AAA–Treasury (`kvj2012aggregate`) | Treasury vs. AAA/CP/deposit | **Rel** | ~73bp avg | 5 |
| On/off-the-run (`krishnamurthy2002bond`) | new vs. old 30Y Treasury | **Rel** (liq only) | ~12→3bp | 3 |
| Treasury–Refcorp (`longstaff2004flight`) | Treasury vs. guaranteed agency | **Rel** (liq, credit-matched) | up to ~15% of value | 4 |
| TIPS–Treasury (`fleckensteinlongstafflustig2014tips`) | nominal Treasury vs. swapped TIPS | **Rel** (within issuer) | ~54.5bp (to >200) | 4 |
| Negative swap spread (`klinglersundaresan2019negative`) | swap rate vs. Treasury yield | **Rel** (contaminated) | 30Y <0 post-2008 | 3 |
| **Box risk-free rate (`vanbinsbergendiamondgrotteria2022riskfree`)** | Treasury vs. option-implied riskless | **ABS** | **~40bp**; ~4× crisis | 5 |
| **Option-implied, per currency (`diamondvantassel2024riskfree`)** | safe bond vs. option-implied riskless, by ccy | **ABS (per ccy)** | rises ~linearly w/ rate; US 5th of G11 | 5 |
| **Zero-beta rate (`ditella2025zerobeta`)** | safe bond vs. zero-beta equity portfolio | **ABS** | **~2–4%** | 5 |
| `acharyalaarits2023when` | Treasury conv. via stock-bond covariance | Abs (conditional) | low when Treasury-equity cov. high | 4 |
| `nenova2025global` | demand elasticities / substitution | (segmentation) | US Treasuries substituted *globally*; Bunds *regionally* | 4 |
| `duschreger2016local`, `engel2016exchange`, `avdjievdukochshin2019dollar`, `andersonduschlusche2021arbitrage` | local-ccy sovereign / UIP wedge / $-leverage / arb capital | Rel | — | 3/4/3/2 |
| `duhebertli2023intermediary` *(added)* | Treasury curve vs. intermediary balance sheet | Rel (interm.-priced) | dealer net-short→net-long post-GFC | 3 |
| `englewu2023liquidity` *(added)* | $ liquidity yield in UIP | Rel | liquidity/convenience priced into FX | 4 |

---

## The frontier: where the cells are empty

Reading the table down the **demand-origin × measured-object** plane exposes the
gap big_cy occupies:

```
                    EXOGENOUS demand              ENDOGENOUS demand
                    (BiU / preference)            (collateral/interm./search)
  RELATIVE       │ KL2024, JKL2021/24,         │ Maggiori2017, GRV2022,
  convenience    │ Gabaix-Maggiori,            │ Kekre-Lenel2025
  (basis-cal.)   │ GHSS2023                    │
  ───────────────┼─────────────────────────────┼──────────────────────────────
  ABSOLUTE       │  ►►  big_cy's cell  ◄◄       │ Devereux-Engel-Wu,
  convenience    │ Choi-Kirpalani-Perez,       │ Bianchi-Bigio-Engel,
  (level)        │ Koijen-Yogo (estimated),    │ He-Krishnamurthy-Milbradt,
                 │ Farhi-Maggiori (premium)    │ Coppola-Krishnamurthy-Xu
```

Three observations the map makes precise:

1. **The structural convenience models calibrate to a *relative* basis.** KL2024
   and JKL discipline the wedge with the Treasury basis (axis-1 = relative),
   even though the object is conceptually a level. No structural open-economy
   model has asked whether a **DiTella-magnitude *absolute* wedge** (axis-1 =
   absolute, axis-2 = large) can carry the weight currently borne by
   **risk-bearing heterogeneity** in KL2024. That cell — exogenous-BiU,
   absolute-large, open-economy — is essentially empty.

2. **The "absolute level" papers are mostly closed-economy or
   single-currency.** DiTella (one aggregate wedge), vBDG (US Treasury),
   Koijen-Yogo (a US-asset characteristic). Diamond-Van Tassel is the exception
   (per-currency absolute levels) and is therefore the natural disciplining
   anchor for a two-country, two-currency absolute-convenience model.

3. **Quantity and level have not been combined.** The "quantity of safe debt"
   channel (JKL2024, GHSS2023, Coppola-Krishnamurthy-Xu, Jiang-Richmond-Zhang)
   and the "absolute level" channel (DiTella, Koijen-Yogo) live in separate
   literatures. big_cy's "not enough dollar debt to satisfy US diversification
   demand" reinterpretation sits at their intersection.

See `synthesis.md` for the narrative tying these together, the points of
agreement/disagreement, and the open questions.
