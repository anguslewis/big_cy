# Strand 1 — Frontier Notes: The Absolute Level of Safe-Asset Convenience Yields

## The central tension

There are two fundamentally different objects in this literature, and the field has not
reconciled them:

1. **Relative convenience yields** — the gap between two near-substitute safe assets. These
   are well-measured, robust, and *small to moderate*:
   - Treasury vs. AAA/BAA corporates, CDs (Krishnamurthy-Vissing-Jorgensen 2012): **~73 bp**
     average 1926–2008; strongly supply-elastic.
   - Treasury vs. box-spread (genuinely risk-free, no-convenience) rate (van Binsbergen-
     Diamond-Grotteria 2022): **~40 bp** average, larger at short maturities, ~4x in crises.
   - Dollar vs. currency-hedged foreign govt bonds, the "Treasury basis" (Jiang-
     Krishnamurthy-Lustig 2021): a few tens of bp, time-varying; drives up to 41% of dollar.
   - Treasury vs. CDS-cleaned corporate (Mota): a sizable, time-varying *safety* component.

2. **The absolute wedge** — the gap between safe bond yields and the *true intertemporal
   price of consumption* (the zero-beta rate), or the true return to capital m. This is
   *large*:
   - Di Tella-Hébert-Kurlat-Wang (2025, JPE): real zero-beta rate **~8.3%** (8.55% ex-2020),
     SD ~9.3%, implying a **~7.6%/yr** spread over the real T-bill. Equivalently a ~6%
     equity-side "convenience"/wedge (Sec 7, their AAA-vs-equity comparison).
   - Reis (2021): the m–r gap (several %) as fiscal-relevant safe-asset/bubble premium.
   - Black-Jensen-Scholes (1972) / Black (1972): the flat SML with a high zero-beta
     intercept — the original average-level evidence the keystone makes time-varying.

The big_cy thesis lives in the *ratio* of these: the absolute wedge (~6–8%) is an order of
magnitude larger than the relative spreads (~0.4–0.75%).

## Side-by-side: competing estimates of "the level"

| Source | Object measured | Magnitude | Absolute or relative? | Identification basis |
|---|---|---|---|---|
| DiTella-HKW 2025 (JPE) | Zero-beta rate vs. T-bill | ~7.6%/yr (ZB ~8.3% real) | **Absolute** (intertemporal price) | Equity-only min-var zero-beta portfolio + linear-factor SDF + GMM |
| van Binsbergen-Diamond-Grotteria 2022 (JFE) | Treasury vs. box risk-free | **~40 bp** avg | **Absolute** (no-default risk-free) | Options put-call parity / box spreads |
| Krishnamurthy-Vissing-Jorgensen 2012 (JPE) | Treasury vs. corporates/CDs | **~73 bp** avg | Relative | Supply regressions on yield spreads |
| Mota (WP) | Treasury *safety* premium | large; R² 43–90% of spread var. | Relative | CDS-cleaned corporate bond |
| Jiang-Krishnamurthy-Lustig 2021 (JF) | Dollar vs. foreign safe (basis) | tens of bp; 41% of dollar | Relative | Treasury basis |
| Reis 2021 (BIS WP) | r vs. m (safe-asset/bubble premium) | several % | **Absolute-ish** (vs. MPK) | Calibrated GE bubble premium |
| Black-Jensen-Scholes 1972 | Zero-beta return level | high, > T-bill | **Absolute** (avg only) | Beta-sorted SML intercept |

The striking fact: the two genuinely *absolute* market-based estimates disagree by ~15x.
VBDG's box-vs-Treasury gap (~40 bp) measures convenience relative to a *no-default
collateralized derivative*; DiTella-HKW's ~6–8% measures the gap to the *zero-beta
intertemporal price*. These are different notions of "the true risk-free rate" — the box
rate is itself well below the zero-beta rate. Reconciling them is an open question.

## What is settled

- **There is a wedge.** Safe bond yields are not the intertemporal price; the SML is flat
  with a high intercept (BJS 1972 → DiTella-HKW 2025), and Treasuries are special relative
  to near substitutes (KVJ 2012; VBDG 2022; Mota). This is not seriously contested.
- **Relative convenience yields are small and supply-elastic.** KVJ's downward-sloping
  demand for safety is robust and widely replicated.
- **The wedge is time-varying and countercyclical.** Convenience/zero-beta wedge widens in
  crises (VBDG: ~4x; DiTella-HKW: zero-beta rate falls when recession likely / EBP high).
- **A high zero-beta return survives rich factor controls** (Lopez-Lira-Roussanov 2020;
  DiTella-HKW with 5-factor FF + bond factors).

## What is contested / open

1. **How big is the *absolute* level, really?** The 15x gap between VBDG (~40 bp) and
   DiTella-HKW (~6–8%) is unresolved. Are they measuring different wedges (no-default vs.
   full intertemporal-price), or is the zero-beta level overstated by an omitted high-priced
   factor (DiTella-HKW's own stated caveat; their level estimate is "less robust" than their
   Euler-equation test)?
2. **Interpretation of the zero-beta rate.** Is it the intertemporal price (DiTella-HKW), or
   the SDF of leverage-constrained investors (Frazzini-Pedersen 2014; DiTella-HKW Sec 6.4),
   or an artifact of multiple SDFs in incomplete markets (DiTella-HKW Sec 6.4, Chamberlain-
   Rothschild)? The keystone is careful that its *facts* hold regardless, but the *level*
   interpretation does not.
3. **Safety vs. liquidity decomposition of the absolute wedge.** KVJ and Mota decompose the
   *relative* spread; no one has cleanly decomposed the *absolute* ~6–8% wedge into safety,
   liquidity, segmentation, and risk-premium components.
4. **Source of the safe-rate puzzle.** Liquidity/incomplete-markets (DiTella-HKW 2024;
   Kaplan-Violante; Kocherlakota 2025; Bansal-Coleman; Lagos; Herrenbrueck) vs. segmentation
   of bond and stock SDFs (Fama-French; Collin-Dufresne et al.). Unsettled.
5. **Sign anomaly.** Occasionally the estimated zero-beta rate falls *below* the T-bill
   (also Campbell-Shiller 2001 for the market; Greenwood-Hanson 2013 for high-yield). Is
   this overfitting/linear extrapolation, or genuine episodes where "riskier/less-liquid"
   assets earn less than Treasuries? Open (DiTella-HKW call it an open question).
6. **The "different safe assets, different yields" problem.** Any absolute-level theory must
   also explain why distinct safe assets (T-bills, AAA, box, OIS) have different yields —
   i.e. why the convenience wedge is asset-specific (Diamond 2020; Angeletos et al. 2023).

## Implications flagged for the absolute-CY level (relevant downstream)

- **Seigniorage / fiscal capacity** is "much higher than traditionally thought" if the wedge
  is the m–r or zero-beta–r gap rather than the AAA-Treasury gap (Reis 2021; Kocherlakota
  2025; DiTella-HKW Sec 7).
- **Cost of capital / capital structure**: banks' reluctance to fund with equity is rational
  if the zero-beta rate (true cost of riskless equity-funded investment) is ~6% over
  Treasuries (Baker-Wurgler 2015; DiTella-HKW Sec 7). Lower-bound estimates from leverable
  arbitrages (CIP, box) give a bank-perceived riskless equity return of 4%+ over safe rates
  (Boyarchenko et al. 2018; Du-Tepper-Verdelhan 2018; VBDG 2019).
- **International finance**: a large absolute dollar-convenience yield (vs. the small
  *relative* Treasury basis of Jiang-Krishnamurthy-Lustig) is the bridge to big_cy's
  home-bias / non-carry-trader puzzles — but the *absolute* dollar wedge has not been
  estimated directly (only the relative basis).
