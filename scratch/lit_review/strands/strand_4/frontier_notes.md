# Strand 4 — Frontier notes: what is relative, what is absolute, and the bridge to Strand 1

## The taxonomy in one table

| Measure | Two legs | Relative or Absolute | What it identifies | What it does NOT identify | Reported magnitude |
|---|---|---|---|---|---|
| Cross-currency basis (DTV; Du-Schreger; Avdjiev et al.) | USD cash rate vs. FX-swap-implied USD rate | **Relative** (cross-currency) | USD-vs-foreign funding/convenience gap | Absolute USD convenience | ~24 bps (3M), ~27 bps (5Y); JPY 5Y ~−90 bps |
| U.S. Treasury Premium (Du-Im-Schreger) | US Treasury vs. hedged foreign govt bond | **Relative** (US vs. foreign sovereign) | US-minus-foreign safe convenience | Common (absolute) safe convenience | ~21 bps pre-GFC; ~−8 bps post |
| Treasury basis (JKL 2021) | US govt vs. hedged foreign govt | **Relative** | Foreign-investor convenience on Treasuries | Absolute level | Explains up to 41% of $ variation |
| AAA–Treasury spread (KVJ 2012) | Treasury vs. AAA corporate/CP/deposit | **Relative** (Treasury vs. near-safe) | Treasury-over-AAA convenience | Convenience shared by AAA + Treasury | ~73 bps avg (1926–2008) |
| On/off-the-run (Krishnamurthy 2002) | New vs. old 30Y Treasury | **Relative** (liquidity only) | Pure liquidity convenience | Safety; absolute level | ~12→3 bps over cycle |
| Treasury–Refcorp (Longstaff 2004) | Treasury vs. Treasury-guaranteed agency | **Relative** (liquidity, credit-matched) | Flight-to-liquidity premium | Safety; absolute level | up to ~15% of bond value |
| TIPS–Treasury (FLL 2014) | Nominal Treasury vs. inflation-swapped TIPS | **Relative** (within issuer) | Nominal-over-TIPS convenience | Common Treasury convenience; absolute | ~54.5 bps avg, up to >200 bps |
| Negative swap spread (Klingler-Sundaresan) | Swap rate vs. Treasury yield | **Relative** (contaminated) | Treasury richness + pension/dealer frictions | Clean convenience level | 30Y spread persistently <0 post-2008 |
| **Box-spread risk-free rate (vBDG 2022)** | Treasury vs. option-implied riskless rate | **ABSOLUTE (level)** | Level of Treasury convenience | Long maturities (>~3Y) | **~40 bps**; ~4x in crisis |
| **Option-implied rates across currencies (Diamond-Van Tassel 2024)** | Safe bond vs. option-implied riskless, per currency | **ABSOLUTE (level), per currency** | Per-currency convenience level | Macro/structural intertemporal price | rises ~linearly w/ rate level; US 5th |
| **Zero-beta rate (DiTella et al. 2025)** | Safe bond vs. zero-beta equity portfolio | **ABSOLUTE (level)** | Absolute wedge / "safe-rate puzzle" | Which friction; instrument split | **~2–4% (200–400 bps)** |

## The core conceptual point

**Relative measures net out the convenience common to both legs.** A basis or a
spread is a difference. If both legs share a large common convenience component,
the spread can be small — even zero or negative — while the *absolute* level on
each leg is large. Three facts make this concrete and non-hypothetical:

1. The **US Treasury Premium** (Treasury vs. hedged foreign sovereign) is ~21 bps
   pre-GFC and **negative (~−8 bps) post-GFC** — yet no one believes US Treasury
   convenience went to zero. Foreign safe govt bonds simply carry comparable
   convenience, so the *difference* is tiny.
2. **Diamond-Van Tassel (2024)** confirm this directly with absolute per-currency
   measures: convenience *levels* rise with each country's rate level and are
   substantial everywhere, while the **US-minus-foreign difference does not grow
   in crises** even as levels do. Absolute and relative move differently.
3. The **AAA–Treasury spread** measures Treasury-over-AAA, not
   Treasury-over-(convenience-free rate). If AAA corporates themselves are
   money-like, KVJ's ~73 bps understates the absolute level.

## What CAN recover an absolute level

Only measures that price the safe asset against a **convenience-free benchmark
constructed from risky assets**:

- **Option box spreads** (van Binsbergen-Diamond-Grotteria; Diamond-Van Tassel):
  a synthetic riskless rate from put-call parity carries no Treasury-style
  convenience. ⇒ absolute Treasury convenience ≈ **40 bps** (US), maturities ≤ ~3Y.
- **Zero-beta equity portfolio** (DiTella-Hébert-Kurlat-Wang): an equity-only
  intertemporal price, model-free about the wedge. ⇒ absolute wedge ≈ **2–4%**.

## The bridge to Strand 1 (and the open tension)

Strand 1's anchor (the zero-beta rate) is the *equity-based absolute* measure.
This strand supplies the *menu of alternatives* and the key unresolved gap:

> **The two absolute measures disagree by an order of magnitude.**
> Option-box methods say ~40 bps; the zero-beta rate says 200–400 bps.

Possible reconciliations big_cy should engage:
- **Maturity / horizon:** box spreads are short (≤3Y) financial rates; the
  zero-beta rate is a long-horizon intertemporal price. Convenience may rise with
  horizon and the equity wedge may include term/risk components beyond pure
  convenience.
- **Object measured:** box spreads price a *tradable* riskless payoff; the
  zero-beta rate prices the *consumption-Euler* intertemporal margin. The wedge
  between them is exactly what DiTella et al. leave model-free, and what
  Acharya-Laarits decompose state-by-state (hedging/negative-beta channel).
- **Composition:** the zero-beta wedge may bundle convenience with a broader
  "missing-risk-compensation" or preference component; the box-spread convenience
  is narrowly the Treasury-vs-riskless gap.

**Implication for big_cy's calibration.** When big_cy "cranks the absolute
convenience to DiTella magnitudes (2–4%)," it is choosing the *equity-based*
absolute benchmark over the *option-based* one. The relative measures (CIP basis,
Treasury premium, AAA spread) are then *consistent with* a large common absolute
level precisely because they difference it out — the small/negative measured gaps
are not evidence against a large absolute level. Diamond-Van Tassel's finding
(absolute levels large and rate-linked; cross-country gap small and stable) is the
single most useful empirical anchor: it operationalizes the relative-vs-absolute
wedge in international data and is the natural place to discipline the model's
US-vs-RoW convenience parameters.

**Caveat to flag for the critic/strategist.** big_cy's strongest empirical claim
(absolute >> relative) leans on the equity-based number being the "right"
intertemporal price. The option-based literature is the principal countervailing
evidence (the level is modest, ~40 bps). A neutral survey must present both; the
model can be framed as exploring the *upper-bound* (equity-based) scenario.
