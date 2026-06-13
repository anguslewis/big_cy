# Literature Review — librarian-critic
**Date:** 2026-06-13
**Score:** 84/100

## Summary judgment

This is a strong, genuinely broad and neutral survey — exactly what was
commissioned. The four-strand structure is coherent, the relative-vs-absolute
organizing axis is applied consistently and correctly, the frontier map's
demand-origin × measured-object plane is a real contribution (not decoration), and
the "inclusion against interest" (vBDG ~40 bp as the principal countervailing
evidence to DiTella) is handled with unusual intellectual honesty. The
categorization (Rel/Abs, exo/endo, RF/microfounded, Prox) is accurate. The bib is
clean and well-deduplicated.

The deductions are almost entirely **coverage gaps**: a deliberately broad, neutral
map has a higher bar for completeness than a focused paper's lit review, and a
handful of papers that the strand frontier notes *themselves name* never made it
into the master bib, plus several adjacent literatures the brief flagged as
candidates are genuinely load-bearing and missing. I do **not** penalize breadth —
I penalize the specific holes below.

## Issues Found

### A. Coverage gaps (the main deductions)

**A1. Valchev (bond convenience yields & exchange rates) — MISSING. Severity: high.**
Rosen Valchev, "Bond Convenience Yields and Exchange Rate Dynamics" (AEJ:
Macroeconomics, 2020) is arguably the single most on-point omission for big_cy. It
is a structural model where a time-varying *bond convenience yield* drives
exchange-rate dynamics and produces the UIP/forward-premium reversal — the exact
mechanism class (convenience-in-the-Euler-equation → FX) the project sits in, and a
direct intellectual neighbor of JKL2021 and Engel (2016), both included. **-20.**

**A2. Gourinchas-Rey exorbitant-privilege accounting — MISSING. Severity: high.**
The synthesis invokes "exorbitant privilege" repeatedly and the project's three
target facts are NFA-composition facts, yet the foundational accounting —
Gourinchas & Rey, "From World Banker to World Venture Capitalist" (2007) — is
absent. The seminal empirical statement of the very facts big_cy and `rier` aim to
explain. **-15.**

**A3. HQLA/LCR & reserve-demand plumbing — UNDER-COVERED. Severity: medium.**
The annotated bib explicitly flags this as a "coverage extension" and frontier-map
axis 3 names "Regulatory (HQLA) demand … a third, institutional, source" — but no
HQLA/regulatory-demand paper is actually in the bib. At minimum d'Avernas &
Vandeweyer (reserves/repo) and Acharya-Rajan belong; Diamond-Jiang-Schrimpf is
optional. **-10.**

**A4. Other genuinely-relevant omissions (cumulative). Severity: medium.**
Gorton-Lewellen-Metrick, "The Safe-Asset Share" (AER P&P 2012); Caballero-Krishnamurthy,
"Global Imbalances and Financial Fragility" (AER P&P 2009); Bocola-Lorenzoni,
"Financial Crises, Dollarization, and Lending of Last Resort" (AER 2020);
Engel-Wu, "Liquidity and Exchange Rates" (ReStud 2023; distinct from the included
DEW2023 "Collateral Advantage"); Du-Hébert-Li, "Intermediary Balance Sheets and the
Treasury Yield Curve" (JFE 2023); Bansal-Coleman (1996, JPE) — *cited by name in
strand_1 frontier notes yet absent from the bib*. Lenel-Piazzesi-Schneider /
Piazzesi-Rogers-Schneider optional. **-10** (one cumulative deduction).

**Genuinely out of scope (do NOT add):** Liao (credit-currency basis; `lx1`
territory); Diamond-Van Tassel already in; further Koijen-Yogo follow-ons;
Gertler-Karadi / He-Krishnamurthy (HKM2019 already represents that line).

### B. Recency (2023–2026)
Well-covered to 2025 (Mian-Straub-Sufi 2025, Kekre-Lenel GFC Nov 2025, Nenova 2025,
Kocherlakota 2025, Clayton-Maggiori-Schreger 2024, Diamond-Van Tassel 2025). No
additional deduction beyond A4.

### C. Categorization quality
Largely correct and consistent. Rel/Abs, exo/endo, and proximity tags spot-check
out. **One soft inconsistency:** `farhimaggiori2017ims` is tagged "Abs (safety
premium)" — an *endogenous monopoly rent*, arguably a price/preference object,
slightly overloading the "Abs" column alongside Koijen-Yogo's estimated level. **-5.**

### D. Journal quality / published-vs-WP balance
Healthy. ~66 entries (initial sweep estimate), ~30% working papers — well under the
50% threshold. Strong top-5 + top-finance coverage + two Handbook chapters. The
"verify before use" header is the right hygiene. No deduction.

### E. BibTeX correctness
Clean merge; no duplicate keys survived. Minor field items (not scored): fix
`farhimaggiori2017ims` year-vs-volume (QJE 133(1) is **2018**); verify
`michaillatsaez2022economical` venue; `koijenyogo2020globaldemand` now published
(update venue); `diamondvantassel2024riskfree` 2024-vs-2025 already self-flagged.

## Score Breakdown
- Starting: 100
- A1 Valchev missing: **-20**
- A2 Gourinchas-Rey privilege accounting missing: **-15**
- A3 HQLA/LCR/reserve-demand named but unpopulated: **-10**
- A4 Cumulative adjacent omissions: **-10**
- C Single soft categorization inconsistency: **-5**
- D / E: 0
- **Final: 84/100 — strike 1.**

## Top fixes (highest marginal value)
1. Add Valchev (2020, AEJ:Macro) to Strand 2.
2. Add Gourinchas-Rey privilege accounting (2007).
3. Populate the HQLA/regulatory source (d'Avernas-Vandeweyer at minimum).
4. Reconcile Bansal-Coleman (cited but absent); fix `farhimaggiori2017ims`
   year-vs-volume and verify `michaillatsaez2022economical` venue.

No three-strikes escalation (strike 1). This is a scope/coverage matter the
Librarian can address by adding the named items; no breadth-vs-depth ruling needed
from the user.

*(This report was returned in full by the librarian-critic; reproduced here as the
written record. The resolution below documents the fixes applied.)*

---

## Resolution (big_cy-lit, 2026-06-13)

All four top fixes addressed; 9 web-verified entries added (initial 79 → **88**
unique references). De-duplication re-checked: **0 duplicate citekeys, 0 orphan
citations**.

- **A1 Valchev (2020, AEJ:Macro 12(2):124–166)** — added to Strand 2 as the
  closest existing structural convenience-FX model. ✅
- **A2 Gourinchas-Rey (2007, "World Banker to World Venture Capitalist," in
  Clarida ed., pp. 11–66)** — added to Strand 2 as the privilege/NFA-accounting
  foundation of the target facts. ✅
- **A3 HQLA/regulatory cell populated** — d'Avernas-Vandeweyer (2024, JF
  79(6):4083–4141) added to Strand 3; frontier-map axis-3 regulatory source now
  cites it. ✅
- **A4 adjacent omissions** — added Bocola-Lorenzoni (2020, AER 110(8)),
  Engel-Wu "Liquidity and Exchange Rates" (2023, ReStud 90(5); *distinct* from
  the already-included Devereux-Engel-Wu "Collateral Advantage"),
  Du-Hébert-Li (2023, JFE 150(3)), Gorton-Lewellen-Metrick (2012, AER P&P
  102(3):101–106), Caballero-Krishnamurthy (2009, AER P&P 99(2):584–588), and
  Bansal-Coleman (1996, JPE 104(6):1135–1171; previously cited in the Strand-1
  frontier notes but absent from the bib). ✅
- **BibTeX corrections** — `farhimaggiori2017ims` year → **2018** (QJE 133(1));
  `koijenyogo2020globaldemand` updated to its publication (JPE 2026, 134(8)).
  `michaillatsaez2022economical` re-verified as **correct** (OEP 74(2):382–411,
  2022) — no change needed. ✅
- **Categorization soft-fix** — `farhimaggiori2017ims` retagged "Abs level *as an
  endogenous monopoly rent*" to distinguish it from Koijen-Yogo's estimated level. ✅
- **Deliberately NOT added** (per critic's own "out of scope" guidance): Liao
  (credit-currency basis; `lx1` territory), further Koijen-Yogo follow-ons,
  Gertler-Karadi (HKM2019 already represents the He-Krishnamurthy intermediary line).
- **Still flagged to verify before manuscript use** (header of `references.bib`):
  `andersonduschlusche2021arbitrage` (NBER WP 28658 vs. a JF version),
  `diamondvantassel2024riskfree` (JF forthcoming/2025, cited as 2024), and the
  WP-status items (kocherlakota, mota×2, diamond2020, lopezlira).

Resulting state: the score's coverage deductions (A1–A4 = −55, plus the soft
categorization −5) are addressed; remaining open items are verify-before-use
citations, not coverage gaps.
