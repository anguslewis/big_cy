# big_cy — project instructions

These project-wide rules apply to **every** agent working in this repo, on top
of the global conventions in `~/code/claude/`.

## Research vision

**Core thesis.** Standard convenience-yield estimates are *relative* — dollar
vs. euro via the CIP basis, Treasuries vs. AAA corporates via spreads — so they
measure the *gap* between near-substitutes, not the *absolute* level of the
convenience yield on safe assets. Di Tella, Hébert, Kurlat & Wang (2025, JPE)
suggest the absolute level is large: their equity-based **zero-beta rate** sits
~2–4% above safe bond yields, implying a big absolute convenience/wedge on safe
bonds. **big_cy asks whether a convenience yield that large can jointly explain
several puzzles in international finance** — theoretically, empirically, and
quantitatively.

**The model we want to build.** A Kekre-Lenel (2024, AER)–style two-country
(US / RoW) open-economy model with **bonds in the utility function**:
- US households value safe **dollar** bonds (with a potentially time-varying
  convenience/safety shock).
- Foreign households value safe bonds in **both** dollars **and** their own
  currency.

Crank the *absolute* convenience yield to DiTella magnitudes (2–4%) and ask
whether the model can **microfound home currency bias** and thereby resolve
Kekre-Lenel's central tension: KL need *greater US risk-bearing capacity* to
stop the US from running a giant carry portfolio, but high US risk tolerance is
counterfactual. Large convenience yields may deliver KL's attractive results
(dollar's negative beta; carry-return patterns) **without** implying the US
should run the carry trade.

**Three target facts (international finance):**
1. The US is **not** a carry trader (it is net *short* the carry trade), even
   though carry-return patterns exist.
2. Strong **home currency bias** — the US holds dollars; foreigners hold mostly
   own-currency debt, dollars second.
3. **Equity bias of US external assets** — foreign holdings of US assets mirror
   the US market debt/equity split, while US foreign assets are heavily
   equity-tilted. Reinterpretation: not US risk tolerance, but that there isn't
   enough dollar debt outstanding worldwide to satisfy US diversification
   demand, and foreign-currency debt yields ~2–3% below the true (convenience-
   adjusted) risk-free rate, making it unattractive for the US to hold.

**Anchor papers** (in `…/Dropbox/Angus Lewis/papers/`):
- DiTella, Hébert, Kurlat & Wang (2025, JPE), *The Zero-Beta Interest Rate* —
  the large-absolute-convenience-yield evidence and the "safe-rate puzzle."
- Kekre & Lenel (2024, AER), *The Flight to Safety and International Risk
  Sharing* — the structural template (safety shocks, convenience yield ω,
  nominal rigidity) and the carry/risk-tolerance tension we want to resolve.

**In-house antecedents** (Angus's other projects; read-only references):
- **`rier`** ("U.S. exorbitant privilege: Expected vs. realized", Lewis) —
  the empirical backbone. Granular NFA expected-return analysis establishing all
  three target facts (US net short carry; foreigners ≈ US market portfolio; US
  equity-biased vs all-foreign but debt-biased vs foreign-equity+foreign-dollar-
  debt; CY sensitivity machinery). Latest draft: `draft_2026_04_27.tex` in the
  rier Overleaf folder (note: that folder's `draft_current.tex` is currently a
  stray copy of the lx1 paper — use the dated draft).
- **`lx1`** (Lewis & Xie, cross-currency corporate issuance) — same
  absolute-vs-relative logic on the financing side: the CIP basis *understates*
  the gains from foreign-currency issuance because it misses inframarginal
  diversification savings. Intellectual cousin, not a dependency.

## Standing policies

1. **Local by default; Sherlock allowed for specific approved uses.** Work runs
   **locally by default** — do not submit batch jobs, run `sherlock-agent-com`,
   or take `/oak` / cluster actions **unless the task has been explicitly
   approved for Sherlock**. Approved Sherlock uses to date:
   - **klrep KL-model solves** (approved 2026-06-13): the full-grid solve is
     ~hours/benchmark and ~a day for all 9 calibrations on local CPU but minutes
     on GPU, so the solves run on the cluster as a GPU array-per-specification
     job (`batch_submit.sh` + `sherlock-agent-com`; see `management/proposals/
     04b_assignment_klrep_sherlock.md`).

   Everything else stays local until Angus approves it for the cluster. The
   `batch_submit.sh` / `batch_controller.{sh,do}` files and the `$sherlock==1`
   path branches support these cluster runs. Submission workflow: `cc-merge` →
   `sherlock-agent-com big_cy sync_git` → `sherlock-agent-com big_cy submit`.

2. **Tensor-native quantitative model code.** Any quantitative / structural
   model work (estimation, simulation, solving) must be written
   **tensor-native** — vectorized tensor operations rather than scalar loops,
   device-agnostic where practical — so the code can move to GPU hardware
   (H100s) later with minimal rewriting. Prefer a framework with a GPU backend
   (PyTorch / JAX / NumPy). Reduced-form / descriptive analysis is exempt
   unless it involves heavy numerical model computation.

## Layout

- **Code + management** live in this repo (`~/code/big_cy`).
- **Data and outputs** live in Dropbox: `…/Angus Lewis/big_cy/`
  (`data/{raw,temp,output}`, `figtab`, `logs`). The path globals in
  `project_globals.do` / `project_strings.py` resolve there locally; `data/`,
  `logs/`, and `figtab/` are git-ignored in this repo.
