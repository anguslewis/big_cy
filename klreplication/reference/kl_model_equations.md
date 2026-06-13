# Kekre & Lenel (2024, AER) — "The Flight to Safety and International Risk Sharing"
## Model Equation Reference

**Purpose.** Equation reference for a faithful Python (tensor-native) re-implementation
of the KL (2024) quantitative model. Transcribed from primary sources.

**Sources used (all read-only):**
- **Main text** (model section 2, eqs (1)-(21)): working paper
  `Kekre_Lenel_2024_WP_The_flight_to_safety_and_international_risk_sharing.pdf`, pp. 7-14.
- **Online Appendix A** (full equilibrium, FOCs, and the *re-scaled economy* equations
  (33)-(72)): downloaded from AEA, `https://assets.aeaweb.org/asset-server/files/20818.pdf`
  (55 pp.). **Equations (33)-(72) live in Appendix A.6 "Re-scaled economy", PDF pp. 9-13
  (= appendix-numbered pp. 63-66).** Market clearing (23)-(32) is A.2 (p. 57); equilibrium
  definition is A.3 (p. 58); the un-scaled FOCs are A.5 (pp. 59-62); Foreign optimization
  problem is A.1 (pp. 55-56).
- **Solution algorithm** (defines the scaling and (N1)-(N3), references (33)-(72)):
  `…/big_cy/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/SolutionAlgorithm.pdf`.
- **Parameter name map**: Fortran `mod_param.f90` / `mod_calc.f90` in that package.

> **Important provenance note.** The two PDFs the task pointed at (the WP and the AER
> published file in the papers folder) contain only the *main text* — they do **not**
> contain Appendix A. Equations (33)-(72) are in the **separate online appendix**, which
> repeatedly is the document that `SolutionAlgorithm.pdf` means by "appendix A.6". I
> located and downloaded that online appendix; all (33)-(72) below are transcribed from
> readable rendered PDF pages, not from lossy text extraction. Confidence is high; a few
> spots are flagged `⚠ VERIFY`.

**Notation conventions** (used throughout this file):
- `*` = Foreign country variable (e.g. `c*` is Foreign consumption).
- `H`,`F` subscripts = Home-produced / Foreign-produced good (or Home/Foreign-currency bond).
- `s` subscript on bonds = *safe* dollar bonds; `o` subscript = *other* (non-safe) dollar bonds.
- Tilde `x̃` = re-scaled variable (divided by `z_t`, the global TFP trend — see §6).
  Double-tilde `x̃̃` appears on capital/wage *state* variables that are further deflated by
  the disaster term (see the re-scaling definitions). I write `b̃`, `c̃`, `k̃`, `k̄̃` etc.
- Bar `x̄` = aggregate taken as given by the atomistic household, OR a steady-state level
  (`b̄g`, `ī`); `k̄` = global aggregate capital.
- Hat `x̂` = value/consumption *scaled by agent wealth* (`v̂ ≡ µ̃⁻¹ ṽ`), a numerical device
  from SolutionAlgorithm.pdf (see §11).
- `φ_t` (script phi, written `φ` / `varphi`) = the **disaster** shock term (scales capital
  by `exp(φ_t)`, with `φ_t < 0` in a disaster). Do **not** confuse with `φ` (`phi`) the
  **Taylor-rule inflation coefficient** in (10). The paper uses `φ_t`/`varphi_t` for the
  disaster and a plain `φ` for the Taylor coefficient. In the Fortran the Taylor coefficient
  is `phi_h`/`phi_f`.

---

## 1. Symbol ↔ name (parameter) map

Inferred from the main text and confirmed against `mod_param.f90` (the order is the
`param_input_vec` order in `init_setup`).

| Math symbol | Meaning | Fortran name |
|---|---|---|
| `ζ*` (zeta) | relative size of Foreign (RoW vs US); Home mass = 1, Foreign mass = `ζ*` | `zeta` |
| `β`, `β*` | discount factors (Home/Foreign differ to match NFA level) | `bbeta_vec(1)`,`bbeta_vec(2)` (`bbeta_h`,`bbeta_f`) |
| `γ`, `γ*` | risk aversion (Home/Foreign differ — Home has *more* risk-bearing capacity) | `gma_vec(1)`,`gma_vec(2)` (`gma_h`,`gma_f`) |
| `ψ` | IES (jointly controls intertemporal subst. + consumption-labor complementarity) | `ies_vec(1)`,`ies_vec(2)` (`ies_h`,`ies_f`) |
| `χ^W` | Rotemberg wage-adjustment-cost scale (same in both countries, ≈ 400 → Calvo 5 qtr) | `chi` |
| `σ` | trade elasticity (Armington elasticity across Home/Foreign goods) | `sigma` |
| `ς` (varsigma) | home-bias parameter in the CES consumption aggregator | `varsigma_vec(1)`,`varsigma_vec(2)` |
| `δ` | depreciation rate of capital | `ddelta` |
| `α` | capital share in production (labor share = 1−α) | `aalpha` |
| `χ^x` | investment / capital-producer adjustment-cost curvature | `chiX` |
| `φ_w`, `ε_w` (varepsilon_w) | wage Phillips-curve objects (see note) | `phi_w`, `vareps_w` |
| `φ` (Taylor) | Taylor-rule coefficient on inflation, Home / Foreign | `phi_h`, `phi_f` |
| `φ_yh`,`φ_yf` | Taylor-rule output-gap coefficients (alt. rule / robustness) | `phi_yh`,`phi_yf` |
| `ρ_i` | interest-rate smoothing in Taylor rule | `rho_i` |
| `b̄g` | fiscal rule: safe-dollar govt-debt-to-global-consumption ratio (eq (12)) | `b_lmbd` (shift), `bg_yss` (ss level) |
| `ν̄`, `ν̄*` | disutility-of-labor level (Foreign differs only to normalize ℓ*=1) | (in `chi0_vec`) |
| `ν` | controls Frisch elasticity of labor | (in Φ; see (3)) |
| `ε` (epsilon) | elasticity of substitution across labor varieties | `vareps_w`/related |
| `ε^d` (eps^d) | elasticity of bond demand to non-pecuniary value (in Ω) | (in convenience-yield block) |
| `ω^d`, `ω^d*` | exogenous private demand for safe dollar bonds (the **safety shock** driver) | `omg_*` (`omg_mean`,`omg_std`,`omg_shift`,`omg_rho`) |
| `Δ^ω` | shift so that mean(ω_t)=0 | `omg_shift` |
| `σ^z`, `ρ^p`,`σ^p`, `ρ^ω`,`σ^ω`, `ρ^F`,`σ^F`, `ρ^{pω}` | shock SDs / AR(1) / correlation | `sig_z`, `disast_rho`/`disast_std`, `omg_rho`/`sig_omg`, `rho_zf`/`zf_std`, `corr_omg_dis` |

`n_I = 2` agents (Home, Foreign). `idx_*` and `sidx_*` in the Fortran index the 7 state
variables and the 5 shocks respectively (z, zF, p, ω, disaster).

---

## 2. Households — Home problem (main text §2.1, eqs (1)-(5))

**(1) Epstein-Zin recursion with bonds in utility.** The representative Home household has
recursive preferences over consumption `c_t`, labor `ℓ_t`, and the real value of *safe*
dollar bonds `B_{Ht,s}/P_t`:

```
v_t = ( (1−β) ( c_t · Φ(ℓ_t) · Ω_t( B_{Ht,s}/P_t ) )^{1−1/ψ}
        + β · E_t[ (v_{t+1})^{1−γ} ]^{(1−1/ψ)/(1−γ)} )^{1/(1−1/ψ)}
```

The certainty equivalent is `CE_t = E_t[(v_{t+1})^{1−γ}]^{1/(1−γ)}` (named `ce_t` later).

**(2) CES consumption aggregator** (Home-/Foreign-produced goods, home bias via `ς`):

```
c_t = ( (1/(1+ζ*) + ς)^{1/σ} · c_{Ht}^{(σ−1)/σ}
        + (ζ*/(1+ζ*) − ς)^{1/σ} · c_{Ft}^{(σ−1)/σ} )^{σ/(σ−1)}
```

**(3) Disutility of labor Φ(ℓ)** (Shimer 2010 / Trabandt-Uhlig 2011 form; multiplies `c` in
utility, giving consumption-labor complementarity governed by ψ):

```
Φ(ℓ_t) = ( 1 + (1/ψ − 1) · ν̄ · ℓ_t^{1+1/ν} / (1+1/ν) )^{(1/ψ)/(1−1/ψ)}
```
⚠ VERIFY the outer exponent: rendered as `1/ψ / (1−1/ψ)`. (Foreign uses `ν̄*`.)

**(4) Budget / resource constraint** (nominal; household supplies labor, chooses
consumption + portfolio of: safe dollar bonds `B_{Ht,s}`, other dollar bonds `B_{Ht,o}`,
Foreign-currency bonds `B_{Ft}`, capital `k_t`):

```
P_{Ht} c_{Ht} + E_t^{-1} P*_{Ft} c_{Ft} + B_{Ht,s} + B_{Ht,o} + E_t^{-1} B_{Ft} + Q_t^k k_t
  ≤ (1+i_{t−1}) B_{Ht−1,s} + (1+ι_{t−1}) B_{Ht−1,o} + E_t^{-1}(1+i*_{t−1}) B_{Ft−1}
    + (Π_t + (1−δ)Q_t^k) k_{t−1} exp(φ_t)
    + ∫_0^1 W_t(j) ℓ_t(j) dj − ∫_0^1 AC_t^W(j) dj + T_t
```
Here `E_t` = nominal exchange rate (Foreign units of account per dollar); `i_t` = safe
dollar rate, `ι_t` = other-dollar rate, `i*_t` = Foreign rate. `exp(φ_t)` = disaster scaling
of capital. `T_t` = government transfer.

**(5) Rotemberg wage adjustment cost** (the model uses **Rotemberg**, NOT Calvo; disaster
term `exp(φ_t)` is in the denominator for computational simplicity):

```
AC_t^W(j) = (χ^W / 2) · W_t ℓ_t · ( W_t(j) / (W_{t−1}(j) exp(φ_t)) − 1 )^2
```

**Foreign household** (Appendix A.1, pp. 55-56) is analogous, with its own `Ω*`, risk
aversion `γ*`, discount `β*`, labor disutility `ν̄*`, *same* `ψ, ς, σ, ν, χ^W`. Foreign value
function:
```
v*_t = ( (1−β*)( c*_t Φ*(ℓ*_t) Ω*_t( B*_{Ht,s}/(E_t^{-1} P*_t) ) )^{1−1/ψ}
         + β* E_t[(v*_{t+1})^{1−γ*}]^{(1−1/ψ)/(1−γ*)} )^{1/(1−1/ψ)}
```
Foreign also values **safe dollar bonds** `B*_{Ht,s}` (real value `B*_{Ht,s}/(E_t^{-1}P*_t)`)
— this is the channel by which RoW demands dollars. Foreign budget constraint (A.1):
```
E_t P*_{Ht} c*_{Ht} + P*_{Ft} c*_{Ft} + E_t B*_{Ht,s} + E_t B*_{Ht,o} + B*_{Ft} + E_t Q_t^k k*_t
  ≤ E_t(1+i_{t−1})B*_{Ht−1,s} + E_t(1+ι_{t−1})B*_{Ht−1,o} + (1+i*_{t−1})B*_{Ft−1}
    + E_t(Π_t+(1−δ)Q_t^k)k*_{t−1}exp(φ_t) + ∫W*_t(j*)ℓ*_t(j*)dj* − ∫ AC_t^{W*}(j*)dj* + T*_t
```

---

## 3. THE CENTRAL FEATURE — bonds-in-utility / convenience yield (Ω, ω) — main text §2.4

This is the mechanism big_cy cares about. Document carefully.

**Convenience-yield wedge / investor indifference (14).** Because safe and "other" dollar
bonds both pay in dollars and are risk-free, indifference requires:
```
(1 + i_t) / ( 1 − c_t · Ω'_t(B_{Ht,s}/P_t)/Ω_t(B_{Ht,s}/P_t) )  =  1 + ι_t
```
LHS = effective (convenience-inclusive) return on *safe* dollar bonds; RHS = return on
*other* dollar bonds. The non-pecuniary marginal value, equated across agents, is **`ω_t`**.

**Definition of ω_t (15)** (the marginal convenience yield, equated Home = Foreign):
```
ω_t ≡ c_t · Ω'_t(B_{Ht,s}/P_t) / Ω_t(B_{Ht,s}/P_t)
    = c*_t · Ω*'_t(B*_{Ht,s}/(E_t^{-1}P*_t)) / Ω*_t(B*_{Ht,s}/(E_t^{-1}P*_t))
```

**Functional form for Ω_t** (the "particularly convenient" form; bars = aggregates the
household takes as given):
```
Ω_t( B_{Ht,s}/P_t ) =
  exp(  ω^d_t · (B_{Ht,s}/(P_t c̄_t))  −  (1/2)(1/ε^d) (B_{Ht,s}/(P_t c̄_t))^2
      − [ ω^d_t · (B̄_{Ht,s}/(P_t c̄_t)) − (1/2)(1/ε^d)(B̄_{Ht,s}/(P_t c̄_t))^2 ]  )
```
The first two terms give a time-varying non-pecuniary value diminishing in the position;
the bracketed (aggregate) terms ensure **`Ω_t(B_{Ht,s}/P_t) = 1` in equilibrium**, so a
time-varying convenience yield does **not** mechanically perturb the stochastic discount
factor. (Foreign `Ω*` is analogous with `ω^d*_t`, `ε^d`.)

**Equilibrium ω_t (16)-(17).** The Ω-form + market clearing in safe dollar bonds give the
real safe-bond-to-consumption ratios equalized across countries (16):
```
B_{Ht,s}/(P_t c_t) = B*_{Ht,s}/(E_t^{-1}P*_t c*_t) = (−B^g_{Ht,s}) / (P_t c_t + ζ* E_t^{-1} P*_t c*_t)      (16)
```
and the **equilibrium convenience yield (17):**
```
ω_t = ω^d_t − (1/ε^d) · (−B^g_{Ht,s}) / ( P_t c_t + ζ* E_t^{-1} P*_t c*_t )                                 (17)
```
**Intuition:** ω rises with private demand `ω^d_t` (the safety shock) and falls with public
supply `−B^g_{Ht,s}`; sensitivity governed by `ε^d`. Because govt supply (12) is a fixed
ratio of global consumption, `ω_t` inherits the stochastic properties of `ω^d_t` and is
effectively exogenous in the baseline.

**How ω enters returns / budget (recurring throughout the scaled system below):**
- The *effective gross return* on safe dollar bonds is `(1+i_t)/(1−ω_t)` — see the FOC (38)
  and the term `(1+r_{t+1})/(1−ω_t)` in the wealth transition (51) and budget constraints.
- **Total dollar-bond position** (composition irrelevant except for seigniorage):
  ```
  B_{Ht}  ≡ (1−ω_t)(B_{Ht,s} + B^g_{Ht,s}) + B_{Ht,o}
  B*_{Ht} ≡ (1−ω_t) B*_{Ht,s} + B*_{Ht,o}
  ```
  each earning `1 + ι_t = (1+i_t)/(1−ω_t)`.
- **Seigniorage / transfer.** Home earns seigniorage on safe dollar debt held by Foreign.
  As derived in the appendix, Foreign transfers to Home each period
  ```
  ω_t · [ (ζ* q_t^{-1} c*_t) / (c_t + ζ* q_t^{-1} c*_t) ] · (−B^g_{Ht,s})    (units of Home bundle)
  ```
  (this is the `seigniorage_nxt` block in `mod_calc.f90 :: calc_unexpected_transition`,
  computed as `cf_spending/(ch_spending+cf_spending) * (bg + shock)*(...) * omg`).

---

## 4. Government / fiscal block (main text §2.3, eqs (12),(13))

**(12) Safe-dollar government debt rule** (constant ratio of safe dollar debt `−B^g_{Ht,s}`
to global consumption; empirically relevant case `b̄g > 0` = govt *borrows* in T-bills):
```
−B^g_{Ht,s} = b̄g · ( P_t c_t + ζ* E_t^{-1} P*_t c*_t )
```

**(13) Home government transfer** (rebates wage-adjustment costs + rolls the debt):
```
T_t = ∫_0^1 AC_t^W(j) dj + (1+i_{t−1}) B^g_{Ht−1,s} − B^g_{Ht,s}
```
Foreign government provides wage subsidies + lump-sum transfers but does **not** participate
in asset markets (assumed unable to create safe-dollar liquidity), so Ricardian equivalence
holds there. Foreign fiscal: `T*_t = ∫_0^1 AC_t^{W*}(j*) dj*` (A.1).

---

## 5. Supply side (main text §2.2, eqs (6)-(9)) and Taylor/price (10)-(11)

**(6) Labor packer** (CES across labor varieties, elasticity `ε`):
```
W_t [ ∫_0^1 ℓ_t(j)^{(ε−1)/ε} ]^{ε/(ε−1)} − ∫_0^1 W_t(j) ℓ_t(j) dj
```
**(7) Producer** (Home; productivity `z_t`, labor share `1−α`):
```
P_{Ht} (z_t ℓ_t)^{1−α} (κ_t)^α − W_t ℓ_t − Π_t κ_t
```
(`κ_t` = capital rented; `Π_t` = dollar return per unit capital, equalized across countries
since capital is freely deployable.)
**(8) Global capital producer aggregator** (no home bias in investment):
```
x_t = ( (1/(1+ζ*))^{1/σ}(x_{Ht})^{(σ−1)/σ} + (ζ*/(1+ζ*))^{1/σ}(x_{Ft})^{(σ−1)/σ} )^{σ/(σ−1)}
```
**(9) Capital producer profit** (zero in equilibrium; `χ^x` adjustment curvature):
```
Q_t^k x_t − ( k̄_t / (k̄_{t−1} exp(φ_t)) )^{χ^x} ( P_{Ht} x_{Ht} + E_t^{-1} P*_{Ft} x_{Ft} )
```

**(10) Taylor rule (Home)** — CPI-targeting; Foreign analogous with same `φ`:
```
1 + i_t = (1 + ī) · ( P_t / P_{t−1} )^φ
```
(Fortran adds smoothing `ρ_i` and optional output term `φ_yh`; baseline is this form.)

**(11) Ideal price index** (Home):
```
P_t = [ (1/(1+ζ*) + ς) P_{Ht}^{1−σ} + (ζ*/(1+ζ*) − ς)(E_t^{-1} P*_{Ft})^{1−σ} ]^{1/(1−σ)}
```
Foreign ideal price index (A.1):
```
P*_t = [ (1/(1+ζ*) − ς/ζ*)(E_t P_{Ht})^{1−σ} + (ζ*/ζ* + ς/ζ*) P*_{Ft}^{1−σ} ]^{1/(1−σ)}
```
⚠ VERIFY the Foreign weights: rendered `(1/(1+ζ*) − ς/ζ*)` and `(ζ*/(1+ζ*) + ς/ζ*)`.

---

## 6. Driving forces (main text §2.5, eqs (18)-(21))

Four shocks; innovations `{ε^z, ε^p, ε^ω, ε^F}` are iid N(0,1); disaster/ω may be
correlated `ρ^{pω}`, all others zero.

**(18) Global TFP** — unit root + rare disaster:
```
log(z_t) = log(z_{t−1}) + σ^z ε^z_t + φ_t
```
where `φ_t = 0` w.p. `1−p_t`, and `φ_t = φ̲ < 0` (a fixed negative draw) w.p. `p_t`.

**(19) Disaster probability** — AR(1) in logs (captures skewness):
```
log p_t − log p̄ = ρ^p ( log p_{t−1} − log p̄ ) + σ^p ε^p_t
```

**(20) Safety / convenience demand** — `ω^d_t = Δ^ω + ω̃^d_t`, with:
```
log ω̃^d_t − log ω̄^d = ρ^ω ( log ω̃^d_{t−1} − log ω̄^d ) + σ^ω ε^ω_t
```
(`Δ^ω` shifts so mean(ω_t)=0; insights rely only on *time-variation* in ω.)

**(21) Relative Foreign productivity** — AR(1):
```
log z_{Ft} = ρ^F log z_{Ft−1} + σ^F ε^F_t
```

---

## 7. Market clearing (Appendix A.2, eqs (23)-(32)) — un-scaled

```
(23) Home good:    c_{Ht} + ζ* c*_{Ht} + (k̄_t/(k̄_{t−1}exp(φ_t)))^{χ^x} x_{Ht} = (z_t ℓ_t)^{1−α}(κ_t)^α
(24) Foreign good: c_{Ft} + ζ* c*_{Ft} + (k̄_t/(k̄_{t−1}exp(φ_t)))^{χ^x} x_{Ft} = (z_t z_{Ft} ζ* ℓ*_t)^{1−α}(κ*_t)^α
(25) Home labor:   [ ∫_0^1 ℓ_t(j)^{(ε−1)/ε} dj ]^{ε/(ε−1)} = ℓ_t
(26) Foreign labor:[ ∫_0^1 ℓ*_t(j*)^{(ε−1)/ε} dj* ]^{ε/(ε−1)} = ℓ*_t
(27) Capital rental: κ_t + κ*_t = k̄_{t−1} exp(φ_t)
(28) Capital:        k_{t−1} + ζ* k*_{t−1} = k̄_{t−1}
(29) Capital LoM:    (1−δ) k̄_{t−1} exp(φ_t) + x_t = k̄_t
(30) Safe dollar bonds:  B_{Ht,s} + ζ* B*_{Ht,s} + B^g_{Ht,s} = 0
(31) Other dollar bonds: B_{Ht,o} + ζ* B*_{Ht,o} = 0
(32) Foreign bonds:      B_{Ft} + ζ* B*_{Ft} = 0
```

---

## 8. Definition of equilibrium (Appendix A.3, Definition 1)

> An equilibrium is a sequence of prices and policies such that:
> - each Home representative household chooses `{c_{Ht}, c_{Ft}, B_{Ht,s}, B_{Ht,o}, B_{Ft}, k_t}`
>   to maximize (1) subject to (2)-(5), and analogously in Foreign;
> - each Home union `j` chooses `{W_t(j), ℓ_t(j)}` to maximize utilitarian social welfare of
>   its members subject to (5), and analogously in Foreign;
> - the representative Home labor packer chooses `{ℓ_t(j)}` to maximize profits (6), analogously
>   in Foreign;
> - the representative Home producer chooses `{ℓ_t, κ_t}` to maximize profits (7), analogously
>   in Foreign;
> - the representative global capital producer chooses `{x_{Ht}, x_{Ft}, x_t}` to maximize
>   profits (9) subject to (8);
> - the Home government sets `B^g_{Ht,s}` per (12) and `{i_t, {T_t}}` per (10) and (13); the
>   Foreign government analogously sets the latter;
> - goods, factor, and asset markets clear per (23)-(32).

**Recursive representation (main text §2.6).** All real variables (except labor) divided by
`z_t`, nominal by `P_t z_t`. The aggregate **state in period t (7 variables)**:
```
{ p_t, ω_t, z_{Ft}, θ_t, k̄̃_{t−1}, w̃̃_{t−1}, w̃̃*_{t−1} }
```
= disaster probability, convenience yield, relative Foreign productivity, **Home financial
wealth share θ_t**, scaled lagged aggregate capital, and scaled lagged Home & Foreign real
wages. (Confirmed verbatim at end of A.6, p. 66.) In the Fortran the 9 state slots add the
two lagged nominal rates `i_{t−1}`, `i*_{t−1}` for the smoothed Taylor rule.

---

## 9. Additional variables (Appendix A.4) — used by (33)-(72)

```
Real exchange rate:        q_t ≡ E_t P_t / P*_t          (increase = Home appreciation)
Home real rate:            1 + r_{t+1} ≡ (1+i_t) P_t/P_{t+1}
Foreign real rate:         1 + r*_{t+1} ≡ (1+i*_t) P*_t/P*_{t+1}
Real return to capital:    1 + r^k_{t+1} ≡ ((Π_{t+1}+(1−δ)Q^k_{t+1})/Q^k_t)(P_t/P_{t+1}) exp(φ_{t+1})
Terms of trade:            s_t ≡ E_t P_{Ht} / P*_{Ft}
Output (Home):             y_t ≡ (z_t ℓ_t)^{1−α}(κ_t)^α
Real aggregate saving:     a_t ≡ (1/P_t)( B_{Ht} + E_t^{-1} B_{Ft} + Q^k_t k_t )
Real net foreign assets:   nfa_t ≡ a_t − (Q^k_t/P_t) κ_{t+1} exp(−φ_{t+1})
```

---

## 10. ★ EQUATIONS (33)-(72) — Appendix A.6 "Re-scaled economy" (PDF pp. 9-13) ★

**Re-scaling definitions** (so the system is stationary). For Home (Foreign analogous with `*`):
```
c̃_t ≡ c_t/z_t,  c̃_{Ht} ≡ c_{Ht}/z_t,  c̃_{Ft} ≡ c_{Ft}/z_t,  CE_t̃ ≡ CE_t/z_t,
m̃_{t,t+1} ≡ m_{t,t+1}(z_{t+1}/z_t)^γ,
b̃_{Ht} ≡ B_{Ht,s}/(P_t z_t),  b̃^g_{Ht,s} ≡ B^g_{Ht,s}/(P_t z_t),  b̃_{Ft} ≡ B_{Ft}/(P_t z_t),
b̃̃_{Ht−1} ≡ b̃_{Ht−1}/exp(σ^z ε^z_t + φ_t),   b̃̃_{Ft−1} ≡ b̃_{Ft−1}/exp(σ^z ε^z_t + φ_t),
k̃_t ≡ k_t/z_t,  κ̃_t ≡ κ_t/z_t,  k̃̃_{t−1} ≡ k̃_{t−1}/exp(σ^z ε^z_t),
w̃_t ≡ w_t/z_t,  w̃̃_{t−1} ≡ w̃_{t−1}/exp(σ^z ε^z_t),
x̃_{Ht} ≡ x_{Ht}/z_t,  x̃_{Ft} ≡ x_{Ft}/z_t,  x̃_t ≡ x_t/z_t,  k̄̃_t ≡ k̄_t/z_t,  k̄̃̃_{t−1} ≡ k̄̃_{t−1}/exp(σ^z ε^z_t).
```
(`b̃_{Ht}` here is the household's *total scaled dollar bond*; the s/o split only matters for
seigniorage.) Re-scaled real rates: `r`, `r*`, `r^k` defined by the §9 formulas with scaled
prices; `φ_{t+1}` written as `φ_{t+1}` below.

### Home FOCs and constraint

**(33) Re-scaled value function** (Epstein-Zin recursion, scaled):
```
ṽ_t = ( (1−β)(c̃_t Φ(ℓ_t))^{1−1/ψ} + β (CE_t̃)^{1−1/ψ} )^{1/(1−1/ψ)}
```

**(34) Re-scaled certainty equivalent:**
```
CE_t̃ = E_t[ exp((1−γ)(σ^z ε^z_{t+1} + φ_{t+1})) (ṽ_{t+1})^{1−γ} ]^{1/(1−γ)}
```

**(35) Consumption aggregator (scaled):**
```
c̃_t = ( (1/(1+ζ*) + ς)^{1/σ}(c̃_{Ht})^{(σ−1)/σ} + (ζ*/(1+ζ*) − ς)^{1/σ}(c̃_{Ft})^{(σ−1)/σ} )^{σ/(σ−1)}
```

**(36) Home/Foreign goods demand (intratemporal, Home household):**
```
c̃_{Ht}/c̃_{Ft} = [ (1/(1+ζ*) + ς) / (ζ*/(1+ζ*) − ς) ] · s_t^{−σ}
```

**(37) Pricing kernel / Euler (consumption-savings):**
```
m̃_{t,t+1} = β (c̃_t/c̃_{t+1}) ( c̃_{t+1} Φ(ℓ_{t+1}) / (c̃_t Φ(ℓ_t)) )^{1−1/ψ} ( ṽ_{t+1}/CE_t̃ )^{1/ψ−γ}
```
⚠ The leading factor renders as `β (c̃_t/c̃_{t+1})·(…)` — i.e. β times the consumption-growth
and complementarity ratios. (In the un-scaled A.5 form it is `β (c_t/c_{t+1})(c_{t+1}Φ_{t+1}/(c_tΦ_t))^{1−1/ψ}(v_{t+1}/ce_t)^{1/ψ−γ}`. The SolutionAlgorithm (N1) writes this as
`m̃ = β µ̃_{t+1}^{−γ} (v̂_{t+1}/CE)^{1/ψ−γ} (ĉ_{t+1})^{−1/ψ}Φ(ℓ_{t+1})^{1−1/ψ} / [(c̃_t)^{−1/ψ}Φ(ℓ_t)^{1−1/ψ}]`.)

**(38) Portfolio FOC — safe dollar bond** (note the `(1+r)/(1−ω)` convenience-inclusive return):
```
1 = E_t m̃_{t,t+1} exp(−γ(σ^z ε^z_{t+1} + φ_{t+1})) · (1 + r_{t+1})/(1 − ω_t)
```

**(39) Portfolio FOC — Foreign bond:**
```
1 = E_t m̃_{t,t+1} exp(−γ(σ^z ε^z_{t+1} + φ_{t+1})) · (q_t/q_{t+1}) (1 + r*_{t+1})
```

**(40) Portfolio FOC — capital:**
```
1 = E_t m̃_{t,t+1} exp(−γ(σ^z ε^z_{t+1} + φ_{t+1})) · (1 + r^k_{t+1})
```

**(41) Home budget constraint (scaled):**
```
c̃_t + b̃_{Ht} + q_t^{-1} b̃_{Ft} + q_t^k k̃_t = w̃_t ℓ_t + θ_t (π_t + (1−δ)q_t^k) k̄̃̃_{t−1}
```
⚠ VERIFY: the dividend/wealth term uses the double-tilde lagged aggregate capital `k̄̃̃_{t−1}`
and `π_t` = scaled `Π_t/P_t`. (The full A.5 un-scaled version carries the `(1+r)/(1−ω)` on
lagged bonds and a separate `b^g` term; (41) is the wealth-share-collapsed form using θ_t.)

### Foreign FOCs and constraint

**(42) Foreign value function:**
```
ṽ*_t = ( (1−β*)(c̃*_t Φ*(ℓ*_t))^{1−1/ψ} + β* (CE*_t̃)^{1−1/ψ} )^{1/(1−1/ψ)}
```

**(43) Foreign certainty equivalent:**
```
CE*_t̃ = E_t[ exp((1−γ*)(σ^z ε^z_{t+1} + φ_{t+1})) (ṽ*_{t+1})^{1−γ*} ]^{1/(1−γ*)}
```

**(44) Foreign consumption aggregator (scaled):**
```
c̃*_t = ( (1/(1+ζ*) − ς/ζ*)^{1/σ}(c̃*_{Ht})^{(σ−1)/σ} + (ζ*/(1+ζ*) + ς/ζ*)^{1/σ}(c̃*_{Ft})^{(σ−1)/σ} )^{σ/(σ−1)}
```

**(45) Foreign goods demand:**
```
c̃*_{Ht}/c̃*_{Ft} = [ (1/(1+ζ*) − ς/ζ*) / (ζ*/(1+ζ*) + ς/ζ*) ] · s_t^{−σ}
```

**(46) Foreign pricing kernel:**
```
m̃*_{t,t+1} = β* (c̃*_t/c̃*_{t+1}) ( c̃*_{t+1}Φ*(ℓ*_{t+1}) / (c̃*_t Φ*(ℓ*_t)) )^{1−1/ψ} ( ṽ*_{t+1}/CE*_t̃ )^{1/ψ−γ*}
```

**(47) Foreign FOC — safe dollar bond:**
```
1 = E_t m̃*_{t,t+1} exp(−γ*(σ^z ε^z_{t+1} + φ_{t+1})) · (q_{t+1}/q_t)(1 + r_{t+1})/(1 − ω_t)
```

**(48) Foreign FOC — Foreign bond:**
```
1 = E_t m̃*_{t,t+1} exp(−γ*(σ^z ε^z_{t+1} + φ_{t+1})) · (1 + r*_{t+1})
```

**(49) Foreign FOC — capital:**
```
1 = E_t m̃*_{t,t+1} exp(−γ*(σ^z ε^z_{t+1} + φ_{t+1})) · (q_{t+1}/q_t)(1 + r^k_{t+1})
```

**(50) Foreign budget constraint (scaled):**
```
q_t^{-1} c̃*_t + b̃*_{Ht} + q_t^{-1} b̃*_{Ft} + q_t^k k̃*_t = w̃*_t ℓ*_t + (1/ζ*)(1−θ_t)(π_t + (1−δ)q_t^k) k̄̃̃_{t−1}
```

### Wealth share, unions, producers

**(51) Global wealth share of Home (inclusive of seigniorage) — θ transition:**
```
θ_{t+1} = ( 1 / ((π_{t+1}+(1−δ)q^k_{t+1}) k̄̃̃_t) ) ·
          [ ((1+r_{t+1})/(1−ω_t)) b̃̃_{Ht}
            + (1/q_{t+1})(1+r*_{t+1}) b̃̃_{Ft}
            + (π_{t+1}+(1−δ)q^k_{t+1}) k̃̃_t
            − ω_{t+1} · ( ζ* q^{-1}_{t+1} c̃*_{t+1} / (c̃_{t+1} + ζ* q^{-1}_{t+1} c̃*_{t+1}) ) · b̃^g_{Ht+1,s} ]
```
The last term is the **seigniorage** Home owes/earns on the safe dollar debt held by Foreign
(matches the `mod_calc.f90` seigniorage block). ⚠ VERIFY signs/`k̄̃̃_t` vs `k̃̃_t` in the
denominator vs numerator (denominator uses aggregate `k̄̃̃_t`, third numerator term uses
household `k̃̃_t`).

**(52) Home union optimal labor supply** (Rotemberg wage Phillips curve; supply-side optimality):
```
w̃_t − c̃_t Φ'(ℓ_t)/Φ(ℓ_t)
  + w̃_t (χ^W/ε) [ (w̃_t/w̃̃_{t−1})(P_t/P_{t−1}) ( (w̃_t/w̃̃_{t−1})(P_t/P_{t−1}) − 1 )
      − E_t m_{t,t+1} exp(−γ(σ^z ε^z_{t+1}+φ_{t+1})) ×
        (w̃_{t+1}/w̃_t)^2 (P_{t+1}/P_t)(ℓ_{t+1}/ℓ_t)( (w̃_{t+1}/w̃_t)(P_{t+1}/P_t) − 1 ) ] = 0
```

**(53) Foreign union optimal labor supply** (analogous; note `E_t P_t` and `q`/`m` Foreign):
```
w̃*_t − q_t^{-1} c̃*_t Φ*'(ℓ*_t)/Φ*(ℓ*_t)
  + w̃*_t (χ^W/ε)[ (w̃*_t/w̃̃*_{t−1})(E_t P_t/(E_{t−1}P_{t−1}))( (w̃*_t/w̃̃*_{t−1})(E_t P_t/(E_{t−1}P_{t−1})) − 1 )
      − E_t m_{t,t+1} exp(−γ*(σ^z ε^z_{t+1}+φ_{t+1})) ×
        (w̃*_{t+1}/w̃*_t)^2 (q_{t+1}/q_t)(E_{t+1}P_{t+1}/(E_tP_t))(ℓ*_{t+1}/ℓ*_t)( (w̃*_{t+1}/w̃*_t)(E_{t+1}P_{t+1}/(E_tP_t)) − 1 ) ] = 0
```
⚠ VERIFY the Foreign-union deflators: the wage inflation is in `E_t P_t` (the Foreign nominal
unit price level `E_t P_t`); rendered with `q_{t+1}/q_t` multiplying the continuation term.

**(54)-(57) Wages & profits / capital allocation across countries:**
```
(54) Home real wage:    w̃_t = (P_{Ht}/P_t)(1−α) ℓ_t^{−α} κ̃_t^{α}
(55) Foreign real wage:  w̃*_t = q_t^{-1}(P*_{Ft}/P*_t)(1−α) z_{Ft}^{1−α}(ζ* ℓ*_t)^{−α} κ̃*_t^{α}
(56) Home profit:        π_t = (P_{Ht}/P_t) α ℓ_t^{1−α} κ̃_t^{α−1}
(57) Foreign profit:     π_t = q_t^{-1}(P*_{Ft}/P*_t) α (z_{Ft} ζ* ℓ*_t)^{1−α} κ̃*_t^{α−1}
```
(56)=(57) jointly pin down the capital allocation `κ̃_t`, `κ̃*_t` (returns equalized across
countries). ⚠ VERIFY exponent placement in (55),(57): `z_{Ft}^{1−α}` and the `(ζ*ℓ*_t)`
grouping render as above.

### Goods-market / terms-of-trade / investment block

**(58) Investment aggregator (scaled):**
```
x̃_t = ( (1/(1+ζ*))^{1/σ}(x̃_{Ht})^{(σ−1)/σ} + (ζ*/(1+ζ*))^{1/σ}(x̃_{Ft})^{(σ−1)/σ} )^{σ/(σ−1)}
```

**(59) Investment good demand:**
```
x̃_{Ht}/x̃_{Ft} = (1/ζ*) s_t^{−σ}
```

**(60) Price of capital q_t^k:**
```
q_t^k = ( k̄̃̃_t / k̄̃̃_{t−1} )^{χ^x} ( (1/(1+ζ*))(P_{Ht}/P_t)^{1−σ} + (ζ*/(1+ζ*))(P_{Ft}/P_t)^{1−σ} )^{1/(1−σ)}
```
⚠ VERIFY: second `P_{Ft}` should be the dollar price of the Foreign good `E_t^{-1}P*_{Ft}/P_t`;
rendered as `P_{Ft}/P_t`. The ratio `k̄̃̃_t/k̄̃̃_{t−1}` is the scaled investment-adjustment term.

**(61) Home-good market clearing (scaled):**
```
c̃_{Ht} + ζ* c̃*_{Ht} + ( k̄̃̃_t / k̄̃̃_{t−1} )^{χ^x} x̃_{Ht} = (ℓ_t)^{1−α}(κ̃_t)^α
```

**(62) Foreign-good market clearing (scaled):**
```
c̃_{Ft} + ζ* c̃*_{Ft} + ( k̄̃̃_t / k̄̃̃_{t−1} )^{χ^x} x̃_{Ft} = (z_{Ft} ζ* ℓ*_t)^{1−α}(κ̃*_t)^α
```

**(63) Capital rental clearing (scaled):**
```
κ̃_t + κ̃*_t = k̄̃̃_{t−1}
```

**(64) Aggregate capital savings clearing (scaled):**
```
k̃_t + ζ* k̃*_t = k̄̃_t
```

**(65) Capital law of motion (scaled):**
```
(1−δ) k̄̃̃_{t−1} + x̃_t = k̄̃_t
```

**(66) Home bond market clearing (scaled):**
```
b̃_{Ht} + ζ* b̃*_{Ht} = 0
```
(By Walras' Law the Foreign-bond market then clears.)

### Return definitions and price identities

**(67) Home real return:**
```
1 + r_{t+1} = (1 + i_t) · P_t/P_{t+1}
```
**(68) Foreign real return:**
```
1 + r*_{t+1} = (1 + i*_t) · P*_t/P*_{t+1}
```
**(69) Capital real return:**
```
1 + r^k_{t+1} = ( (π_{t+1} + (1−δ) q^k_{t+1}) / q^k_t ) · exp(φ_{t+1})
```
⚠ Note: in (69) the `P_t/P_{t+1}` of the §9 definition is absorbed because π and q^k are
already real (scaled). The disaster term `exp(φ_{t+1})` is explicit.

**(70) Price identity — P_t/P_{Ht}:**
```
P_t/P_{Ht} = [ (1/(1+ζ*) + ς) + (ζ*/(1+ζ*) − ς) s_t^{σ−1} ]^{1/(1−σ)}
```
**(71) Price identity — P*_t/P*_{Ft}:**
```
P*_t/P*_{Ft} = [ (1/(1+ζ*) − ς/ζ*) s_t^{1−σ} + (ζ*/(1+ζ*) + ς/ζ*) ]^{1/(1−σ)}
```
**(72) Real exchange rate ↔ terms of trade:**
```
q_t = (E_t P_{Ht}/P*_{Ft}) · (P_t/P_{Ht}) / (P*_t/P*_{Ft})
    = s_t · [ ( (1/(1+ζ*)+ς) + (ζ*/(1+ζ*)−ς) s_t^{σ−1} ) / ( (1/(1+ζ*)−ς/ζ*) s_t^{1−σ} + (ζ*/(1+ζ*)+ς/ζ*) ) ]^{1/(1−σ)}
```

> "Together with the Taylor and fiscal rules and specification of driving forces, (33)-(72)
> define the equilibrium. By Walras' Law, the Foreign bond market clears as well." (A.6, p.66)

---

## 11. Numerical scaling devices (SolutionAlgorithm.pdf) — needed to implement (33)-(72)

The solver stores value/consumption **scaled by agent wealth**:
```
v̂_t ≡ µ̃_t^{-1} ṽ_t,   ĉ_t ≡ µ̃_t^{-1} c̃_t
```
where `µ̃_t` is an affine transform of the representative Home agent's current-state wealth:
```
µ̃_t = ( θ_t (π_t + (1−δ) q^k_t) k̄̃̃_{t−1} + ā ) / b̄
```
(`ā`, `b̄` chosen so `µ̃_t` stays positive and ≈ 1). Analogously `µ̃*_t, v̂*_t, ĉ*_t`.

**(N1) Pricing kernel in scaled form:**
```
m̃_{t,t+1} = β (µ̃_{t+1})^{−γ} (v̂_{t+1}/CE_t̃)^{1/ψ−γ} · (ĉ_{t+1})^{−1/ψ}Φ(ℓ_{t+1})^{1−1/ψ} / [ (c̃_t)^{−1/ψ}Φ(ℓ_t)^{1−1/ψ} ]
```
**(N2) Certainty equivalent in scaled form:**
```
CE_t̃ = E_t[ exp((1−γ)(σ^z ε^z_{t+1} + φ_{t+1})) (µ̃_{t+1} v̂_{t+1})^{1−γ} ]^{1/(1−γ)}
```
**(N3) Next-period µ̃ as a function of t savings & portfolio choice** (this is where ω and
the `(1+r)/(1−ω)`, `1/q` and seigniorage terms enter the transition):
```
µ̃_{t+1} = (1/b̄) [ (1/(1−ω_t))(1+r_{t+1}) b̃̃_{Ht}
                  + (1/q_{t+1})(1+r*_{t+1}) b̃̃_{Ft}
                  + (π_{t+1}+(1−δ)q^k_{t+1}) k̃̃_t
                  − ω_{t+1} ( ζ* q^{-1}_{t+1} c̃*_{t+1} / (c̃_{t+1} + ζ* q^{-1}_{t+1} c̃*_{t+1}) ) b̃^g_{Ht+1,s}
                  + ā ]
```
**Implication for the implementation:** at each iteration savings/portfolio choices at `t`
affect the pricing kernel and CE between `t` and `t+1` *only* through `µ̃_{t+1}` and `c̃_t`
(and analogously for Foreign). The solver (subroutine `calc_equilibrium_and_update` in
`mod_calc.f90`) loops the 10 steps in §A.7 / SolutionAlgorithm.pdf, mapping each step to the
equation groups above:
1. current production → (54)-(57),(63),(70)-(72); RHS of (41),(50)
2. next-period production → (69)
3. price of capital → (60)
4. Home bond market clearing → (66) with (38),(40),(47),(49); returns (67),(68)
5. Foreign bond market clearing → (38),(39),(47),(48)
6. value functions + θ transition + aggregate savings → (33),(34),(42),(43),(51),(64)
7. terms of trade → (35),(36),(44),(45),(58),(59),(61),(62),(65)
8. consumption-savings → (37),(40),(41),(46),(49),(50) [(37) rewritten as (N1)]
9. labor supply → (52),(53)
10. inflation rates → Taylor rules (10)

---

## 12. Equation-by-equation coverage checklist

All of (1)-(21), (23)-(72) located and transcribed. (Eq (22) is not used by the referenced
set; equations between (21) and (23) are absorbed in the main-text/appendix transition — no
"(22)" appears in the referenced algorithm.)

| Eq | Status | Eq | Status | Eq | Status |
|----|--------|----|--------|----|--------|
| (1)-(11) main | ✅ | (33)-(41) Home scaled | ✅ | (58)-(62) goods/inv | ✅ |
| (12),(13) fiscal | ✅ | (42)-(50) Foreign scaled | ✅ | (63)-(66) clearing | ✅ |
| (14)-(17) conv. yield | ✅ | (51) θ transition | ✅ | (67)-(69) returns | ✅ |
| (18)-(21) shocks | ✅ | (52),(53) unions | ✅ | (70)-(72) prices/RER | ✅ |
| (23)-(32) clearing | ✅ | (54)-(57) wage/profit/κ | ✅ | (N1)-(N3) scaling | ✅ |

Outstanding `⚠ VERIFY` items (low risk; re-check against rendered PDF pp. 9-13 / 59-62 if a
discrepancy surfaces during implementation): Φ outer exponent (3); Foreign price-index
weights (11); (41) lagged-wealth term form; (51) `k̄̃̃` vs `k̃̃` placement and seigniorage sign;
(52),(53) Foreign-union deflators; (55),(57) exponent grouping; (60) the second good's dollar
price; (37)/(46) leading factor.
