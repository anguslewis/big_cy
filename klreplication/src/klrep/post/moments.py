"""KL Table-2 moments — port of the Table-2 block of `calc_moments.m` (lines
33-47) plus the cross-sim average of `collect_moments.m`.

Delivers the 12 moments that need only the deterministic Table-2 series:
{1,2,3,4,5,6,7,11,12,13,14,15}. Moments 8/9/10 require the levered-equity return
`rA`, which blends the 20-period bond realized returns — i.e. the `calc_bond_prices`
recursive bond ladder — and are deferred (PLAN sect 13.6).

Each moment is computed per simulation path (rows of the (n_sims, T) series) and
then averaged across paths, matching `collect_moments.m`. MATLAB `std`/`corr` use
the sample (N-1) covariance — `torch.std(unbiased=True)` and the Pearson formula
below match that.
"""
import torch

# KL Table-2 "Model" column targets, by moment number (from table_2.tex).
KL_TARGETS = {
    1: 1.6, 2: 0.5, 3: 0.8, 4: 0.5, 5: 1.6, 6: 2.0, 7: -23.0,
    8: 0.5, 9: 5.2, 10: 0.6, 11: 3.8, 12: 1.0, 13: 1.0, 14: 0.0, 15: 0.0,
}

# Short labels (the LaTeX moment expressions, plain-text).
LABELS = {
    1: "y*/(s y)",
    2: "100 sd(dlog c)",
    3: "100 sd(dlog y*)",
    4: "rho(y*/y, lag)",
    5: "100 sd(dlog x)",
    6: "4 E r",
    7: "nfa/(4y) %",
    8: "rho(r^e, carry)",
    9: "4 E[r^e - r] %",
    10: "beta(dnfa/y, r^e-r)",
    11: "b*_Hs/(4y) %",
    12: "l",
    13: "l*",
    14: "100 E log P/P_-1",
    15: "100 E log P*/P*_-1",
}

# The 12 moments deliverable without the bond ladder.
DELIVERABLE = [1, 2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15]


def _std(x):
    """Sample (N-1) std along time (dim 1)."""
    return x.std(dim=1, unbiased=True)


def _dlog(x):
    """diff(log(x)) along time."""
    return torch.log(x[:, 1:]) - torch.log(x[:, :-1])


def _corr(a, b):
    """Per-row Pearson correlation of (n_sims, L) tensors a, b."""
    ac = a - a.mean(dim=1, keepdim=True)
    bc = b - b.mean(dim=1, keepdim=True)
    num = (ac * bc).sum(dim=1)
    den = torch.sqrt((ac ** 2).sum(dim=1) * (bc ** 2).sum(dim=1))
    return num / den


def _ols(cols_list, y):
    """Batched per-sim OLS. `cols_list` = regressors (each (n_sims,L)); an intercept
    is prepended. Returns (beta (n_sims, 1+k), resid (n_sims, L))."""
    n_sims, L = y.shape
    X = torch.stack([torch.ones_like(y)] + list(cols_list), dim=2)   # (n_sims, L, 1+k)
    sol = torch.linalg.lstsq(X, y.unsqueeze(-1)).solution            # (n_sims, 1+k, 1)
    resid = (y.unsqueeze(-1) - X @ sol).squeeze(-1)
    return sol.squeeze(-1), resid


def _bond_moments_per_sim(s):
    """Moments 8/9/10 + their regressions (calc_moments.m:15-31,40-42). Requires the
    bond-dependent series (exc_retA, exc_rf, rA, uip_pvt, ...). Returns dict 8/9/10."""
    # R1: excess equity returns on lagged smoothed dividend-price -> resid_dp.
    y1 = 4.0 * 100.0 * s["exc_retA"][:, 4:]
    x1 = s["div_price_smoothed_1"][:, 3:-1]
    _, resid_dp = _ols([x1], y1)
    # R2: UIP residual on lagged yield spread + lagged output growth -> resid_carry.
    xdep = s["uip_pvt"][:, 4:]
    _, resid_carry = _ols([s["fama_yield_pvt"][:, 3:-1], s["y_growth"][:, 3:-1]], xdep)
    # R3: nfa growth on excess equity + excess foreign-bond returns -> slope on equity.
    y3 = s["nfa_rel_growth"][:, 2:-1]
    beta3, _ = _ols([s["exc_retA"][:, 2:-1], s["exc_rf"][:, 2:-1]], y3)
    return {
        8: _corr(resid_dp, resid_carry),
        9: 4.0 * 100.0 * (s["rA"][:, 1:-1].mean(dim=1) - s["rfh"][:, 1:-1].mean(dim=1)),
        10: beta3[:, 1],     # loading on excess equity return
    }


def compute_table2_moments_per_sim(s, *, bg_yss):
    """Per-simulation Table-2 moments (each value a (n_sims,) tensor) — the moment
    computed on each path before the cross-sim average of `collect_moments.m`.
    Moments 8/9/10 are included only if the bond-dependent series are present."""
    ratio = s["yf"] / s["yh"]                                   # y*/y
    out = {
        1: (s["yf"] / (s["s"] * s["yh"])).mean(dim=1),
        2: 100.0 * _std(_dlog(s["ch"])),
        3: 100.0 * _std(_dlog(s["yf"])),
        4: _corr(ratio[:, :-1], ratio[:, 1:]),
        5: 100.0 * _std(_dlog(s["inv"])),
        6: 4.0 * 100.0 * s["rfh"].mean(dim=1),
        7: 100.0 * (s["nfa_rel"] / 4.0).mean(dim=1),
        11: 100.0 * bg_yss * ((s["cf"] / s["qx"]) / (4.0 * s["yh"])).mean(dim=1),
        12: s["lh"].mean(dim=1),
        13: s["lf"].mean(dim=1),
        14: 100.0 * torch.log(s["infl_h"]).mean(dim=1),
        15: 100.0 * torch.log(s["infl_f"]).mean(dim=1),
    }
    if "exc_retA" in s:
        out.update(_bond_moments_per_sim(s))
    return out


def compute_table2_moments(s, *, bg_yss):
    """Compute the 12 deliverable Table-2 moments, averaged across sims.

    `s` is the dict of named (n_sims, T) series from `build_table2_series`.
    Returns dict[moment_number -> float].
    """
    return {k: float(v.mean()) for k, v in
            compute_table2_moments_per_sim(s, bg_yss=bg_yss).items()}


# KL NFA-decomposition targets, BENCHMARK (spec-1) columns.
#   Table-3 memo (% of 4y): "Model" col = benchmark → (k-kappa),b_H,b_F (sum≈nfa/4y=-23).
#   Table-10 portfolio shares (% of a=h_sav): columns are specs [7,2,8,9,1], so the
#   benchmark (spec 1) is the LAST column → k/a,b_H/a,b_F/a (sum 100).
KL_NFA_DECOMP = {
    "t3_kmk": 59.8, "t3_bH": -102.5, "t3_bF": 19.8,
    "t10_k": 142.41, "t10_bH": -51.94, "t10_bF": 9.53,
}


def nfa_decomposition(s):
    """Localize the NFA level: Table-3 memo (k-kappa, b_H, b_F as % of 4y) and
    Table-10 portfolio shares (k, b_H, b_F as % of total savings a). Per-sim then
    averaged. Matches calc_moments.m table_3 rows 4-6 and table_10 rows 1-3."""
    h_bf_sav = (s["h_sav"] - s["h_ksav"]) - s["h_bh_sav"]
    fy = 4.0 * s["yh"]
    per_sim = {
        "t3_kmk": 100.0 * ((s["h_ksav"] - s["h_kap"]) / fy).mean(dim=1),
        "t3_bH": 100.0 * (s["h_bh_sav"] / fy).mean(dim=1),
        "t3_bF": 100.0 * (h_bf_sav / fy).mean(dim=1),
        "t10_k": 100.0 * (s["h_ksav"] / s["h_sav"]).mean(dim=1),
        "t10_bH": 100.0 * (s["h_bh_sav"] / s["h_sav"]).mean(dim=1),
        "t10_bF": 100.0 * (h_bf_sav / s["h_sav"]).mean(dim=1),
    }
    return {k: float(v.mean()) for k, v in per_sim.items()}


def format_nfa_decomposition(dec):
    """Render the NFA decomposition vs KL (localizes the nfa/4y gap)."""
    rows = [
        ("Table-3 (k-kappa)/4y %", "t3_kmk"),
        ("Table-3 b_H/4y %", "t3_bH"),
        ("Table-3 b_F/4y %", "t3_bF"),
        ("  sum = nfa/4y %", None),
        ("Table-10 k/a %", "t10_k"),
        ("Table-10 b_H/a %", "t10_bH"),
        ("Table-10 b_F/a %", "t10_bF"),
    ]
    lines = ["", "NFA decomposition vs KL (localizes the nfa/4y gap):",
             f"  {'component':<24} {'model':>10} {'KL':>10}"]
    nfa_sum = dec["t3_kmk"] + dec["t3_bH"] + dec["t3_bF"]
    kl_sum = KL_NFA_DECOMP["t3_kmk"] + KL_NFA_DECOMP["t3_bH"] + KL_NFA_DECOMP["t3_bF"]
    for label, key in rows:
        if key is None:
            lines.append(f"  {label:<24} {nfa_sum:>10.2f} {kl_sum:>10.2f}")
        else:
            lines.append(f"  {label:<24} {dec[key]:>10.2f} {KL_NFA_DECOMP[key]:>10.2f}")
    return "\n".join(lines)


def format_table2(moments, per_sim=None):
    """Render a comparison table (Model vs KL target) to a string. If `per_sim`
    (dict[int -> (n_sims,) tensor]) is given, also report the cross-sim standard
    deviation and the standard error of the mean (sd/sqrt(n_sims)) — so a gap vs
    KL can be read against sampling error (NFA, the persistent level, has the
    largest SE; the RNG is not bit-reproducible vs KL's NAG stream)."""
    show_se = per_sim is not None
    head = f"  {'#':>3}  {'moment':<22} {'model':>10} {'KL':>8}"
    if show_se:
        head += f" {'sd_sim':>9} {'se_mean':>9} {'(KL-mdl)/se':>12}"
    lines = ["", "KL Table-2 moment comparison (all 15; 8/9/10 use the bond ladder):", head]
    for k in range(1, 16):
        kl = KL_TARGETS[k]
        if k in moments:
            row = f"  {k:>3}  {LABELS[k]:<22} {moments[k]:>10.3f} {kl:>8.2f}"
            if show_se and k in per_sim:
                v = per_sim[k]
                n = v.shape[0]
                sd = float(v.std(unbiased=True))
                se = sd / (n ** 0.5)
                z = (kl - moments[k]) / se if se > 0 else float("nan")
                row += f" {sd:>9.3f} {se:>9.3f} {z:>12.2f}"
            lines.append(row)
        else:
            lines.append(f"  {k:>3}  {LABELS[k]:<22} {'(deferred)':>10} {kl:>8.2f}")
    return "\n".join(lines)
