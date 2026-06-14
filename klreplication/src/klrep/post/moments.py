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


def compute_table2_moments(s, *, bg_yss):
    """Compute the 12 deliverable Table-2 moments, averaged across sims.

    `s` is the dict of named (n_sims, T) series from `build_table2_series`.
    Returns dict[moment_number -> float].
    """
    ratio = s["yf"] / s["yh"]                                   # y*/y

    per_sim = {
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
    return {k: float(v.mean()) for k, v in per_sim.items()}


def format_table2(moments):
    """Render a comparison table (Model vs KL target) to a string."""
    lines = ["", "KL Table-2 moment comparison (12 deliverable; 8/9/10 deferred):",
             f"  {'#':>3}  {'moment':<22} {'model':>10} {'KL':>8}"]
    for k in range(1, 16):
        kl = KL_TARGETS[k]
        if k in moments:
            lines.append(f"  {k:>3}  {LABELS[k]:<22} {moments[k]:>10.3f} {kl:>8.2f}")
        else:
            lines.append(f"  {k:>3}  {LABELS[k]:<22} {'(deferred)':>10} {kl:>8.2f}")
    return "\n".join(lines)
