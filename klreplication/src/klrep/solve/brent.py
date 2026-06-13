"""Batched Brent root-finder — vectorized version of the inline Brent solves in
`mod_calc.f90` (inverse-quadratic / secant / bisection with the mflag guard).

Solves a batch of independent scalar problems f(x)=0 simultaneously: `f` maps a
(batch,) tensor to a (batch,) tensor, and brackets `a`, `b` are (batch,). This is
how the per-grid-point nom_i and portfolio-share solves are batched in the
tensor-native solver. Same algorithm as the Fortran, so converged roots match to
tolerance (we do not need bit-identical iterates — the outer backward iteration is
itself a fixed point).
"""
import torch

_DELTA = 5e-16


def batched_brent(f, a, b, *, xtol=1e-14, max_iter=200):
    """Return roots (batch,) of f on the per-element brackets [a, b].

    Assumes f(a) and f(b) have opposite signs elementwise (caller handles the
    constraint-binding / non-bracketed cases). Iterates a fixed max_iter with a
    convergence mask so converged elements stop moving."""
    a = a.clone()
    b = b.clone()
    fa = f(a)
    fb = f(b)

    # Ensure b is the better iterate (|fb| <= |fa|).
    swap = fa.abs() < fb.abs()
    a, b = torch.where(swap, b, a), torch.where(swap, a, b)
    fa, fb = torch.where(swap, fb, fa), torch.where(swap, fa, fb)

    c = a.clone()
    fc = fa.clone()
    d = c.clone()
    mflag = torch.ones_like(a, dtype=torch.bool)

    for _ in range(max_iter):
        # Inverse-quadratic interpolation where possible, else secant.
        denom_iqi = (fa - fb) * (fa - fc)
        use_iqi = (fa != fc) & (fc != fb)
        s_iqi = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
        # guard division by zero in secant
        df = fb - fa
        df = torch.where(df == 0, torch.full_like(df, 1e-300), df)
        s_sec = b - fb * (b - a) / df
        s = torch.where(use_iqi, s_iqi, s_sec)

        lo = (3.0 * a + b) / 4.0
        not_between = ~(((s > lo) & (s < b)) | ((s < lo) & (s > b)))
        cond2 = mflag & ((s - b).abs() >= (b - c).abs() / 2.0)
        cond3 = (~mflag) & ((s - b).abs() >= (b - d).abs() / 2.0)
        cond4 = mflag & ((b - c).abs() < _DELTA)
        cond5 = (~mflag) & ((b - d).abs() < _DELTA)
        bisect = not_between | cond2 | cond3 | cond4 | cond5
        s = torch.where(bisect, (a + b) / 2.0, s)
        mflag = bisect

        fs = f(s)

        d = c
        c = b
        fc = fb

        left = (fa * fs) < 0  # root in [a, s] -> keep a, set b=s; else set a=s
        b = torch.where(left, s, b)
        fb = torch.where(left, fs, fb)
        a = torch.where(left, a, s)
        fa = torch.where(left, fa, fs)

        swap = fa.abs() < fb.abs()
        a, b = torch.where(swap, b, a), torch.where(swap, a, b)
        fa, fb = torch.where(swap, fb, fa), torch.where(swap, fa, fb)

        if bool(((a - b).abs() <= xtol).all()) or bool((fb == 0).all()):
            break

    return b


def batched_bracket_expand(f, x0, *, step, max_iter=5000):
    """Expand a bracket around x0 until f changes sign, elementwise.

    Steps x in the direction that reduces |f| toward a sign change: where f(x0)>0
    we move by -step, where f(x0)<0 by +step (matching the Fortran nom_i search).
    Returns (a, b, fa, fb) brackets straddling the root. Elements already at a
    root (f==0) get a degenerate bracket [x0, x0]."""
    a = x0.clone()
    fa = f(a)
    direction = torch.where(fa > 0, torch.full_like(fa, -1.0), torch.full_like(fa, 1.0))
    b = a.clone()
    fb = fa.clone()
    active = fa != 0
    for _ in range(max_iter):
        if not bool(active.any()):
            break
        b_try = torch.where(active, b + direction * step, b)
        fb_try = f(b_try)
        # keep previous bracket endpoint a = last b before crossing
        crossed = active & ((fa * fb_try) <= 0)
        # advance b for active, not-yet-crossed elements
        a = torch.where(active & ~crossed, b, a)
        fa = torch.where(active & ~crossed, fb, fa)
        b = torch.where(active, b_try, b)
        fb = torch.where(active, fb_try, fb)
        active = active & ~crossed
    return a, b, fa, fb
