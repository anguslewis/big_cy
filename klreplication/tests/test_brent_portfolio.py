"""Validation — batched Brent root-finder and the portfolio kernel."""
import torch

from klrep.solve.brent import batched_brent, batched_bracket_expand
from klrep.model.portfolio import portfolio_return, portfolio_foc


def test_brent_quadratic_batch():
    targets = torch.tensor([0.25, 1.0, 4.0, 9.0, 2.0], dtype=torch.float64)
    f = lambda x: x ** 2 - targets
    a = torch.zeros_like(targets)
    b = torch.full_like(targets, 5.0)
    root = batched_brent(f, a, b)
    assert torch.allclose(root, targets.sqrt(), atol=1e-12)


def test_brent_transcendental_batch():
    targets = torch.tensor([1.0, 2.0, 5.0, 10.0], dtype=torch.float64)
    f = lambda x: x.exp() - targets
    a = torch.full_like(targets, -2.0)
    b = torch.full_like(targets, 3.0)
    root = batched_brent(f, a, b)
    assert torch.allclose(root, targets.log(), atol=1e-10)


def test_bracket_expand_finds_sign_change():
    target = torch.tensor([2.0, -1.0, 0.5], dtype=torch.float64)
    f = lambda x: x - target  # root at x=target
    x0 = torch.zeros_like(target)
    a, b, fa, fb = batched_bracket_expand(f, x0, step=0.1)
    assert torch.all(fa * fb <= 0)  # straddles the root
    root = batched_brent(f, a, b)
    assert torch.allclose(root, target, atol=1e-10)


def test_portfolio_return_formula():
    rh = torch.tensor([[1.01, 1.02]], dtype=torch.float64)
    rf = torch.tensor([[1.00, 1.03]], dtype=torch.float64)
    rk = torch.tensor([[1.05, 0.95]], dtype=torch.float64)
    bond_share = torch.tensor([[0.7]], dtype=torch.float64)
    foreign_share = torch.tensor([[0.2]], dtype=torch.float64)
    r = portfolio_return(bond_share, foreign_share, rh, rf, rk)
    expect = (0.7 - 0.2) * rh + 0.2 * rf + (1 - 0.7) * rk
    assert torch.allclose(r, expect)


def test_portfolio_foc_sign():
    # If r1 > r2 uniformly and SDF, weights > 0, the FOC numerator is positive.
    n_quad = 5
    v = torch.ones(1, n_quad, dtype=torch.float64)
    mc = torch.ones(1, n_quad, dtype=torch.float64)
    w = torch.full((1, n_quad), 1.0 / n_quad, dtype=torch.float64)
    nw = torch.ones(1, n_quad, dtype=torch.float64)
    ccg = torch.ones(1, n_quad, dtype=torch.float64)
    r1 = torch.full((1, n_quad), 1.05, dtype=torch.float64)
    r2 = torch.full((1, n_quad), 1.00, dtype=torch.float64)
    foc = portfolio_foc(nw, v, mc, w, r1, r2, ccg, gma=10.0, ies=0.75)
    assert float(foc) > 0
    # Reversing r1, r2 flips the sign.
    foc2 = portfolio_foc(nw, v, mc, w, r2, r1, ccg, gma=10.0, ies=0.75)
    assert float(foc2) < 0
