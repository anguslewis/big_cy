"""Validation ladder #1 — Smolyak grid, basis, and interpolation.

Self-validating (no Fortran needed): element-count checks against known JMMV
values, the nested Chebyshev-extrema node ordering, first-kind Chebyshev basis
identity vs numpy, anisotropic pruning, and — the keystone — interpolation
exactness for any function in the Smolyak polynomial space (fit on grid, evaluate
off grid, recover exactly).
"""
import numpy as np
import torch

from klrep.grid.smolyak_elem_isotrop import smolyak_elem_isotrop
from klrep.grid.smolyak_elem_anisotrop import smolyak_elem_anisotrop
from klrep.grid.smolyak_grid import smolyak_grid, build_points_1d
from klrep.grid.smolyak_polynomial import smolyak_polynomial
from klrep.grid.interp_fit import interp_factor, interp_solve, basis_condition_number
from klrep.grid.interp_eval import interp_eval


def test_isotrop_counts_1d():
    # d=1: number of points = m(mu+1) = 1 if mu==0 else 2^mu + 1.
    for mu in range(5):
        elem = smolyak_elem_isotrop(1, mu)
        expected = 1 if mu == 0 else 2 ** mu + 1
        assert elem.shape == (expected, 1)


def test_isotrop_counts_2d_known_jmmv():
    # Known Smolyak (Chebyshev-extrema) point counts in 2D.
    assert smolyak_elem_isotrop(2, 1).shape[0] == 5
    assert smolyak_elem_isotrop(2, 2).shape[0] == 13
    assert smolyak_elem_isotrop(2, 3).shape[0] == 29


def test_elements_are_unique():
    for d, mu in [(2, 3), (3, 2), (4, 2)]:
        elem = smolyak_elem_isotrop(d, mu)
        uniq = np.unique(elem, axis=0)
        assert uniq.shape[0] == elem.shape[0]


def test_points_1d_nested_ordering():
    # Subindex -> node: 1->0, 2->-1, 3->1, 4->-cos(pi/4), 5->+cos(pi/4), ...
    p = build_points_1d(2)  # 2^2+1 = 5 nodes
    expected = np.array([0.0, -1.0, 1.0, -np.cos(np.pi / 4), np.cos(np.pi / 4)])
    assert np.allclose(p, expected, atol=1e-12)
    # Endpoints snapped exactly.
    assert p[1] == -1.0 and p[2] == 1.0 and p[0] == 0.0


def test_grid_1d_matches_points_1d():
    mu = 3
    elem = smolyak_elem_isotrop(1, mu)
    grid = smolyak_grid(elem, mu)
    p = build_points_1d(mu)
    assert np.allclose(grid[:, 0], p[elem[:, 0] - 1], atol=1e-14)


def test_basis_is_first_kind_chebyshev_1d():
    # For d=1, column t of B is T_{degree}(x), degree = subindex-1.
    mu = 3
    elem = smolyak_elem_isotrop(1, mu)
    x = np.linspace(-1, 1, 11)
    B = smolyak_polynomial(torch.tensor(x).unsqueeze(1), elem, mu).numpy()
    degrees = elem[:, 0] - 1
    for t, deg in enumerate(degrees):
        coef = np.zeros(deg + 1)
        coef[deg] = 1.0
        T_deg = np.polynomial.chebyshev.chebval(x, coef)
        assert np.allclose(B[:, t], T_deg, atol=1e-12)


def test_basis_square_and_well_conditioned():
    for d, mu in [(2, 3), (3, 2), (4, 2)]:
        elem = smolyak_elem_isotrop(d, mu)
        grid = torch.tensor(smolyak_grid(elem, mu))
        B = smolyak_polynomial(grid, elem, mu)
        assert B.shape[0] == B.shape[1] == elem.shape[0]  # square
        cond = basis_condition_number(B)
        assert np.isfinite(cond) and cond < 1e6  # well conditioned


def test_anisotropic_pruning():
    d, mu = 2, 2
    elem_iso = smolyak_elem_isotrop(d, mu)
    # mus=[2,0]: dim 2 collapses to subindex 1 only; dim 1 keeps up to 2^2+1=5.
    elem_ani = smolyak_elem_anisotrop(elem_iso, [2, 0])
    assert (elem_ani[:, 1] == 1).all()
    assert elem_ani[:, 0].max() <= 5
    # mus=[2,2] reproduces the isotropic set (caps not binding within mu).
    elem_full = smolyak_elem_anisotrop(elem_iso, [2, 2])
    assert elem_full.shape[0] == elem_iso.shape[0]


def test_interpolation_exact_on_polynomial_space():
    # Any function spanned by the Smolyak basis is recovered exactly: build
    # f = B @ c_true, fit c_hat from grid values, check c_hat == c_true and that
    # evaluation at random off-grid points matches the truth.
    torch.manual_seed(0)
    for d, mu in [(2, 3), (3, 2)]:
        elem = smolyak_elem_isotrop(d, mu)
        T = elem.shape[0]
        grid = torch.tensor(smolyak_grid(elem, mu))
        B_grid = smolyak_polynomial(grid, elem, mu)

        c_true = torch.randn(T, dtype=torch.float64)
        v_grid = B_grid @ c_true  # function values at grid points

        lu = interp_factor(B_grid)
        c_hat = interp_solve(lu, v_grid)
        assert torch.allclose(c_hat, c_true, atol=1e-9)

        # Off-grid evaluation matches the true function.
        pts = (torch.rand(50, d, dtype=torch.float64) * 2 - 1)
        B_new = smolyak_polynomial(pts, elem, mu)
        f_hat = interp_eval(B_new, c_hat)
        f_true = B_new @ c_true
        assert torch.allclose(f_hat, f_true, atol=1e-9)


def test_interp_solve_multiple_rhs():
    # Stacked RHS (k columns) solved together, mirroring the solver's reuse.
    d, mu, k = 3, 2, 4
    elem = smolyak_elem_isotrop(d, mu)
    grid = torch.tensor(smolyak_grid(elem, mu))
    B = smolyak_polynomial(grid, elem, mu)
    torch.manual_seed(1)
    C = torch.randn(elem.shape[0], k, dtype=torch.float64)
    V = B @ C
    lu = interp_factor(B)
    C_hat = interp_solve(lu, V)
    assert torch.allclose(C_hat, C, atol=1e-9)
