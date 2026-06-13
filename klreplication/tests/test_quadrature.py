"""Validation ladder #2 — Gauss-Hermite quadrature.

Self-validating (no Fortran needed): checks the probabilists'-Hermite convention
(weights sum to 1; nodes integrate standard-normal moments exactly), and asserts
the nodes equal the eigenvalues of the sqrt(i)-off-diagonal Jacobi matrix — the
guard against the classic hermegauss-vs-hermgauss sqrt(2) footgun (critic S3).
"""
import math

import numpy as np
import torch

from klrep.quad.gauss_hermite import gauss_hermite_nb


def test_weights_sum_to_one():
    for n in [1, 2, 3, 5, 7, 10]:
        _, w = gauss_hermite_nb(n)
        assert math.isclose(float(w.sum()), 1.0, abs_tol=1e-12)


def test_standard_normal_moments():
    # n-point Gauss-Hermite integrates polynomials up to degree 2n-1 exactly.
    nodes, w = gauss_hermite_nb(5)
    x = nodes
    e = lambda p: float((w * x ** p).sum())
    assert math.isclose(e(0), 1.0, abs_tol=1e-12)   # E[1]
    assert math.isclose(e(1), 0.0, abs_tol=1e-12)   # E[X]
    assert math.isclose(e(2), 1.0, abs_tol=1e-12)   # E[X^2]
    assert math.isclose(e(3), 0.0, abs_tol=1e-12)   # E[X^3]
    assert math.isclose(e(4), 3.0, abs_tol=1e-12)   # E[X^4]
    assert math.isclose(e(6), 15.0, abs_tol=1e-10)  # E[X^6], within deg 2n-1=9


def test_nodes_match_jacobi_eigenvalues():
    # The defining identity: nodes are the eigenvalues of the probabilists'
    # Hermite Jacobi matrix (zero diagonal, off-diagonal sqrt(i)).
    n = 7
    nodes, _ = gauss_hermite_nb(n)
    off = np.sqrt(np.arange(1, n))
    jac = np.diag(off, 1) + np.diag(off, -1)
    evals = np.sort(np.linalg.eigvalsh(jac))
    assert np.allclose(nodes.numpy(), evals, atol=1e-12)


def test_matches_numpy_hermite_e_convention():
    # numpy probabilists' Gauss-Hermite, weights normalized by sqrt(2*pi).
    n = 6
    nodes, w = gauss_hermite_nb(n)
    np_nodes, np_w = np.polynomial.hermite_e.hermegauss(n)
    assert np.allclose(nodes.numpy(), np_nodes, atol=1e-12)
    assert np.allclose(w.numpy(), np_w / math.sqrt(2 * math.pi), atol=1e-12)


def test_dtype_is_float64():
    nodes, w = gauss_hermite_nb(3)
    assert nodes.dtype == torch.float64
    assert w.dtype == torch.float64
