"""Standard-normal Gauss-Hermite quadrature nodes and weights.

Faithful port of `mod_normal.f90 :: get_quadrature_points_nb` (the only quadrature
variant the Safety model calls). Golub-Welsch on the **probabilists'** Hermite
Jacobi matrix: symmetric tridiagonal with zero diagonal and off-diagonal
`sqrt(i)` (i = 1..n-1). Eigenvalues are the nodes; the squared first components of
the eigenvectors are the weights. Weights are **probability-normalized** (sum to
1), so for X ~ N(0,1):  E[f(X)] = sum_k w_k f(x_k).

Computed once on the host in float64 (an n x n eig at setup), then returned as
device tensors — never recomputed in the solver hot loop.
"""
import numpy as np
import torch

from ..config import DTYPE


def gauss_hermite_nb(n: int, *, dtype=DTYPE, device=None):
    """Return (nodes, weights) for the standard-normal expectation.

    Parameters
    ----------
    n : int
        Number of quadrature nodes.
    dtype, device : torch dtype / device
        Output tensor dtype (default float64) and device.

    Returns
    -------
    nodes, weights : torch.Tensor, each shape (n,)
        Nodes x_k and probability weights w_k with sum(w_k) == 1.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if n == 1:
        nodes_np = np.array([0.0], dtype=np.float64)
        weights_np = np.array([1.0], dtype=np.float64)
    else:
        # Probabilists' Hermite Jacobi matrix: zero diagonal, off-diag sqrt(i).
        off = np.sqrt(np.arange(1, n, dtype=np.float64))
        jacobi = np.diag(off, 1) + np.diag(off, -1)
        evals, evecs = np.linalg.eigh(jacobi)  # ascending eigenvalues
        nodes_np = evals
        weights_np = evecs[0, :] ** 2  # rows orthonormal => weights sum to 1

    # Sanity: probability weights.
    if not np.isclose(weights_np.sum(), 1.0, atol=1e-12):
        raise AssertionError(f"GH weights sum to {weights_np.sum()}, expected 1")

    return (
        torch.tensor(nodes_np, dtype=dtype, device=device),
        torch.tensor(weights_np, dtype=dtype, device=device),
    )
