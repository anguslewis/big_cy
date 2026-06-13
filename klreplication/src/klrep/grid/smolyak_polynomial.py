"""Smolyak interpolation basis — first-kind Chebyshev tensor products.

Port of `mod_smolyak.f90 :: Smolyak_Polynomial` / `Smolyak_Polynomial2` (one
function here; the Fortran split is just allocation hygiene). Builds the basis
matrix B (n_points x n_terms): each term (a row of the element matrix) is the
tensor product over dimensions of the 1-D first-kind Chebyshev polynomial
T_{n-1}, where n is that dimension's subindex (so column subindex n <-> degree
n-1). T is built by the recurrence T_0=1, T_1=x, T_j = 2x T_{j-1} - T_{j-2}.

This is the solver hot path: vectorized as a batched Chebyshev recurrence plus an
advanced-indexing gather and a product reduction — no Python loop over points,
terms, or dimensions (only over polynomial degree, <= ~17). Pure tensor op, so it
moves to GPU unchanged.
"""
import torch

from ..config import DTYPE


def smolyak_polynomial(points, elem, mu: int, *, dtype=DTYPE) -> torch.Tensor:
    """Evaluate the Smolyak basis at `points`.

    Parameters
    ----------
    points : torch.Tensor, shape (N, d), values in [-1, 1].
    elem : torch.Tensor (int) or array-like, shape (T, d), 1-based subindices.
    mu : int
        Max active per-dim level; the max Chebyshev degree is m_max-1 = 2^mu.

    Returns
    -------
    B : torch.Tensor, shape (N, T) — basis matrix (B[i, t] = term t at point i).
    """
    if not torch.is_tensor(points):
        points = torch.as_tensor(points, dtype=dtype)
    points = points.to(dtype)
    N, d = points.shape

    elem_t = torch.as_tensor(elem, dtype=torch.long, device=points.device)
    if elem_t.shape[1] != d:
        raise ValueError(
            f"elem has {elem_t.shape[1]} cols but points have d={d}"
        )

    m_max = 2 ** mu + 1  # number of 1-D Chebyshev polynomials T_0..T_{m_max-1}
    if int(elem_t.max()) > m_max:
        raise ValueError(
            f"elem max subindex {int(elem_t.max())} exceeds m_max={m_max} (mu={mu})"
        )

    # 1-D Chebyshev T_0..T_{m_max-1} at every point/dim -> phi (N, d, m_max).
    phi = torch.empty((N, d, m_max), dtype=dtype, device=points.device)
    phi[:, :, 0] = 1.0
    if m_max > 1:
        phi[:, :, 1] = points
    for j in range(2, m_max):
        phi[:, :, j] = 2.0 * points * phi[:, :, j - 1] - phi[:, :, j - 2]

    # Gather the per-term, per-dim Chebyshev factor, then product over dims.
    dim_idx = torch.arange(d, device=points.device).unsqueeze(0)  # (1, d)
    factors = phi[:, dim_idx, elem_t - 1]  # (N, T, d) via advanced indexing
    return factors.prod(dim=-1)  # (N, T)
