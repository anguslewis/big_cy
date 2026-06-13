"""Smolyak interpolation — evaluate at new points.

Port of the `dgemm` step in `mod_calc.f90`: given coefficients c (from interp_fit)
and the basis B_new evaluated at off-grid next-period states (from
smolyak_polynomial), the interpolated values are the single batched matmul
B_new @ c. This is the operation that dominates the solver, so it stays a dense
matmul — GPU-friendly, no Python loops.
"""
import torch


def interp_eval(B_new: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    """Interpolated values at the new points: B_new (M x T) @ c (T x k) -> (M x k).

    `c` may be (T,) for a single function or (T, k) for k stacked functions.
    """
    return B_new @ c
