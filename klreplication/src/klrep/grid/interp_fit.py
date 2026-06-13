"""Smolyak interpolation — fit coefficients (factor once, solve many RHS).

Port of the `F07ABF` (LAPACK dgesvx) usage in `mod_calc.f90`: factor the square
basis matrix B (evaluated on the grid) once per outer iteration, then solve
B c = v for the coefficient columns c against multiple right-hand sides (policies,
prices, values). The Fortran uses the expert driver (equilibration + condition
estimation); we use LU (lu_factor/lu_solve) in float64 and surface the condition
number so ill-conditioning at high anisotropic levels is caught, not silent.
"""
import torch


def interp_factor(B: torch.Tensor):
    """LU-factor the square basis matrix B for reuse across RHS. Returns lu_piv."""
    return torch.linalg.lu_factor(B)


def interp_solve(lu_piv, v: torch.Tensor) -> torch.Tensor:
    """Solve B c = v using the factorization from `interp_factor`.

    `v` may be (n,) or (n, k); the return matches (coeffs c, same trailing shape).
    """
    LU, piv = lu_piv
    squeeze = v.ndim == 1
    rhs = v.unsqueeze(1) if squeeze else v
    c = torch.linalg.lu_solve(LU, piv, rhs)
    return c.squeeze(1) if squeeze else c


def basis_condition_number(B: torch.Tensor) -> float:
    """2-norm condition number of B (diagnostic for the coefficient solve)."""
    return float(torch.linalg.cond(B))
