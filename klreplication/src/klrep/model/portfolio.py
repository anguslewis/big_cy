"""Portfolio return and the portfolio FOC residual — port of
`mod_calc.f90 :: portfolio_return` and `portfolio_foc` (4388-4410).

These are the kernel of the bond-clearing / portfolio-share Brent solves. Batched
(tensor-native): the quadrature axis is the LAST dim; the leading dims batch over
grid points (and, where relevant, agents). `bond_share`/`foreign_share`/`gma`/
`ies` broadcast against the (..., n_quad) return tensors.
"""
import torch


def portfolio_return(bond_share, foreign_share, rh, rf, rk):
    """r_alpha = (bond_share - foreign_share)*rh + foreign_share*rf
                 + (1 - bond_share)*rk.   (mod_calc.f90:4402-4410)

    `rh` is the *home* return already including the convenience wedge
    (rf_home = rf_vec[...,0] / (1 - omg)); `rf` is the foreign return; `rk` the
    capital return. Shapes broadcast; returns (..., n_quad)."""
    return (bond_share - foreign_share) * rh + foreign_share * rf + (1.0 - bond_share) * rk


def portfolio_foc(next_period_wealth, v, mc, big_weight_vec, r1, r2,
                  c_cost_growth, gma, ies):
    """Scaled portfolio FOC (mod_calc.f90:4388-4400):

        sdf = v^(1/ies - gma) * mc
        foc = sum_q w * sdf * nw^(-gma) * (r1 - r2) / c_cost_growth
              / | sum_q w * sdf |

    The `/|sum w*sdf|` is just numerical scaling. All tensors are (..., n_quad);
    the sum is over the last (quadrature) axis. Returns (...,)."""
    sdf = v ** (1.0 / ies - gma) * mc
    num = (big_weight_vec * sdf * next_period_wealth ** (-gma)
           * (r1 - r2) / c_cost_growth).sum(dim=-1)
    den = (big_weight_vec * sdf).sum(dim=-1).abs()
    return num / den
