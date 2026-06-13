"""CES consumption bundle — port of `mod_calc.f90 :: calc_bundle`.

    c = ( varsigma^(1/sigma) c_h^((s-1)/s) + (1-varsigma)^(1/sigma) c_f^((s-1)/s) )^(s/(s-1))

Elementwise: works on Python floats or torch tensors (used both at steady state
and, batched over grid points, in the equilibrium step).
"""


def calc_bundle(c_h, c_f, varsigma, sigma):
    e = (sigma - 1.0) / sigma
    return (
        varsigma ** (1.0 / sigma) * c_h ** e
        + (1.0 - varsigma) ** (1.0 / sigma) * c_f ** e
    ) ** (sigma / (sigma - 1.0))
