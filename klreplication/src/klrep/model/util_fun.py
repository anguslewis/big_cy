"""Felicity, its consumption derivative, and the labor part — port of
`mod_calc.f90 :: util_fun`.

Trabandt-Uhlig / Shimer consumption-labor form. For IES != 1:

    labor_part   = ( 1 + (1/ies - 1) * chi0 * chi/(1+chi) * l^((1+chi)/chi) )^(1/ies)
    util         = c^(1-1/ies) * labor_part
    util_c_deriv = util / c

For IES == 1 (log case): labor_part=1, util = log c - chi0*chi/(1+chi) l^((1+chi)/chi),
util_c_deriv = 1/c. `chi0` is the per-country disutility level from calc_steady;
`chi` is the global Frisch parameter. Elementwise (floats or torch tensors).
"""
import math

SQRT_EPS = math.sqrt(2.220446049250313e-16)


def util_fun(consumption, labor, chi0, ies, chi):
    """Return (util, util_c_deriv, labor_part). `ies`/`chi0` are per-country
    scalars (the log vs. non-log branch is chosen on `ies`); `consumption` and
    `labor` may be floats or torch tensors."""
    labor_pow = labor ** ((1.0 + chi) / chi)
    if 1.0 - SQRT_EPS < ies < 1.0 + SQRT_EPS:  # log utility
        log_c = consumption.log() if hasattr(consumption, "log") else math.log(consumption)
        util = log_c - chi0 * chi / (1.0 + chi) * labor_pow
        labor_part = labor * 0.0 + 1.0 if hasattr(labor, "shape") else 1.0
        return util, 1.0 / consumption, labor_part
    labor_part = (
        1.0 + (1.0 / ies - 1.0) * chi0 * chi / (1.0 + chi) * labor_pow
    ) ** (1.0 / ies)
    util = consumption ** (1.0 - 1.0 / ies) * labor_part
    util_c_deriv = util / consumption
    return util, util_c_deriv, labor_part
