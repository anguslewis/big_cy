"""Isotropic Smolyak element (multi-index) set — JMMV (2014) construction.

Port of `mod_smolyak.f90 :: Smolyak_Elem_Isotrop`. Returns the integer matrix of
multidimensional subindices (n_elem x d), 1-based, for the isotropic level-`mu`
Smolyak rule on [-1,1]^d. The SAME matrix indexes both the grid points (via
smolyak_grid) and the basis terms (via smolyak_polynomial), so the basis matrix
evaluated on the grid is square.

Construction (disjoint-set form): enumerate level vectors i = (i_1,...,i_d),
i_k >= 1, with d <= sum(i) <= d+mu; for each, expand every dimension's level i_k
into the *new* 1-D subindices A(i_k) it introduces, and take the Cartesian product
across dimensions. The disjoint sets guarantee each multi-index appears once.

This is a setup-time host computation (NumPy/itertools), not a hot path.
"""
from itertools import product

import numpy as np


def _disjoint_indices(level: int):
    """A(i): the NEW 1-D subindices (1-based) introduced at level `level`.

    level=1 -> {1}; level=2 -> {2,3}; level>=3 -> {2^(l-2)+2, ..., 2^(l-1)+1}.
    Matches the 1-D node growth m(i) = 2^(i-1)+1 (1,3,5,9,...).
    """
    if level == 1:
        return [1]
    if level == 2:
        return [2, 3]
    lo = 2 ** (level - 2) + 2
    hi = 2 ** (level - 1) + 1
    return list(range(lo, hi + 1))


def _compositions(total: int, d: int):
    """Yield all length-d non-negative integer tuples summing to `total`."""
    if d == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, d - 1):
            yield (first,) + rest


def smolyak_elem_isotrop(d: int, mu: int) -> np.ndarray:
    """Return the (n_elem x d) int64 matrix of 1-based subindices, isotropic mu."""
    if d < 1 or mu < 0:
        raise ValueError(f"need d>=1, mu>=0; got d={d}, mu={mu}")
    rows = []
    # Level vector i_k = a_k + 1 with sum(a_k) = total, total = 0..mu.
    for total in range(mu + 1):
        for a in _compositions(total, d):
            sets = [_disjoint_indices(ak + 1) for ak in a]
            rows.extend(product(*sets))
    return np.array(rows, dtype=np.int64)
