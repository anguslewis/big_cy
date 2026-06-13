"""Anisotropic Smolyak element set — prune the isotropic set per-dimension.

Port of `mod_smolyak.f90 :: Smolyak_Elem_Anisotrop`. Given per-dimension levels
mu_k (`vector_mus_dimensions` restricted to active dims), keep only the rows of
the isotropic element matrix whose subindex in every dimension is within that
dimension's cap. The cap is the max 1-D subindex allowed at level mu_k:

    cap_k = 1                if mu_k == 0   (dimension collapses to the center node)
    cap_k = 2^{mu_k} + 1     otherwise      (= m(mu_k + 1))

Fully vectorized boolean mask. Setup-time host computation.
"""
import numpy as np


def smolyak_elem_anisotrop(elem_iso: np.ndarray, mus) -> np.ndarray:
    """Return the anisotropic subset of `elem_iso` (n x d) for per-dim levels `mus`.

    `mus` has length d (the active dimensions); `elem_iso` has d columns.
    """
    mus = np.asarray(mus, dtype=np.int64)
    if elem_iso.shape[1] != mus.shape[0]:
        raise ValueError(
            f"elem_iso has {elem_iso.shape[1]} cols but mus has length {mus.shape[0]}"
        )
    caps = np.where(mus == 0, 1, 2 ** mus + 1)  # (d,) max subindex per dimension
    mask = (elem_iso <= caps).all(axis=1)
    return elem_iso[mask]
