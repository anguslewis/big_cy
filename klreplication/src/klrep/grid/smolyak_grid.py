"""Smolyak collocation grid — Chebyshev extrema (Lobatto) nodes.

Port of `mod_smolyak.f90 :: Smolyak_Grid`. Builds the nested 1-D node sequence
`points_1d` indexed by subindex, then gathers the grid points by the element
matrix's subindices. Nodes are Chebyshev *extrema* on [-1,1] (endpoints
included): for level i with m_i = m(i) nodes,

    x_j = -cos( pi * (j-1) / (m_i - 1) ),   j = 1..m_i,    m(1) = 1 -> {0}.

with m(i) = 2^(i-1)+1 for i>1. Levels are concatenated and order-preserving
deduped so subindex n maps to the n-th distinct node in the nested ordering
(matching the disjoint-set indexing in smolyak_elem_isotrop). Setup-time host
computation; returns a NumPy array (model code rescales states into [-1,1]).
"""
import numpy as np

_SNAP_TOL = 1e-12


def _extrema_1d(m: int) -> np.ndarray:
    """Chebyshev-extrema nodes of order m on [-1,1] (m=1 -> [0]), snapped."""
    if m == 1:
        return np.array([0.0])
    j = np.arange(1, m + 1, dtype=np.float64)
    x = -np.cos(np.pi * (j - 1) / (m - 1))
    x[np.abs(x) < _SNAP_TOL] = 0.0
    x[np.abs(x - 1.0) < _SNAP_TOL] = 1.0
    x[np.abs(x + 1.0) < _SNAP_TOL] = -1.0
    return x


def build_points_1d(mu: int) -> np.ndarray:
    """Nested 1-D node sequence; index (1-based) is the Smolyak subindex."""
    seen = []
    for i in range(1, mu + 2):  # level i = 1..mu+1
        m_i = 1 if i == 1 else 2 ** (i - 1) + 1
        for x in _extrema_1d(m_i):
            if not any(abs(x - s) < _SNAP_TOL for s in seen):
                seen.append(float(x))
    return np.array(seen, dtype=np.float64)


def smolyak_grid(elem: np.ndarray, mu: int) -> np.ndarray:
    """Return collocation points (n_elem x d) in [-1,1] for element matrix `elem`.

    `elem` holds 1-based subindices; `mu = max(active per-dim levels)`.
    """
    points_1d = build_points_1d(mu)
    if elem.max() > points_1d.shape[0]:
        raise ValueError(
            f"elem max subindex {elem.max()} exceeds points_1d length "
            f"{points_1d.shape[0]} for mu={mu}"
        )
    return points_1d[elem - 1]  # gather; subindices are 1-based
