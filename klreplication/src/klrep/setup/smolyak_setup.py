"""Smolyak grid wiring — port of the grid-build block of init_setup.

`mod_param.f90:329-347`: restrict the 9 per-dimension levels to the active dims
(mu>0), build the isotropic element set at the max active level, prune to the
anisotropic set, evaluate the collocation grid and the square basis matrix on it.
The basis matrix `smol_polynom` is what the solver factorizes each iteration.

The mapping from the [-1,1] collocation grid to economic state values
(`grid_setup`) is deferred — it needs the steady-state grid centers/half-widths
(capital, wages) computed by calc_steady. See state_grid.py.
"""
from dataclasses import dataclass
from typing import List

import numpy as np
import torch

from ..config import DTYPE
from ..params import Params
from ..grid.smolyak_elem_isotrop import smolyak_elem_isotrop
from ..grid.smolyak_elem_anisotrop import smolyak_elem_anisotrop
from ..grid.smolyak_grid import smolyak_grid
from ..grid.smolyak_polynomial import smolyak_polynomial


@dataclass
class SmolyakSetup:
    elem: np.ndarray          # (n_states, d) anisotropic element matrix (1-based)
    smol_grid: np.ndarray     # (n_states, d) collocation points in [-1,1]
    smol_polynom: torch.Tensor  # (n_states, n_states) basis on the grid (square)
    n_states: int
    mu: int                   # max active per-dim level
    active_dims: List[int]    # 0-based state dims with mu>0, in order
    mus_redux: List[int]      # per-dim levels restricted to active dims


def build_smolyak(p: Params, *, dtype=DTYPE, device=None) -> SmolyakSetup:
    active_dims = [d for d, m in enumerate(p.mus_dims) if m > 0]
    mus_redux = [p.mus_dims[d] for d in active_dims]
    d = len(active_dims)
    if d == 0:
        raise ValueError("no active Smolyak dimensions (all mus == 0)")
    mu = max(mus_redux)

    elem_iso = smolyak_elem_isotrop(d, mu)
    elem = smolyak_elem_anisotrop(elem_iso, mus_redux)
    smol_grid = smolyak_grid(elem, mu)
    grid_t = torch.tensor(smol_grid, dtype=dtype, device=device)
    smol_polynom = smolyak_polynomial(grid_t, elem, mu, dtype=dtype)

    return SmolyakSetup(
        elem=elem,
        smol_grid=smol_grid,
        smol_polynom=smol_polynom,
        n_states=elem.shape[0],
        mu=mu,
        active_dims=active_dims,
        mus_redux=mus_redux,
    )
