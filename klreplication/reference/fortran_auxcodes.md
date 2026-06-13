# Kekre–Lenel (2024 AER) — AuxCodes numerical-library digest

Scope: the **four** AuxCodes files that `runme.sh` actually compiles. The
other AuxCodes files (`sandia_rules`, `sgmga`, `sparse_grid_*`, `mod_bspline`,
`mod_sparsequad`, `*_test`) are **not compiled** and are ignored here.

Source root (read-only):
`…/big_cy/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/Code/src/fortran/AuxCodes/`

All line references below are `file:line`. NAG dependencies are flagged because
the Python port must replace them (NAG = proprietary). `dp = nag_wp` is just
IEEE double (float64); see §1.

---

## 0. How these are used by the model (call-site map)

From `mod_param.f90` and `mod_calc.f90` (outside AuxCodes):

- **Grid build** (`mod_param.f90:343-346`), one-time at setup:
  ```fortran
  smolyak_elem_iso = Smolyak_Elem_Isotrop(n_active_dims, maxval(mus_dimensions_redux))
  smol_elem_ani    = Smolyak_Elem_Anisotrop(smolyak_elem_iso, n_active_dims, mus_dimensions_redux)
  smol_grid        = Smolyak_Grid(n_active_dims, maxval(mus_dimensions_redux), smol_elem_ani)
  smol_polynom     = Smolyak_Polynomial(smol_grid, n_active_dims, maxval(mus_dimensions_redux), smol_elem_ani)
  ```
  - `n_active_dims = count(vector_mus_dimensions > 0)`; only dims with µ>0 enter
    the Smolyak machinery. `mus_dimensions_redux` is `vector_mus_dimensions`
    restricted to the active dims (`mod_param.f90:329-343`).
  - `vector_mus_dimensions` has length `smolyak_d` (max 9), read from the
    param file (`mod_param.f90:224-233`). So Smolyak dimension `d = n_active_dims ≤ 9`,
    and the level is **per-dimension** (anisotropic).
  - `smol_polynom` is the **square** basis matrix `B` evaluated at the grid
    itself (`n_states × n_states`, where `n_states = #collocation points = #basis terms`).

- **Quadrature** (`mod_param.f90:358,366,374,381`), one-time:
  ```fortran
  call get_quadrature_points(0.0_dp, 1.0_dp, n_uni_quad_z, quad_vec_z, uni_weight_vec_z)
  ```
  Four scalar shocks (z, zf, disaster p, omega). The 4-arg call `(mu=0, sigma=1, n, points, weights)`
  resolves to **`get_quadrature_points_nb`** (no bounds) = Gauss–Hermite for a
  standard-normal expectation. See §2.

- **Interpolation in the solver hot loop** (`mod_calc.f90`), repeated every iteration:
  1. Solve for coefficients: `F07ABF` (LAPACK `dgesvx`) solves
     `smol_polynom · smol_coeffs = values` for `smol_coeffs`
     (`mod_calc.f90:465, 725, 810, 886, 919, 1044`).
  2. Evaluate at next-period states: `Smolyak_Polynomial2(next_state, …)` builds
     basis at new points, then `dgemm`/matmul `polyn_points · smol_coeffs`
     gives interpolated values (`mod_calc.f90:493-499, 753-763, …`).
  `mod_results.f90` uses `Smolyak_Polynomial` (the allocatable variant) the same way.

- **Markov:** `make_markov_2d/3d` / `calc_rouwen` are **public** in `mod_markov`
  but ⚠ I found **no call sites** in `mod_param`/`mod_calc`/`mod_results`/`main`.
  All four model shocks go through Gauss–Hermite `get_quadrature_points`, not the
  Markov discretization. The Markov module appears to be a dormant/legacy path
  for this Safety model. **Lower port priority** — port only if a later config
  switches a shock to a discrete chain. (Documented in §3 for completeness.)

---

## 1. `base_lib.f90` — primitives

Module `base_lib` (`base_lib.f90:1`). Public list at `:6-7`.

### Kind & constants
| Entity | Def | Value | Python |
|---|---|---|---|
| `dp` | `:9` `dp = nag_wp` | IEEE double | `float64` / `torch.float64` |
| `pi` | `:11` `4*atan(1)` | π | `math.pi` / `torch.pi` |
| `eps` | `:12` `epsilon(eps)` | ≈2.22e-16 | `np.finfo(np.float64).eps` |
| `sqrt_eps` | `:13` `sqrt(eps)` | ≈1.49e-8 | `sqrt(eps)` |

### Utility routines (public)
- **`Fill_linspace_dp(xmin, xmax, x)`** `:18`. In-place fill of `x(:)` with
  `n=size(x)` evenly spaced points, endpoints exact. n=1 requires xmin==xmax.
  → `np.linspace(xmin, xmax, n)` / `torch.linspace`.
- **`m_choose_r(m, r)`** `:41` (+ helper `m_choose_r_second :60`). Integer
  binomial C(m,r) via recursion with `r=min(m-r, r)` symmetry. Used only inside
  `mod_normal` truncated-moment recursions. → `math.comb(m, r)`.
- **`choldc(a, n, p, ifail)`** `:71`. Cholesky (lower) à la Numerical Recipes.
  `a(n,n)` in/out (subdiagonal of L written into lower triangle of `a`), diagonal
  of L returned in `p(n)`. ⚠ On non-PD pivot it sets `ifail=1` and **fudges**
  `sum = sqrt(sum**2) = |sum|` rather than failing — silent regularization.
  Used by `mod_normal`'s `make_points_and_weights`. → `torch.linalg.cholesky`
  (but replicate the fudge only if matching exact non-PD behavior; normally PD).
- **`QsortC(A)` / `Partition`** `:98,:109`. Recursive quicksort, real 1-D,
  in-place ascending. → `torch.sort` / `np.sort`.
- **`calc_var(vec, n)`** `:156`, **`calc_cov(vec1, vec2, n)`** `:168` (and private
  `calc_mean :147`). Sample variance/cov with **(n−1)** denominator. → `np.var(ddof=1)`,
  `np.cov`.
- **`remove_dups_2d(input_mat)`** `:182`. Removes duplicate **rows** of an
  integer matrix, **preserving first-seen order** (no sort). Critical for Smolyak
  index construction. → stable unique on rows: e.g.
  `np.unique(a, axis=0)` **does not preserve order** — instead use a seen-set or
  `pd`-style first-occurrence dedup. ⚠ Order matters downstream (grid/point order).
- **`remove_dups_1d_real(input_mat)`** `:212`. Same, for a real vector
  (exact `==` match), first-seen order. Used in `Smolyak_Grid` to build the 1-D
  node set. → order-preserving unique on floats (exact equality — the rounding in
  `Smolyak_Grid` makes exact dedup safe).
- **`initBrent(...)`** `:242`, **`updateBrent(...)`** `:279`. Brent root-finder
  (inverse-quadratic / secant / bisection) operating on a 4-slot history. Used by
  the model's scalar root solves (not the numerical-library hot path). → SciPy
  `brentq` for a drop-in; reimplement only if you need bit-identical iterates.

---

## 2. `mod_normal.f90` — Gauss quadrature + truncated-normal moments

Module `mod_normal` (`mod_normal.f90:1`). Public: `trunc_normal_moment_{ub,lb,lb_ub}`,
`normal_pdf`, **`get_quadrature_points`** (`:8-9`).
NAG deps: `S15ABF` (normal CDF), `F08FAF` (= LAPACK `dsyev`, symmetric eig),
`F08JEF` (= LAPACK `dstev`, symmetric tridiagonal eig).

### `get_quadrature_points` — generic interface (`:11-14`)
Resolves by argument signature to one of four procedures:
| Variant | Args | Meaning |
|---|---|---|
| `_nb` | `(mu, sigma, n, points, weights)` | **no bounds** — plain Gauss–Hermite for N(µ,σ²) |
| `_ub` | `(mu, sigma, lb=-99, ub, n, …)` | upper-truncated normal |
| `_lb` | `(mu, sigma, lb, ub=-99, n, …)` | lower-truncated normal |
| `_lb_ub` | `(mu, sigma, lb, ub, n, …)` | doubly-truncated normal |

**The model only ever calls the 5-argument `_nb` form** (`get_quadrature_points(0.0_dp, 1.0_dp, n, pts, wts)`),
i.e. standard-normal Gauss–Hermite. The truncated variants are compiled but unused
by the Safety model paths I checked. ⚠ Port `_nb` first; truncated variants are
optional.

#### `get_quadrature_points_nb(mu, sigma, n, points, weights)` `:19-43` — **port this exactly**
Golub–Welsch via the **physicists'-vs-probabilists'** Hermite Jacobi matrix:
- Builds symmetric tridiagonal `J` with **zero diagonal** (`diag=0`, `:28`) and
  **off-diagonal `off_diag(i) = sqrt(i)`** for `i=1..n−1` (`:31-33`).
- Eigen-decomposes `J` via `F08JEF('I', n, diag, off_diag, J_mat, …)` (`:34`);
  eigenvalues → `points`, first row of eigenvectors → `weights = J_mat(1,1:n)**2`
  (`:38-40`).
- **Normalization:** weights sum to 1 (probability weights for an **expectation**,
  not the √π convention). Because off-diag = √i (not √(i/2)), the nodes are the
  **probabilists' Hermite** roots → these compute `E[f(X)]` for `X ~ N(0,1)`
  directly, i.e. `∑ w_k f(x_k) ≈ ∫ f(x) φ(x) dx`.
- Then shifted/scaled: `points = mu + sigma*points` (`:39`). With the model's
  `mu=0, sigma=1` the nodes are unshifted.

→ **Python equivalent that reproduces IDENTICAL nodes/weights:** use the
probabilists' convention. `numpy.polynomial.hermite_e.hermegauss(n)` returns nodes
`x_k` and weights `w_k` with `∑ w_k = √(2π)`; normalize `w_k ← w_k / √(2π)` to get
the probability weights this routine returns, and the nodes match directly. (Equivalently
`hermgauss` (physicists') with node rescale by √2 and weight `/√π` — but `hermegauss`
is the cleaner match.) For tensor code, precompute once on CPU and store as
`float64` tensors; do **not** recompute per call. Verify against the Fortran by
checking off-diag=√i Jacobi-matrix eigen-decomposition rather than trusting library
conventions blindly.

#### Truncated-normal machinery (used only by `_ub/_lb/_lb_ub`)
- `make_points_and_weights(n, moments_vec, points, weights)` `:109`. Given raw
  moments `m_0..m_{2n}`, builds the Hankel/Gram moment matrix (`matrix_maker :266`),
  Cholesky-factorizes it (`choldc`), forms the Jacobi recurrence coefficients
  (`alpha_vec`, `beta_vec`), assembles symmetric tridiagonal `J`, eigen-decomposes
  via `F08FAF` (`dsyev`), and returns nodes=eigenvalues, weights=first-eigvec². This
  is generic Gauss quadrature from moments. → `scipy` orthogonal-poly-from-moments,
  or replicate the Golub–Welsch-from-moments steps.
- `trunc_normal_moment_ub/lb/lb_ub(mu, sigma, [a,]b, m, moment)` `:155,:188,:221`.
  Recursions for moments of a (one/two-sided) truncated normal (Dhrymes notes),
  using `normal_pdf` and `S15ABF` (normal CDF). → `scipy.stats.truncnorm.moment`
  or port the recursion; needs `norm.cdf` in place of `S15ABF`.
- `normal_pdf(x)` `:260`. Standard-normal pdf `(1/√(2π))exp(−x²/2)`. → closed form.

---

## 3. `mod_markov.f90` — Rouwenhorst / Gospodinov–Lkhagvasuren (⚠ no call sites found)

Module `mod_markov` (`mod_markov.f90:3`). Implements the **MM method of
Gospodinov & Lkhagvasuren (2013)** for VAR(1) discretization. Public:
`calc_rouwen`, `make_markov_2d`, `make_markov_3d` (`:9`). `dp = nag_wp` (`:11`).
**No invocation found** in the compiled model driver files — document for
completeness, port lazily.

- **`calc_rouwen(rho, mu_uncond, sig_uncond, n_R, P_Rouw, z_Rouw)`** `:339`.
  Standard **Rouwenhorst** discretization of a univariate AR(1):
  - Grid `z_Rouw`: `n_R` evenly spaced on `[mu_uncond − step, mu_uncond + step]`
    with `step = sig_uncond·√(n_R−1)` (`:362-367`).
  - Transition built recursively from the 2×2 base `[[p,1−p],[1−p,p]]`,
    `p=(rho+1)/2`, the usual 4-corner padding recursion (`:382-391`), then
    **transposed** and column-normalized (`:393-397`) so **columns sum to 1**
    (column = "from" state). ⚠ Note the column-stochastic orientation (transpose
    of the usual row-stochastic convention).
  → `quantecon.markov.rouwenhorst` (but transpose to match orientation).
- **`make_markov_2d(A, Omg, n_stoch, n_grid, n_tune, nodes, probs)`** `:184` and
  **`make_markov_3d(...)`** `:19`. Discretize a 2-/3-variable VAR(1)
  `y' = A y + ε`, `Var(ε)=Omg`, into a tensor-product Rouwenhorst grid with the
  GL "tuning" of conditional variances. Helpers: `norm_var` `:403` (solves the
  discrete Lyapunov eq by fixed-point iteration for the unconditional variance,
  then normalizes), `calc_cond_var` `:462` (mixture conditional variance, picks
  bracketing grid nodes `n_a,n_b` and mixing weight `p`).
  - Output `nodes`: `(n_grid**n_stoch) × n_stoch`; `probs`: column-stochastic
    `(n_grid**n_stoch) × (n_grid**n_stoch)`.
  → If ever needed: port `calc_rouwen` first, then `norm_var` (Lyapunov fixed
    point → `scipy.linalg.solve_discrete_lyapunov`), then the tuning loop. All
    small (n_grid, n_stoch≤3) — vectorization gains negligible; correctness first.

---

## 4. `mod_smolyak.f90` — anisotropic Smolyak (Judd–Maliar–Maliar–Valero 2014) — **MOST IMPORTANT**

Module `mod_smolyak` (`mod_smolyak.f90:14`), a Fortran port (by Lenel) of the
JMMV (2014) MATLAB code. Public: `Smolyak_Elem_Isotrop`, `Smolyak_Elem_Anisotrop`,
`Smolyak_Grid`, `Smolyak_Polynomial`, `Smolyak_Polynomial2` (`:22`).

### Key conventions (get these exactly right)
- Domain is the hypercube **[−1, 1]^d**. Model states are rescaled into [−1,1]
  before calling (state centering/scaling lives in `mod_param`/`mod_calc`, not here).
- 1-D growth rule **m(i)**: `m(1)=1`, `m(i)=2^(i−1)+1` for `i>1`
  (`Smolyak_Grid:442-446`, `Smolyak_Polynomial:606-611`). So levels give
  1, 3, 5, 9, 17, … nodes. The level index `i` runs `1..mu+1`.
- A per-dimension **subindex** `n ∈ {1,…,m_i_max}` selects a 1-D node / a 1-D
  Chebyshev degree (degree = n−1). A multidimensional element is a length-`d`
  row of such subindices. The same index matrix `smol_elem` is used for both
  **grid points** (via `Smolyak_Grid`) and **basis terms** (via `Smolyak_Polynomial`)
  — points and basis are in one-to-one correspondence, so the basis matrix on the
  grid is square.

### `Smolyak_Elem_Isotrop(d, mu)` `:44-300`
Returns the integer matrix of multidimensional subindices (`n_elem × d`) for the
**isotropic** level-`mu` Smolyak rule.
- **Construction rule:** build the set of level-vectors `i=(i_1,…,i_d)`,
  `i_k ≥ 1`, satisfying the Smolyak rule `d ≤ |i| ≤ d+mu` (`|i|=∑i_k`). Built
  incrementally by `j=0..mu`: at each `j`, generate vectors with `|i| = d+j` by
  adding 1 in each dimension to the previous shell and **deduping rows**
  (`remove_dups_2d`, order-preserving) (`:76-156`).
- For each level-vector ("`one_comb`"), expand each dimension's level `i_k` into
  the **disjoint set of new 1-D indices** A_{i_k} (`:208-239`):
  - `i_k=1` → index `{1}` (the single point 0).
  - `i_k=2` → indices `{2, 3}` (= `2 .. 2^(i_k−1)+1`).
  - `i_k≥3` → indices `{2^(i_k−2)+2 , … , 2^(i_k−1)+1}` (the *new* extrema added
    at this level).
- Take the **Cartesian product** across dimensions of these per-dimension index
  sets (`:241-280`, manual product into `z`), concatenate across all level-vectors
  (`:282-298`). Result: the disjoint-set construction so every multi-index appears
  exactly once. → Python: `itertools.product` per level-vector or vectorized
  `meshgrid`; **preserve row order** to match Fortran (downstream solves rely on
  consistent point/term ordering — but since coeffs are solved from the same matrix,
  any consistent permutation is internally valid; match the Fortran only if you need
  cross-checking against KL outputs). ⚠ Recommend matching order for validation.

### `Smolyak_Elem_Anisotrop(smol_elem_iso, d, vector_mus_dimensions)` `:331-379`
Prunes the isotropic index set to the **anisotropic** one given per-dimension
levels `µ_k = vector_mus_dimensions(k)`.
- `points_dimensions(k) = 1` if `µ_k=0`, else `2^{µ_k}+1` (`:343-354`).
  ⚠ Note exponent is `2^{µ_k}` here, **not** `2^{µ_k−1}+1`: this is the max 1-D
  subindex allowed in dim k (i.e. m(µ_k+1) using the level→index shift, since the
  level cap is µ_k so the index cap is m(µ_k+1)=2^{µ_k}+1).
- Keep a row of `smol_elem_iso` iff **every** dimension's subindex
  `≤ points_dimensions(k)`; rows with any over-cap entry are dropped (`:358-377`).
  → Python: boolean mask `(iso <= caps).all(axis=1)`; `iso[mask]`. Fully vectorized.

### `Smolyak_Grid(d, mu, smol_elem)` `:412-537`
Returns the collocation points `smol_grid` (`n_elem × d`), real, in [−1,1].
- Builds the ordered 1-D node vector `points_1d` by, for `i=1..mu+1`, computing
  the **extrema of Chebyshev polynomials** (Chebyshev–Gauss–**Lobatto** /
  "extrema" nodes):
  `extrem_Cheb_1d(j) = −cos(π·(j−1)/(m_i−1))`, `j=1..m_i` (`:458-460`),
  i.e. nodes on [−1,1] including endpoints, with a leading minus sign so they run
  −1 → +1. m(1)=1 special-cases to `{0}` (`:452-453`). Values within 1e-12 of
  0/±1 are snapped exactly (`:464-472`). Appended across levels and **order-preserving
  deduped** (`remove_dups_1d_real`) so `points_1d` is the nested sequence
  `0, −1, 1, …` indexed by subindex `n` (`:480-498`).
- Each grid point: `smol_grid(jp, jd) = points_1d( smol_elem(jp, jd) )` — gather
  by subindex (`:519-533`). → Python: precompute `points_1d` once; grid is a
  fancy-index gather `points_1d[smol_elem - 1]` (−1 for 0-based). Fully vectorized.
- ⚠ Node type: **Chebyshev extrema (Lobatto), not Gauss (roots).** Endpoints
  ±1 are included. Matches `−cos(π(j−1)/(m−1))`.

### `Smolyak_Polynomial(points, d, mu, smol_elem)` `:573-698` — **the interpolation basis**
Builds the basis matrix `smol_bases` (`numb_pts × numb_terms`), where
`numb_terms = #rows of smol_elem`.
- 1-D Chebyshev basis `T_{n-1}` (first-kind), built by recurrence into
  `phi(point, dim, n)`, `n=1..m_i_max`, `m_i_max = 2^{mu}+1` (`:606-640`):
  - `phi(:,:,1) = 1` (T_0),
  - `phi(:,:,2) = points` (T_1 = x),
  - `phi(:,:,j) = 2·points·phi(:,:,j−1) − phi(:,:,j−2)` (Chebyshev recurrence,
    `:636`).
  So **column subindex `n` ↔ Chebyshev degree `n−1`** in that dimension.
- Each multidimensional term `jt` (a row of `smol_elem`) is the **tensor product**
  `∏_{jd=1}^{d} phi(:, jd, n_{jd})` where `n_{jd}=smol_elem(jt,jd)`; terms with
  `n=1` contribute T_0=1 and are skipped in the product (`:657-696`). Column `jt`
  of the basis matrix = that product evaluated at all points.
- So a basis column maps to the multi-index of Chebyshev **degrees**
  `(n_1−1, …, n_d−1)` read off row `jt` of `smol_elem`.

**Interpolation procedure (this is used everywhere in the solver):**
1. `B = Smolyak_Polynomial(smol_grid, d, mu, smol_elem)` — square `n×n` basis at
   the grid (`smol_polynom` in `mod_param:346`).
2. Solve `B · c = v` for coefficients `c`, where `v` are function values at the
   grid points (`F07ABF`/`dgesvx`, `mod_calc.f90:465` etc.). `B` is **reused/factored
   once per outer iteration** and solved against multiple RHS (policies, bond prices,
   values).
3. Evaluate at new points: `B_new = Smolyak_Polynomial(new_points, d, mu, smol_elem)`,
   then `f̂(new_points) = B_new · c` (`dgemm`/matmul, `mod_calc.f90:499` etc.).

### `Smolyak_Polynomial2(points, d, numb_pts, numb_terms, mu, smol_elem)` `:700-826`
Functionally **identical** basis to `Smolyak_Polynomial`, but with sizes passed
explicitly (fixed-shape arrays, no `allocate`) and `phi` dimensioned `2^mu+1`.
Written for speed inside the OpenMP hot loop (`mod_calc.f90` uses this variant;
`mod_results.f90` uses the allocatable `Smolyak_Polynomial`). **Port as one
function** — the math is the same; the split is just Fortran allocation hygiene.

---

## 5. Vectorization notes (PyTorch/JAX, tensor-native, GPU-ready)

- **`base_lib`:** linspace/var/cov/sort → direct tensor ops. `m_choose_r` → host
  scalars (`math.comb`), used only in index math. `choldc` → `torch.linalg.cholesky`.
  `remove_dups_*` → order-preserving unique (seen-set on CPU at setup time; not in
  hot path). Brent → SciPy or a small batched custom solver if the model needs many
  roots at once (then vectorize over the batch).

- **`mod_normal` (`get_quadrature_points_nb`):** compute the `(n,)` Gauss–Hermite
  nodes/weights **once** at setup (probabilists' convention, `hermegauss` then
  normalize weights by √(2π)); store as `float64` device tensors. Expectations
  become a **batched contraction**: `E[f] = einsum('q,...q->...', w, f_at_nodes)`,
  batching over states/grid points along the leading dims and quadrature nodes on
  the last. No per-call eigen-decomposition.

- **`mod_markov`:** small fixed sizes; port for correctness, not speed. `norm_var`
  Lyapunov fixed point → `scipy.linalg.solve_discrete_lyapunov` (host, setup time).
  Build transition matrices once at setup; the solver just indexes them.

- **`mod_smolyak` (hot path — make this a batched matmul):**
  - `Smolyak_Elem_*` and `Smolyak_Grid`: **setup-time, host**. Build the integer
    index matrix `smol_elem` (shape `(T, d)`) and `points_1d` once; grid =
    `points_1d[smol_elem-1]` gather.
  - **Basis evaluation** = the critical kernel. Vectorize as:
    1. Evaluate 1-D Chebyshev up to max degree at all points via the recurrence,
       producing `phi` of shape `(N_pts, d, m_max)` — recurrence is a loop over
       degree only (≤ ~17 iters), each step a tensor op; **no Python loop over
       points or dims**.
    2. Gather per-term per-dim factors with the precomputed index tensor
       `smol_elem` (shape `(T, d)`): `factors = phi[:, arange(d), smol_elem-1]`
       → shape `(N_pts, T, d)` (advanced indexing, batched).
    3. Term values = `factors.prod(dim=-1)` → basis `B` of shape `(N_pts, T)`.
  - **Coefficient solve:** factor the square grid basis once
    (`torch.linalg.lu_factor`), then `torch.linalg.lu_solve` against the
    stacked RHS (policies/prices/values as columns) — mirrors the Fortran
    `dgesvx` reuse.
  - **Evaluate at new points:** `f̂ = B_new @ c` — a single batched matmul,
    GPU-friendly. This `B_new @ c` is the operation that dominates the solver, so
    it must be a dense `matmul`, never a Python loop over terms or quad nodes.

---

## 6. The single most important thing to get right

**The Smolyak interpolation identity must be bit-faithful: basis = first-kind
Chebyshev `T_{n-1}` (recurrence `T_j = 2x·T_{j-1} − T_{j-2}`, `T_0=1`, `T_1=x`)
evaluated on Chebyshev *extrema/Lobatto* nodes `x = −cos(π(j−1)/(m−1))` over
[−1,1], with the 1-D growth rule m(i)=2^{i−1}+1; the column↔degree map is
`smol_elem[term, dim] − 1`; and interpolation is solve `B·c = v` then evaluate
`B_new·c`.** If the node type (extrema vs roots), the degree-offset (n−1), or the
anisotropic index cap (`2^{µ_k}+1`) is wrong, every policy/value interpolation in
the solver is silently wrong. Validate the Python `Smolyak_Polynomial` against the
Fortran on a fixed `smol_elem` and random `points` before trusting any model output.
