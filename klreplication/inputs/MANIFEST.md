# Provenance manifest — KL (2024 AER) quantitative-model inputs

Copied for the big_cy tensorized-Python replication (task Q0, agent big_cy-klrep)
on 2026-06-13.

**Source root:**
`…/Angus Lewis/big_cy/kekre_lenel_2024_aer_replication_package/Safety/Quantitative Model/Code/src/`

The Python replication is intended to run **from just the files under `params/`
and `data/` below** (plus, eventually, a Python port of the generators). The
`generators/` files are copied for reference only — they *produce* the
`param_file_*.csv`; they are not themselves model inputs.

## `params/` — true numeric inputs (read directly by the Fortran)

| File | Source path | Raw or generated | Notes |
|---|---|---|---|
| `param_file_1.csv` … `param_file_9.csv` | `src/params/param_file_N.csv` | **generated** by `create_param_files.m` | One header row (65 param names) + one row of 65 values. One file per calibration. `param_file_1` = benchmark. Read by `mod_param.f90:init_setup` (`./main N`). Copied as the concrete inputs the solver consumes; also to be reproduced by a Python port of the generator. |
| `sample1_mat.txt`, `sample2_mat.txt`, `sample3_mat.txt` | `src/params/sampleK_mat.txt` | **raw** (pre-drawn) | Single-line comma-separated pre-generated standardized shock draws used by the simulation in `mod_results.f90`. Layout: blocks of `n_shocks=4` per period; in the observed sample the first two entries per block are 0. Exact read order to be confirmed from `mod_results.f90`. Three files → three simulation samples. |
| `additional_params.mat` | `src/params/additional_params.mat` | **raw** (small, 229 B) | MATLAB v? binary; a handful of extra scalars consumed by `create_param_files.m`. Read in Python via `scipy.io.loadmat`. |

## `data/` — data series for model-vs-data comparison plots

| File | Source path | Raw or generated | Notes |
|---|---|---|---|
| `fig_10.csv` | `src/data/fig_10.csv` | raw (external) | Data series plotted against the model in Figure 10 (`create_safety_fig.m`). |
| `modelvsdata.csv` | `src/data/modelvsdata.csv` | raw (external) | Data column(s) for model-vs-data figures/tables. |

## `generators/` — reference only (NOT model inputs)

| File | Source path | Role |
|---|---|---|
| `create_param_files.m` | `src/params/create_param_files.m` | Builds all `param_file_*.csv` for the 9 calibrations; writes `output/tmp/n_comp.txt` (number of calibrations). To be ported to Python. |
| `get_moments.m` | `src/params/get_moments.m` | Supplies calibration moments/targets used by `create_param_files.m`. |
| `sim_moments.m` | `src/params/sim_moments.m` | Simulation-moment helper used in calibration. |

## The 65 parameters (header order of `param_file_*.csv`)

```
zeta, bbeta_h, bbeta_f, gma_h, gma_f, ies_h, ies_f, chi, sigma, varsigma_h,
varsigma_f, l_target, delta, aalpha, chiX, inv_share_h, varphi_w, vareps_w,
phi_pi_h, phi_pi_f, phi_y_h, phi_y_f, tayl_ic_h, tayl_ic_f, rho_i, b_lmbd,
bg_yss, sig_z, std_zf, rho_zf, disast_p, disast_shock, p_rho, p_std, omg_mean,
omg_std, omg_shift, omg_rho, corr_omg_dis, tht_trgt_h, theta_h_dev, k_dev_param,
k_grid_adj, w_dev_param, w_grid_adj, ih_grid_mean, if_grid_mean, ih_grid_dev,
if_grid_dev, mu_k, mu_tht, mu_zf, mu_wh, mu_wf, mu_p, mu_ih, mu_if, mu_omg,
n_quad_z, n_quad_zf, n_quad_p, n_quad_omg, foreign_trading, run_bg, run_samp
```

(Maps 1:1 to the `param_input_vec(1..65)` reads in `mod_param.f90:init_setup`.
`mu_*` are the per-dimension Smolyak approximation levels
`vector_mus_dimensions`; `n_quad_*` the per-shock Gauss-Hermite node counts.)

## NOT copied (deliberately)

- The Fortran/MATLAB source itself (read in place; this is a language refactor).
- The `Empirical/` folder (separate task; produces the data moments/IRFs).
- Reference outputs (`Code/output/figures/*.pdf`, `tables/*.tex`) — read in place
  for comparison, not inputs.
