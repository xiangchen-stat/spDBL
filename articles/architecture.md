# Package architecture

This page is for people who want to read, modify, or extend spDBL. It
maps the workflow functions to the code that implements them. For the
statistical model, see
[`vignette("calibration-details")`](https://xiangchen-stat.github.io/spDBL/articles/calibration-details.md)
and [Banerjee et al. (2025)](https://jmlr.org/papers/v26/22-0896.html).

## Layers

| Layer | Files | What lives there |
|----|----|----|
| Workflow API | `R/func_wrapper.R`, `R/calibration_prepare.R`, `R/calibration_fit.R` | [`prepare_data()`](https://xiangchen-stat.github.io/spDBL/reference/prepare_data.md), [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md), [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md), [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md), [`calibration_validate()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_validate.md), [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md), [`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md), and the S3 methods for `spDBL_calibration` |
| Dynamic linear model kernels | `R/func_FF.R`, `R/func_BS.R`, `R/func_FFBS.R` | Forward filter, backward sampler, posterior sampling and prediction for each covariance model |
| Out of core | `R/func_FF.R` ([`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md)), `R/func_bigdata.R` | Filter, smoother and exact predictor that stream CSV blocks from disk |
| Compiled kernels | `src/FF.cpp`, `src/BS.cpp`, `src/rmn.cpp` (RcppEigen) | One forward-filter step, one backward-smoothing step, matrix-normal draws |
| Calibration internals | `R/calibration_recover.R`, `inst/stan/*.stan` | Stan data assembly, the two Stan programs, R reference log densities, recovery of discrepancy and field draws |
| Distributions | `R/func_MNIW.R`, `R/func_MNIW_sampler.R`, `R/func_lppd_ana.R` | MNIW/MNIG updates and samplers, matrix-T density, one-step log predictive densities |
| Data and design construction | `R/func_gen_data.R`, `R/func_generate_grid_Dan.R` | Kernels, autoregressive design matrices, episode partitioning, grid traversal, synthetic data |
| Simulator and legacy helpers | `R/func_calibrate.R` | The SIR reaction–diffusion simulator ([`SIR()`](https://xiangchen-stat.github.io/spDBL/reference/SIR.md), [`gen_pde()`](https://xiangchen-stat.github.io/spDBL/reference/gen_pde.md)) and helpers from the original research code |
| Plotting | `R/func_plot.R` | Heatmap panels and interval summaries |

## Covariance models and naming

The emulator is a matrix-variate dynamic linear model. Its right
(spatial) covariance has three forms, and kernel functions carry a
matching suffix:

| Model | Suffix | Example functions | Available in [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md) |
|----|----|----|----|
| Inverse-Wishart Σ | none | [`FF()`](https://xiangchen-stat.github.io/spDBL/reference/FF.md), [`FFBS()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md), [`FFBS_sampling()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md) | `covariance = "iw"` (default) |
| σ² R, R fixed | `_sigma2R` | [`FF_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_sigma2R.md), [`FFBS_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md), [`FFBS_sampling_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_sigma2R.md) | `covariance = "sigma2R"` (required for calibration) |
| Identity | `_I` | [`FF_I()`](https://xiangchen-stat.github.io/spDBL/reference/FF_I.md), [`FFBS_I()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md), [`FFBS_sampling_I()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_I.md) | no; kernel functions only |

## Emulation: what `emulator_learn()` does

1.  Validates the grid (`Nx`, `Ny`), the blocking (`episode_window` must
    divide `Nx`), and the kernel settings.
2.  Builds two exponential kernels with range
    `3 / (gp_tune * max distance)`: `V` over simulator inputs and
    `R_episode` over the sites of one spatial block. The full-grid
    `R_field` (S × S) is not formed unless supplied; only calibration
    uses it, and `.cal_R_field()` builds it at the observed sites.
3.  Splits each time step’s field into `B = Nx / episode_window` spatial
    blocks
    ([`generate.grid.rowsnake()`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.rowsnake.md)).
    Every (time, block) pair is one *episode* of the dynamic model; this
    is the dynamic transfer-learning partition.
4.  Builds the autoregressive design for every episode with
    [`gen_F_ls_AR1_EP()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1_EP.md)
    or
    [`gen_F_ls_AR2_EP()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2_EP.md),
    unless `F_ls_train` is given.
5.  Runs forward filtering and backward smoothing:
    [`FFBS()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md)
    or
    [`FFBS_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md).
    The inverse-Wishart forward step is in C++ (`FF_1step_cpp`); the σ²
    R step is in R
    ([`FF_1step_R_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_sigma2R.md));
    every backward step calls `BS_1step_cpp`.
6.  Draws `nsam` posterior samples with
    [`FFBS_sampling()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md)
    or
    [`FFBS_sampling_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_sigma2R.md).
    `sampling = "marginal"` draws each episode’s smoothed marginal given
    a shared covariance draw; `sampling = "joint"` draws whole
    coefficient paths.
7.  Returns an object of class `spDBL_emulator` with `para_ffbs`
    (filtering and smoothing moments), `res_ffbs` (posterior draws), and
    `setup` (every setting and design needed for prediction).

[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
calls
[`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md)
(default) or
[`FFBS_predict_MC()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_MC.md)
(`MC = TRUE`), then maps episodes back to the time × space layout with
[`recover_from_EP_exact()`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_exact.md)
or
[`recover_from_EP_MC()`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_MC.md).
Given the output of
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md),
it predicts with GP-interpolated design rows instead, so no new
simulator runs are needed.

## Calibration: from emulator to posterior

1.  [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)
    requires a `"sigma2R"` emulator. It reuses the emulator’s draws or
    resamples them (`sample_calibration_emulator()`), then
    `make_offline_calibration_input()` fits the design-row GP
    (`fit_offline_ar2()`; by default a Matern 5/2 kernel with one length
    scale per input and separate hyperparameters for each time step) and
    caches Cholesky factors. Result: class `spDBL_calibration_input`.
    `predict_calibration()` turns one coefficient draw into the field
    mean and covariance at an input. With `design_uncertainty = TRUE`
    (the default) the covariance includes the design GP’s error, as a
    rank-N factor per time (`.cal_design_factors()`).
2.  [`calibration_validate()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_validate.md)
    compares these PDE-free predictions with held-out simulations.
3.  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md)
    validates the settings, compiles or loads the Stan program
    (`.cal_stan_model()`), and calls `fit_calibration_rstan()`. That
    function runs one conditional Stan fit per emulator draw and pools
    an equal number of samples from each fit. Result: class
    `spDBL_calibration`, with
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`print()`](https://rdrr.io/r/base/print.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
4.  [`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md)
    recovers discrepancy paths (`recover_calibration_discrepancy()`) and
    field draws (`recover_calibration_predictions()`).

There are two Stan programs in `inst/stan/`:

- `calibration_marginal.stan`, the default
  (`marginal_discrepancy = TRUE`), integrates the discrepancy states out
  with a Kalman filter.
- `calibration_offline.stan` samples the discrepancy innovations
  explicitly.

Both integrate out the emulator field analytically. R reference
implementations of the same log densities (`calibration_log_density()`,
`calibration_marginal_log_density()`) let the tests check the Stan code.

## Large grids out of core

For grids too large to hold every episode in memory, three functions
stream CSV blocks from disk (`R/func_FF.R`, `R/func_bigdata.R`):

1.  [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md)
    filters forward and writes each time step’s filtered moments, all
    blocks side by side, to one set of CSV files.
2.  [`BS_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/BS_bigdata_R.md)
    reads those files backwards, one time step at a time, and writes the
    smoothed moments.
3.  [`FFBS_predict_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_bigdata_R.md)
    streams each block’s responses, designs and smoothed means, and
    writes one prediction file per time step.

Memory holds one time step of moments and one data block at a time. The
tests check all three against the in-memory
[`FF()`](https://xiangchen-stat.github.io/spDBL/reference/FF.md),
[`BS()`](https://xiangchen-stat.github.io/spDBL/reference/BS.md) and
[`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md).
Prediction needs every training run in each block (`bnrow = fnrow`), as
[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md)
arranges them.

## Reproducibility mechanisms

- Every function that takes a `seed` sets it locally and restores the
  caller’s random-number state on exit (`.cal_seed()`).
- Compiled Stan models are cached in memory and on disk
  (`tools::R_user_dir("spDBL", "cache")`, or `cache_dir`). The cache key
  combines the Stan file’s MD5 hash and the R, rstan, StanHeaders and
  platform versions, so editing a Stan file triggers a recompile.
- `calibration_fit(checkpoint_dir = ...)` saves each conditional fit to
  disk as it finishes. Nothing reloads these files automatically; a
  rerun starts again from the first emulator draw.
- [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md)
  records `execution_only = TRUE` when there is no warmup or fewer than
  two post-warmup iterations, as with `iter = 1`. `iter` has no default,
  so every analysis states its sampler settings.

## Extending spDBL

- **A new covariance model** needs a one-step forward filter, a sampler
  for the covariance, an `FFBS_*()`/`FFBS_sampling_*()` pair, and a new
  `covariance` option in
  [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).
  The `_sigma2R` functions are the closest template.
- **A new input kernel**: the exponential kernel and its range heuristic
  are defined inside
  [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md);
  [`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md)
  and the calibration code read the stored `phi_para`, `gp_sigma2`, and
  `gp_tau2`.
- **A new sampling backend**: `calibration_fit(backend = )` accepts only
  `"rstan"` today. The backend-specific code is `.cal_stan_model()` and
  `fit_calibration_rstan()`.

## Tests

The tests are in `tests/testthat/`. The calibration tests compare the
Gaussian conditioning, the Stan log densities, and the pooled sampler
with independent R calculations. Tests that compile Stan are skipped
unless `SPDBL_RUN_STAN_TESTS=true`; the coverage workflow on GitHub
Actions turns them on. See the [contributing
guide](https://github.com/xiangchen-stat/spDBL/blob/main/.github/CONTRIBUTING.md)
for how to run tests and coverage locally.

## Known limitations

- The out-of-core path covers filtering, smoothing and exact (mean)
  prediction. Posterior sampling and calibration work in memory.
- The `delta` (discount) argument of the one-step filters is accepted
  but not used.
- Calibration supports the σ² R emulator only, with a fixed discrepancy
  covariance scale `U`, complete observations, and the training grid and
  times.
- The design-row error is treated as independent across observation
  times, although it is in fact persistent; see
  [`vignette("calibration-details")`](https://xiangchen-stat.github.io/spDBL/articles/calibration-details.md).
