# spDBL 1.2.0

## Changes that affect existing code or results

* The default autoregressive designs are now causal and deterministic. A lag
  that does not exist yet is a column of ones (a per-site intercept), so the
  AR(1) design at t = 1 is no longer the response itself, and the AR(2)
  designs at t = 1, 2 no longer use random noise or the future value Y_3.
  `emulator_learn()` and `emulator_predict()` results change, and
  predictions no longer draw from the session's random numbers. AR(2)
  emulators no longer need three time steps.
* `calibration_fit()` requires `iter`. For inference use several chains and
  enough iterations (for example `chains = 4, iter = 2000`); `iter = 1` is a
  quick execution check.
* `emulator_learn()` no longer builds the dense S x S `R_field` unless it is
  supplied, so `setup$R_field` is `NULL` by default. Calibration builds the
  same kernel at the observed sites when it needs it, and results are
  unchanged. This removes an O(S^3) step and O(S^2) memory from every fit.
* The nine helpers from the original research code (`update_y_eta()`,
  `update_y_eta_one()`, `update_muSigma_eta_one()`, `sample_y_eta_one()`,
  `cal_Bt_bt()`, `scale_uniform()`, `scale_back_uniform()`, `expit()`,
  `cal_jacobian_logit_uniform()`) are deprecated. They still work and warn
  once per session; `options(spDBL.warn_deprecated = FALSE)` silences them.
* `gen_ran_matrix()` now scales by the largest absolute entry, so its entries
  stay in [-1, 1] as documented.

## Calibration uncertainty

* Calibration intervals were too narrow because the predictive variance of
  the design-row GP was never used. `calibration_prepare()` now propagates it
  by default (`design_uncertainty = TRUE`): at each time the field covariance
  gains the design GP's error through each coefficient draw, including the
  correlation between design columns and between spatial blocks. Both Stan
  programs, the R reference densities, `calibration_fit()` and
  `calibration_predict()` use it. Predicting held-out runs of a 20 x 20 SIR
  model, pointwise 95% intervals covered 54% of outputs before and 94% after.
  `design_uncertainty = FALSE` restores the 1.1.1 behaviour, and objects
  prepared by 1.1.1 keep it. The recorded result in `vignette("calibration")`
  was regenerated. Design-row errors are still treated as independent across
  times (see `vignette("calibration-details")`).

## New features

* Out-of-core smoothing and prediction: `BS_bigdata_R()` runs the backward
  pass over the files written by `FF_bigdata_R()`, and
  `FFBS_predict_bigdata_R()` predicts at new inputs block by block. Both
  match the in-memory `BS()` and `FFBS_predict_exact()`. `FF_bigdata_R()`
  now checks the design file's width to find lag-2 columns instead of
  silently ignoring read errors.
* `print()`, `summary()` and `predict()` methods for fitted emulators. The
  summary reports the in-sample fit of the smoothed coefficients and the
  posterior mean residual variance. Prepared calibration objects gain
  `print()` and `predict()`.
* `calibration_fit(init = )` sets the chains' starting values: `"midpoint"`
  (the default, unchanged), `"random"` (uniform within the bounds,
  reproducible from `seed`), or a vector of parameter values.
* Every exported function now has a runnable example. The calibration
  examples that compile Stan run only interactively.

## Bug fixes

* `prepare_data()` keeps a one-row matrix when only one run goes to the
  training or test set, and stops with a clear message if `prop_train`
  leaves either set empty.
* The discount factor `delta` was accepted but ignored by the filters,
  smoothers and predictors. Values other than 1 now give a warning; results
  are unchanged. (The samplers already required `delta = 1`.)
* `generate.grid.rowsnake()` and `generate.grid.lr()` returned overlapping
  blocks when the grid had a single block, and left rows uncovered when the
  blocks spanned the full width. They now compute each block directly; the
  output is unchanged wherever it was a valid partition, including the calls
  made by `emulator_learn()` and `emulator_predict()`.
* `verbose = TRUE` no longer fails on short series: `FF()` with one time
  step, `FF_sigma2R()` and `FF_I()` with fewer than five, and
  `FFBS_sampling_I()` with fewer than three.
* `FFBS_sampling_I()` no longer reads past the end of a time-varying `G_ls`
  list. It never used `G_ls`, `F_ls` or `delta`; they are kept for
  compatibility. It also handles a single coefficient row or column.
* `gen_pde()` sizes the ODE solver's work array by grid size instead of
  requesting a fixed 1.6 GB. Results are identical.
* `gen_Jt()` no longer adds the nugget to the cross-covariance when a new
  input equals a training input.
* `update_y_eta()`, `update_y_eta_one()` and `update_muSigma_eta_one()` use
  the input-GP prior variance `gp_sigma2 + gp_tau2`, not 1, so the predictive
  variance cannot go negative near a training input.
  `update_muSigma_eta_one()` uses the selected variance draw, like
  `update_y_eta_one()`.
* `cal_jacobian_logit_uniform()` skips fixed coordinates (equal bounds)
  rather than coordinates equal to zero.

## Other changes

* `emulator_predict()` on a prepared calibration object and
  `calibration_validate()` loop over coefficient draws first, which reuses
  each draw's cached residual weights. Results are identical.

## Infrastructure

* Test coverage rises from 31% to 99.9% of R and C++ lines (`covr`). New
  tests check the filters, smoothers, samplers, log predictive densities,
  out-of-core filtering, grid utilities, plotting, legacy helpers, and the
  calibration back end against independent calculations. The RStan back-end
  logic is also tested with a fake rstan, so it runs without compiling Stan.
  Suggests now requires testthat >= 3.3.2.
* Every hand-written source file now carries an MIT license header
  (`SPDX-License-Identifier: MIT`), checked in CI.
* GitHub Actions: `R CMD check` on macOS, Windows, and Ubuntu (R devel,
  release, previous release); test coverage including the RStan integration
  tests; license-header check; pkgdown website.
* New pkgdown website with a grouped function reference and an architecture
  article for contributors.
* New contributing guide, code of conduct, and issue and pull request
  templates.
* README: CRAN installation, a runnable quick example, and documentation
  links.

# spDBL 1.1.1

* Calibration accepts a positive-definite discrepancy covariance scale `U`
  with non-unit diagonal, preserving the original research experiment's
  discrepancy-to-observation variance ratio.

* Default emulator coefficient sampling to smoothed marginals conditional on
  a shared covariance draw. Add sampling = "joint" for coherent backward paths
  under both scalar-variance and inverse-Wishart covariance models.
* Let calibration preparation inherit or explicitly change the fitted sampling
  method, with sampling metadata retained in saved objects and fit summaries.
* Default calibration to one chain, one iteration, zero warmup, and one retained
  sample per emulator draw. Longer iterations and other controls remain options;
  warmup, retained-sample counts, and default controls adapt to the settings.
* Default the installed reproduction script to the marginal oneiter500 profile;
  update the vignettes and recorded example to these defaults.

# spDBL 1.1.0

* Add optional RStan calibration with Gaussian discrepancy marginalization,
  coherent posterior recovery, equal pooling of emulator conditionals, and
  conditional-fit diagnostics.
* Add a sigma2R option and complete metadata to the emulator wrapper.
* Add an offline GP for predicting AR design rows without further PDE solves.
  Its interpolation uncertainty remains excluded from calibration.
* Extend exact FFBS prediction with fixed coefficient trajectories, cached
  factors, and conditional covariance; preserve the default mean interface.
* Correct sigma2R trajectory sampling to condition backward on sampled states.
* Add compact calibration data, two package vignettes, and installed
  reproduction scripts. One-iteration runs are labeled execution experiments.
* Show 500 one-transition calibration samples in the end-to-end vignette,
  with eta histograms and a variance/discrepancy/field comparison panel.
