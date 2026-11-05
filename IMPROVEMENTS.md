# Improvements and extensions in spDBL

This document lists what spDBL adds to the method of Banerjee, Chen,
Frankenburg and Zhou (2025), *Dynamic Bayesian Learning for Spatiotemporal
Mechanistic Models*, Journal of Machine Learning Research 26(146):1–43, and to
the research code written for that paper
(<https://github.com/xiangchen-stat/Bayesian_Modeling_Mechanistic_Systems>).
It accompanies the submission of spDBL 1.2.0 to the JMLR Machine Learning Open
Source Software track. The paper does not present spDBL, and the package has
not been published elsewhere. Changes by version are in [NEWS.md](NEWS.md).

## Starting point

* **The paper** sets out the method:
  * a matrix-variate dynamic linear model emulator with exact forward
    filtering and backward sampling (FFBS), under inverse-Wishart and
    σ²R covariance models;
  * Gaussian-process interpolation over simulator inputs;
  * episode-wise processing of spatial blocks;
  * analytic log predictive densities;
  * Bayesian calibration with a dynamic spatial discrepancy, fitted by
    Metropolis-within-Gibbs.
* **The research code** is a set of R and C++ scripts for the paper's
  experiments, without a package structure or tests. Its calibration runs the
  simulator for every proposed parameter value. In
  `dev/calibrate_EP_sigma2R.R`, each Metropolis proposal calls
  `update_muSigma_eta_one()`, which solves the PDE with `gen_pde()`. The paper
  reports 59 minutes for one calibration of a 6 × 6 SIR model (30,000 draws).
* **spDBL 1.0.2**, the first CRAN release (August 2026), packaged the
  emulator: `emulator_learn()`, `emulator_predict()`, the filters and
  samplers, out-of-core forward filtering, and the paper's SIR simulator. It
  had no calibration interface.

## Extensions of the method

1. **Calibration without simulator runs.** `calibration_prepare()` fits a
   Gaussian process over the inputs that predicts the autoregressive design
   rows at untried parameter values, so the sampler never calls the simulator
   and its cost no longer depends on the simulator's.
   * Since 1.2.0 the GP uses a Matérn 5/2 kernel with separate
     hyperparameters for each time step.
   * Its predictive variance enters the field covariance
     (`design_uncertainty = TRUE`, the default).
   * `calibration_validate()` checks it on held-out runs.
2. **Marginalised calibration with NUTS.** `calibration_fit()` conditions on
   each emulator posterior draw and then:
   * integrates the emulated field out in closed form and the dynamic
     discrepancy by a Kalman filter;
   * samples only (η, τ²) with Stan's NUTS;
   * recovers the fields and the discrepancy afterwards by exact Gaussian
     backward sampling (`calibration_predict()`);
   * pools the conditional posteriors over emulator draws with equal weights,
     reporting R̂, effective sample sizes and divergences for each.

   Two Stan programs, one with the discrepancy marginalised and one with it
   explicit, have R reference implementations of the same log densities.
3. **More general discrepancy:** a positive-definite discrepancy covariance
   scale `U` with a non-unit diagonal (1.1.1).
4. **Emulator coefficient sampling** either jointly along backward paths or
   from smoothed marginals (`sampling = "joint"` or `"marginal"`), recorded in
   saved objects (1.1.1).
5. **Exact prediction** with fixed coefficient trajectories, cached factors
   and the conditional covariance (`FFBS_predict_exact()`, 1.1.0).
6. **Causal default autoregressive designs** (1.2.0). A lag that does not
   exist yet is an intercept column, so early-time designs no longer use the
   response itself, random noise or future values.

## Software

* **Interface.** Seven workflow functions run the analysis:
  * `prepare_data()`, `emulator_learn()`, `emulator_predict()`;
  * `calibration_prepare()`, `calibration_validate()`, `calibration_fit()`,
    `calibration_predict()`.

  They return classed, serialisable objects (`spDBL_emulator`,
  `spDBL_calibration_input`, `spDBL_calibration`) with `print()`,
  `summary()`, `predict()` and `plot()` methods. The filters, smoothers,
  samplers and predictors remain exported as building blocks for new
  emulators.
* **Large grids.**
  * Emulators no longer build the dense S × S field kernel unless it is
    supplied, which removes an O(S³) step and O(S²) memory from every fit.
  * Out-of-core backward smoothing (`BS_bigdata_R()`) and prediction
    (`FFBS_predict_bigdata_R()`) join forward filtering (`FF_bigdata_R()`),
    streaming CSV blocks. They match the in-memory versions.
* **Reproducibility.**
  * Functions that take a seed restore the caller's random-number state.
  * Compiled Stan models are cached on disk, keyed by the program and the R,
    rstan and platform versions.
  * Long calibrations can checkpoint each conditional fit.
  * Bundled datasets record their provenance.
  * An installed script reproduces the calibration example with recorded
    settings.
* **Model comparison.** Analytical and draw-based WAIC for multivariate
  regression with inverse-Wishart covariance (1.1.2).
* **Fixes found by the new tests** (1.2.0), among them:
  * the discount factor `delta` was accepted but ignored;
  * grid partitions overlapped or left rows uncovered in edge cases;
  * the paper-era calibration helpers could give a negative predictive
    variance near a training input.

## Quality assurance and documentation

* **Tests.**
  * Coverage is 99.9% of R and C++ lines (`covr`), up from 31% in 1.1.1.
  * Numerical components are checked against independent calculations: the
    filter against exact Gaussian conditioning, the Stan log densities
    against R code, and the samplers by their moments.
  * The Stan back end is also tested without compiling Stan.
* **Continuous integration.** GitHub Actions run:
  * `R CMD check` on Ubuntu (R devel, release and the previous release),
    macOS and Windows;
  * test coverage with the Stan tests;
  * a license-header check;
  * the pkgdown build.
* **Documentation.**
  * Every exported function has a runnable example.
  * Vignettes cover emulation, calibration, the calibration model in detail,
    and model comparison.
  * The website (<https://xiangchen-stat.github.io/spDBL/>) adds a grouped
    function reference and an architecture article for contributors.
  * A contributing guide, code of conduct, and issue and pull-request
    templates are in `.github/`.
* **Licensing.** MIT; every hand-written source file carries an SPDX license
  header.

## Benchmarks

Scripts, Slurm job files and raw results are in [`benchmarks/`](benchmarks/).
The paper-scale results are in `benchmarks/results/*/paper/`.

* **Emulator scaling (E1).** Training is closed form; with 10⁴ sites a fit
  takes 1.5 to 15.5 minutes, depending on the episode window.
* **Emulator accuracy (E2).** Predictions at untried inputs, made without
  simulator runs, were compared with RobustGaSP (`ppgasp`) and a
  principal-component GP on 50 held-out runs.
  * spDBL's 95% intervals covered 95–96% of held-out values on every grid.
  * Its RMSE was the lowest on the 20 × 20 grid and within 4–6% of the best
    otherwise.
* **Calibration (E3).** In 30 replicate calibrations against
  RobustCalibration:
  * accuracy matched (relative RMSE 0.103 against 0.102–0.105);
  * spDBL was the most accurate for the hardest parameter;
  * its interval coverage was the highest (82% against 76–80%);
  * no simulator runs were needed.
