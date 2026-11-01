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
