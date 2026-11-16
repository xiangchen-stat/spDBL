# Fit modular calibration conditional on emulator posterior draws

Targets one conditional posterior per fixed emulator coefficient set and
pools equal numbers of samples. The emulator is not updated by field
observations. Gaussian emulator fields are integrated out. By default
discrepancy states are also integrated out using a Kalman filter.

## Usage

``` r
calibration_fit(
  object,
  z,
  U,
  locations = seq_len(object$S),
  backend = "rstan",
  draws = seq_len(object$M),
  chains = 1L,
  iter,
  warmup = iter%/%2L,
  draws_per_emulator = min(200L, chains * (iter - warmup)),
  n0 = 30,
  d0 = 2,
  discount = 1,
  m0 = rep(0, length(locations)),
  M0 = diag(length(locations)),
  seed = 1234L,
  cores = 1L,
  marginal_discrepancy = TRUE,
  keep_predictions = FALSE,
  keep_discrepancy = TRUE,
  control = if (warmup == 0L) list(adapt_engaged = FALSE, stepsize = 0.001, max_treedepth
    = 12L) else list(adapt_delta = 0.95, max_treedepth = 12L),
  refresh = 100L,
  checkpoint_dir = NULL,
  cache = TRUE,
  cache_dir = NULL,
  compiled_model = NULL,
  init = "mode",
  mode_starts = 20L,
  parallel_draws = 1L
)
```

## Arguments

- object:

  Prepared emulator from
  [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md).

- z:

  Finite time-by-observed-location matrix of field observations.

- U:

  Full-grid symmetric positive-definite discrepancy covariance scale.
  Correlation matrices are supported, but the diagonal need not equal
  one. Discrepancy increments have covariance
  tau2[t](https://rdrr.io/r/base/t.html) \* U; observation errors retain
  covariance tau2[t](https://rdrr.io/r/base/t.html) \* I.

- locations:

  Ordered original-grid indices corresponding to columns of z.

- backend:

  Currently "rstan".

- draws:

  Emulator draw indices to use, without duplicates.

- chains, iter, warmup:

  Number of chains, total iterations per chain, and warmup iterations
  for each emulator-conditional fit. `iter` has no default and must be
  given. For inference use several chains and enough iterations, for
  example `chains = 4, iter = 2000`; `iter = 1` (with no warmup) is a
  quick execution check, not inference. `chains` defaults to one and
  `warmup` to `floor(iter / 2)`.

- draws_per_emulator:

  Equal number of retained samples per conditional. Defaults to the
  smaller of 200 and the available post-warmup samples.

- n0, d0:

  Positive inverse-gamma shape and scale hyperparameters.

- discount:

  Calibration variance discount factor in (0, 1\].

- m0, M0:

  Initial discrepancy mean and covariance scale at observed sites.

- seed:

  Nonnegative integer RNG seed, local to the function call.

- cores:

  Positive number of parallel chains. Defaults to one.

- marginal_discrepancy:

  Whether to integrate out discrepancy states.

- keep_predictions, keep_discrepancy:

  Whether to recover and retain latent fields, replicated observations,
  and discrepancy states.

- control:

  Named list passed to rstan sampling. With zero warmup the defaults
  disable adaptation and use stepsize = 0.001. With warmup they use
  adapt_delta = 0.95. Both use max_treedepth = 12.

- refresh:

  Stan progress reporting interval.

- checkpoint_dir:

  Optional directory for completed conditional fits. Each fit is saved
  as it completes; a fit already saved there with the same data,
  settings, seed and Stan program is reused instead of refitted, so an
  interrupted calibration resumes where it stopped.

- cache, cache_dir:

  Whether to reuse compiled models and optional user cache directory.
  Installed package files are never modified.

- compiled_model:

  Optional compatible compiled Stan model for advanced use.

- init:

  Starting values of the free parameters for every chain. `"mode"` (the
  default) runs L-BFGS
  ([`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html))
  from `mode_starts` random starts in each conditional fit and starts
  every chain, slightly jittered, at the optimum with the highest log
  posterior density, with the precisions of that optimum. Chains started
  far from the main mode can settle in local modes of negligible
  posterior mass, which this avoids; the optima found are recorded in
  `diagnostics[[j]]$mode_search`. `"random"` draws each chain's start
  uniformly within the bounds, reproducibly from `seed`, so that \\\hat
  R\\ can reveal chains that settle in different modes; `"midpoint"`
  starts every chain at the centre of the bounds; or give a numeric
  vector of parameter values in training-column order (named vectors are
  matched by name), strictly inside the bounds of the free parameters.
  Values given for fixed parameters are ignored. Each conditional fit
  starts afresh; except for `"mode"`, precisions start at n0 / d0.

- mode_starts:

  Number of random starts of the mode search for `init = "mode"`.

- parallel_draws:

  Number of conditional fits (one per emulator draw) run at the same
  time, in forked processes. Each forked fit runs its chains one after
  another, so `cores` is then not used and `parallel_draws` processes
  run at once. With a multithreaded BLAS (for example oneMKL or
  OpenBLAS), start R with `OMP_NUM_THREADS=1`: forking after the BLAS
  has used several threads can leave the forked fits hung, and
  `calibration_fit()` warns when it detects that setup. Results do not
  depend on it: every fit draws its random numbers from its own seeds,
  and rstan's draws do not depend on the number of cores. Forking is
  unavailable on Windows, where fits run one after another. Defaults to
  one.

## Value

An ordinary serializable list of class spDBL_calibration containing
samples, per-conditional diagnostics, settings, and prediction context.
No compiled Stan model is stored in the result.

## Details

Requires the suggested package rstan and a C++ toolchain on first
compilation. One iteration per emulator draw is an execution experiment,
not converged posterior inference; its samples stay near `init`. When
the midpoint of the bounds is close to the true value, as in
`dt_calibration`, use `init = "random"` to check that results do not
depend on the start. Missing observations and new time/grid locations
are not supported in this initial interface.

## See also

[`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md),
[`calibration_validate()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_validate.md)

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("rstan", quietly = TRUE)
data("dt_calibration")
d <- dt_calibration
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
                     covariance = "sigma2R", episode_window = 2, nsam = 3,
                     seed = 1)
cal <- calibration_prepare(em, d$lower, d$upper)
fit <- calibration_fit(cal, d$field$z, d$field$U,
                       locations = d$field$locations, draws = 1:2,
                       chains = 1, iter = 200, seed = 1)
summary(fit)
}
```
