# Prepare a fitted emulator for PDE-free calibration

Fits a Gaussian process (the design GP) to the recorded autoregressive
design matrices, so that design rows can be predicted at untried inputs.
Reuses the fitted emulator's coefficient samples or draws additional
samples from saved filtering results, without refitting the dynamic
emulator.

## Usage

``` r
calibration_prepare(
  emulator,
  lower,
  upper,
  n_draws = NULL,
  seed = 1234,
  sampling = NULL,
  design_uncertainty = TRUE,
  design_kernel = c("matern52", "sqexp"),
  design_groups = c("time", "shared")
)
```

## Arguments

- emulator:

  A sigma2R object from
  [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).

- lower, upper:

  Finite parameter bounds in training-column order. Equal bounds specify
  fixed coordinates. Named bounds are matched to the training parameter
  names.

- n_draws:

  Number of fresh coefficient sample sets. NULL reuses existing samples
  when the selected sampling method agrees with the fitted emulator.

- seed:

  Optional local seed used when drawing fresh samples.

- sampling:

  NULL inherits the fitted emulator's method (normally "marginal").
  "marginal" draws independent smoothed episode coefficients conditional
  on a shared variance; "joint" draws coherent backward paths. Changing
  the method resamples the existing number of draws unless n_draws is
  supplied.

- design_uncertainty:

  If TRUE (the default), the predictive covariance of the AR design GP
  is added to the field covariance used by
  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md)
  and
  [`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md).
  FALSE uses the GP mean only, as spDBL 1.1.1 did; its intervals are too
  narrow.

- design_kernel:

  Kernel of the design GP: `"matern52"` (Matern 5/2, the default) or
  `"sqexp"` (squared exponential). Both are anisotropic, with one
  lengthscale per free parameter on the unit scale.

- design_groups:

  `"time"` (the default) fits separate kernel hyperparameters, by
  maximum likelihood, for the design columns of each time step;
  `"shared"` fits one set for all times.
  `design_kernel = "sqexp", design_groups = "shared"` reproduces the
  design GP of earlier spDBL versions.

## Value

A spDBL_calibration_input object containing the training data,
coefficient samples, offline design predictor, and shared GP factors. At
most one draw's residual weights are cached at a time.

## Details

The AR design GP supplies the design rows at a new input. At each time
step, all design columns share one kernel (with
`design_groups = "time"`; one kernel for all times with `"shared"`), so
their prediction errors at an input \\\eta\\ are jointly Gaussian with
covariance \\C\_{F,t}(\eta) = v_t(\eta)\hat\Sigma\_{F,t}\\:
\\v_t(\eta)\\ is the GP's latent variance factor and
\\\hat\Sigma\_{F,t}\\ is the maximum-likelihood covariance between the
design columns of time \\t\\. With `design_uncertainty = TRUE` this is
propagated through each coefficient draw: at each time, the field
covariance gains \\\Theta^\top C\_{F,t}(\eta) \Theta\\ over the episodes
of that time, including covariance between spatial blocks. Design errors
are treated as independent across times and independent of the
emulator's own prediction error. The emulator's conditional covariance
and coefficient uncertainty are always retained. Calibration support is
limited to the scalar sigma2R model on the original time/grid. Objects
prepared by spDBL 1.1.1 or earlier have no design variance.

## See also

[`calibration_validate()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_validate.md),
[`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md),
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)

## Examples

``` r
data("dt_calibration")
d <- dt_calibration
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
                     covariance = "sigma2R", episode_window = 2, nsam = 3,
                     seed = 1)
cal <- calibration_prepare(em, d$lower, d$upper)
class(cal)
#> [1] "spDBL_calibration_input"      "offline_calibration_emulator"
c(times = cal$T, locations = cal$S, draws = cal$M)
#>     times locations     draws 
#>        26        36         3 
cal$parameter_names
#> [1] "eta1"   "eta2"   "alpha1" "alpha2" "alpha3"
```
