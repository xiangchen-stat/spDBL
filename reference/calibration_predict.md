# Recover fields at the calibrated observation times and locations

Recovers discrepancy trajectories when needed, then draws latent
emulator fields conditional on the observations and fresh replicated
observations. This is not forecasting at new times or new spatial sites.

## Usage

``` r
calibration_predict(object, seed = 1234L)
```

## Arguments

- object:

  Result from
  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md).

- seed:

  Local RNG seed for Gaussian recovery.

## Value

A list containing y, u, and z_rep arrays, each with dimensions posterior
draws by observation times by observed locations.

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
pred <- calibration_predict(fit)
dim(pred$y)
}
```
