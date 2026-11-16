# Summarize modular calibration samples

Summarize modular calibration samples

## Usage

``` r
# S3 method for class 'spDBL_calibration'
summary(object, ...)
```

## Arguments

- object:

  A calibration result.

- ...:

  Unused.

## Value

A list of parameter summaries, conditional diagnostics, and settings.

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
s <- summary(fit)
s$parameters
}
```
