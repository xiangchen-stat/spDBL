# Print and predict from a prepared calibration object

Methods for the result of
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md).
[`print()`](https://rdrr.io/r/base/print.html) shows the parameter
bounds, the design GP and the emulator draws.
[`predict()`](https://rdrr.io/r/stats/predict.html) calls
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md),
which needs no new simulator runs, or with `interval = "prediction"`
returns pointwise predictive distributions.

## Usage

``` r
# S3 method for class 'spDBL_calibration_input'
print(x, ...)

# S3 method for class 'spDBL_calibration_input'
predict(
  object,
  newdata,
  MC = FALSE,
  interval = c("none", "prediction"),
  level = 0.95,
  ...
)
```

## Arguments

- x, object:

  An `spDBL_calibration_input` from
  [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md).

- ...:

  Unused.

- newdata:

  Matrix (or vector, for one input) of mechanistic inputs.

- MC:

  Passed to
  [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md).

- interval:

  `"none"` (the default) returns what
  [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
  returns. `"prediction"` returns the pointwise predictive distribution
  of the field: given each emulator draw it is Gaussian, with the
  design-row uncertainty when
  [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)
  propagates it; the equal-weight mixture over draws is summarised by
  its mean and variance, and the intervals are Gaussian with those
  moments.

- level:

  Coverage of the central prediction intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`predict()`](https://rdrr.io/r/stats/predict.html) returns what
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
returns, or, with `interval = "prediction"`, a list with time-indexed
lists `mean`, `sd`, `lower` and `upper` of input by site matrices, plus
`level`, `distribution` and `n_draws`.

## See also

[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md),
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)

## Examples

``` r
data("dt_calibration")
d <- dt_calibration
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
  covariance = "sigma2R", episode_window = 2, nsam = 3, seed = 1)
cal <- calibration_prepare(em, d$lower, d$upper)
cal
#> spDBL calibration input: 3 emulator draws (marginal sampling), 36 sites, 26 times
#>  parameter lower upper status lengthscale
#>       eta1   2.0   4.0   free       0.828
#>       eta2   0.2   0.4   free       1.100
#>     alpha1   0.0   0.0  fixed          NA
#>     alpha2   0.0   0.2   free       0.503
#>     alpha3   0.0   0.0  fixed          NA
#> Design GP: Matern 5/2 kernel, hyperparameters per time step (26) (lengthscales on the unit scale; medians over time steps)
#> Design GP: leave-one-out RMSE 0.431, nugget 1e-08
#> Design-row GP variance is propagated within each time, independently across times 
band <- predict(cal, d$validation$input[1:2, ], interval = "prediction")
round(band$sd[[10]][, 1:4], 3)
#>       [,1]  [,2]  [,3]  [,4]
#> [1,] 0.183 0.173 0.191 0.181
#> [2,] 0.147 0.140 0.156 0.147
pred <- predict(cal, d$validation$input[1:2, ])
dim(pred[[5]])
#> [1]  2 36
```
