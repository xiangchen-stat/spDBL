# Print, summarize and predict from a fitted emulator

Methods for the result of
[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).
[`print()`](https://rdrr.io/r/base/print.html) shows the model settings.
[`summary()`](https://rdrr.io/r/base/summary.html) adds the in-sample
fit of the smoothed coefficients and the posterior mean of the residual
variance. [`predict()`](https://rdrr.io/r/stats/predict.html) calls
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md),
or with `interval = "prediction"` returns pointwise predictive
distributions.

## Usage

``` r
# S3 method for class 'spDBL_emulator'
print(x, ...)

# S3 method for class 'spDBL_emulator'
summary(object, ...)

# S3 method for class 'spDBL_emulator'
predict(
  object,
  newdata,
  dt_pde_test = NULL,
  F_new_ls = "default",
  MC = FALSE,
  interval = c("none", "prediction"),
  level = 0.95,
  ...
)
```

## Arguments

- x, object:

  An `spDBL_emulator` from
  [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).

- ...:

  Unused.

- newdata:

  Matrix (or vector, for one run) of new mechanistic inputs.

- dt_pde_test, F_new_ls, MC:

  Passed to
  [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md).

- interval:

  `"none"` (the default) returns what
  [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
  returns. `"prediction"` returns the pointwise predictive distribution
  given the autoregressive rows: for each site and time, a Student-t
  distribution from integrating the coefficients, the residual
  covariance and the input Gaussian process analytically.

- level:

  Coverage of the central prediction intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.spDBL_emulator`, a list with

- settings:

  One-row data frame of model settings.

- rmse_by_time:

  Root mean squared error of the smoothed fit \\F_e \hat\Theta_e\\ to
  the training responses, at each time.

- rmse:

  The same over all times.

- response_sd:

  Standard deviation of the training responses.

- residual_variance:

  Posterior mean of the residual variance, averaged over the sites of a
  block: \\D/(n-C-1)\\ for the inverse-Wishart model and \\D/(n-1)\\
  times the mean diagonal of `R_episode` for the \\\sigma^2 R\\ model,
  from the final filter.

[`predict()`](https://rdrr.io/r/stats/predict.html) returns what
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
returns, or, with `interval = "prediction"`, a list with time-indexed
lists `mean`, `sd`, `lower` and `upper` of new-run by site matrices (the
bounds of central intervals with coverage `level`), plus `level`,
`distribution` and the degrees of freedom `df`.

## See also

[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md),
[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)

## Examples

``` r
data("dt_emulation")
d <- dt_emulation
em <- emulator_learn(d$pde_para_train, d$dt_pde_train, Nx = 10, Ny = 10,
                     nsam = 5, seed = 1)
em
#> spDBL emulator (inverse-Wishart covariance, marginal coefficient sampling)
#>   Training runs: 40 with 5 inputs
#>   Grid: 10 x 10 (100 sites) over 26 times
#>   Episodes: 260 = 26 times x 10 blocks of 10 sites
#>   Design: AR(2), 20 columns
#>   Stored posterior draws: 5
summary(em)
#> spDBL emulator summary
#>   Covariance: iw; sampling: marginal; AR order: 2
#>   40 runs x 100 sites x 26 times; 260 episodes; 5 stored draws
#>   In-sample RMSE of the smoothed fit: 0.05785 (response SD 3.137)
#>   Posterior mean residual variance per site: 0.003574
pred <- predict(em, d$pde_para_test, dt_pde_test = d$dt_pde_test)
dim(pred[[3]])
#> [1]  10 100
band <- predict(em, d$pde_para_test, dt_pde_test = d$dt_pde_test,
                interval = "prediction")
# Share of held-out outputs at time 3 inside the 95% prediction intervals
mean(d$dt_pde_test[[3]] >= band$lower[[3]] & d$dt_pde_test[[3]] <= band$upper[[3]])
#> [1] 0.918
```
