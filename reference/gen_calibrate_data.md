# Generate synthetic calibration data with correlated discrepancy

Simulates observed calibration data \\z_t\\ by adding a spatially
correlated random discrepancy \\u_t\\ and an independent noise term to
the computer model output \\y_t\\. The variance \\\tau^2_t\\ decays over
time via the discount factor \\b\\.

## Usage

``` r
gen_calibrate_data(ycal_mat, para_gen_cal, U_gen)
```

## Arguments

- ycal_mat:

  Numeric matrix of dimension `c(t_cal, s_cal)`. Computer model output
  at calibration locations and time steps.

- para_gen_cal:

  Named list with elements:

  b

  :   Numeric scalar. Discount factor in \\(0, 1)\\.

  n0

  :   Numeric scalar. Prior shape of \\\tau^2_0\\.

  d0

  :   Numeric scalar. Prior scale of \\\tau^2_0\\.

  u0

  :   Numeric vector of length `s_cal`. Initial discrepancy.

- U_gen:

  Numeric matrix. Spatial covariance matrix for the discrepancy
  (`s_cal x s_cal`).

## Value

A named list with:

- u:

  Matrix of dimension `c(t_cal, s_cal)` of discrepancy realisations.

- z:

  Matrix of dimension `c(t_cal, s_cal)` of simulated observations.

- tau2:

  Numeric vector of length `t_cal` with sampled variances.

## See also

[`gen_prior_u_tau2`](https://xiangchen-stat.github.io/spDBL/reference/gen_prior_u_tau2.md)

## Examples

``` r
set.seed(1)
ycal <- matrix(rnorm(3 * 4), nrow = 3, ncol = 4)
para <- list(b = 0.9, n0 = 20, d0 = 2, u0 = rep(0, 4))
sim <- gen_calibrate_data(ycal, para, U_gen = 0.5 * diag(4))
names(sim)
#> [1] "u"    "z"    "tau2"
dim(sim$z)
#> [1] 3 4
sim$tau2
#> [1] 0.12001881 0.09979815 0.12033094
```
