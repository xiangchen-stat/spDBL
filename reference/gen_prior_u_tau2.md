# Sample prior discrepancy trajectory and variance sequence

Simulates a prior trajectory \\u_1, \ldots, u_T\\ of the spatial
discrepancy field and the associated inverse-gamma variance sequence
\\\tau^2_1, \ldots, \tau^2_T\\ from the dynamic calibration model.

## Usage

``` r
gen_prior_u_tau2(n0, d0, b, m0, M0, U, nT_cal)
```

## Arguments

- n0:

  Numeric scalar. Prior shape of \\\tau^2_0\\.

- d0:

  Numeric scalar. Prior scale of \\\tau^2_0\\.

- b:

  Numeric scalar. Discount factor in \\(0, 1)\\.

- m0:

  Numeric vector. Prior mean of the initial discrepancy \\u_0\\.

- M0:

  Numeric matrix. Prior covariance of \\u_0\\.

- U:

  Numeric matrix. Spatial covariance matrix for the discrepancy
  increments.

- nT_cal:

  Integer. Number of calibration time steps.

## Value

A named list with:

- u0:

  Numeric vector. Sampled initial discrepancy.

- u:

  Matrix of dimension `c(nT_cal, s_cal)` of discrepancy realisations.

- tau2_0:

  Numeric scalar. Sampled \\\tau^2_0\\.

- tau2:

  Numeric vector of length `nT_cal`.

## See also

[`gen_calibrate_data`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data.md)

## Examples

``` r
set.seed(1)
prior <- gen_prior_u_tau2(n0 = 20, d0 = 2, b = 0.95, m0 = rep(0, 4),
                          M0 = diag(4), U = 0.5 * diag(4), nT_cal = 3)
names(prior)
#> [1] "u0"     "u"      "tau2_0" "tau2"  
dim(prior$u)
#> [1] 3 4
prior$tau2
#> [1] 0.12906170 0.09390695 0.10341340
```
