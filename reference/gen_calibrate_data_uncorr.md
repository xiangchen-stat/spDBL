# Generate synthetic calibration data with uncorrelated discrepancy

Like
[`gen_calibrate_data`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data.md)
but uses independent (uncorrelated) discrepancy increments drawn from a
univariate normal with variance \\\tau^2_t \cdot U\_{11}\\.

## Usage

``` r
gen_calibrate_data_uncorr(ycal_mat, para_gen_cal, U_gen)
```

## Arguments

- ycal_mat:

  Numeric matrix. Computer model output (`t_cal x s_cal`).

- para_gen_cal:

  Named list (same structure as in
  [`gen_calibrate_data`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data.md)).

- U_gen:

  Numeric matrix. Only the `[1,1]` entry is used as the marginal
  variance scale.

## Value

A named list with:

- z:

  Matrix of dimension `c(t_cal, s_cal)` of simulated observations.

- tau2:

  Numeric vector of length `t_cal`.

## See also

[`gen_calibrate_data`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data.md)

## Examples

``` r
set.seed(1)
ycal <- matrix(rnorm(3 * 4), nrow = 3, ncol = 4)
para <- list(b = 0.9, n0 = 20, d0 = 2, u0 = rep(0, 4))
sim <- gen_calibrate_data_uncorr(ycal, para, U_gen = 0.5 * diag(4))
names(sim)
#> [1] "z"    "tau2"
dim(sim$z)
#> [1] 3 4
sim$tau2
#> [1] 0.12001881 0.09979815 0.12033094
```
