# Forward Filter Backward Sampler (scalar sigma-squared-times-R model)

Runs the complete FFBS algorithm under the model where the right
covariance is \\\sigma^2 R\\ with \\\sigma^2 \sim IG(n_0, d_0)\\.

## Usage

``` r
FFBS_sigma2R(
  Y,
  F_ls,
  G_ls,
  W_ls,
  V_ls,
  m0,
  M0,
  n0,
  D0,
  nT,
  R,
  delta = 1,
  verbose = FALSE
)
```

## Arguments

- Y:

  List of length `nT`. Data matrices.

- F_ls:

  Covariate matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- G_ls:

  State transition matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- W_ls:

  State noise left-covariance matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- V_ls:

  Observation noise left-covariance matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- m0:

  Numeric matrix. Prior state mean.

- M0:

  Numeric matrix. Prior state left-covariance.

- n0:

  Numeric scalar. Prior shape of \\\sigma^2\\.

- D0:

  Numeric scalar. Prior rate of \\\sigma^2\\.

- nT:

  Integer. Number of time steps.

- R:

  Numeric matrix. Fixed spatial correlation matrix.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list with `ff` (output of
[`FF_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FF_sigma2R.md))
and `bs` (output of
[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md)).

## See also

[`FF_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FF_sigma2R.md),
[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md),
[`FFBS_sampling_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_sigma2R.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
R <- 0.5^abs(outer(1:3, 1:3, "-"))
res <- FFBS_sigma2R(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
                    n0 = 3, D0 = 1, nT = 4, R = R)
res$bs$st[, , 1]
#>             [,1]       [,2]       [,3]
#> [1,] -0.65045166 -0.6325265  0.2084102
#> [2,] -0.05196911  1.0530088 -0.5105561
all.equal(res$ff, FF_sigma2R(dat$Y, par$F0, par$G0, par$W0, par$V0,
                             par$m0, par$M0, 3, 1, nT = 4, R = R))
#> [1] TRUE
```
