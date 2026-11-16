# Forward Filter Backward Sampler (identity right-covariance)

Runs the complete FFBS algorithm with the right covariance fixed at the
identity matrix (no inverse-Wishart update for \\\Sigma\\).

## Usage

``` r
FFBS_I(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1, verbose = FALSE)
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

- nT:

  Integer. Number of time steps.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list with `ff` (output of
[`FF_I`](https://xiangchen-stat.github.io/spDBL/reference/FF_I.md)) and
`bs` (output of
[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md)).

## See also

[`FF_I`](https://xiangchen-stat.github.io/spDBL/reference/FF_I.md),
[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md),
[`FFBS_sampling_I`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_I.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS_I(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
              nT = 4)
res$bs$st[, , 1]
#>             [,1]       [,2]       [,3]
#> [1,] -0.65045166 -0.6325265  0.2084102
#> [2,] -0.05196911  1.0530088 -0.5105561
all.equal(res$ff, FF_I(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0,
                       par$M0, nT = 4))
#> [1] TRUE
```
