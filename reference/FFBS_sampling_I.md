# Draw posterior samples from FFBS output (identity right-covariance)

Draws `nsam` posterior samples of the state \\\Theta_t\\ from the
smoothed distributions produced by
[`FFBS_I`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md),
assuming an identity right-covariance matrix.

## Usage

``` r
FFBS_sampling_I(nsam, para_ffbs, F_ls, G_ls, nT, delta = 1, verbose = FALSE)
```

## Arguments

- nsam:

  Integer. Number of posterior samples.

- para_ffbs:

  List. Output of
  [`FFBS_I`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md).

- F_ls, G_ls, delta:

  Not used: each draw comes from the smoothed marginal at its own time.
  Kept for compatibility with the other samplers.

- nT:

  Integer. Number of time steps.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list of length `nT`. Each element `"T<t>"` is an array of
dimension `c(p, q, nsam)` with posterior samples of the state at time
\\t\\.

## See also

[`FFBS_I`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS_I(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
              nT = 4)
sam <- FFBS_sampling_I(500, res, par$F0, par$G0, nT = 4)
dim(sam$T2)
#> [1]   2   3 500
round(apply(sam$T2, 1:2, mean) - res$bs$st[, , 2], 2)
#>      [,1]  [,2]  [,3]
#> [1,] 0.03 -0.06 -0.03
#> [2,] 0.02 -0.03 -0.01
```
