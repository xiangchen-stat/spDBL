# One-step log posterior predictive density (identity right-covariance model)

Computes the log posterior predictive density of \\Y_t\\ at a single
time step under the model with an identity right-covariance matrix,
using a matrix-normal density.

## Usage

``` r
lppd_id_1t(Yt, Ft, Vt, st, St)
```

## Arguments

- Yt:

  Numeric matrix. Observed data at time \\t\\ (\\N \times S\\).

- Ft:

  Numeric matrix. Covariate matrix at time \\t\\ (\\N \times p\\).

- Vt:

  Numeric matrix. Observation noise left-covariance (\\N \times N\\).

- st:

  Numeric matrix. Smoothed state mean at time \\t\\ (\\p \times S\\).

- St:

  Numeric matrix. Smoothed state left-covariance at time \\t\\ (\\p
  \times p\\).

## Value

Numeric scalar. Log posterior predictive density.

## See also

[`lppd_IW_1t`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IW_1t.md),
[`lppd_IG_1t`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IG_1t.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS_I(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
              nT = 4)
lppd_id_1t(dat$Y[[2]], par$F0, par$V0, st = res$bs$st[, , 2],
           St = res$bs$St[, , 2])
#> [1] -9.465089
```
