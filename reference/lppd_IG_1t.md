# One-step log posterior predictive density (scalar sigma-squared-times-R model)

Computes the log posterior predictive density of \\Y_t\\ at a single
time step under the scalar-sigma model (\\\sigma^2 \sim IG\\), using
[`dMTig`](https://xiangchen-stat.github.io/spDBL/reference/dMTig.md).

## Usage

``` r
lppd_IG_1t(Yt, Ft, Vt, st, St, nt, Dt, R)
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

- nt:

  Numeric scalar. Filtered shape parameter at time \\t\\.

- Dt:

  Numeric scalar. Filtered rate parameter at time \\t\\.

- R:

  Numeric matrix. Fixed right correlation matrix (\\S \times S\\).

## Value

Numeric scalar. Log posterior predictive density.

## See also

[`dMTig`](https://xiangchen-stat.github.io/spDBL/reference/dMTig.md),
[`lppd_IW_1t`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IW_1t.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
R <- 0.5^abs(outer(1:3, 1:3, "-"))
res <- FFBS_sigma2R(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
                    n0 = 3, D0 = 1, nT = 4, R = R)
lppd_IG_1t(dat$Y[[2]], par$F0, par$V0, st = res$bs$st[, , 2],
           St = res$bs$St[, , 2], nt = res$ff$T2$nt, Dt = res$ff$T2$Dt,
           R = R)
#> [1] -2.155358
```
