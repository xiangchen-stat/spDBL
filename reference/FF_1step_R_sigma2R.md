# Single forward filter step (scalar right-covariance, sigma-squared times R)

Performs one step of the forward filter under the model where the right
covariance of \\Y\\ is \\\sigma^2 R\\ with \\\sigma^2 \sim IG(n, d)\\.
The inverse of \\R\\ (`Rinv`) must be precomputed and passed in.

## Usage

``` r
FF_1step_R_sigma2R(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, Rinv, delta = 1)
```

## Arguments

- Yt:

  Numeric matrix. Data matrix at time \\t\\ (\\N \times S\\).

- Ft:

  Numeric matrix. Covariate matrix at time \\t\\ (\\N \times p\\).

- Gt:

  Numeric matrix. State transition matrix at time \\t\\ (\\p \times
  p\\).

- Wt:

  Numeric matrix. State noise left-covariance at time \\t\\ (\\p \times
  p\\).

- Vt:

  Numeric matrix. Observation noise left-covariance at time \\t\\ (\\N
  \times N\\).

- mt_1:

  Numeric matrix. Filtered state mean at \\t-1\\ (\\p \times S\\).

- Mt_1:

  Numeric matrix. Filtered state left-covariance at \\t-1\\ (\\p \times
  p\\).

- nt_1:

  Numeric scalar. Shape parameter of \\\sigma^2 \| D\_{t-1}\\.

- Dt_1:

  Numeric scalar. Rate parameter of \\\sigma^2 \| D\_{t-1}\\.

- Rinv:

  Numeric matrix. Precomputed inverse of the spatial correlation matrix
  \\R\\ (\\S \times S\\).

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

## Value

A named list with updated filtering parameters:

- nt:

  Updated shape parameter.

- Dt:

  Updated rate parameter.

- at:

  One-step-ahead state mean (\\p \times S\\).

- At:

  One-step-ahead state left-covariance (\\p \times p\\).

- mt:

  Filtered state mean (\\p \times S\\).

- Mt:

  Filtered state left-covariance (\\p \times p\\).

- delta:

  Always 1; discounting is not implemented.

## See also

[`FF_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FF_sigma2R.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 1)
par <- dat$para
R <- 0.5^abs(outer(1:3, 1:3, "-"))
one <- FF_1step_R_sigma2R(dat$Y[[1]], par$F0, par$G0, par$W0, par$V0,
                          par$m0, par$M0, nt_1 = 3, Dt_1 = 1,
                          Rinv = solve(R))
c(nt = one$nt, Dt = one$Dt)
#>        nt        Dt 
#> 10.500000  1.523951 
one$mt
#>            [,1]       [,2]       [,3]
#> [1,] -0.8310570 -0.3986266  0.1869505
#> [2,] -0.1131768  1.0532239 -0.5331634
```
