# Single forward filter step (identity right-covariance)

Performs one step of the forward filter under the model where the right
covariance of the state and observation is the identity matrix (no
inverse-Wishart update for \\\Sigma\\).

## Usage

``` r
FF_1step_R_I(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, delta = 1)
```

## Arguments

- Yt:

  Numeric matrix. Data matrix at time \\t\\ (\\N \times S\\).

- Ft:

  Numeric matrix. Covariate matrix at time \\t\\ (\\N \times p\\).

- Gt:

  Numeric matrix. State transition matrix (\\p \times p\\).

- Wt:

  Numeric matrix. State noise left-covariance (\\p \times p\\).

- Vt:

  Numeric matrix. Observation noise left-covariance (\\N \times N\\).

- mt_1:

  Numeric matrix. Filtered state mean at \\t-1\\ (\\p \times S\\).

- Mt_1:

  Numeric matrix. Filtered state left-covariance at \\t-1\\ (\\p \times
  p\\).

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

## Value

A named list with:

- at:

  One-step-ahead state mean.

- At:

  One-step-ahead state left-covariance.

- mt:

  Filtered state mean.

- Mt:

  Filtered state left-covariance.

- delta:

  Always 1; discounting is not implemented.

## See also

[`FF_I`](https://xiangchen-stat.github.io/spDBL/reference/FF_I.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 1)
par <- dat$para
one <- FF_1step_R_I(dat$Y[[1]], par$F0, par$G0, par$W0, par$V0,
                    par$m0, par$M0)
one$mt
#>            [,1]       [,2]       [,3]
#> [1,] -0.8310570 -0.3986266  0.1869505
#> [2,] -0.1131768  1.0532239 -0.5331634
one$Mt
#>           [,1]      [,2]
#> [1,] 0.5375842 0.1300801
#> [2,] 0.1300801 0.2427817
```
