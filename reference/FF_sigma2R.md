# Forward Filter for the scalar-sigma-squared-times-R model

Runs the full forward filtering pass over `nT` time steps under the
model where the right covariance of \\Y\\ is \\\sigma^2 R\\ with
\\\sigma^2 \sim IG(n_0, d_0)\\.

## Usage

``` r
FF_sigma2R(
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

  List of length `nT`. Each element is the \\N \times S\\ data matrix at
  the corresponding time step.

- F_ls:

  Covariate matrix or list of matrices (see
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

  Numeric matrix. Prior state mean (\\p \times S\\).

- M0:

  Numeric matrix. Prior state left-covariance (\\p \times p\\).

- n0:

  Numeric scalar. Prior shape of \\\sigma^2\\.

- D0:

  Numeric scalar. Prior rate of \\\sigma^2\\.

- nT:

  Integer. Number of time steps.

- R:

  Numeric matrix. Fixed spatial correlation matrix (\\S \times S\\). Its
  Cholesky inverse is used internally.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list of length `nT + 1` (same structure as
[`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)) with an
additional scalar `Dt` in place of a matrix.

## See also

[`FF_1step_R_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_sigma2R.md),
[`FFBS_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
R <- 0.5^abs(outer(1:3, 1:3, "-"))
res_ff <- FF_sigma2R(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
                     n0 = 3, D0 = 1, nT = 4, R = R)
c(nt = res_ff$T4$nt, Dt = res_ff$T4$Dt)  # IG posterior of sigma^2
#>        nt        Dt 
#> 33.000000  6.012184 
```
