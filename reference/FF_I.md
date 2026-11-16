# Forward Filter with identity right-covariance

Runs the full forward filtering pass over `nT` time steps under the
model where the right covariance is fixed at the identity (no
inverse-Wishart update). Calls
[`FF_1step_R_I`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_I.md)
at each step.

## Usage

``` r
FF_I(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1, verbose = FALSE)
```

## Arguments

- Y:

  List of length `nT`. Each element is the data matrix at time \\t\\.

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

A named list of length `nT + 1`. Elements `"T1"` through `"T<nT>"`
contain `at`, `At`, `mt`, `Mt`. The element `prior` stores `m0` and
`M0`.

## See also

[`FF_1step_R_I`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_I.md),
[`FFBS_I`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res_ff <- FF_I(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0, nT = 4)
# The first step equals a single call of FF_1step_R_I
one <- FF_1step_R_I(dat$Y[[1]], par$F0, par$G0, par$W0, par$V0,
                    par$m0, par$M0)
all.equal(res_ff$T1, one)
#> [1] TRUE
```
