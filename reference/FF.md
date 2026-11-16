# Forward Filter for the MNIW dynamic linear model

Runs the full forward filtering pass over `nT` time steps under the
Matrix Normal Inverse Wishart (MNIW) model. At each step, calls
`FF_1step_cpp` to update the state mean and covariance as well as the
inverse-Wishart parameters for the right-covariance matrix \\\Sigma\\.

## Usage

``` r
FF(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1, verbose = FALSE)
```

## Arguments

- Y:

  List of length `nT`. Each element `Y[[t]]` is the \\N \times q\\ data
  matrix at time \\t\\.

- F_ls:

  Either a single \\N \times p\\ covariate matrix (constant over time)
  or a list of `nT` such matrices (time-varying).

- G_ls:

  Either a single \\p \times p\\ state transition matrix or a list of
  `nT` such matrices.

- W_ls:

  Either a single \\p \times p\\ state noise left-covariance matrix or a
  list of `nT` such matrices.

- V_ls:

  Either a single \\N \times N\\ observation noise left-covariance
  matrix or a list of `nT` such matrices.

- m0:

  Numeric matrix. Prior mean of the state \\\beta_0 \| D_0\\ (\\p \times
  q\\).

- M0:

  Numeric matrix. Prior left-covariance of \\\beta_0 \| D_0\\ (\\p
  \times p\\).

- n0:

  Numeric scalar. Prior degrees of freedom of the inverse-Wishart on
  \\\Sigma \| D_0\\.

- D0:

  Numeric matrix. Prior scale matrix of the inverse-Wishart on \\\Sigma
  \| D_0\\ (\\q \times q\\).

- nT:

  Integer. Number of time steps.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list of length `nT + 1`. Elements `"T1"` through `"T<nT>"` each
contain a list with filtered parameters `nt`, `Dt`, `at`, `At`, `mt`,
`Mt`. The additional element `prior` stores `m0`, `M0`, `n0`, `D0`.

## See also

[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md),
[`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res_ff <- FF(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
             par$n0, par$D0, nT = 4)
names(res_ff)
#> [1] "T1"    "T2"    "T3"    "T4"    "prior"
res_ff$T4$mt  # filtered state mean at the last time (p x S)
#>             [,1]       [,2]       [,3]
#> [1,] -0.31333924 -0.9620115  0.1916426
#> [2,]  0.04449267  0.9590536 -0.5153927
```
