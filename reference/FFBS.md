# Forward Filter Backward Sampler (MNIW model)

Runs the complete FFBS algorithm under the Matrix Normal Inverse Wishart
(MNIW) model: first applies the forward filter
([`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)) and
then the backward sampler
([`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md)).

## Usage

``` r
FFBS(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1, verbose = FALSE)
```

## Arguments

- Y:

  List of length `nT`. Each element is the \\N \times q\\ data matrix at
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

  Numeric matrix. Prior state mean (\\p \times q\\).

- M0:

  Numeric matrix. Prior state left-covariance (\\p \times p\\).

- n0:

  Numeric scalar. Prior degrees of freedom of the inverse-Wishart.

- D0:

  Numeric matrix. Prior scale matrix of the inverse-Wishart.

- nT:

  Integer. Number of time steps.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list with:

- ff:

  Output of
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md):
  filtered distributions for each time step.

- bs:

  Output of
  [`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md):
  smoothed state means (`st`) and left-covariances (`St`).

## See also

[`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md),
[`BS`](https://xiangchen-stat.github.io/spDBL/reference/BS.md),
[`FFBS_sampling`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
            par$n0, par$D0, nT = 4)
res$bs$st[, , 1]  # smoothed state mean at time 1
#>             [,1]       [,2]       [,3]
#> [1,] -0.65045166 -0.6325265  0.2084102
#> [2,] -0.05196911  1.0530088 -0.5105561
# The forward part equals a separate call of FF
all.equal(res$ff, FF(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0,
                     par$M0, par$n0, par$D0, nT = 4))
#> [1] TRUE
```
