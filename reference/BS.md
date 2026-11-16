# Backward sampler for the Forward Filter Backward Sampler (FFBS)

Performs the backward sampling pass of the FFBS algorithm using the
filtered distributions produced by
[`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md). Iterates
from time \\T\\ back to \\1\\, drawing smoothed state samples at each
step via `BS_1step_cpp`.

## Usage

``` r
BS(res_ff, G_ls, nT, delta = 1, verbose = FALSE)
```

## Arguments

- res_ff:

  List of length `nT` (plus a `prior` element) returned by
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md). Each
  element `res_ff[[t]]` must contain `mt`, `Mt`, `at`, and `At`.

- G_ls:

  Either a single numeric matrix (state transition matrix \\G\\,
  constant over time) or a list of `nT` matrices (time-varying).

- nT:

  Integer. Number of time steps.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

## Value

A named list with:

- st:

  Array of dimension `c(p, q, nT)` containing the smoothed state means.

- St:

  Array of dimension `c(p, p, nT)` containing the smoothed state
  left-covariance matrices.

## See also

[`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md),
[`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res_ff <- FF(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
             par$n0, par$D0, nT = 4)
res_bs <- BS(res_ff, par$G0, nT = 4)
dim(res_bs$st)
#> [1] 2 3 4
# At the last time the smoothed and filtered means agree
all.equal(res_bs$st[, , 4], res_ff$T4$mt)
#> [1] TRUE
```
