# Draw posterior samples from FFBS output (MNIW model)

Given the filtered and smoothed distributions from
[`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md),
draws `nsam` posterior sample sets of the state matrices \\\Theta_1,
\ldots, \Theta_T\\ and the covariance \\\Sigma\\ using the MNIW sampler.

## Usage

``` r
FFBS_sampling(
  nsam,
  para_ffbs,
  F_ls,
  G_ls,
  nT,
  delta = 1,
  verbose = FALSE,
  sampling = c("marginal", "joint"),
  seed = NULL,
  keep_draws = TRUE
)
```

## Arguments

- nsam:

  Integer. Number of posterior samples to draw.

- para_ffbs:

  List. Output of
  [`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md),
  containing elements `ff` and `bs`.

- F_ls:

  Covariate matrix or list of matrices (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- G_ls:

  State transition matrix or list (see
  [`FF`](https://xiangchen-stat.github.io/spDBL/reference/FF.md)).

- nT:

  Integer. Number of time steps.

- delta:

  Discount factor. Must be 1: discounting is not implemented. Defaults
  to `1`.

- verbose:

  Logical. If `TRUE`, emits progress messages. Defaults to `FALSE`.

- sampling:

  "marginal" (default) samples the smoothed marginal at each episode
  conditional on the same covariance draw. "joint" samples backward
  conditionals to obtain coherent trajectories.

- seed:

  Optional local RNG seed; NULL follows the session RNG.

- keep_draws:

  If TRUE (default), store every episode's coefficient draws. If FALSE
  (marginal sampling only), store the random-number state at the start
  of each episode's draws instead; spDBL's prediction and calibration
  functions then regenerate identical draws on demand, which cuts the
  memory of large emulators. Indexing the returned draws by episode
  (`[[e]]`) then fails.

## Value

A named list of length `nT + 1`. Elements `"T1"` through `"T<nT>"` are
arrays of dimension `c(p, q, nsam)` containing posterior samples of
\\\Theta_t\\. The element `Sigma` is an array of dimension
`c(q, q, nsam)` containing posterior samples of \\\Sigma\\. Attributes
sampling_method and joint_trajectories identify the chosen law.

## See also

[`FFBS`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md),
[`MNIW_sampler`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_sampler.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
            par$n0, par$D0, nT = 4)
sam <- FFBS_sampling(500, res, par$F0, par$G0, nT = 4, seed = 1)
dim(sam$T2)
#> [1]   2   3 500
dim(sam$Sigma)
#> [1]   3   3 500
# Sample mean versus smoothed mean at time 2
round(apply(sam$T2, 1:2, mean) - res$bs$st[, , 2], 2)
#>      [,1] [,2]  [,3]
#> [1,]    0 0.01 -0.01
#> [2,]    0 0.00  0.00
```
