# Monte Carlo prediction using FFBS output (MNIW model)

Estimates posterior predictive means at new spatial locations using
Monte Carlo integration over the posterior samples of the state
\\\Theta_t\\. Uses an exponential GP kernel to compute the
cross-covariance between observed and new locations.

## Usage

``` r
FFBS_predict_MC(
  nsam,
  Y,
  res_ffbs,
  F_ls,
  F_new_ls,
  input,
  input_new,
  nT,
  phi_para,
  gp_sigma2 = 1.1,
  gp_tau2 = 1e-04,
  delta = 1
)
```

## Arguments

- nsam:

  Integer. Number of posterior samples to average over.

- Y:

  List of length `nT`. Observed data matrices.

- res_ffbs:

  List. Posterior samples of \\\Theta_t\\, as returned by
  [`FFBS_sampling`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md).

- F_ls:

  Covariate matrix or list for observed locations.

- F_new_ls:

  Covariate matrix or list for new (prediction) locations.

- input:

  Numeric matrix. Coordinates of observed locations (\\N \times d\\).

- input_new:

  Numeric matrix or vector. Coordinates of new locations.

- nT:

  Integer. Number of time steps.

- phi_para:

  Numeric scalar. Range parameter of the exponential kernel.

- gp_sigma2:

  Numeric scalar. Variance parameter of the GP kernel. Defaults to
  `1.1`.

- gp_tau2:

  Numeric scalar. Nugget variance of the GP kernel. Defaults to `1e-4`.

- delta:

  Discount factor. Not implemented: values other than 1 give a warning
  and are ignored. Defaults to `1`.

## Value

A named list of length `nT`. Each element `"T<t>"` is an array of
dimension `c(N_new, q, nsam)` with posterior predictive mean samples at
the new locations.

## See also

[`FFBS_predict_exact`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md),
[`gen_exp_kernel`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md)

## Examples

``` r
set.seed(1)
dat <- gen_ffbs_data(N = 5, S = 3, p = 2, nT = 4)
par <- dat$para
res <- FFBS(dat$Y, par$F0, par$G0, par$W0, par$V0, par$m0, par$M0,
            par$n0, par$D0, nT = 4)
sam <- FFBS_sampling(50, res, par$F0, par$G0, nT = 4, seed = 1)
x_new <- matrix(c(-0.5, 0.5))  # two new inputs
pred <- FFBS_predict_MC(50, dat$Y, sam, F_ls = par$F0,
                        F_new_ls = cbind(1, x_new), input = par$F0[, 2],
                        input_new = x_new, nT = 4, phi_para = 1)
dim(pred$T4)
#> [1]  2  3 50
apply(pred$T4, 1:2, mean)  # Monte Carlo predictive mean at time 4
#>            [,1]       [,2]       [,3]
#> [1,] -0.4608624 -1.5027306  0.2486866
#> [2,] -0.4236159 -0.3317425 -0.1456163
```
