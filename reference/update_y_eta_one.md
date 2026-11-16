# Update the likelihood of observations given PDE parameters (single sample)

Like
[`update_y_eta`](https://xiangchen-stat.github.io/spDBL/reference/update_y_eta.md)
but uses only one FFBS sample (indexed by `ind_sam`) instead of
averaging over all samples, making it faster for MCMC proposals.

## Usage

``` r
update_y_eta_one(
  eta,
  ind_sam,
  Nx,
  Ny,
  N_people,
  nT_ori,
  nT,
  nT_cal,
  N_sp,
  N_sp_train,
  AR_choice,
  bncol,
  res_ffbs,
  n_b_cal,
  ind_block_cal,
  t_VinvY_FTheta_ls,
  Vinv,
  R_cal_train,
  input,
  input_cal,
  phi_para,
  gp_sigma2,
  gp_tau2,
  ind_train_cal
)
```

## Arguments

- eta:

  Numeric vector. PDE parameters (length 5).

- ind_sam:

  Integer. Index of the FFBS sample to use.

- Nx:

  Integer. Grid cells in x-direction.

- Ny:

  Integer. Grid cells in y-direction.

- N_people:

  Numeric. Total population size.

- nT_ori:

  Integer. Number of original time steps.

- nT:

  Integer. Number of episode time steps.

- nT_cal:

  Integer. Number of calibration time steps.

- N_sp:

  Integer. Total number of spatial locations.

- N_sp_train:

  Integer. Number of training spatial locations.

- AR_choice:

  Integer. AR order (1 or 2).

- bncol:

  Integer. Number of columns per block.

- res_ffbs:

  List. FFBS posterior samples.

- n_b_cal:

  Integer. Number of calibration blocks.

- ind_block_cal:

  Data frame. Block indices for calibration.

- t_VinvY_FTheta_ls:

  List. Precomputed \\V^{-1}(Y - F\Theta)\\ terms.

- Vinv:

  Numeric matrix. Inverse of the observation covariance.

- R_cal_train:

  Numeric matrix. Correlation matrix for training locations.

- input:

  Numeric matrix. Observed location coordinates.

- input_cal:

  Numeric matrix. Calibration location coordinates.

- phi_para:

  Numeric scalar. GP range parameter.

- gp_sigma2:

  Numeric scalar. GP marginal variance.

- gp_tau2:

  Numeric scalar. GP nugget variance.

- ind_train_cal:

  Integer vector. Indices of training locations.

## Value

A named list with `y`, `mu`, and `Sigma`.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## See also

[`update_y_eta`](https://xiangchen-stat.github.io/spDBL/reference/update_y_eta.md)

## Examples

``` r
set.seed(1)
input <- matrix(runif(10), nrow = 5, ncol = 2)
res_ffbs <- lapply(1:12, function(e) array(rnorm(8), c(2, 2, 2)))
res_ffbs$Sigma <- c(0.2, 0.3)
args <- list(eta = c(0.3, 0.1, 0.01, 0.01, 0.01), Nx = 2, Ny = 3, N_people = 100,
  nT_ori = 4, nT = 12, nT_cal = 3, N_sp = 6, N_sp_train = 4, AR_choice = 1,
  bncol = 2, res_ffbs = res_ffbs, n_b_cal = 3,
  ind_block_cal = generate_grid(1, 6, 1, 2),
  t_VinvY_FTheta_ls = lapply(1:3, function(t) array(rnorm(60), c(6, 5, 2))),
  Vinv = solve(gen_exp_kernel(input, phi = 2, sigma2 = 1.1, tau2 = 1e-4)),
  R_cal_train = diag(4), input = input, input_cal = c(0.4, 0.6), phi_para = 2,
  gp_sigma2 = 1.1, gp_tau2 = 1e-4, ind_train_cal = c(1, 3, 4, 6))
out <- do.call(update_y_eta_one, c(args, ind_sam = 2))
#> Warning: update_y_eta_one() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
dim(out$y)  # nT_cal x N_sp_train
#> [1] 3 4
dim(out$Sigma)
#> [1] 4 4
```
