# Update the posterior mean and covariance of the discrepancy field

Computes the posterior mean vector \\B_t b_t\\ and covariance matrix
\\B_t\\ of the discrepancy \\u_t\\ given the current predictive mean
\\\mu_t\\, observations \\z_t - u_t\\, and variance \\\tau^2_t\\.

## Usage

``` r
cal_Bt_bt(mut, invSigma, tau2t, zt_ut, I_Stilde)
```

## Arguments

- mut:

  Numeric vector. Predictive mean at calibration locations.

- invSigma:

  Numeric matrix. Inverse of the predictive spatial covariance.

- tau2t:

  Numeric scalar. Current variance \\\tau^2_t\\.

- zt_ut:

  Numeric vector. Observation residuals \\z_t - u_t\\.

- I_Stilde:

  Numeric matrix. Identity-like matrix for the discrepancy prior
  precision (`N_sp_train x N_sp_train`).

## Value

A named list with:

- Bt:

  Numeric matrix. Posterior covariance of the discrepancy.

- Btbt:

  Numeric vector. Posterior mean of the discrepancy.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## Examples

``` r
set.seed(1)
Sigma <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.3, 0.1, 0.3, 1), 3)
mu <- rnorm(3)
zu <- rnorm(3)
post <- cal_Bt_bt(mu, invSigma = solve(Sigma), tau2t = 0.4, zt_ut = zu, I_Stilde = diag(3))
#> Warning: cal_Bt_bt() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
K <- Sigma %*% solve(Sigma + 0.4 * diag(3))  # Kalman gain
all.equal(post$Bt, Sigma - K %*% Sigma)
#> [1] TRUE
all.equal(as.vector(post$Btbt), as.vector(mu + K %*% (zu - mu)))
#> [1] TRUE
```
