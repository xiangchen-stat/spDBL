# Draw predictive samples from a precomputed mean and covariance

Draws one sample from the predictive distribution \\y_t \sim N(\mu_t,
\Sigma)\\ at each calibration time step, given precomputed mean and
covariance (e.g., from
[`update_muSigma_eta_one`](https://xiangchen-stat.github.io/spDBL/reference/update_muSigma_eta_one.md)).

## Usage

``` r
sample_y_eta_one(ind_sam, nT_cal, N_sp_train, mu_arry_train, Sigma_eta_train)
```

## Arguments

- ind_sam:

  Integer. Sample index (used only for loop iteration; single sample is
  drawn).

- nT_cal:

  Integer. Number of calibration time steps.

- N_sp_train:

  Integer. Number of training spatial locations.

- mu_arry_train:

  Numeric array of dimension `c(nT_cal, N_sp_train)`. Predictive means.

- Sigma_eta_train:

  Numeric matrix. Predictive covariance (`N_sp_train x N_sp_train`).

## Value

Numeric matrix of dimension `c(nT_cal, N_sp_train)` with one draw per
time step.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## See also

[`update_muSigma_eta_one`](https://xiangchen-stat.github.io/spDBL/reference/update_muSigma_eta_one.md)

## Examples

``` r
set.seed(1)
mu <- matrix(1:6, nrow = 2, ncol = 3)
y <- sample_y_eta_one(ind_sam = 1, nT_cal = 2, N_sp_train = 3,
                      mu_arry_train = mu, Sigma_eta_train = 0.01 * diag(3))
#> Warning: sample_y_eta_one() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
dim(y)
#> [1] 2 3
round(y - mu, 2)
#>       [,1] [,2]  [,3]
#> [1,] -0.06 0.02 -0.08
#> [2,]  0.16 0.03 -0.08
```
