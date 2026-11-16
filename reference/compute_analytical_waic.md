# Analytical WAIC for an MNIW multivariate regression posterior

Computes the log pointwise predictive density and the posterior-variance
form of WAIC without posterior simulation. Each row of the response
matrix is one multivariate observation.

## Usage

``` r
compute_analytical_waic(Y, X, Mn, mn, Psi_n, nu_n)
```

## Arguments

- Y:

  Finite numeric response matrix with dimensions n by q.

- X:

  Finite numeric design matrix with dimensions n by p. Include an
  intercept column explicitly when required.

- Mn:

  Symmetric positive-definite p by p posterior row covariance of the
  regression coefficients, conditional on the response covariance.

- mn:

  Finite numeric p by q posterior coefficient mean. Supply the actual
  mean, not its precision-weighted parameterization.

- Psi_n:

  Symmetric positive-definite q by q inverse-Wishart posterior scale
  matrix.

- nu_n:

  Finite scalar inverse-Wishart posterior degrees of freedom, strictly
  greater than q - 1.

## Value

A list containing the original scalar summaries `WAIC`, `lppd`, and
`p_waic`; their difference `elpd_waic`; and an n by 4 `pointwise` matrix
with columns `lppd`, `p_waic`, `elpd_waic`, and `WAIC`.

## Details

The likelihood and posterior are \$\$Y_i \mid B, \Sigma \sim N_q(x_i B,
\Sigma),\$\$ independently across rows, and \$\$B \mid \Sigma, D \sim
MN(m_n, M_n, \Sigma), \quad \Sigma \mid D \sim IW(\nu_n, \Psi_n).\$\$
The inverse-Wishart density is proportional to
\\\|\Sigma\|^{-(\nu_n+q+1)/2}
\exp\\-\mathrm{tr}(\Psi_n\Sigma^{-1})/2\\\\. In particular, `mn` is the
coefficient mean itself. When using
[`MNIW_R()`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R.md)
with `H = diag(n)`, pass `Vbnew`, `Cnew`, `Snew`, and `vnew` as `Mn`,
`mn`, `Psi_n`, and `nu_n`, respectively.

Let \\h_i=x_i M_n x_i^T\\ and \\d_i=(y_i-x_i m_n)\Psi_n^{-1}(y_i-x_i
m_n)^T\\. The row predictive density is multivariate Student t with
\\\nu_n-q+1\\ degrees of freedom, mean \\x_i m_n\\, and scale
\\(1+h_i)\Psi_n/(\nu_n-q+1)\\. The WAIC penalty for row i is \$\$p_i =
\tfrac14\left\\ \sum\_{j=1}^q\psi_1((\nu_n-j+1)/2) + 2q h_i^2 +
4h_i\nu_n d_i + 2\nu_n d_i^2 - 4d_i\right\\.\$\$ Here \\\psi_1\\ is the
trigamma function. WAIC is \\-2\sum_i(\mathrm{lppd}\_i-p_i)\\.

Use the observations that fitted the posterior to compute training-data
WAIC. Scoring other observations evaluates the supplied posterior at
those observations; it does not refit the model. These functions assume
independent likelihood rows and do not directly compute WAIC for
correlated spatial GP rows or a full state-space emulator. Compare
models using the same response data and the same definition of an
observation. Lower WAIC indicates higher estimated predictive accuracy
among the compared models.

## References

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
Journal of Machine Learning Research, 11, 3571–3594.

## See also

[`compute_sample_waic()`](https://xiangchen-stat.github.io/spDBL/reference/compute_sample_waic.md),
[`MNIW_R()`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R.md),
[`lppd_IW_1t()`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IW_1t.md)

## Author

Sudipto Banerjee (original derivation and development implementation).

## Examples

``` r
X <- cbind(1, c(-1, 0, 1, 2))
Y <- cbind(c(0.1, 0.4, 1.2, 1.7), c(1.2, 0.8, 0.5, 0.1))
post <- MNIW_R(X, Y, H = diag(nrow(Y)), v = 4, S = diag(2),
               C = matrix(0, 2, 2), Vb = diag(2))
score <- compute_analytical_waic(Y, X, post$Vbnew, post$Cnew,
                                 post$Snew, post$vnew)
score$WAIC
#> [1] 10.1313
score$pointwise
#>            lppd    p_waic  elpd_waic     WAIC
#> [1,] -0.9815339 0.6702952 -1.6518291 3.303658
#> [2,] -0.6549797 0.2292402 -0.8842199 1.768440
#> [3,] -0.7079442 0.2629221 -0.9708664 1.941733
#> [4,] -0.9120960 0.6466405 -1.5587365 3.117473
```
