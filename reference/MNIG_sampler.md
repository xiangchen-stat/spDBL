# Sample from the Matrix Normal Inverse Gamma (MNIG) distribution

Draws `nsam` posterior samples of the regression coefficient matrix
\\B\\ and the scalar variance \\\sigma^2\\ under the MNIG model where
the right covariance of \\B \| \sigma^2\\ is \\\sigma^2 R\\:
\$\$\sigma^2 \sim IG(v, S), \quad B \| \sigma^2 \sim MN(C, V_b, \sigma^2
R).\$\$

## Usage

``` r
MNIG_sampler(nsam, X, v, S, C, Vb, R)
```

## Arguments

- nsam:

  Integer. Number of samples to draw.

- X:

  Numeric matrix. Covariate (design) matrix (\\n \times p\\).

- v:

  Numeric scalar. Shape parameter of the inverse-gamma prior on
  \\\sigma^2\\.

- S:

  Numeric scalar. Rate parameter of the inverse-gamma prior on
  \\\sigma^2\\.

- C:

  Numeric matrix. Prior mean of \\B \| \sigma^2\\ (\\p \times q\\).

- Vb:

  Numeric matrix. Prior row covariance of \\B \| \sigma^2\\ (\\p \times
  p\\, positive definite).

- R:

  Numeric matrix. Fixed right correlation matrix (\\q \times q\\,
  positive definite).

## Value

A named list with:

- sigma2:

  Numeric vector of length `nsam` with samples of \\\sigma^2\\.

- B:

  Array of dimension `c(p, q, nsam)` with samples of \\B\\.

## See also

[`MNIW_sampler`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_sampler.md),
[`FFBS_sampling_sigma2R`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_sigma2R.md)

## Examples

``` r
set.seed(1)
draws <- MNIG_sampler(nsam = 2000, X = matrix(1, 5, 1), v = 5, S = 2,
                      C = matrix(c(1, -1), nrow = 1), Vb = matrix(1), R = diag(2))
dim(draws$B)
#> [1]    1    2 2000
mean(draws$sigma2)  # close to S / (v - 1) = 0.5
#> [1] 0.5044806
round(rowMeans(draws$B[1, , ]), 2)  # close to C
#> [1]  1.00 -0.99
```
