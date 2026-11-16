# Sample from the Matrix Normal Inverse Wishart (MNIW) distribution

Draws `nsam` posterior samples of the regression coefficient matrix
\\B\\ and the covariance matrix \\\Sigma\\ under the MNIW model: \$\$Y =
X B + E, \quad E \sim MN(0, H, \Sigma)\$\$ \$\$\Sigma \sim IW(v, S),
\quad B \| \Sigma \sim MN(C, V_b, \Sigma).\$\$

## Usage

``` r
MNIW_sampler(nsam, X, v, S, C, Vb)
```

## Arguments

- nsam:

  Integer. Number of samples to draw.

- X:

  Numeric matrix. Covariate (design) matrix (\\n \times p\\).

- v:

  Numeric scalar. Degrees of freedom of the inverse-Wishart prior on
  \\\Sigma\\.

- S:

  Numeric matrix. Scale matrix of the inverse-Wishart prior (\\q \times
  q\\, positive definite).

- C:

  Numeric matrix. Prior mean of \\B \| \Sigma\\ (\\p \times q\\).

- Vb:

  Numeric matrix. Prior row covariance of \\B \| \Sigma\\ (\\p \times
  p\\, positive definite).

## Value

A named list with:

- sigma:

  Array of dimension `c(q, q, nsam)` with samples of \\\Sigma\\.

- B:

  Array of dimension `c(p, q, nsam)` with samples of \\B\\.

## See also

[`MNIG_sampler`](https://xiangchen-stat.github.io/spDBL/reference/MNIG_sampler.md),
[`FFBS_sampling`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md)

## Examples

``` r
set.seed(1)
C <- matrix(c(1, -1), nrow = 1)
draws <- MNIW_sampler(nsam = 2000, X = matrix(1, 5, 1), v = 10, S = diag(2),
                      C = C, Vb = matrix(1))
dim(draws$B)
#> [1]    1    2 2000
# E[Sigma] = S / (v - q - 1) = diag(2) / 7
round(apply(draws$sigma, 1:2, mean), 2)
#>      [,1] [,2]
#> [1,] 0.14 0.00
#> [2,] 0.00 0.14
round(rowMeans(draws$B[1, , ]), 2)  # close to C
#> [1]  1 -1
```
