# Draw multiple samples from a matrix-normal distribution (Cholesky parameterisation)

Draws `nsam` independent samples from \\MN(m, U, \Sigma)\\ using the
Cholesky factors of \\U\\ and \\\Sigma\\.

## Usage

``` r
rmn_chol_more(nsam, m, RM, RSigma)
```

## Arguments

- nsam:

  Integer. Number of samples to draw.

- m:

  Numeric matrix. Mean matrix (\\p \times S\\).

- RM:

  Numeric matrix. Upper-triangular Cholesky factor of the row covariance
  \\U\\ (\\p \times p\\).

- RSigma:

  Numeric matrix. Upper-triangular Cholesky factor of the column
  covariance \\\Sigma\\ (\\S \times S\\).

## Value

An array of dimension `c(p, S, nsam)`.

## See also

[`rmn_chol`](https://xiangchen-stat.github.io/spDBL/reference/rmn_chol.md)

## Examples

``` r
set.seed(1)
U <- matrix(c(1, 0.5, 0.5, 1), 2)
Sigma <- diag(c(1, 2, 3))
draws <- rmn_chol_more(2000, m = matrix(0, 2, 3), RM = chol(U), RSigma = chol(Sigma))
dim(draws)
#> [1]    2    3 2000
var(draws[1, 3, ])                # close to U[1, 1] * Sigma[3, 3] = 3
#> [1] 2.933414
cor(draws[1, 1, ], draws[2, 1, ]) # close to U[1, 2] = 0.5
#> [1] 0.5231691
```
