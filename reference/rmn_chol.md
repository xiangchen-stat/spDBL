# Draw one sample from a matrix-normal distribution (Cholesky parameterisation)

Samples \\\Theta \sim MN(m, U, \Sigma)\\ using the upper-triangular
Cholesky factors `RM` and `RSigma` of \\U\\ and \\\Sigma\\ respectively:
\\\Theta = m + \text{RM}^\top Z \\ \text{RSigma}\\, where \\Z\\ has
i.i.d. standard normal entries.

## Usage

``` r
rmn_chol(m, RM, RSigma)
```

## Arguments

- m:

  Numeric matrix. Mean matrix (\\p \times S\\).

- RM:

  Numeric matrix. Upper-triangular Cholesky factor of the row covariance
  \\U\\ (\\p \times p\\).

- RSigma:

  Numeric matrix. Upper-triangular Cholesky factor of the column
  covariance \\\Sigma\\ (\\S \times S\\).

## Value

A numeric matrix of dimension `c(p, S)`.

## See also

[`rmn_chol_more`](https://xiangchen-stat.github.io/spDBL/reference/rmn_chol_more.md)

## Examples

``` r
set.seed(1)
U <- matrix(c(1, 0.5, 0.5, 1), 2)
Sigma <- diag(c(1, 2, 3))
draw <- rmn_chol(m = matrix(0, 2, 3), RM = chol(U), RSigma = chol(Sigma))
dim(draw)
#> [1] 2 3
draw
#>            [,1]      [,2]       [,3]
#> [1,] -0.6264538 -1.181757  0.5707242
#> [2,] -0.1541871  1.362933 -0.9453405
```
