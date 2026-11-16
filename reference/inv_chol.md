# Invert a matrix via its Cholesky factorisation

Computes \\X^{-1}\\ using `chol` followed by `chol2inv` for numerical
stability.

## Usage

``` r
inv_chol(X)
```

## Arguments

- X:

  Numeric symmetric positive definite matrix.

## Value

The inverse of `X`.

## Examples

``` r
A <- matrix(c(4, 1, 1, 3), 2)
inv_chol(A)
#>             [,1]        [,2]
#> [1,]  0.27272727 -0.09090909
#> [2,] -0.09090909  0.36363636
all.equal(inv_chol(A), solve(A))
#> [1] TRUE
```
