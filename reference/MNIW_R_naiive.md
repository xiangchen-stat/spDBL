# Naive MNIW posterior update

Computes the posterior parameters of the Matrix Normal Inverse Wishart
(MNIW) model using explicit matrix transposes. Equivalent to
[`MNIW_R`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R.md)
but without `crossprod` optimisations. Intended for reference and
testing.

## Usage

``` r
MNIW_R_naiive(X, Y, newX = NULL, newY = NULL, H, v, S, C, Vb)
```

## Arguments

- X:

  Numeric matrix. Covariate (design) matrix of dimension \\n \times p\\.

- Y:

  Numeric matrix. Response matrix of dimension \\n \times q\\.

- newX:

  Unused legacy parameter. Defaults to `NULL`.

- newY:

  Unused legacy parameter. Defaults to `NULL`.

- H:

  Numeric matrix. Row covariance matrix of the error term (\\n \times
  n\\).

- v:

  Numeric scalar. Prior degrees of freedom of the inverse-Wishart on
  \\\Sigma\\.

- S:

  Numeric matrix. Prior scale matrix of the inverse-Wishart on
  \\\Sigma\\.

- C:

  Numeric matrix. Prior mean matrix of \\B \| \Sigma\\ (\\p \times q\\).

- Vb:

  Numeric matrix. Prior row covariance matrix of \\B \| \Sigma\\ (\\p
  \times p\\).

## Value

A named list with posterior parameters:

- vnew:

  Posterior degrees of freedom.

- Snew:

  Posterior scale matrix of the inverse-Wishart.

- Cnew:

  Posterior mean matrix of \\B \| \Sigma\\.

- Vbnew:

  Posterior row covariance matrix of \\B \| \Sigma\\.

## Details

The model is \\Y = X B + E\\, where \\E \sim MN(0, H, \Sigma)\\,
\\\Sigma \sim IW(v, S)\\, and \\B \| \Sigma \sim MN(C, V_b, \Sigma)\\.

## See also

[`MNIW_R`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R.md)

## Examples

``` r
set.seed(1)
n <- 30
X <- cbind(1, rnorm(n))
Y <- X %*% matrix(c(1, -1, 0.5, 2), 2) + matrix(rnorm(n * 2, sd = 0.1), n, 2)
prior <- list(H = diag(n), v = 4, S = diag(2), C = matrix(0, 2, 2), Vb = diag(10, 2))
post <- do.call(MNIW_R_naiive, c(list(X = X, Y = Y), prior))
round(post$Cnew, 2)
#>       [,1] [,2]
#> [1,]  1.01 0.51
#> [2,] -0.99 1.99
# Same result as the crossprod version
all.equal(post, do.call(MNIW_R, c(list(X = X, Y = Y), prior)))
#> [1] TRUE
```
