# Log density of the matrix-T distribution with inverse-gamma right covariance

Evaluates the log (or raw) density of the matrix-T distribution that
arises when the right covariance is \\\sigma^2 R\\ and \\\sigma^2 \sim
IG(\nu, d)\\.

## Usage

``` r
dMTig(Y, m, M, nu, d, R, log = TRUE)
```

## Arguments

- Y:

  Numeric matrix. Observation matrix (\\p \times S\\).

- m:

  Numeric matrix. Mean matrix (\\p \times S\\).

- M:

  Numeric matrix. Left covariance matrix (\\p \times p\\, positive
  definite).

- nu:

  Numeric scalar. Degrees of freedom (shape parameter of the
  inverse-gamma on \\\sigma^2\\).

- d:

  Numeric scalar. Rate parameter of the inverse-gamma on \\\sigma^2\\).

- R:

  Numeric matrix. Right correlation matrix (\\S \times S\\, positive
  definite).

- log:

  Logical. If `TRUE` (default), returns the log density.

## Value

Numeric scalar. Log density (or density if `log = FALSE`).

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(6), 2, 3)
R <- 0.5^abs(outer(1:3, 1:3, "-"))
dMTig(Y, m = matrix(0, 2, 3), M = diag(2), nu = 3, d = 2, R = R)
#> [1] -9.534796
dMTig(Y, m = matrix(0, 2, 3), M = diag(2), nu = 3, d = 2, R = R, log = FALSE)
#> [1] 7.229208e-05
```
