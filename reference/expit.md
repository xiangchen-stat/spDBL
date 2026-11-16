# Logistic (expit) function

Computes the logistic sigmoid \\1 / (1 + e^{-x})\\.

## Usage

``` r
expit(x)
```

## Arguments

- x:

  Numeric vector or matrix.

## Value

Numeric vector or matrix with values in \\(0, 1)\\.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## Examples

``` r
expit(c(-2, 0, 2))
#> Warning: expit() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
#> [1] 0.1192029 0.5000000 0.8807971
all.equal(expit(1.5), stats::plogis(1.5))
#> [1] TRUE
```
