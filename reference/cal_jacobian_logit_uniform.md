# Log absolute Jacobian of the logit-uniform transformation

Computes the log (or raw) absolute Jacobian of the transformation from
the constrained parameter `eta` (bounded in
`[eta_limit_low, eta_limit_high]`) to the unconstrained logit scale.
Fixed coordinates (equal lower and upper bounds) are excluded from the
product. Versions up to 1.1.1 excluded zero entries of `eta` instead,
which wrongly dropped a free coordinate equal to zero.

## Usage

``` r
cal_jacobian_logit_uniform(eta, eta_limit_low, eta_limit_high, log = FALSE)
```

## Arguments

- eta:

  Numeric vector. Constrained parameter values.

- eta_limit_low:

  Numeric vector. Lower bounds for each element of `eta`.

- eta_limit_high:

  Numeric vector. Upper bounds for each element of `eta`.

- log:

  Logical. If `TRUE`, returns the log Jacobian. Defaults to `FALSE`.

## Value

Numeric scalar (log absolute Jacobian or absolute Jacobian).

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## Examples

``` r
eta <- c(0.3, 2.5)
low <- c(0, 1)
high <- c(1, 4)
cal_jacobian_logit_uniform(eta, low, high)
#> Warning: cal_jacobian_logit_uniform() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
#> [1] 6.349206
cal_jacobian_logit_uniform(eta, low, high, log = TRUE)
#> [1] 1.84833
# A fixed coordinate (low == high) does not contribute
all.equal(cal_jacobian_logit_uniform(c(0.3, 2), c(0, 2), c(1, 2)), 1 / (0.3 * 0.7))
#> [1] TRUE
```
