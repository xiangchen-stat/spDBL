# Invert a uniform scaling transformation

Maps `x` from the unit interval back to the original range
\\\[\text{low}, \text{high}\]\\: \\x \cdot (\text{high} - \text{low}) +
\text{low}\\.

## Usage

``` r
scale_back_uniform(x, low, high)
```

## Arguments

- x:

  Numeric vector or matrix. Values in \\\[0, 1\]\\.

- low:

  Numeric. Lower bound of the target range.

- high:

  Numeric. Upper bound of the target range.

## Value

Numeric vector or matrix rescaled to \\\[\text{low}, \text{high}\]\\.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## See also

[`scale_uniform`](https://xiangchen-stat.github.io/spDBL/reference/scale_uniform.md)

## Examples

``` r
u <- c(0, 0.25, 1)
x <- scale_back_uniform(u, low = -1, high = 3)
#> Warning: scale_back_uniform() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
x
#> [1] -1  0  3
all.equal(scale_uniform(x, low = -1, high = 3), u)
#> Warning: scale_uniform() is deprecated since spDBL 1.2.0 and will be removed in a future release. It comes from the original research code; new work should use calibration_prepare() and calibration_fit().
#> [1] TRUE
```
