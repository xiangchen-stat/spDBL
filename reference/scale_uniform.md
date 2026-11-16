# Scale a vector to the unit interval via a uniform transformation

Maps `x` to \\(x - \text{low}) / (\text{high} - \text{low})\\. Raises an
error if `x` is outside `[low, high]`. Sets `NaN` results to 0.

## Usage

``` r
scale_uniform(x, low, high)
```

## Arguments

- x:

  Numeric vector or matrix. Values to scale.

- low:

  Numeric. Lower bound of the original range.

- high:

  Numeric. Upper bound of the original range.

## Value

Numeric vector or matrix with values in \\\[0, 1\]\\.

## Lifecycle

Deprecated since spDBL 1.2.0 and will be removed in a future release.
The first call in a session warns; set
`options(spDBL.warn_deprecated = FALSE)` to silence it.

## See also

[`scale_back_uniform`](https://xiangchen-stat.github.io/spDBL/reference/scale_back_uniform.md)

## Examples

``` r
x <- c(2, 3.5, 5)
u <- scale_uniform(x, low = 2, high = 5)
u
#> [1] 0.0 0.5 1.0
all.equal(scale_back_uniform(u, low = 2, high = 5), x)
#> [1] TRUE
```
