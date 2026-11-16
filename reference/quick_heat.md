# Quick raster heatmap

Creates a `ggplot2` raster heatmap using the `col_bgr` colour palette
with a squished colour scale capped at `max_y`.

## Usage

``` r
quick_heat(dt, max_y, col_bgr = spdbl_default_col_bgr())
```

## Arguments

- dt:

  Data frame with columns `x` (horizontal position), `y` (vertical
  position), and `fill` (value to display).

- max_y:

  Numeric scalar. Upper limit of the colour scale. Values above this are
  squished to the maximum colour.

- col_bgr:

  Character vector of colours for the heatmap gradient.

## Value

A `ggplot` object.

## See also

[`quick_save`](https://xiangchen-stat.github.io/spDBL/reference/quick_save.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggpubr", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE)) {
  data("dt_calibration")
  d <- dt_calibration
  dt <- data.frame(x = d$grid$col, y = d$grid$row, fill = d$field$z[10, ])
  p <- quick_heat(dt, max_y = max(d$field$z))
  class(p)
  p
}
```
