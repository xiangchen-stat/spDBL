# Save a ggplot to a timestamped PNG file

Saves `plot` to a PNG file in `path_fig` whose name is `filename`
followed by the current Unix timestamp, ensuring unique file names
across repeated calls.

## Usage

``` r
quick_save(filename, path_fig, plot)
```

## Arguments

- filename:

  Character. Base name for the output file (no extension).

- path_fig:

  Character. Directory path where the file is saved.

- plot:

  A `ggplot` object to save.

## Value

The path of the saved file, invisibly (the value of
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)).
Called for its side effect.

## See also

[`quick_heat`](https://xiangchen-stat.github.io/spDBL/reference/quick_heat.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggpubr", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE)) {
  data("dt_calibration")
  d <- dt_calibration
  dt <- data.frame(x = d$grid$col, y = d$grid$row, fill = d$field$z[10, ])
  p <- quick_heat(dt, max_y = max(d$field$z))
  path_fig <- tempfile("spdbl_fig_")
  dir.create(path_fig)
  quick_save("heat", path_fig, p)
  list.files(path_fig)
  unlink(path_fig, recursive = TRUE)
}
```
