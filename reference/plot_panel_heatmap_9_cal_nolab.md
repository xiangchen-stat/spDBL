# Plot a 3-by-3 panel of calibration heatmaps without axis labels

Like
[`plot_panel_heatmap_9_cal`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9_cal.md)
but suppresses x- and y-axis labels and supports a custom lower
colour-scale limit (`min_y`).

## Usage

``` r
plot_panel_heatmap_9_cal_nolab(
  dat,
  tstamp,
  min_y = 0,
  max_y,
  loc_cal,
  Nx,
  Ny,
  filename = "plot_panel",
  savei = TRUE,
  col_bgr = spdbl_default_col_bgr(),
  path_fig = "."
)
```

## Arguments

- dat:

  Numeric matrix of dimension `c(nT, n_locations)`.

- tstamp:

  Integer vector of length 9. Time-step indices to display.

- min_y:

  Numeric scalar. Lower limit of the colour scale. Defaults to `0`.

- max_y:

  Numeric scalar. Upper limit of the colour scale.

- loc_cal:

  Data frame with columns `row` and `col`.

- Nx:

  Integer. Number of grid points in the x-direction.

- Ny:

  Integer. Number of grid points in the y-direction.

- filename:

  Character. Base name for the saved PNG. Defaults to `"plot_panel"`.

- savei:

  Logical. Whether to save the panel to disk. Defaults to `TRUE`.

- col_bgr:

  Character vector of colours for the heatmap gradient.

- path_fig:

  Character. Directory path where the file is saved when `savei = TRUE`.
  Defaults to the current working directory.

## Value

A list of nine `ggplot` objects.

## See also

[`plot_panel_heatmap_9_cal`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9_cal.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggpubr", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE)) {
  data("dt_calibration")
  d <- dt_calibration
  u <- d$truth$u
  p <- plot_panel_heatmap_9_cal_nolab(u, tstamp = seq(2, 26, by = 3),
                                      min_y = min(u), max_y = max(u),
                                      loc_cal = d$grid, Nx = d$Nx, Ny = d$Ny,
                                      savei = FALSE)
  length(p)
  p[[1]]
}
```
