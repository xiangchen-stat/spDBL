# Plot a 3-by-3 panel of calibration heatmaps

Produces nine raster heatmaps for calibration data arranged in a 3-by-3
grid, one for each time stamp in `tstamp`. Locations are specified via
`loc_cal` rather than a regular grid.

## Usage

``` r
plot_panel_heatmap_9_cal(
  dat,
  tstamp,
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

  Numeric matrix of dimension `c(nT, n_locations)`. Rows are time steps;
  columns are calibration locations.

- tstamp:

  Integer vector of length 9. Time-step indices (rows of `dat`) to
  display.

- max_y:

  Numeric scalar. Upper limit of the colour scale.

- loc_cal:

  Data frame with columns `row` and `col` giving the spatial coordinates
  of calibration locations.

- Nx:

  Integer. Number of grid points in the x-direction (used for reference
  only).

- Ny:

  Integer. Number of grid points in the y-direction (used for reference
  only).

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

[`plot_panel_heatmap_9_cal_nolab`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9_cal_nolab.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggpubr", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE)) {
  data("dt_calibration")
  d <- dt_calibration
  p <- plot_panel_heatmap_9_cal(d$field$z, tstamp = seq(2, 26, by = 3),
                                max_y = max(d$field$z), loc_cal = d$grid,
                                Nx = d$Nx, Ny = d$Ny, savei = FALSE)
  length(p)
  p[[1]]
}
```
