# Validate offline calibration predictions against held-out simulations

Held-out trajectories are used only for this accuracy assessment. They
are not added to either GP's training data.

## Usage

``` r
calibration_validate(object, input_new, responses, seed = 1234)
```

## Arguments

- object:

  A prepared object from
  [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md).

- input_new:

  Matrix of held-out mechanistic inputs.

- responses:

  Time-indexed list of held-out response matrices.

- seed:

  Not used since spDBL 1.2.0, whose AR designs are deterministic; kept
  for compatibility.

## Value

A list with episode design errors and per-input field RMSE using known
versus GP-predicted autoregressive rows.

## Examples

``` r
data("dt_calibration")
d <- dt_calibration
steps <- 1:6  # the first 6 of 26 time steps, to keep the example fast
em <- emulator_learn(d$training$input, d$training$responses[steps], d$Nx, d$Ny,
                     covariance = "sigma2R", episode_window = 2, nsam = 3,
                     seed = 1)
cal <- calibration_prepare(em, d$lower, d$upper)
held_out <- lapply(d$validation$responses[steps], function(y) y[1:3, , drop = FALSE])
val <- calibration_validate(cal, d$validation$input[1:3, ], held_out)
val$field
#>   input rmse_known_design rmse_offline_design
#> 1     1        0.13372704          0.18827710
#> 2     2        0.06139812          0.06267718
#> 3     3        0.07369800          0.07945046
```
