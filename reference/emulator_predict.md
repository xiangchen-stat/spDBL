# Predict simulator output from a fitted emulator

With an ordinary
[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md)
result, supply either new design matrices or saved simulator outputs
used to construct them. With a
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)
result, the offline GP supplies those rows and dt_pde_test is
unnecessary.

## Usage

``` r
emulator_predict(
  emulator,
  input_new,
  dt_pde_test = NULL,
  F_new_ls = "default",
  MC = FALSE
)
```

## Arguments

- emulator:

  An emulator or prepared calibration-input object.

- input_new:

  Matrix or vector of new mechanistic parameter inputs.

- dt_pde_test:

  Optional saved responses at new inputs, used only to construct default
  AR designs for an ordinary emulator.

- F_new_ls:

  "default" or explicit episode design matrices.

- MC:

  If TRUE, return coefficient-conditioned mean samples, not full draws
  of the emulator field. If FALSE, return smoothed means for an ordinary
  emulator, or draw-averaged means for a prepared emulator.

## Value

A time-indexed list of matrices, or arrays with a final posterior draw
dimension when MC = TRUE.

## See also

[`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md),
[`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md)

## Examples

``` r
data("dt_emulation")
d <- dt_emulation
em <- emulator_learn(d$pde_para_train, d$dt_pde_train, Nx = 10, Ny = 10,
                     nsam = 5, seed = 1)
pred <- emulator_predict(em, d$pde_para_test, dt_pde_test = d$dt_pde_test)
dim(pred[[2]])
#> [1]  10 100
sqrt(mean((unlist(pred) - unlist(d$dt_pde_test))^2))
#> [1] 0.1241679
```
