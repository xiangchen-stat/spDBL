# Prepare PDE emulator training and testing data from CSV files

Reads a matrix of PDE input parameters and a set of PDE output CSV
files, splits the inputs into training and testing sets, and reshapes
the simulator outputs into time-indexed lists for emulator fitting and
evaluation.

## Usage

``` r
prepare_data(file_para, file_data, prop_train = 0.8)
```

## Arguments

- file_para:

  Character scalar. Path to a CSV file containing the PDE input
  parameter matrix. Rows correspond to simulator runs.

- file_data:

  Character vector. Paths to CSV files containing PDE outputs, one file
  per simulator run. The length should match the number of rows in
  `file_para`.

- prop_train:

  Numeric scalar in \\(0, 1)\\. Proportion of simulator runs used for
  training. Defaults to `0.8`.

## Value

A named list with:

- pde_para_train:

  Numeric matrix of training input parameters.

- pde_para_test:

  Numeric matrix of testing input parameters.

- dt_pde_train:

  Named list of training PDE outputs by time step.

- dt_pde_test:

  Named list of testing PDE outputs by time step.

- n_train:

  Integer. Number of training simulator runs.

- n_test:

  Integer. Number of testing simulator runs.

## Details

The function assumes that each file in `file_data` contains a matrix
whose rows correspond to time steps and whose columns correspond to
spatial locations. The number of time steps is inferred from the first
PDE output file.

## See also

[`emulator_learn`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md),
[`emulator_predict`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)

## Examples

``` r
data("dt_emulation")
d <- dt_emulation
dir <- tempfile("spdbl_")
dir.create(dir)
file_para <- file.path(dir, "para.csv")
utils::write.csv(d$pde_para_train[1:5, ], file_para, row.names = FALSE)
file_data <- file.path(dir, paste0("run_", 1:5, ".csv"))
for (i in 1:5) utils::write.csv(t(sapply(d$dt_pde_train, function(y) y[i, ])),
                                file_data[i], row.names = FALSE)
dat <- prepare_data(file_para, file_data, prop_train = 0.6)
c(n_train = dat$n_train, n_test = dat$n_test)
#> n_train  n_test 
#>       3       2 
dim(dat$dt_pde_train$T1)
#> [1]   3 100
unlink(dir, recursive = TRUE)
```
