# Example emulation dataset

A list of example training and testing objects used to illustrate the
emulation workflow in the package.

## Usage

``` r
dt_emulation
```

## Format

A list with 4 elements:

- pde_para_train:

  Training set of PDE parameter values.

- pde_para_test:

  Test set of PDE parameter values.

- dt_pde_train:

  Training set of PDE-based outputs corresponding to `pde_para_train`.

- dt_pde_test:

  Test set of PDE-based outputs corresponding to `pde_para_test`.

## Source

Simulated data / package example data

## Examples

``` r
data("dt_emulation")
d <- dt_emulation
dim(d$pde_para_train)
#> [1] 40  5
length(d$dt_pde_train)
#> [1] 26
dim(d$dt_pde_train$T1)
#> [1]  40 100
```
