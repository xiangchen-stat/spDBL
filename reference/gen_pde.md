# Simulate a spatially extended SIR PDE model

Numerically integrates the spatial SIR ODE
([`SIR`](https://xiangchen-stat.github.io/spDBL/reference/SIR.md)) on a
two-dimensional grid using
[`deSolve::ode.2D`](https://rdrr.io/pkg/deSolve/man/ode.2D.html) and
returns the log(1 + I) values of the infected compartment.

## Usage

``` r
gen_pde(eta, Nx, Ny, N, nT_ori)
```

## Arguments

- eta:

  Numeric vector of length 5. PDE parameters: transmission rate,
  recovery rate, and diffusion coefficients for S, I, R.

- Nx:

  Integer. Number of grid cells in the x-direction (rows).

- Ny:

  Integer. Number of grid cells in the y-direction (columns).

- N:

  Numeric. Total population size.

- nT_ori:

  Integer. Number of time steps to simulate (times \\0, 1, \ldots,
  \text{nT\\ori} - 1\\).

## Value

Numeric matrix of dimension `c(nT_ori, Nx * Ny)` containing \\\log(1 +
I_t)\\ for each time step and grid cell.

## See also

[`SIR`](https://xiangchen-stat.github.io/spDBL/reference/SIR.md),
[`update_y_eta`](https://xiangchen-stat.github.io/spDBL/reference/update_y_eta.md)

## Examples

``` r
if (requireNamespace("deSolve", quietly = TRUE) &&
    requireNamespace("ReacTran", quietly = TRUE)) {
y <- gen_pde(eta = c(0.5, 0.1, 0.05, 0.05, 0.05), Nx = 4, Ny = 4, N = 100,
             nT_ori = 3)
dim(y)  # nT_ori x (Nx * Ny)
round(rowSums(y), 3)  # total log(1 + I) grows over time
}
#> [1] 0.693 1.048 1.521
```
