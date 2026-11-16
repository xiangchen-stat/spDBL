# Right-hand side of the spatially extended SIR ODE

Defines the right-hand side of a spatial SIR
(Susceptible-Infected-Recovered) partial differential equation suitable
for use with
[`deSolve::ode.2D`](https://rdrr.io/pkg/deSolve/man/ode.2D.html).
Diffusion is implemented via
[`ReacTran::tran.2D`](https://rdrr.io/pkg/ReacTran/man/tran.2D.html)
with zero-flux boundary conditions.

## Usage

``` r
SIR(t, y, parms)
```

## Arguments

- t:

  Numeric. Current time (passed by the ODE solver).

- y:

  Numeric vector. Flattened state vector of length \\3 N_x N_y\\: \\S\\
  values, then \\I\\ values, then \\R\\ values.

- parms:

  Named list of parameters:

  Nx

  :   Integer. Number of grid cells in the x-direction.

  Ny

  :   Integer. Number of grid cells in the y-direction.

  N

  :   Numeric. Total population size.

  eta

  :   Numeric vector of length 5: `eta[1]` transmission rate, `eta[2]`
      recovery rate, `eta[3]` diffusion coefficient for S, `eta[4]`
      diffusion coefficient for I, `eta[5]` diffusion coefficient for R.

  Gridx

  :   1D grid object for x (from
      [`ReacTran::setup.grid.1D`](https://rdrr.io/pkg/ReacTran/man/setup.grid.1D.html)).

  Gridy

  :   1D grid object for y.

## Value

A list with one element: a numeric vector of derivatives of the same
length as `y`, as required by `deSolve`.

## See also

[`gen_pde`](https://xiangchen-stat.github.io/spDBL/reference/gen_pde.md)

## Examples

``` r
if (requireNamespace("deSolve", quietly = TRUE) &&
    requireNamespace("ReacTran", quietly = TRUE)) {
set.seed(1)
parms <- list(eta = c(0.3, 0.1, 0, 0, 0), Nx = 2, Ny = 2, N = 100,
              Gridx = ReacTran::setup.grid.1D(x.up = 0, x.down = 2, N = 2),
              Gridy = ReacTran::setup.grid.1D(x.up = 0, x.down = 2, N = 2))
y <- c(runif(4, 50, 100), runif(4, 0, 10), rep(0, 4))
d <- SIR(t = 0, y = y, parms = parms)[[1]]
length(d)
all.equal(sum(d), 0)  # without diffusion, S + I + R is conserved
}
#> [1] TRUE
```
