# Recover fields at the observed sites and predict them at new locations

Recovers discrepancy trajectories when needed, then draws latent
emulator fields conditional on the observations and fresh replicated
observations at the observed locations. With `sites` (grid sites without
observations) or `coords` (points anywhere in the grid's cells), it also
predicts the emulator field, the discrepancy and new field observations
there, at the observation times. It does not forecast new times.

## Usage

``` r
calibration_predict(
  object,
  seed = 1234L,
  sites = NULL,
  joint = FALSE,
  m0_sites = NULL,
  M0_sites = NULL,
  coords = NULL,
  discrepancy_kernel = NULL,
  field = c("cell", "bilinear")
)
```

## Arguments

- object:

  Result from
  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md).

- seed:

  Local RNG seed for Gaussian recovery.

- sites:

  Optional grid indices, numbered as `locations` in
  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md),
  of sites without observations to predict at.

- joint:

  If `FALSE` (default), each site or point is drawn from its own
  predictive distribution, so draws at different ones are independent;
  this avoids site-by-site matrices. If `TRUE`, they are drawn jointly,
  as needed for sums or maps over several sites or points.

- m0_sites, M0_sites:

  Prior mean vector and covariance scale matrix of the starting
  discrepancy at `sites` or `coords`. The defaults, zero and the
  identity, match those of
  [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md).

- coords:

  Optional two-column matrix of points to predict at, in the grid's
  coordinates: the units of the grid's row and column indices, with cell
  centres at integers (as `grid` in
  [dt_calibration](https://xiangchen-stat.github.io/spDBL/reference/dt_calibration.md)).
  Each point must lie inside the grid's cells. Give `sites` or `coords`,
  not both.

- discrepancy_kernel:

  With `coords`: the discrepancy covariance function, called as
  `discrepancy_kernel(a, b)` with two coordinate matrices and returning
  the `nrow(a)` by `nrow(b)` matrix of its values.

- field:

  With `coords`: `"cell"` (default) or `"bilinear"`.

## Value

A list containing y, u, and z_rep arrays, each with dimensions posterior
draws by observation times by observed locations. With `sites`, also
`sites`, a list with the site indices `locations` and arrays `y`, `u`
and `z` of dimensions posterior draws by observation times by sites.
With `coords`, also `points`, a list with `coords`, the grid `cell`
containing each point, `field`, and arrays `y`, `u` and `z` of
dimensions posterior draws by observation times by points.

## Details

Predictions are composition draws: each one conditions on the same
posterior draw of the field and the discrepancy at the observed sites.
At each time, the emulator field at the observed sites and at grid sites
is Gaussian with the mean and covariance that calibration uses,
design-row term included; at `sites` it is drawn given its values at the
observed sites. The emulator is defined only on its grid, so at `coords`
the field is that of the cell containing the point (`field = "cell"`,
exact for simulators whose outputs are cell averages, such as
[`gen_pde()`](https://xiangchen-stat.github.io/spDBL/reference/gen_pde.md))
or the bilinear interpolation of the four surrounding cell centres
(`field = "bilinear"`, which assumes the field is smooth within a cell);
the difference from the true field at a point is left to the
discrepancy.

The discrepancy's increments are drawn given their values at the
observed sites. At `sites` their covariance scale comes from `U`, the
full-grid matrix given to
[`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md);
at `coords` it comes from `discrepancy_kernel`, the covariance function
behind `U`, evaluated at the points themselves (it must reproduce `U` at
the observed sites). The starting discrepancy at the new locations has
prior mean `m0_sites` and covariance \\\tau_0^2\\ `M0_sites`,
independent of its starting value at the observed sites. This extends
the model without changing the calibration posterior. New observations
add independent noise of variance \\\tau_t^2\\. The discrepancy at new
locations is informed only through `U` or the kernel, which is fixed,
not estimated.

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("rstan", quietly = TRUE)
data("dt_calibration")
d <- dt_calibration
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
                     covariance = "sigma2R", episode_window = 2, nsam = 3,
                     seed = 1)
cal <- calibration_prepare(em, d$lower, d$upper)
# Leave two grid sites out of the calibration, then predict them.
held_out <- c(8, 15)
fit <- calibration_fit(cal, d$field$z[, -held_out], d$field$U,
                       locations = d$field$locations[-held_out],
                       draws = 1:2, chains = 1, iter = 200, seed = 1)
pred <- calibration_predict(fit, sites = held_out)
dim(pred$y)
dim(pred$sites$z)
# Predictive medians against the held-out observations at the first time
rbind(median = apply(pred$sites$z[, 1, ], 2, stats::median),
      observed = d$field$z[1, held_out])
# Predict at two points between cell centres, with the kernel behind U.
phi <- 3 / (0.5 * max(dist(d$grid)))
kernel <- function(a, b) {
  h <- sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
  (1.1 * exp(-phi * h) + 1e-4 * (h == 0)) / 1.1001
}
points <- rbind(c(2.5, 3.25), c(4.6, 1.2))
off <- calibration_predict(fit, coords = points, discrepancy_kernel = kernel,
                           field = "bilinear")
dim(off$points$z)
}
```
