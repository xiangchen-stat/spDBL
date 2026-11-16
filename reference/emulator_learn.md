# Fit an FFBS-based dynamic emulator

Fits an episode-based emulator to saved simulator runs. The default
inverse-Wishart model is retained for compatibility. Use covariance =
"sigma2R" for scalar-variance coefficient samples that can be passed to
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md).

## Usage

``` r
emulator_learn(
  pde_para_train,
  dt_pde_train,
  Nx,
  Ny,
  F_ls_train = "default",
  F_ls_test = "default",
  N_people = 10000,
  nsam = 10,
  AR_choice = 2,
  episode_window = 1,
  gp_tune = 0.5,
  gp_sigma2 = 1.1,
  gp_tau2 = 1e-04,
  covariance = c("iw", "sigma2R"),
  R_episode = NULL,
  R_field = NULL,
  seed = NULL,
  verbose = FALSE,
  sampling = c("marginal", "joint"),
  keep_draws = TRUE
)
```

## Arguments

- pde_para_train:

  Finite numeric matrix of training parameter inputs.

- dt_pde_train:

  Time-indexed list of training response matrices, with simulator runs
  in rows and original-grid locations in columns.

- Nx, Ny:

  Positive integer grid dimensions.

- F_ls_train:

  "default" or a list of episode design matrices.

- F_ls_test:

  Unused compatibility argument.

- N_people:

  Population scale recorded in the emulator metadata.

- nsam:

  Positive integer number of posterior coefficient draws.

- AR_choice:

  Autoregressive order, 1 or 2.

- episode_window:

  Number of spatial strips per block; must divide Nx.

- gp_tune:

  Positive range tuning factor.

- gp_sigma2, gp_tau2:

  Positive input/spatial kernel amplitude and nugget.

- covariance:

  "iw" (default) or "sigma2R".

- R_episode, R_field:

  Optional spatial covariance matrices for one block and for the full
  grid. `R_episode` defaults to the exponential kernel over the sites of
  a block. `R_field` is used only by calibration. Unless it is supplied
  it is not formed here, since an S x S matrix costs \\O(S^2)\\ memory
  and its check \\O(S^3)\\ time; calibration builds the same exponential
  kernel at the observed sites when it needs it.

- seed:

  Optional local RNG seed. NULL follows the session RNG.

- verbose:

  Whether to report progress.

- sampling:

  "marginal" (default) uses smoothed episode marginals conditional on a
  shared covariance draw, following Eq. (9) and Algorithm 4. "joint"
  uses backward conditionals to preserve dependence across episodes.

- keep_draws:

  If TRUE (default), store every coefficient draw. If FALSE (marginal
  sampling only), store only where each episode's draws start in the
  random-number stream and regenerate identical draws when prediction or
  calibration needs them. This removes the largest part of a large
  emulator's memory at the cost of re-drawing (about the original
  sampling time per full pass over the draws). Each episode then costs
  about 2.5 KB instead of its draws, so it pays off once an episode's
  draws exceed that.

## Value

An object of class spDBL_emulator, retaining para_ffbs, res_ffbs, and a
complete setup list for prediction and calibration.

## See also

[`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md),
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)

## Examples

``` r
data("dt_emulation")
d <- dt_emulation
em <- emulator_learn(d$pde_para_train, d$dt_pde_train, Nx = 10, Ny = 10,
                     nsam = 5, seed = 1)
class(em)
#> [1] "spDBL_emulator"
c(n_train = em$setup$n_train, episodes = em$setup$nT, locations = em$setup$N_sp)
#>   n_train  episodes locations 
#>        40       260       100 
```
