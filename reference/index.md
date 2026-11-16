# Package index

## Emulation

Fit a dynamic Bayesian emulator to simulator runs and predict the
spatiotemporal output at new inputs. Exact (MCMC-free) posterior
updates.

- [`prepare_data()`](https://xiangchen-stat.github.io/spDBL/reference/prepare_data.md)
  : Prepare PDE emulator training and testing data from CSV files
- [`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md)
  : Fit an FFBS-based dynamic emulator
- [`emulator_predict()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_predict.md)
  : Predict simulator output from a fitted emulator
- [`print(`*`<spDBL_emulator>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/spDBL_emulator-methods.md)
  [`summary(`*`<spDBL_emulator>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/spDBL_emulator-methods.md)
  [`predict(`*`<spDBL_emulator>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/spDBL_emulator-methods.md)
  : Print, summarize and predict from a fitted emulator

## Calibration

Infer mechanistic parameters from field data without new PDE solves.
Requires an emulator fitted with `covariance = "sigma2R"`; the sampler
uses RStan (suggested package).

- [`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)
  : Prepare a fitted emulator for PDE-free calibration
- [`print(`*`<spDBL_calibration_input>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/spDBL_calibration_input-methods.md)
  [`predict(`*`<spDBL_calibration_input>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/spDBL_calibration_input-methods.md)
  : Print and predict from a prepared calibration object
- [`calibration_validate()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_validate.md)
  : Validate offline calibration predictions against held-out
  simulations
- [`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md)
  : Fit modular calibration conditional on emulator posterior draws
- [`calibration_predict()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_predict.md)
  : Recover fields at the calibrated observation times and locations
- [`summary(`*`<spDBL_calibration>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/summary.spDBL_calibration.md)
  : Summarize modular calibration samples
- [`plot(`*`<spDBL_calibration>`*`)`](https://xiangchen-stat.github.io/spDBL/reference/plot.spDBL_calibration.md)
  : Plot calibrated mechanistic parameter samples

## Example data

- [`dt_emulation`](https://xiangchen-stat.github.io/spDBL/reference/dt_emulation.md)
  : Example emulation dataset
- [`dt_calibration`](https://xiangchen-stat.github.io/spDBL/reference/dt_calibration.md)
  : Saved SIR simulations and synthetic field observations for
  calibration

## Forward filtering and backward sampling

Advanced. The matrix-variate dynamic linear model kernels behind
[`emulator_learn()`](https://xiangchen-stat.github.io/spDBL/reference/emulator_learn.md).
Three covariance models: inverse-Wishart (no suffix), scalar variance
times a fixed correlation (`_sigma2R`), and identity (`_I`).

- [`FF()`](https://xiangchen-stat.github.io/spDBL/reference/FF.md) :
  Forward Filter for the MNIW dynamic linear model
- [`FF_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_sigma2R.md)
  : Forward Filter for the scalar-sigma-squared-times-R model
- [`FF_I()`](https://xiangchen-stat.github.io/spDBL/reference/FF_I.md) :
  Forward Filter with identity right-covariance
- [`FF_1step_R_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_sigma2R.md)
  : Single forward filter step (scalar right-covariance, sigma-squared
  times R)
- [`FF_1step_R_I()`](https://xiangchen-stat.github.io/spDBL/reference/FF_1step_R_I.md)
  : Single forward filter step (identity right-covariance)
- [`BS()`](https://xiangchen-stat.github.io/spDBL/reference/BS.md) :
  Backward sampler for the Forward Filter Backward Sampler (FFBS)
- [`FFBS()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS.md) :
  Forward Filter Backward Sampler (MNIW model)
- [`FFBS_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sigma2R.md)
  : Forward Filter Backward Sampler (scalar sigma-squared-times-R model)
- [`FFBS_I()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_I.md)
  : Forward Filter Backward Sampler (identity right-covariance)
- [`FFBS_sampling()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling.md)
  : Draw posterior samples from FFBS output (MNIW model)
- [`FFBS_sampling_sigma2R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_sigma2R.md)
  : Draw posterior samples from FFBS output (scalar
  sigma-squared-times-R model)
- [`FFBS_sampling_I()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_sampling_I.md)
  : Draw posterior samples from FFBS output (identity right-covariance)
- [`FFBS_predict_exact()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_exact.md)
  : Exact posterior predictive mean using FFBS smoothed states (MNIW
  model)
- [`FFBS_predict_MC()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_MC.md)
  : Monte Carlo prediction using FFBS output (MNIW model)

## Large data and episode partitioning

Out-of-core filtering, smoothing and prediction from CSV blocks, grid
traversal, autoregressive design construction, and mapping episode
blocks back to time.

- [`FF_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FF_bigdata_R.md)
  : Forward Filter for big data stored in CSV files (MNIW model)
- [`BS_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/BS_bigdata_R.md)
  : Out-of-core backward smoothing
- [`FFBS_predict_bigdata_R()`](https://xiangchen-stat.github.io/spDBL/reference/FFBS_predict_bigdata_R.md)
  : Out-of-core exact prediction
- [`read_big_csv_quick()`](https://xiangchen-stat.github.io/spDBL/reference/read_big_csv_quick.md)
  : Read a rectangular block from a large CSV file
- [`generate_grid()`](https://xiangchen-stat.github.io/spDBL/reference/generate_grid.md)
  : Generate block indices for big data grid traversal
- [`generate.grid.exact()`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.exact.md)
  : Generate an exact block grid analytically
- [`generate.grid.rowsnake()`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.rowsnake.md)
  : Generate a flexible block grid with snake traversal
- [`generate.grid.lr()`](https://xiangchen-stat.github.io/spDBL/reference/generate.grid.lr.md)
  : Generate a flexible block grid with left-to-right traversal
- [`gen_F_ls_AR1()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1.md)
  : Build AR(1) covariate list from a list of response matrices
- [`gen_F_ls_AR2()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2.md)
  : Build AR(2) covariate list from a list of response matrices
- [`gen_F_ls_AR1_EP()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR1_EP.md)
  : Build AR(1) covariate list for the episode-block model
- [`gen_F_ls_AR2_EP()`](https://xiangchen-stat.github.io/spDBL/reference/gen_F_ls_AR2_EP.md)
  : Build AR(2) covariate list for the episode-block model
- [`recover_from_EP_exact()`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_exact.md)
  : Recover episode-partitioned data to original time dimension (exact)
- [`recover_from_EP_MC()`](https://xiangchen-stat.github.io/spDBL/reference/recover_from_EP_MC.md)
  : Recover episode-partitioned posterior samples to original time
  dimension

## Distributions and conjugate updates

Advanced. Matrix-normal, MNIW and MNIG sampling and updates.

- [`MNIW_R()`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R.md)
  : MNIW posterior update
- [`MNIW_R_naiive()`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_R_naiive.md)
  : Naive MNIW posterior update
- [`MNIW_sampler()`](https://xiangchen-stat.github.io/spDBL/reference/MNIW_sampler.md)
  : Sample from the Matrix Normal Inverse Wishart (MNIW) distribution
- [`MNIG_sampler()`](https://xiangchen-stat.github.io/spDBL/reference/MNIG_sampler.md)
  : Sample from the Matrix Normal Inverse Gamma (MNIG) distribution
- [`rmn_chol()`](https://xiangchen-stat.github.io/spDBL/reference/rmn_chol.md)
  : Draw one sample from a matrix-normal distribution (Cholesky
  parameterisation)
- [`rmn_chol_more()`](https://xiangchen-stat.github.io/spDBL/reference/rmn_chol_more.md)
  : Draw multiple samples from a matrix-normal distribution (Cholesky
  parameterisation)
- [`dMTig()`](https://xiangchen-stat.github.io/spDBL/reference/dMTig.md)
  : Log density of the matrix-T distribution with inverse-gamma right
  covariance
- [`check_pds()`](https://xiangchen-stat.github.io/spDBL/reference/check_pds.md)
  : Check and repair a matrix to be positive definite and symmetric
- [`make_pds()`](https://xiangchen-stat.github.io/spDBL/reference/make_pds.md)
  : Force a matrix to be positive definite and symmetric
- [`inv_chol()`](https://xiangchen-stat.github.io/spDBL/reference/inv_chol.md)
  : Invert a matrix via its Cholesky factorisation
- [`gen_pd_matrix()`](https://xiangchen-stat.github.io/spDBL/reference/gen_pd_matrix.md)
  : Generate a random positive definite matrix
- [`gen_ran_matrix()`](https://xiangchen-stat.github.io/spDBL/reference/gen_ran_matrix.md)
  : Generate a random matrix with entries scaled to \\\[-1, 1\]\\

## Model assessment

One-step log posterior predictive densities for each covariance model,
and WAIC for MNIW multivariate regression.

- [`lppd_IW_1t()`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IW_1t.md)
  : One-step log posterior predictive density (MNIW / inverse-Wishart
  model)
- [`lppd_IG_1t()`](https://xiangchen-stat.github.io/spDBL/reference/lppd_IG_1t.md)
  : One-step log posterior predictive density (scalar
  sigma-squared-times-R model)
- [`lppd_id_1t()`](https://xiangchen-stat.github.io/spDBL/reference/lppd_id_1t.md)
  : One-step log posterior predictive density (identity right-covariance
  model)
- [`compute_analytical_waic()`](https://xiangchen-stat.github.io/spDBL/reference/compute_analytical_waic.md)
  : Analytical WAIC for an MNIW multivariate regression posterior
- [`compute_sample_waic()`](https://xiangchen-stat.github.io/spDBL/reference/compute_sample_waic.md)
  : WAIC from paired coefficient and covariance draws

## Gaussian process kernels

- [`gen_gp_kernel()`](https://xiangchen-stat.github.io/spDBL/reference/gen_gp_kernel.md)
  : Compute a Gaussian Process covariance kernel matrix
- [`gen_exp_kernel()`](https://xiangchen-stat.github.io/spDBL/reference/gen_exp_kernel.md)
  : Compute an exponential GP kernel matrix
- [`gen_expsq_kernel()`](https://xiangchen-stat.github.io/spDBL/reference/gen_expsq_kernel.md)
  : Compute a squared-exponential (Gaussian) GP kernel matrix
- [`gen_Jt()`](https://xiangchen-stat.github.io/spDBL/reference/gen_Jt.md)
  : Compute the cross-covariance matrix between observed and new
  locations

## Simulators and synthetic data

- [`SIR()`](https://xiangchen-stat.github.io/spDBL/reference/SIR.md) :
  Right-hand side of the spatially extended SIR ODE
- [`gen_pde()`](https://xiangchen-stat.github.io/spDBL/reference/gen_pde.md)
  : Simulate a spatially extended SIR PDE model
- [`gen_ffbs_data()`](https://xiangchen-stat.github.io/spDBL/reference/gen_ffbs_data.md)
  : Generate synthetic FFBS data in memory
- [`gen_ffbs_csv()`](https://xiangchen-stat.github.io/spDBL/reference/gen_ffbs_csv.md)
  : Generate synthetic FFBS data and write to CSV files
- [`gen_calibrate_data()`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data.md)
  : Generate synthetic calibration data with correlated discrepancy
- [`gen_calibrate_data_uncorr()`](https://xiangchen-stat.github.io/spDBL/reference/gen_calibrate_data_uncorr.md)
  : Generate synthetic calibration data with uncorrelated discrepancy
- [`gen_prior_u_tau2()`](https://xiangchen-stat.github.io/spDBL/reference/gen_prior_u_tau2.md)
  : Sample prior discrepancy trajectory and variance sequence

## Plotting

- [`quick_heat()`](https://xiangchen-stat.github.io/spDBL/reference/quick_heat.md)
  : Quick raster heatmap
- [`quick_save()`](https://xiangchen-stat.github.io/spDBL/reference/quick_save.md)
  : Save a ggplot to a timestamped PNG file
- [`plot_panel_heatmap_9()`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9.md)
  : Plot a 3-by-3 panel of heatmaps across selected time stamps
- [`plot_panel_heatmap_9_cal()`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9_cal.md)
  : Plot a 3-by-3 panel of calibration heatmaps
- [`plot_panel_heatmap_9_cal_nolab()`](https://xiangchen-stat.github.io/spDBL/reference/plot_panel_heatmap_9_cal_nolab.md)
  : Plot a 3-by-3 panel of calibration heatmaps without axis labels
- [`cal_errorbar()`](https://xiangchen-stat.github.io/spDBL/reference/cal_errorbar.md)
  : Compute median and 95% credible interval across rows
- [`cal_errorbar_mean()`](https://xiangchen-stat.github.io/spDBL/reference/cal_errorbar_mean.md)
  : Compute mean and 95% credible interval across rows

## Legacy calibration helpers

Building blocks from the research code of Banerjee et al. (2025), where
calibration solved the PDE at every MCMC iteration. No other package
function calls them. New work should use
[`calibration_prepare()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_prepare.md)
and
[`calibration_fit()`](https://xiangchen-stat.github.io/spDBL/reference/calibration_fit.md).

- [`update_y_eta()`](https://xiangchen-stat.github.io/spDBL/reference/update_y_eta.md)
  : Update the likelihood of observations given PDE parameters (Monte
  Carlo)
- [`update_y_eta_one()`](https://xiangchen-stat.github.io/spDBL/reference/update_y_eta_one.md)
  : Update the likelihood of observations given PDE parameters (single
  sample)
- [`update_muSigma_eta_one()`](https://xiangchen-stat.github.io/spDBL/reference/update_muSigma_eta_one.md)
  : Compute posterior predictive mean and covariance without sampling
  (single sample)
- [`sample_y_eta_one()`](https://xiangchen-stat.github.io/spDBL/reference/sample_y_eta_one.md)
  : Draw predictive samples from a precomputed mean and covariance
- [`cal_Bt_bt()`](https://xiangchen-stat.github.io/spDBL/reference/cal_Bt_bt.md)
  : Update the posterior mean and covariance of the discrepancy field
- [`scale_uniform()`](https://xiangchen-stat.github.io/spDBL/reference/scale_uniform.md)
  : Scale a vector to the unit interval via a uniform transformation
- [`scale_back_uniform()`](https://xiangchen-stat.github.io/spDBL/reference/scale_back_uniform.md)
  : Invert a uniform scaling transformation
- [`expit()`](https://xiangchen-stat.github.io/spDBL/reference/expit.md)
  : Logistic (expit) function
- [`cal_jacobian_logit_uniform()`](https://xiangchen-stat.github.io/spDBL/reference/cal_jacobian_logit_uniform.md)
  : Log absolute Jacobian of the logit-uniform transformation

## Package

- [`spDBL`](https://xiangchen-stat.github.io/spDBL/reference/spDBL-package.md)
  [`spDBL-package`](https://xiangchen-stat.github.io/spDBL/reference/spDBL-package.md)
  : spDBL: Dynamic Bayesian Learning for Spatiotemporal Mechanistic
  Models
