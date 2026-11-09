# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# These tests exercise the R side of the rstan back end (argument checks,
# model caching, draw selection, reshaping, diagnostics, checkpoints and
# post-processing) with rstan's entry points replaced by a deterministic fake.
# The compiled Stan programs themselves are checked in test-calibration-stan.R.

calibration_fit_fixture <- function() {
  f <- calibration_fixture()
  f$U <- 1.1001 * f$U
  f
}

# Draws whose values encode their (iteration, chain) position.
fake_draw_array <- function(data, chains, n_iter, with_u) {
  pars <- c(paste0("eta[", seq_len(data$D_full), "]"),
            paste0("tau2[", seq_len(data$T + 1), "]"))
  if (with_u) pars <- c(pars, paste0("u[", rep(seq_len(data$T + 1), data$S), ",",
                                     rep(seq_len(data$S), each = data$T + 1), "]"))
  pars <- c(pars, "lp__")
  raw <- array(0, c(n_iter, chains, length(pars)),
               dimnames = list(iterations = NULL, chains = paste0("chain:", seq_len(chains)),
                               parameters = pars))
  it <- matrix(seq_len(n_iter), n_iter, chains)
  ch <- matrix(seq_len(chains), n_iter, chains, byrow = TRUE)
  raw[, , "eta[1]"] <- 0.1 + 0.01 * it
  raw[, , "eta[2]"] <- 0.1 + 0.1 * ch
  for (k in seq_len(data$T + 1)) raw[, , paste0("tau2[", k, "]")] <- 0.1 * k + 0.001 * it
  if (with_u) for (t in seq_len(data$T + 1)) for (s in seq_len(data$S))
    raw[, , paste0("u[", t, ",", s, "]")] <- t + s / 10 + it / 1000 + ch / 10000
  raw
}

# Install the fake rstan back end for the calling test; returns a call log.
local_fake_rstan <- function(divergent_call = 0L, rhat = 1.001, env = parent.frame()) {
  log <- new.env()
  log$calls <- list()
  log$compiled <- character()
  log$optimizing <- integer()
  local_mocked_bindings(
    sampling = function(object, data, init, chains, iter, warmup, seed, cores,
                        refresh, pars, control, ...) {
      k <- length(log$calls) + 1L
      log$calls[[k]] <- list(object = object, data = data, init = init(),
                             init_fun = init,
                             chains = chains, iter = iter, warmup = warmup,
                             seed = seed, pars = pars, control = control)
      n_iter <- iter - warmup
      structure(list(raw = fake_draw_array(data, chains, n_iter, "u" %in% pars),
                     chains = chains, n_iter = n_iter, divergent = k == divergent_call),
                class = "fake_stanfit")
    },
    extract = function(object, permuted = TRUE, inc_warmup = FALSE, ...) object$raw,
    get_sampler_params = function(object, inc_warmup = TRUE, ...)
      lapply(seq_len(object$chains), function(c) cbind(
        divergent__ = c(as.numeric(object$divergent), rep(0, object$n_iter - 1)),
        treedepth__ = c(12, rep(3, object$n_iter - 1)))),
    # Two local optima: starts with eta_unit[1] < 0.5 reach (0.3, ...) with
    # log density -10, the others (0.8, ...) with -20; one start fails.
    optimizing = function(object, data, init, seed, ...) {
      log$optimizing <- c(log$optimizing, seed)
      eta <- as.numeric(init$eta_unit)
      low <- eta[1] < 0.5
      list(par = list(eta_unit = array(rep(if (low) 0.3 else 0.8, length(eta)), length(eta)),
                      precision = rep(if (low) 20 else 5, length(init$precision))),
           value = if (low) -10 else -20,
           return_code = if (length(log$optimizing) == 2L) 70L else 0L)
    },
    stan_model = function(file, ...) {
      log$compiled <- c(log$compiled, basename(file))
      methods::new("stanmodel", model_name = basename(file))
    },
    Rhat = function(sims) rhat,
    ess_bulk = function(sims) 51,
    ess_tail = function(sims) 47,
    .package = "rstan", .env = env)
  local_mocked_s3_method("summary", "fake_stanfit", function(object, pars, ...) {
    names <- dimnames(object$raw)[[3]]
    names <- names[names != "lp__"]
    # The classic Rhat is deliberately bad: the package must replace it.
    s <- cbind(mean = rep(0, length(names)), n_eff = seq_along(names) + 50,
               Rhat = rep(9, length(names)))
    rownames(s) <- names
    list(summary = s)
  }, frame = env)
  log
}

skip_without_fake_rstan <- function() {
  skip_if_not_installed("rstan")
  # Mocking an S3 method that does not exist yet needs testthat 3.3.2.
  skip_if_not(utils::packageVersion("testthat") >= "3.3.2",
              "testthat >= 3.3.2 is needed to mock S3 methods")
}

# Empty the package's in-memory Stan model cache, returning its contents so
# the caller can put them back with restore_stan_cache() in on.exit().
clear_stan_cache <- function() {
  cache <- spDBL:::.spdbl_stan_cache
  saved <- as.list(cache, all.names = TRUE)
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  saved
}

restore_stan_cache <- function(saved) {
  cache <- spDBL:::.spdbl_stan_cache
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  list2env(saved, cache)
  invisible(NULL)
}

test_that("calibration_fit validates its inputs before touching Stan", {
  f <- calibration_fit_fixture()
  e <- f$emulator
  fit <- function(..., iter = 1L) calibration_fit(e, f$z, f$U,
    compiled_model = "unused", iter = iter, ...)
  expect_error(calibration_fit(e, f$z, f$U), "Set iter")
  expect_error(calibration_fit(e, f$z, f$U, iter = 1L, seed = NULL), "integer seed")
  expect_error(calibration_fit(list(), f$z, f$U, iter = 1L), "calibration_prepare")
  expect_error(fit(backend = "cmdstanr"), "Only the rstan backend")
  expect_error(fit(iter = 4L, warmup = 4L), "Invalid calibration sampling settings")
  expect_error(fit(draws = c(1L, 1L)), "Invalid calibration sampling settings")
  expect_error(fit(draws = 3L), "Invalid calibration sampling settings")
  expect_error(fit(iter = 4L, warmup = 2L, draws_per_emulator = 3L),
               "Invalid calibration sampling settings")
  expect_error(calibration_fit(e, f$z, f$U[1:3, 1:3], iter = 1L), "full-grid covariance")
  expect_error(fit(n0 = Inf), "check priors")
  expect_error(calibration_fit(e, f$z[, 1:3], f$U, compiled_model = "unused", iter = 1L),
               "Invalid calibration observations")
  expect_error(spDBL:::fit_calibration_rstan(e, f$z, f$U, iter = 4, warmup = 4,
                                             compiled_model = "unused"),
               "Invalid sampling settings")
})

test_that("summaries, plots and predictions work on stored calibration draws", {
  # A fit object assembled by hand, so this runs without rstan.
  f <- calibration_fit_fixture()
  e <- f$emulator
  eta <- cbind(a = c(0.3, 0.5, 0.4, 0.6), b = c(0.2, 0.4, 0.6, 0.5))
  fit <- structure(list(mcmc_eta = eta, mcmc_tau2 = matrix(0.2, 4, 2),
    tau2_0 = rep(0.3, 4), emulator_draw = c(1L, 1L, 2L, 2L), emulator = e,
    diagnostics = list(
      list(emulator_draw = 1L, divergences = 0, max_treedepth_hits = 1,
           max_rhat = 1.004, min_n_eff = 120),
      list(emulator_draw = 2L, divergences = 3, max_treedepth_hits = 0,
           max_rhat = NA_real_, min_n_eff = NA_real_)),
    field = list(z = f$z, U = f$U, locations = 1:4, m0 = rep(0, 4), M0 = diag(4)),
    settings = list(), execution_only = FALSE,
    predictor_uncertainty = e$predictor$uncertainty), class = "spDBL_calibration")
  s <- summary(fit)
  expect_identical(s$parameters$parameter, c("a", "b"))
  expect_equal(s$parameters$lower, unname(apply(eta, 2, stats::quantile, 0.025)))
  expect_equal(s$parameters$upper, unname(apply(eta, 2, stats::quantile, 0.975)))
  expect_identical(s$emulator_sampling, "joint")
  out <- capture_output(print(fit))
  expect_match(out, "Divergences: 3 ; tree-depth hits: 1")
  expect_match(out, "Largest conditional Rhat: 1.004; smallest conditional bulk ESS: 120")
  expect_match(out, "Design-row GP variance is propagated")
  pdf(NULL)  # keep plot() from opening a file device in the working directory
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(fit), fit)
  pred <- calibration_predict(fit, seed = 9L)
  expect_identical(dim(pred$u), c(4L, 2L, 4L))
  expect_equal(pred$u, spDBL:::recover_calibration_discrepancy(eta, fit$mcmc_tau2,
    fit$tau2_0, fit$emulator_draw, e, f$z, f$U, seed = 9L)$u)
  expect_equal(pred$y, spDBL:::recover_calibration_predictions(eta, fit$mcmc_tau2,
    pred$u, fit$emulator_draw, e, f$z, seed = 10L)$y)
})

test_that("Stan models are compiled once and cached in memory and on disk", {
  skip_without_fake_rstan()
  saved_cache <- clear_stan_cache()
  on.exit(restore_stan_cache(saved_cache), add = TRUE)
  cache <- spDBL:::.spdbl_stan_cache
  log <- local_fake_rstan()
  dir <- tempfile("stancache")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  model <- spDBL:::.cal_stan_model(TRUE, cache = TRUE, cache_dir = dir)
  expect_true(methods::is(model, "stanmodel"))
  expect_identical(log$compiled, "calibration_marginal.stan")
  cached <- list.files(dir, full.names = TRUE)
  expect_length(cached, 1L)
  expect_identical(spDBL:::.cal_stan_model(TRUE, cache = TRUE, cache_dir = dir), model)
  rm(list = ls(cache), envir = cache)
  from_disk <- spDBL:::.cal_stan_model(TRUE, cache = TRUE, cache_dir = dir)
  expect_identical(from_disk@model_name, "calibration_marginal.stan")
  expect_identical(log$compiled, "calibration_marginal.stan")
  # An unreadable cache file is replaced by a fresh compilation.
  rm(list = ls(cache), envir = cache)
  writeLines("not an rds file", cached)
  spDBL:::.cal_stan_model(TRUE, cache = TRUE, cache_dir = dir)
  expect_identical(log$compiled, rep("calibration_marginal.stan", 2))
  expect_true(methods::is(readRDS(cached), "stanmodel"))
  # A cache path that cannot be written only warns.
  rm(list = ls(cache), envir = cache)
  blocked <- tempfile("stanblocked")
  on.exit(unlink(blocked, recursive = TRUE), add = TRUE)
  dir.create(file.path(blocked, basename(cached)), recursive = TRUE)
  warnings <- capture_warnings(spDBL:::.cal_stan_model(TRUE, cache = TRUE, cache_dir = blocked))
  expect_true(any(grepl("Could not persist the Stan cache", warnings)))
  # Without caching nothing is stored.
  keys <- ls(cache)
  offline <- spDBL:::.cal_stan_model(FALSE, cache = FALSE)
  expect_identical(offline@model_name, "calibration_offline.stan")
  expect_identical(ls(cache), keys)
})

test_that("conditional fits keep chain, iteration and state bookkeeping", {
  skip_without_fake_rstan()
  f <- calibration_fit_fixture()
  e <- f$emulator
  log <- local_fake_rstan(divergent_call = 2L)
  checkpoints <- tempfile("checkpoints")
  on.exit(unlink(checkpoints, recursive = TRUE), add = TRUE)
  msgs <- capture_messages(expect_warning(
    fit <- calibration_fit(e, f$z, f$U, chains = 2L, iter = 10L, warmup = 5L,
      draws_per_emulator = 4L, marginal_discrepancy = FALSE, keep_predictions = TRUE,
      control = list(adapt_delta = 0.9), refresh = 0L, checkpoint_dir = checkpoints,
      compiled_model = "fake-model", seed = 21L, init = "midpoint"),
    "Emulator draw 2 needs additional"))
  expect_match(msgs[1], "Calibration conditional on emulator draw 1 \\(1/2\\)")
  expect_s3_class(fit, "spDBL_calibration")
  expect_length(log$calls, 2L)
  for (j in 1:2) {
    call <- log$calls[[j]]
    expect_identical(call$object, "fake-model")
    expect_identical(call$seed, 21L + j)
    expect_identical(call$pars, c("eta", "tau2", "u"))
    expect_equal(call$data$sigma2, e$posterior$sigma2[j])
    expect_equal(call$init$eta_unit, array(c(0.5, 0.5), 2))
    expect_equal(call$init$precision, rep(15, 3))
    expect_equal(call$init$u_std, matrix(0, 3, 4))
  }
  expect_identical(colnames(fit$mcmc_eta), c("a", "b"))
  expect_identical(fit$emulator_draw, rep(1:2, each = 4))
  it <- fit$iteration
  ch <- fit$chain
  expect_equal(unname(fit$mcmc_eta[, 1]), 0.1 + 0.01 * it)
  expect_equal(unname(fit$mcmc_eta[, 2]), 0.1 + 0.1 * ch)
  expect_equal(fit$tau2_0, 0.1 + 0.001 * it)
  for (t in 1:2) expect_equal(fit$mcmc_tau2[, t], 0.1 * (t + 1) + 0.001 * it)
  for (t in 1:2) for (s in 1:4)
    expect_equal(fit$mcmc_u[, t, s], (t + 1) + s / 10 + it / 1000 + ch / 10000)
  for (s in 1:4) expect_equal(fit$u0[, s], 1 + s / 10 + it / 1000 + ch / 10000)
  for (m in 1:2) expect_false(anyDuplicated(paste(ch, it)[fit$emulator_draw == m]) > 0)
  expect_identical(dim(fit$mcmc_y), c(8L, 2L, 4L))
  expect_identical(dim(fit$mcmc_z), c(8L, 2L, 4L))
  expect_equal(vapply(fit$diagnostics, `[[`, numeric(1), "divergences"), c(0, 2))
  expect_equal(vapply(fit$diagnostics, `[[`, numeric(1), "max_treedepth_hits"), c(2, 2))
  expect_equal(fit$diagnostics[[1]]$max_rhat, 1.001)
  expect_equal(fit$diagnostics[[1]]$min_n_eff, 51)
  # The rank-normalised diagnostics replace rstan::summary()'s classic Rhat.
  expect_equal(fit$diagnostics[[1]]$min_ess_tail, 47)
  expect_true(all(fit$diagnostics[[1]]$summary[, "Rhat"] < 9))
  saved <- readRDS(file.path(checkpoints, "conditional_002.rds"))
  expect_equal(saved$samples$eta, fit$mcmc_eta[5:8, ], ignore_attr = TRUE)
  expect_false(saved$marginal_discrepancy)
  expect_false(fit$execution_only)
  expect_identical(fit$settings$emulator_sampling, "joint")
  expect_identical(fit$field$z, f$z)

  s <- summary(fit)
  expect_s3_class(s, "summary.spDBL_calibration")
  expect_equal(s$parameters$mean, unname(colMeans(fit$mcmc_eta)))
  expect_equal(s$parameters$median, unname(apply(fit$mcmc_eta, 2, stats::median)))
  expect_identical(s$n_samples, 8L)
  out <- capture_output(print(fit))
  expect_match(out, "Modular calibration")
  expect_match(out, "8 retained samples from 2 emulator conditionals")
  expect_match(out, "Divergences: 2")
  expect_match(out, "Largest conditional Rhat: 1.001")
  pdf(NULL)  # keep plot() from opening a file device in the working directory
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(fit), fit)

  pred <- calibration_predict(fit, seed = 3L)
  expect_identical(pred$u, fit$mcmc_u)
  expect_identical(dim(pred$y), c(8L, 2L, 4L))
  refit <- fit
  refit$mcmc_u <- NULL
  recovered <- calibration_predict(refit, seed = 3L)
  expected_u <- spDBL:::recover_calibration_discrepancy(fit$mcmc_eta, fit$mcmc_tau2,
    fit$tau2_0, fit$emulator_draw, e, f$z, f$U, seed = 3L)$u
  expect_equal(recovered$u, expected_u)
  expect_error(calibration_predict(list()), "Supply a calibration_fit")
})

test_that("conditional fits run in parallel give the same result", {
  skip_without_fake_rstan()
  skip_on_os("windows")
  f <- calibration_fit_fixture()
  e <- f$emulator
  local_fake_rstan()
  fit_with <- function(parallel_draws) suppressMessages(calibration_fit(e, f$z,
    f$U, chains = 3L, iter = 4L, warmup = 2L, draws_per_emulator = 2L,
    refresh = 0L, compiled_model = "fake-model", seed = 7L, init = "random",
    parallel_draws = parallel_draws))
  serial <- fit_with(1L)
  parallel <- fit_with(2L)
  expect_identical(parallel$mcmc_eta, serial$mcmc_eta)
  expect_identical(parallel$mcmc_tau2, serial$mcmc_tau2)
  expect_identical(parallel$diagnostics, serial$diagnostics)
  expect_identical(parallel$settings$parallel_draws, 2L)
  expect_error(fit_with(0L), "Invalid calibration sampling settings")
})

test_that("chains start at the midpoint, at random, or at given values", {
  skip_without_fake_rstan()
  f <- calibration_fit_fixture()
  e <- f$emulator
  log <- local_fake_rstan()
  fit_with <- function(init, seed = 5L) suppressMessages(calibration_fit(e, f$z,
    f$U, chains = 3L, iter = 4L, warmup = 2L, draws_per_emulator = 2L,
    refresh = 0L, compiled_model = "fake-model", seed = seed, init = init))
  starts <- function(call) t(vapply(1:3, function(k) call$init_fun(chain_id = k)$eta_unit,
                                    numeric(2)))
  midpoint <- fit_with("midpoint")
  expect_identical(midpoint$settings$init, "midpoint")
  expect_equal(starts(log$calls[[1]]), matrix(0.5, 3, 2))
  random <- fit_with("random")
  first <- starts(log$calls[[3]])
  expect_true(all(first > 0 & first < 1))
  expect_gt(min(dist(first)), 0)                        # chains differ
  expect_false(isTRUE(all.equal(first, starts(log$calls[[4]]))))  # and fits differ
  fit_with("random")
  expect_identical(starts(log$calls[[5]]), first)       # reproducible from seed
  # Random starts do not change which samples are retained.
  expect_identical(random$iteration, midpoint$iteration)
  expect_identical(random$chain, midpoint$chain)
  # Parameter values, in training order or by name, map to the unit box.
  lower <- e$predictor$lower
  width <- e$predictor$upper - lower
  values <- lower + c(0.25, 0.8) * width
  given <- fit_with(values)
  expect_equal(starts(log$calls[[7]]), matrix(c(0.25, 0.8), 3, 2, byrow = TRUE))
  expect_identical(given$settings$init, values)
  fit_with(stats::setNames(rev(values), rev(e$parameter_names)))
  expect_equal(starts(log$calls[[9]]), matrix(c(0.25, 0.8), 3, 2, byrow = TRUE))
  expect_error(fit_with("centre"), "midpoint")
  # The mode search is the default.
  default <- suppressMessages(calibration_fit(e, f$z, f$U, chains = 3L, iter = 4L,
    warmup = 2L, draws_per_emulator = 2L, refresh = 0L,
    compiled_model = "fake-model", seed = 5L))
  expect_identical(default$settings$init, "mode")
  expect_identical(default$settings$spDBL_version,
                   as.character(utils::packageVersion("spDBL")))
  expect_false(is.null(default$diagnostics[[1]]$mode_search))
  expect_error(fit_with(values[1]), "one value per parameter")
  expect_error(fit_with(c(lower[1], values[2])), "strictly inside")
  expect_error(fit_with(stats::setNames(values, c("a", "z"))), "parameter names")
})

test_that("marginal fits recover discrepancies and flag execution-only runs", {
  skip_without_fake_rstan()
  saved_cache <- clear_stan_cache()
  on.exit(restore_stan_cache(saved_cache), add = TRUE)
  f <- calibration_fit_fixture()
  e <- f$emulator
  log <- local_fake_rstan()
  dir <- tempfile("stancache")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  one <- suppressMessages(calibration_fit(e, f$z, f$U, iter = 1L, refresh = 0L,
                                          cache_dir = dir))
  expect_identical(log$compiled, "calibration_marginal.stan")
  expect_identical(log$calls[[1]]$pars, c("eta", "tau2"))
  expect_false(log$calls[[1]]$control$adapt_engaged)
  expect_null(log$calls[[1]]$init$u_std)
  expect_true(one$execution_only)
  expect_true(all(is.na(summary(one)$diagnostics$max_rhat)))
  expect_identical(dim(one$mcmc_u), c(2L, 2L, 4L))
  expect_equal(one$mcmc_u, spDBL:::recover_calibration_discrepancy(one$mcmc_eta,
    one$mcmc_tau2, one$tau2_0, one$emulator_draw, e, f$z, f$U, seed = 1234L + 1000L)$u)
  out <- capture_output(print(one))
  expect_match(out, "Calibration execution experiment")
  expect_match(out, "do not establish convergence")
  lean <- suppressMessages(spDBL:::fit_calibration_rstan(e, f$z, f$U, chains = 1,
    iter = 3, warmup = 1, draws_per_emulator = 2, keep_predictions = FALSE,
    keep_discrepancy = FALSE, refresh = 0))
  expect_null(lean$mcmc_u)
  expect_identical(log$compiled, c("calibration_marginal.stan", "calibration_marginal.stan"))
  expect_identical(nrow(lean$mcmc_eta), 4L)
})

test_that("an empty Stan result is reported", {
  skip_without_fake_rstan()
  f <- calibration_fit_fixture()
  local_fake_rstan()
  # Without max_treedepth in control, Stan's default depth of 10 is used.
  ok <- suppressMessages(spDBL:::fit_calibration_rstan(f$emulator, f$z, f$U,
    chains = 1, iter = 2, warmup = 1, draws_per_emulator = 1, compiled_model = "m",
    refresh = 0, control = list()))
  expect_equal(ok$diagnostics[[1]]$max_treedepth_hits, 1)
  local_mocked_bindings(extract = function(object, ...) array(0, c(0, 1, 3)),
                        .package = "rstan")
  expect_error(suppressMessages(spDBL:::fit_calibration_rstan(f$emulator, f$z, f$U,
    chains = 1, iter = 2, warmup = 1, draws_per_emulator = 1, compiled_model = "m",
    refresh = 0)), "Stan returned no calibration samples")
})

test_that("checkpointed conditional fits are reused only with identical settings", {
  skip_without_fake_rstan()
  f <- calibration_fit_fixture()
  e <- f$emulator
  log <- local_fake_rstan()
  checkpoints <- tempfile("checkpoints")
  on.exit(unlink(checkpoints, recursive = TRUE), add = TRUE)
  run <- function(seed = 3L) calibration_fit(e, f$z, f$U, chains = 2L, iter = 6L,
    warmup = 3L, draws_per_emulator = 3L, refresh = 0L, checkpoint_dir = checkpoints,
    compiled_model = "fake-model", seed = seed)
  first <- suppressMessages(run())
  expect_length(log$calls, 2L)
  expect_false(is.null(readRDS(file.path(checkpoints, "conditional_001.rds"))$fingerprint))
  # The same call resumes from the saved fits without sampling again.
  msgs <- capture_messages(again <- run())
  expect_length(log$calls, 2L)
  expect_true(any(grepl("reusing the saved fit conditional_001.rds", msgs)))
  expect_identical(again$mcmc_eta, first$mcmc_eta)
  expect_identical(again$diagnostics, first$diagnostics)
  # Other settings invalidate the saved fits.
  msgs <- capture_messages(other <- run(seed = 4L))
  expect_length(log$calls, 4L)
  expect_true(any(grepl("saved with other data or settings; refitting", msgs)))
  # A checkpoint without a fingerprint (an older version) is refitted too.
  saved <- readRDS(file.path(checkpoints, "conditional_002.rds"))
  saved$fingerprint <- NULL
  saveRDS(saved, file.path(checkpoints, "conditional_002.rds"))
  suppressMessages(run(seed = 4L))
  expect_length(log$calls, 5L)
})

test_that("init = 'mode' starts the chains at the best of several optima", {
  skip_without_fake_rstan()
  f <- calibration_fit_fixture()
  e <- f$emulator
  log <- local_fake_rstan()
  fit <- suppressMessages(calibration_fit(e, f$z, f$U, draws = 1L, chains = 3L,
    iter = 4L, warmup = 2L, draws_per_emulator = 2L, refresh = 0L,
    compiled_model = "fake-model", seed = 7L, init = "mode", mode_starts = 6L))
  # Six optimizations, seeded per start; the failed one is skipped.
  expect_length(log$optimizing, 6L)
  expect_identical(log$optimizing, 7L + 1000L + 1:6)
  search <- fit$diagnostics[[1]]$mode_search
  expect_identical(search$log_density[2], -Inf)
  expect_identical(max(search$log_density), -10)
  expect_equal(search$eta_unit, c(0.3, 0.3))
  # Every chain starts near the better optimum, with its precisions.
  start <- t(vapply(1:3, function(k) log$calls[[1]]$init_fun(chain_id = k)$eta_unit, numeric(2)))
  expect_true(all(abs(start - 0.3) < 0.05))
  expect_gt(min(dist(start)), 0)
  precision <- log$calls[[1]]$init_fun(chain_id = 1)$precision
  expect_true(all(abs(precision / 20 - 1) < 0.25))
  # Reproducible from the seed.
  again <- suppressMessages(calibration_fit(e, f$z, f$U, draws = 1L, chains = 3L,
    iter = 4L, warmup = 2L, draws_per_emulator = 2L, refresh = 0L,
    compiled_model = "fake-model", seed = 7L, init = "mode", mode_starts = 6L))
  expect_identical(t(vapply(1:3, function(k) log$calls[[2]]$init_fun(chain_id = k)$eta_unit,
                            numeric(2))), start)
  expect_error(calibration_fit(e, f$z, f$U, chains = 1L, iter = 4L, compiled_model = "fake-model",
                               mode_starts = 0), "mode_starts")
})
