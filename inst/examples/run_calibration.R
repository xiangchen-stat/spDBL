# Run from any directory after installing spDBL:
# Rscript /installed/path/spDBL/examples/run_calibration.R oneiter500 output_dir
# Optional third argument: marginal (default) or joint emulator sampling.
# Profiles: oneiter500 (default), small (six sites), full (all sites), smoke.
library(spDBL)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "oneiter500"
if (!mode %in% c("small", "full", "smoke", "oneiter500")) stop("Unknown profile")
out <- if (length(args) > 1L) args[2] else paste0("spDBL-calibration-", mode)
sampling <- if (length(args) > 2L) match.arg(args[3], c("marginal", "joint")) else "marginal"
if (file.exists(out) && length(list.files(out))) stop("Choose an empty output directory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data("dt_calibration", package = "spDBL")
d <- dt_calibration
config <- switch(mode,
  full = list(n_draws = 20L, chains = 4L, iter = 2000L, warmup = 1000L,
              keep = 200L, sites = seq_len(36L)),
  small = list(n_draws = 2L, chains = 4L, iter = 1000L, warmup = 500L,
               keep = 50L, sites = c(1L, 8L, 15L, 22L, 29L, 36L)),
  smoke = list(n_draws = 2L, chains = 1L, iter = 40L, warmup = 20L,
               keep = 4L, sites = c(1L, 8L, 15L, 22L, 29L, 36L)),
  oneiter500 = list(n_draws = 500L, chains = 1L, iter = 1L, warmup = 0L,
                   keep = 1L, sites = seq_len(36L)))
config$emulator_sampling <- sampling
saveRDS(config, file.path(out, "settings.rds"))
em <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
  covariance = "sigma2R", episode_window = 2L, nsam = config$n_draws,
  seed = 1234L, sampling = sampling)
cal <- calibration_prepare(em, d$lower, d$upper)
validation <- calibration_validate(cal, d$validation$input, d$validation$responses)
saveRDS(validation, file.path(out, "validation.rds"))
control <- if (config$warmup == 0L)
  list(adapt_engaged = FALSE, stepsize = 0.001, max_treedepth = 12L) else
  list(adapt_delta = 0.95, max_treedepth = 12L)
fit <- calibration_fit(cal, d$field$z[, config$sites, drop = FALSE], d$field$U,
  locations = config$sites, chains = config$chains, iter = config$iter,
  warmup = config$warmup, draws_per_emulator = config$keep,
  seed = 1237L, control = control, refresh = 0L,
  checkpoint_dir = file.path(out, "conditionals"),
  cache_dir = file.path(out, "stan-cache"))
fit$execution_only <- mode %in% c("smoke", "oneiter500")
saveRDS(fit, file.path(out, "calibration.rds"), compress = "gzip")
print(summary(fit))
pred <- calibration_predict(fit)
saveRDS(pred, file.path(out, "predictions.rds"), compress = "gzip")
png(file.path(out, "parameters.png"), width = 1200, height = 450)
plot(fit)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Saved calibration outputs in", normalizePath(out), "\n")
