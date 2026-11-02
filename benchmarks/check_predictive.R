# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Self-check of the analytic spDBL predictive variance used in E2
# (R/spdbl_predictive.R) against Monte Carlo built only from the package's
# own samplers: coefficient/covariance draws from emulator_learn(nsam = M),
# conditional means from emulator_predict(MC = TRUE), plus the input-GP
# conditional variance q * Sigma_m[c, c]. Uses the bundled dt_emulation data.
#
#   Rscript benchmarks/check_predictive.R [--nsam=400]

local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  options(spdbl_bench_dir = dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
})
source(file.path(getOption("spdbl_bench_dir"), "R", "common.R"))
source(file.path(getOption("spdbl_bench_dir"), "R", "spdbl_predictive.R"))
args <- bench_parse_args()
bench_setup_lib(if (!is.null(args$lib)) args$lib else Sys.getenv("SPDBL_BENCH_LIB"))
suppressPackageStartupMessages(library(spDBL))
nsam <- if (is.null(args$nsam)) 400L else as.integer(args$nsam)
out_dir <- file.path(bench_dir(), "results", "checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data("dt_emulation", package = "spDBL")
d <- dt_emulation
X_new <- d$pde_para_test[1:3, , drop = FALSE]
test <- lapply(d$dt_pde_test, function(y) y[1:3, , drop = FALSE])
rows <- list()
for (cov in c("iw", "sigma2R")) {
  em <- emulator_learn(d$pde_para_train, d$dt_pde_train, 10, 10, nsam = nsam,
                       episode_window = 2, covariance = cov, seed = 1L)
  F_new <- spdbl_known_designs(em, test, 1L)
  analytic <- spdbl_known_predictive(em, X_new, F_new)
  mc <- emulator_predict(em, X_new, F_new_ls = F_new, MC = TRUE)
  s <- em$setup
  g <- .input_gp_terms(s$pde_para_train, X_new, s$phi_para, s$gp_sigma2,
                       s$gp_tau2, chol(s$V))
  # Diagonal of Sigma_m per draw, mapped to the original grid.
  sig <- em$res_ffbs$Sigma
  diag_draw <- if (cov == "iw") apply(sig, 3, diag) else
    outer(diag(s$R_episode), as.numeric(sig))
  diag_draw <- diag_draw[rep(seq_len(s$bncol), s$n_b), , drop = FALSE]
  diag_draw <- diag_draw[order(as.vector(t(s$block_cols))), , drop = FALSE]
  T <- length(mc)
  mc_var <- array(NA_real_, c(nrow(X_new), s$N_sp, T))
  mc_mean <- mc_var
  for (t in seq_len(T)) {
    a <- mc[[t]]
    mc_mean[, , t] <- apply(a, c(1, 2), mean)
    within <- outer(g$q, rowMeans(diag_draw))
    mc_var[, , t] <- within + apply(a, c(1, 2), function(v) mean((v - mean(v))^2))
  }
  keep <- 3:T
  ratio <- (analytic$sd[, , keep]^2) / mc_var[, , keep]
  rows[[cov]] <- data.frame(covariance = cov, nsam = nsam,
    mean_abs_diff = mean(abs(analytic$mean[, , keep] - mc_mean[, , keep])),
    var_ratio_median = stats::median(ratio),
    var_ratio_q05 = unname(stats::quantile(ratio, 0.05)),
    var_ratio_q95 = unname(stats::quantile(ratio, 0.95)))
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
utils::write.csv(tab, file.path(out_dir, "predictive_check.csv"), row.names = FALSE)
