# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Pooled scores and paired differences across replicates for the prediction
# studies (e3_sites.R, e3_offgrid.R), and for off-grid points the coverage of
# points within half a coarse cell of the domain's edge against the interior.
# Writes <study>_pooled.csv, <study>_paired.csv (and offgrid_edge.csv) to the
# results folder.
#   Rscript benchmarks/paper/prediction_summary.R benchmarks/results/e3/sites sites
#   Rscript benchmarks/paper/prediction_summary.R benchmarks/results/e3/offgrid offgrid
args <- commandArgs(TRUE)
study <- args[2]
tab <- utils::read.csv(file.path(args[1], sprintf("e3_%s_results.csv", study)))
ref <- if (study == "sites") "spdbl" else "spdbl_bilinear"
cat(sprintf("%d of %d tasks OK\n", length(unique(tab$rep[tab$status == "OK"])), length(unique(tab$rep))))
tab <- tab[tab$status == "OK", ]
cat("Settings:", unique(tab$settings), "\n\n")
pooled <- do.call(rbind, lapply(split(tab, list(tab$target, tab$predictor), drop = TRUE), function(x)
  data.frame(target = x$target[1], predictor = x$predictor[1], reps = nrow(x),
             rmse = sqrt(mean(x$rmse^2)), crps = mean(x$crps), bias = mean(x$bias),
             coverage95 = mean(x$coverage95), width95 = mean(x$width95))))
pooled <- pooled[order(pooled$target, pooled$crps), ]
cat("Pooled over replicates (coverage: share of truths inside the 95% interval)\n")
print(pooled, row.names = FALSE, digits = 3)
utils::write.csv(pooled, file.path(args[1], sprintf("%s_pooled.csv", study)), row.names = FALSE)
cat(sprintf("\nPaired differences, %s minus each other predictor (mean over replicates, 95%% t-interval; negative = %s better)\n", ref, ref))
rows <- list()
for (tg in unique(tab$target)) {
  a <- tab[tab$target == tg & tab$predictor == ref, ]
  for (p in setdiff(unique(tab$predictor[tab$target == tg]), ref)) {
    b <- tab[tab$target == tg & tab$predictor == p, ]
    m <- merge(a, b, by = "rep")
    for (metric in c("rmse", "crps")) {
      d <- m[[paste0(metric, ".x")]] - m[[paste0(metric, ".y")]]
      ci <- if (length(d) > 1) stats::t.test(d)$conf.int else c(NA, NA)
      rows[[length(rows) + 1]] <- data.frame(target = tg, versus = p, metric = metric,
        mean_diff = mean(d), lower = ci[1], upper = ci[2], reps_better = sum(d < 0), reps = length(d))
    }
  }
}
paired <- do.call(rbind, rows)
print(paired, row.names = FALSE, digits = 3)
utils::write.csv(paired, file.path(args[1], sprintf("%s_paired.csv", study)), row.names = FALSE)
cat(sprintf("\nMean time per replicate: fit %.0f s, recovery %.1f s\n", mean(tab$fit_s), mean(tab$predict_s)))

if (study == "offgrid") {
  rows <- do.call(rbind, lapply(list.files(file.path(args[1], "raw"), "^offgrid_r.*rds$", full.names = TRUE), function(f) {
    dt <- readRDS(f)$detail
    edge <- apply(dt$coords, 1, function(p) any(p < 1 | p > 6))   # within half a coarse cell of the boundary
    do.call(rbind, lapply(c("bilinear_f", "cell_f"), function(b) {
      band <- dt[[b]]
      covered <- dt$truth_f >= band[1, , ] & dt$truth_f <= band[3, , ]
      err <- band[2, , ] - dt$truth_f
      data.frame(field = sub("_f", "", b), where = ifelse(rep(edge, each = nrow(err)), "edge", "interior"),
                 covered = as.vector(covered), err = as.vector(err))
    }))
  }))
  agg <- aggregate(cbind(covered, sq = err^2, err) ~ field + where, rows, mean)
  agg$rmse <- sqrt(agg$sq); agg$n <- aggregate(covered ~ field + where, rows, length)$covered
  agg <- agg[, c("field", "where", "n", "covered", "rmse", "err")]
  names(agg) <- c("field", "where", "n", "coverage95", "rmse", "mean_error")
  cat("\nOff-grid y + u by position (edge: within half a coarse cell of the boundary)\n")
  print(agg, digits = 3, row.names = FALSE)
  utils::write.csv(agg, file.path(args[1], "offgrid_edge.csv"), row.names = FALSE)
}
