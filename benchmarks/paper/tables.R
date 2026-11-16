# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Appendix tables of the MLOSS paper, generated from the committed benchmark
# results so that no number is typed by hand.
# Run from the paper folder (it writes tables/appendix-tables.tex there, which
# spdbl-mloss.tex inputs):
#   Rscript <spDBL>/benchmarks/paper/tables.R [results directory]
args <- commandArgs(TRUE)
# The results directory defaults to benchmarks/results next to this script.
script_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else "."
})
res <- if (length(args)) args[1] else file.path(script_dir, "..", "results")
rd <- function(...) { x <- utils::read.csv(file.path(res, ...)); if ("status" %in% names(x)) x[x$status == "OK", ] else x }
e1 <- rd("e1", "paper", "e1_results.csv")
e1_lazy <- rd("e1", "paper_lazy", "e1_results.csv")
e2 <- rd("e2", "fixed", "e2_results.csv")
e2_summary <- rd("e2", "fixed", "e2_summary.csv")
e2_paired <- rd("e2", "fixed", "e2_paired.csv")
# Times and memory: design 1 rerun after the 1.2.0 speed-ups (same accuracy).
e2_timing <- rd("e2", "fixed_timing", "e2_results.csv")
# E3: 100 Latin-hypercube training runs (main text) and the 40 bundled runs.
e3_100 <- rd("e3", "n100", "e3_results.csv")
e3_40 <- rd("e3", "mode", "e3_results.csv")
gib <- function(mb) mb / 1024
# round() then + 0 so that tiny negatives print as 0, not -0; negatives get a
# true minus sign.
f <- function(x, d = 2) ifelse(is.na(x), "--",
  sub("^-", "$-$", formatC(round(x, d) + 0, format = "f", digits = d)))
sig <- function(x) ifelse(is.na(x), "--", ifelse(x >= 100, formatC(x, format = "f", digits = 0),
  ifelse(x >= 10, formatC(x, format = "f", digits = 1), formatC(x, format = "f", digits = 2))))
out <- character()
add <- function(...) out <<- c(out, ...)

# Table A1: E1 scaling
add("\\begin{table}[ht]", "\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{@{}rr*{3}{rrr}@{}}", "\\toprule",
    " & & \\multicolumn{3}{c}{$w = 1$} & \\multicolumn{3}{c}{$w = 2$} & \\multicolumn{3}{c}{$w = 5$} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    "Grid & $S$ & fit & predict & peak & fit & predict & peak & fit & predict & peak \\\\",
    "\\midrule")
for (S in sort(unique(e1$S))) {
  cells <- unlist(lapply(c(1, 2, 5), function(w) {
    r <- e1[e1$S == S & e1$episode_window == w, ]
    l <- e1_lazy[e1_lazy$S == S & e1_lazy$episode_window == w, ]
    c(sig(r$fit_s), sig(r$predict_s),
      sprintf("%s/%s", f(gib(r$peak_rss_mb), 1), f(gib(l$peak_rss_mb), 1)))
  }))
  add(sprintf("$%d^2$ & %s & %s \\\\", sqrt(S), format(S, big.mark = ","), paste(cells, collapse = " & ")))
}
add("\\bottomrule", "\\end{tabular}",
    "\\caption{E1: inverse-Wishart emulator with 10 posterior draws, 40 training runs and 26 time",
    "steps; fit and prediction (10 held-out runs) time in seconds and peak memory of the",
    "whole task in GiB, one task per node, with stored coefficient draws and, after the slash,",
    "with draws regenerated on demand (\\code{emulator\\_learn(keep\\_draws = FALSE)}; same",
    "accuracy). Times are with stored draws.}", "\\label{tab:e1}", "\\end{table}", "")

# Table A2: E2 (fixed domain) accuracy over 10 designs, cost from design 1
lab2 <- c(spdbl_s2r_designgp = "\\pkg{spDBL}, design GP", designgp_only = "design GP alone$^{a}$",
  spdbl_s2r_known = "\\pkg{spDBL}, known AR rows$^{b}$", ppgasp_default = "\\code{ppgasp}, default",
  ppgasp_tuned = "\\code{ppgasp}, tuned", ppgasp_per_time = "\\code{ppgasp} per time step",
  mlegp_pca = "PCA-GP (\\pkg{mlegp})",
  spdbl_s2r_designgp_lazy = "\\pkg{spDBL}, design GP$^{c}$")
add("\\begin{table}[ht]", "\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{@{}llrrrrrrr@{}}", "\\toprule",
    "Grid & Method & RMSE & CRPS & Cov. & Width & Fit (s) & Predict (s) & Peak (GiB) \\\\", "\\midrule")
for (g in sort(unique(e2$grid))) {
  first <- TRUE
  for (m in names(lab2)) {
    # Regenerated draws give the same predictions: accuracy of stored draws.
    sm <- e2_summary[e2_summary$grid == g & e2_summary$method ==
                     sub("_lazy$", "", m), ]
    t1 <- e2_timing[e2_timing$grid == g & e2_timing$method == m, ]
    if (!nrow(sm)) next
    add(sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", if (first) sprintf("$%d^2$", g) else "",
      lab2[[m]], f(sm$rmse, 3), f(sm$crps, 3), f(sm$coverage95, 3), f(sm$width95, 2),
      sig(t1$fit_s), sig(t1$predict_s), f(gib(t1$peak_rss_mb), 2)))
    first <- FALSE
  }
  if (g != max(e2$grid)) add("\\addlinespace")
}
add("\\bottomrule", "\\end{tabular}",
    "\\caption{E2 on the fixed $10\\times 10$ domain: held-out prediction on 50 runs at times",
    "\\mbox{$t \\ge 3$}, means over 10 independent training/test designs; Cov.\\ and Width are the coverage",
    "and mean width of central 95\\% intervals (\\pkg{spDBL}'s from \\code{predict(..., interval",
    "= \"prediction\")}). Fit and prediction time and peak memory are from design 1, rerun with",
    "spDBL 1.2.0's batched prediction, one task per node; \\code{ppgasp} per time step fits one \\code{ppgasp} per time, and its fit time includes",
    "prediction. $^{a}$Times $3 \\le t \\le 25$; it has no prediction for the last time. $^{b}$Uses the",
    "simulator's outputs at the test inputs to build the autoregressive rows (one-step prediction), so",
    "it is not simulator-free. $^{c}$With coefficient draws regenerated on demand",
    "(\\code{emulator\\_learn(keep\\_draws = FALSE)}); its predictions are identical.}",
    "\\label{tab:e2}", "\\end{table}", "")

# Table A3: E2 paired differences from spDBL over the 10 designs, common times
lab3 <- c(designgp_only = "design GP alone", ppgasp_default = "\\code{ppgasp}, default",
  ppgasp_tuned = "\\code{ppgasp}, tuned", ppgasp_per_time = "\\code{ppgasp} per time step",
  mlegp_pca = "PCA-GP (\\pkg{mlegp})")
cell <- function(g, m, k, scale = 1, d = 3) {
  r <- e2_paired[e2_paired$grid == g & e2_paired$method == m & e2_paired$metric == k, ]
  if (!nrow(r)) return("--")
  sprintf("%s [%s, %s]", f(scale * r$difference, d), f(scale * r$lo95, d), f(scale * r$hi95, d))
}
add("\\begin{table}[ht]", "\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{@{}llrrr@{}}", "\\toprule",
    "Grid & Method & $\\Delta$ RMSE & $\\Delta$ CRPS & $\\Delta$ coverage (points) \\\\", "\\midrule")
for (g in sort(unique(e2_paired$grid))) {
  first <- TRUE
  for (m in names(lab3)) {
    add(sprintf("%s & %s & %s & %s & %s \\\\", if (first) sprintf("$%d^2$", g) else "", lab3[[m]],
      cell(g, m, "rmse_common"), cell(g, m, "crps_common"), cell(g, m, "coverage_common", 100, 1)))
    first <- FALSE
  }
  if (g != max(e2_paired$grid)) add("\\addlinespace")
}
add("\\bottomrule", "\\end{tabular}",
    "\\caption{E2 on the fixed domain: each method minus \\pkg{spDBL}'s simulator-free prediction,",
    "paired by design, mean and 95\\% $t$ interval over the 10 designs, at the times $3 \\le t \\le 25$",
    "that every method predicts. Positive RMSE and CRPS differences favour \\pkg{spDBL}; negative",
    "coverage differences mean the method's 95\\% intervals cover less than \\pkg{spDBL}'s, which",
    "cover 94.5--95.8\\%.}",
    "\\label{tab:ablation}", "\\end{table}", "")

# Tables A4 and A5: E3 per parameter, replicate-bootstrap 95% intervals, for
# 100 training runs and for the 40 bundled runs.
labm <- c(spdbl = "\\pkg{spDBL}", rc_sgasp = "S-GaSP", rc_gasp = "GaSP")
labp <- c(eta1 = "$\\eta_1$", eta2 = "$\\eta_2$", alpha2 = "$\\alpha_2$", all = "all")
e3_table <- function(e3, label, caption, cpu_note) {
set.seed(1); B <- 4000
reps <- sort(unique(e3$rep))
idx <- replicate(B, sample(seq_along(reps), replace = TRUE))
# Per-replicate sums for one method and parameter (or all parameters).
e3_sums <- function(m, p) {
  d <- e3[e3$method == m & (p == "all" | e3$parameter == p), ]
  list(d = d, se = tapply(d$rel_error^2, factor(d$rep, reps), sum),
       nn = tapply(d$rel_error, factor(d$rep, reps), length),
       cv = tapply(d$covered, factor(d$rep, reps), sum))
}
boot_rmse <- function(x) sqrt(colSums(matrix(x$se[idx], nrow(idx))) / colSums(matrix(x$nn[idx], nrow(idx))))
boot_cov <- function(x) colSums(matrix(x$cv[idx], nrow(idx))) / colSums(matrix(x$nn[idx], nrow(idx)))
add("", "\\begin{table}[ht]", "\\centering\\scriptsize\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{@{}llrrrrr@{}}", "\\toprule",
    "& & & & & \\multicolumn{2}{c}{Method minus \\pkg{spDBL}} \\\\ \\cmidrule(l){6-7}",
    "Method & Parameter & Relative RMSE & Coverage & Rel.\\ width & Relative RMSE & Coverage \\\\", "\\midrule")
for (m in names(labm)) {
  first <- TRUE
  for (p in names(labp)) {
    x <- e3_sums(m, p); d <- x$d
    rm_b <- boot_rmse(x); cv_b <- boot_cov(x)
    q <- function(v, d = 3) sprintf("[%s, %s]", f(stats::quantile(v, 0.025), d), f(stats::quantile(v, 0.975), d))
    # Paired replicate bootstrap: the same resampled replicates for both methods.
    paired <- if (m == "spdbl") c("--", "--") else {
      ref <- e3_sums("spdbl", p)
      c(sprintf("%s %s", f(sqrt(sum(x$se) / sum(x$nn)) - sqrt(sum(ref$se) / sum(ref$nn)), 3),
                q(rm_b - boot_rmse(ref))),
        sprintf("%s %s", f(sum(x$cv) / sum(x$nn) - sum(ref$cv) / sum(ref$nn), 2),
                q(cv_b - boot_cov(ref), 2)))
    }
    add(sprintf("%s & %s & %s %s & %s %s & %s & %s & %s \\\\", if (first) labm[[m]] else "", labp[[p]],
      f(sqrt(sum(x$se) / sum(x$nn)), 3), q(rm_b), f(sum(x$cv) / sum(x$nn), 2), q(cv_b, 2),
      f(mean(d$rel_width), 2), paired[1], paired[2]))
    first <- FALSE
  }
  if (m != "rc_gasp") add("\\addlinespace")
}
conv <- sapply(names(labm), function(m) { r <- tapply(e3$rhat_max[e3$method == m], e3$rep[e3$method == m], max); c(sum(r > 1.01), sum(r > 1.1)) })
cpu <- sapply(names(labm), function(m) mean(tapply(e3$sampling_cpu_s[e3$method == m], e3$rep[e3$method == m], `[`, 1)) / 3600)
add("\\bottomrule", "\\end{tabular}",
    sprintf("\\caption{%s Replicates with rank-normalised $\\hat R$ above 1.01 and above 1.1 for some $\\eta$ component (for \\pkg{spDBL} the worst of its 20 conditional fits, for \\pkg{RobustCalibration} its four chains): %s. Mean CPU hours per calibration, %s: %s.}",
      caption,
      paste(sprintf("%s %d and %d of 30", c("\\pkg{spDBL}", "S-GaSP", "GaSP"), conv[1, ], conv[2, ]), collapse = ", "),
      cpu_note, paste(sprintf("%s %.1f", c("\\pkg{spDBL}", "S-GaSP", "GaSP"), cpu), collapse = ", ")),
    sprintf("\\label{%s}", label), "\\end{table}")
}
e3_setting <- "every chain started at its package's optimum (\\pkg{spDBL}: conditional posterior mode, \\code{init = \"mode\"}, 20 emulator draws with 4 NUTS chains of 2{,}000 iterations each; \\pkg{RobustCalibration}: its 10-start maximum-likelihood estimate, 4 chains of 20{,}000 iterations); errors and widths relative to the prior range, with replicate-bootstrap 95\\% intervals, paired over replicates for the differences."
e3_table(e3_100, "tab:e3",
  paste("E3 with 100 Latin-hypercube training runs (one design for every replicate and method): 30 replicate calibrations of the three free inputs $\\eta = (\\eta_1, \\eta_2, \\alpha_2)$ by \\pkg{spDBL} and by \\pkg{RobustCalibration} with its S-GaSP and GaSP discrepancies,", e3_setting),
  "one Slurm job per calibration (\\pkg{spDBL} on one core; \\pkg{RobustCalibration} running its four chains on four cores, CPU time summed over them)")
e3_table(e3_40, "tab:e3-40",
  paste("E3 with the 40 training runs bundled with the package (\\code{dt\\_calibration}), otherwise as Table~\\ref{tab:e3}:", e3_setting),
  "with 18 one-core tasks per 36-core node")
# Table A6: E3 ablations on replicates 1-10, against the 40-run (mode) run on
# the same replicates, with paired replicate-bootstrap 95% intervals of the
# differences.
abl <- c(abl_joint = "Joint paths", abl_m60 = "60 draws",
         abl_nodu = "No design uncertainty", abl_self = "Data from emulator model",
         abl_n40 = "New 40-run design", abl_n80 = "80 training runs",
         n100 = "100 training runs",
         abl_n160 = "160 training runs", abl_two = "40 + 40 targeted runs")
abl <- abl[file.exists(file.path(res, "e3", names(abl), "e3_results.csv"))]
if (length(abl)) {
  ab_reps <- 1:10
  ref <- e3_40[e3_40$method == "spdbl" & e3_40$rep %in% ab_reps, ]
  per_rep <- function(x) sapply(ab_reps, function(r) { y <- x[x$rep == r, ]
    c(cv = sum(y$covered), w = sum(y$rel_width), se = sum(y$rel_error^2), n = nrow(y),
      rh = max(y$rhat_max)) })
  totals <- function(p, i) c(cov = sum(p["cv", i]) / sum(p["n", i]),
    width = sum(p["w", i]) / sum(p["n", i]), rmse = sqrt(sum(p["se", i]) / sum(p["n", i])))
  set.seed(2); ab_idx <- replicate(4000, sample(seq_along(ab_reps), replace = TRUE))
  pr <- per_rep(ref); t0 <- totals(pr, seq_along(ab_reps))
  add("", "\\begin{table}[ht]", "\\centering\\scriptsize\\setlength{\\tabcolsep}{3pt}",
      "\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
      "& \\multicolumn{2}{c}{Coverage} & Relative & \\multicolumn{2}{c}{Relative RMSE} & $\\hat R > 1.01$ \\\\",
      "\\cmidrule(lr){2-3}\\cmidrule(lr){5-6}",
      "Variant & value & minus baseline & width & value & minus baseline & replicates \\\\", "\\midrule",
      sprintf("Baseline (Table~\\ref{tab:e3-40}) & %s & -- & %s & %s & -- & %d \\\\", f(t0[["cov"]]), f(t0[["width"]], 3),
              f(t0[["rmse"]], 3), sum(pr["rh", ] > 1.01)))
  for (v in names(abl)) {
    a <- rd("e3", v, "e3_results.csv"); a <- a[a$method == "spdbl" & a$rep %in% ab_reps, ]
    pa <- per_rep(a); t1 <- totals(pa, seq_along(ab_reps))
    d <- sapply(seq_len(ncol(ab_idx)), function(k) totals(pa, ab_idx[, k]) - totals(pr, ab_idx[, k]))
    ci <- function(j, dd) sprintf("%s [%s, %s]", f(t1[[j]] - t0[[j]], dd),
      f(stats::quantile(d[j, ], 0.025), dd), f(stats::quantile(d[j, ], 0.975), dd))
    add(sprintf("%s & %s & %s & %s & %s & %s & %d \\\\", abl[[v]], f(t1[["cov"]]), ci("cov", 2),
      f(t1[["width"]], 3), f(t1[["rmse"]], 3), ci("rmse", 3), sum(pa["rh", ] > 1.01)))
  }
  add("\\bottomrule", "\\end{tabular}",
      "\\caption{E3 ablations: \\pkg{spDBL} on replicates 1--10 with one change each from the setting of Table~\\ref{tab:e3-40}, against that run on the same replicates and seeds, with paired replicate-bootstrap 95\\% intervals of the differences. Joint paths: \\code{calibration\\_prepare(sampling = \"joint\")}. 60 draws: 60 emulator draws instead of 20. No design uncertainty: \\code{design\\_uncertainty = FALSE}. New 40-run design, 80, 100 and 160 training runs: new Latin-hypercube designs over the prior box instead of the 40 bundled runs (100 runs: replicates 1--10 of Table~\\ref{tab:e3}). 40 + 40 targeted runs: the bundled runs plus 40 placed by greedy maximin where a short first-stage calibration puts the posterior, its spread tripled. Data from emulator model: the simulator output in each replicate's field data replaced by a draw from the emulator's own predictive model at the true parameters (one random emulator draw; the same discrepancy and noise), so that \\pkg{spDBL} is exactly specified.}",
      "\\label{tab:ablation-e3}", "\\end{table}")
}
dir.create("tables", showWarnings = FALSE)
writeLines(out, file.path("tables", "appendix-tables.tex"))
message("Wrote tables/appendix-tables.tex")
