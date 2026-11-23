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
# Times and memory: design 1 rerun after the 1.2.0 speed-ups (same accuracy).
e2_timing <- rd("e2", "fixed_timing", "e2_results.csv")
# E3: 100 Latin-hypercube training runs (main text) and the 40 bundled runs.
e3_100 <- rd("e3", "n100", "e3_results.csv")
e3_40 <- rd("e3", "mode", "e3_results.csv")
# Tuned settings. E2: each method's setting selected on designs 11-12 by the
# lowest mean CRPS (results/e2/gasp_tune*, tune_sp_pca*), run on designs 1-10
# (fixed_tuned: spDBL, PCA-GP; fixed_gasp: RobustGaSP) and timed on design 1
# one task per node (fixed_tuned_timing). E3: RobustCalibration with 20 MLE
# starts, pilot-tuned proposals and an emulator nugget (n100_rc_tuned), against
# the same spDBL runs.
pick <- function(x, cols) x[, intersect(cols, names(x)), drop = FALSE]
e2_cols <- c("grid", "design", "method", "rmse", "crps", "coverage95", "width95",
             "fit_s", "predict_s", "peak_rss_mb")
e2_tuned <- rbind(pick(rd("e2", "fixed_tuned", "e2_results.csv"), e2_cols),
                  pick(rd("e2", "fixed_gasp", "e2_results.csv"), e2_cols))
e2_tuned_timing <- if (file.exists(file.path(res, "e2", "fixed_tuned_timing", "e2_results.csv")))
  rd("e2", "fixed_tuned_timing", "e2_results.csv") else NULL
e2_tune <- do.call(rbind, lapply(c("gasp_tune", "gasp_tune2", "tune_sp_pca", "tune_sp_pca2"),
  function(p) cbind(pick(rd("e2", p, "e2_results.csv"), e2_cols[1:6]),
                    round = if (grepl("2$", p)) 2L else 1L)))
e3_cols <- c("rep", "method", "parameter", "covered", "rel_error", "rel_width", "rhat_max",
             "ess_bulk_min", "sampling_cpu_s", "ess_per_cpu_s", "accept_rate_theta")
e3_tuned <- rbind(pick(e3_100[e3_100$method == "spdbl", ], e3_cols),
                  pick(rd("e3", "n100_rc_tuned", "e3_results.csv"), e3_cols))
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
lab2 <- c(sp_opt = "\\pkg{spDBL}, tuned$^{d}$", ppgasp_opt = "\\code{ppgasp}, tuned$^{d}$",
  pc_opt = "PCA-GP, tuned$^{d}$",
  spdbl_s2r_designgp = "\\pkg{spDBL}, design GP", designgp_only = "design GP alone$^{a}$",
  spdbl_s2r_known = "\\pkg{spDBL}, known AR rows$^{b}$", ppgasp_default = "\\code{ppgasp}, default",
  ppgasp_tuned = "\\code{ppgasp}, nugget, Nelder--Mead", ppgasp_per_time = "\\code{ppgasp} per time step",
  mlegp_pca = "PCA-GP (\\pkg{mlegp})",
  spdbl_s2r_designgp_lazy = "\\pkg{spDBL}, design GP$^{c}$")
tuned2 <- c("sp_opt", "ppgasp_opt", "pc_opt")
add("\\begin{table}[ht]", "\\centering\\scriptsize\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{@{}llrrrrrrr@{}}", "\\toprule",
    "Grid & Method & RMSE & CRPS & Cov. & Width & Fit (s) & Predict (s) & Peak (GiB) \\\\", "\\midrule")
for (g in sort(unique(e2$grid))) {
  first <- TRUE
  for (m in names(lab2)) {
    # Regenerated draws give the same predictions: accuracy of stored draws.
    if (m %in% tuned2) {
      x <- e2_tuned[e2_tuned$grid == g & e2_tuned$method == m, ]
      sm <- as.data.frame(lapply(x[c("rmse", "crps", "coverage95", "width95")], mean))
      t1 <- if (is.null(e2_tuned_timing)) data.frame(fit_s = NA, predict_s = NA, peak_rss_mb = NA) else
        e2_tuned_timing[e2_tuned_timing$grid == g & e2_tuned_timing$method == m, ]
    } else {
      sm <- e2_summary[e2_summary$grid == g & e2_summary$method ==
                       sub("_lazy$", "", m), ]
      t1 <- e2_timing[e2_timing$grid == g & e2_timing$method == m, ]
    }
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
    "(\\code{emulator\\_learn(keep\\_draws = FALSE)}); its predictions are identical.",
    "$^{d}$Each method's setting selected on designs 11--12 (Table~\\ref{tab:e2-tuning}): \\pkg{spDBL} with episode",
    "window 5, \\code{gp\\_tune = 1} and 50 draws; \\code{ppgasp} fitted per time step with a linear trend and",
    "the Mat\\'ern 3/2 kernel; the PCA-GP with a linear mean, 99.99\\% of variance, a nugget and 20",
    "restarts; times and memory from design 1, one task per node; the tuned \\code{ppgasp}'s fit time includes prediction.}",
    "\\label{tab:e2}", "\\end{table}", "")

# Table: E2 tuned settings, each competitor minus tuned spDBL, paired by
# design, at all evaluated times t >= 3.
paired_tuned <- function(g, m, k, scale = 1, d = 3) {
  a <- e2_tuned[e2_tuned$grid == g & e2_tuned$method == m, c("design", k)]
  b <- e2_tuned[e2_tuned$grid == g & e2_tuned$method == "sp_opt", c("design", k)]
  x <- merge(a, b, by = "design"); dlt <- scale * (x[[2]] - x[[3]])
  ci <- mean(dlt) + c(-1, 1) * stats::qt(0.975, length(dlt) - 1) * stats::sd(dlt) / sqrt(length(dlt))
  sprintf("%s [%s, %s]", f(mean(dlt), d), f(ci[1], d), f(ci[2], d))
}
add("\\begin{table}[ht]", "\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{@{}llrrr@{}}", "\\toprule",
    "Grid & Method & $\\Delta$ RMSE & $\\Delta$ CRPS & $\\Delta$ coverage (points) \\\\", "\\midrule")
for (g in sort(unique(e2_tuned$grid))) {
  for (m in c("ppgasp_opt", "pc_opt"))
    add(sprintf("%s & %s & %s & %s & %s \\\\", if (m == "ppgasp_opt") sprintf("$%d^2$", g) else "",
      c(ppgasp_opt = "\\code{ppgasp}, tuned", pc_opt = "PCA-GP, tuned")[[m]],
      paired_tuned(g, m, "rmse"), paired_tuned(g, m, "crps"), paired_tuned(g, m, "coverage95", 100, 1)))
  if (g != max(e2_tuned$grid)) add("\\addlinespace")
}
add("\\bottomrule", "\\end{tabular}",
    "\\caption{E2 with every method's tuned setting (Table~\\ref{tab:e2}, note d): each competitor minus",
    "tuned \\pkg{spDBL}, paired by design, mean and 95\\% $t$ interval over the 10 designs, at times $t \\ge 3$.",
    "Positive RMSE and CRPS differences favour \\pkg{spDBL}; negative coverage differences mean the",
    "competitor's 95\\% intervals cover less.}",
    "\\label{tab:e2-paired}", "\\end{table}", "")

# Table: E2 tuning study on designs 11-12, every setting tried.
spec_lab <- c(
  sp_base = "E2 setting: window 2, \\code{gp\\_tune} 0.5, 20 draws", sp_ar1 = "AR(1) design",
  sp_w1 = "episode window 1", sp_w5 = "episode window 5", sp_tune025 = "\\code{gp\\_tune} 0.25",
  sp_tune1 = "\\code{gp\\_tune} 1", sp_nugget = "\\code{gp\\_tau2} $10^{-2}$",
  sp_sigma2 = "\\code{gp\\_sigma2} 4", sp_sqexp = "squared-exponential design GP",
  sp_shared = "one design GP for all times", sp_nsam50 = "50 draws",
  sp_w5_tune1 = "window 5, \\code{gp\\_tune} 1", sp_w5_nsam50 = "window 5, 50 draws",
  sp_tune1_nsam50 = "\\code{gp\\_tune} 1, 50 draws",
  sp_w5_tune1_nsam50 = "window 5, \\code{gp\\_tune} 1, 50 draws",
  gasp_base = "nugget, Nelder--Mead (E2's tuned)", gasp_starts = "10 optimiser starts",
  gasp_trend = "linear trend", gasp_m32 = "Mat\\'ern 3/2", gasp_powexp = "power exponential",
  gasp_mmle = "\\code{mmle}", gasp_refgamma = "\\code{ref\\_gamma} prior", gasp_pt = "per time step",
  gasp_pt_trend = "per time step, trend", gasp_pt_m32 = "per time step, Mat\\'ern 3/2",
  gasp_pt_powexp = "per time step, power exp.", gasp_pt_trend_m32 = "per time step, trend, Mat\\'ern 3/2",
  gasp_pt_trend_powexp = "per time step, trend, power exp.", gasp_trend_m32 = "trend, Mat\\'ern 3/2",
  pc_base = "E2 setting: 99.9\\% of variance", pc_99 = "99\\% of variance", pc_9999 = "99.99\\% of variance",
  pc_linear = "linear mean", pc_nugget = "nugget", pc_tries = "20 restarts",
  pc_linear_9999 = "linear mean, 99.99\\%", pc_linear_nugget = "linear mean, nugget",
  pc_9999_nugget = "99.99\\%, nugget", pc_linear_9999_nugget = "linear mean, 99.99\\%, nugget",
  pc_all = "linear mean, 99.99\\%, nugget, 20 restarts")
tune_sum <- aggregate(cbind(crps, rmse, coverage95) ~ method + round, e2_tune, mean)
add("\\begin{table}[ht]", "\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{@{}llrrrr@{}}", "\\toprule",
    "Method & Setting & Round & CRPS & RMSE & Cov. \\\\", "\\midrule")
for (pre in c("sp_", "gasp_", "pc_")) {
  x <- tune_sum[startsWith(tune_sum$method, pre), ]
  x <- x[order(x$crps, x$rmse), ]
  best <- x$method[1]
  for (i in seq_len(nrow(x))) {
    lab <- spec_lab[[x$method[i]]]
    if (x$method[i] == best) lab <- paste0("\\textbf{", lab, "}")
    add(sprintf("%s & %s & %d & %s & %s & %s \\\\", if (i == 1)
      c(sp_ = "\\pkg{spDBL}", gasp_ = "\\code{ppgasp}", pc_ = "PCA-GP")[[pre]] else "", lab,
      x$round[i], f(x$crps[i], 4), f(x$rmse[i], 3), f(x$coverage95[i], 3)))
  }
  if (pre != "pc_") add("\\addlinespace")
}
add("\\bottomrule", "\\end{tabular}",
    "\\caption{E2 tuning study: every setting tried for each method on designs 11 and 12, which E2 does not",
    "evaluate, all three grids, means over the six grid--design pairs. Round 1 changes one setting of the",
    "method's E2 setting (first row of round 1); round 2 combines the round-1 changes that lowered CRPS.",
    "The setting with the lowest mean CRPS (bold; ties broken by RMSE), a rule fixed before any result,",
    "is the tuned setting of Tables~\\ref{tab:e2} and~\\ref{tab:e2-paired}. Settings not listed are the package defaults.}",
    "\\label{tab:e2-tuning}", "\\end{table}", "")

# E3 per parameter, replicate-bootstrap 95% intervals: 100 training runs with
# the tuned competitors, and the 40 bundled runs with their default settings.
# The pilot of the proposal scale and the untuned runs (results/e3/rc_pilot,
# n100, e2/fixed) are described in the paper's text, not tabulated.
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
e3_setting_tuned <- "every chain started at its package's optimum from 20 starts (\\pkg{spDBL}: conditional posterior mode from 20 L-BFGS starts, 20 emulator draws with 4 NUTS chains of 2{,}000 iterations each; \\pkg{RobustCalibration}: maximum-likelihood estimate from 20 starts, an emulator nugget, random-walk proposal scales set from a 4{,}000-iteration pilot chain, and 4 chains of 20{,}000 iterations); errors and widths relative to the prior range, with replicate-bootstrap 95\\% intervals, paired over replicates for the differences."
cpu_100 <- "one Slurm job per calibration (\\pkg{spDBL} on one core; \\pkg{RobustCalibration} running its four chains on four cores, CPU time summed over them and including the mode search and pilot)"
e3_table(e3_tuned, "tab:e3",
  paste("E3 with 100 Latin-hypercube training runs (one design for every replicate and method): 30 replicate calibrations of the three free inputs $\\eta = (\\eta_1, \\eta_2, \\alpha_2)$ by \\pkg{spDBL} and by \\pkg{RobustCalibration} with its S-GaSP and GaSP discrepancies,", e3_setting_tuned),
  cpu_100)
e3_table(e3_40, "tab:e3-40",
  paste("E3 with the 40 training runs bundled with the package (\\code{dt\\_calibration}), with \\pkg{RobustCalibration}'s default proposal scales, no emulator nugget and 10 maximum-likelihood starts:", e3_setting),
  "with 18 one-core tasks per 36-core node")

# E3 ablations on replicates 1-10, against the 40-run (mode) run on
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
