#!/bin/bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# End-to-end check of the harness on Hoffman2 before the paper jobs: tiny
# versions of E1, E2 and E3 (about 10-15 min with 4 cores), written to
# $SPDBL_WORK/smoke, never to benchmarks/results. Run after setup.sh, from
# the repository root, in the same kind of interactive session.
set -eo pipefail
. benchmarks/hoffman2/env.sh
OUT=$SPDBL_WORK/smoke
rm -rf "$OUT"

Rscript benchmarks/e1_scaling.R --profile=validation --grids=10 --windows=1,2 \
  --out="$OUT/e1"
Rscript benchmarks/e2_accuracy.R --profile=validation --grids=10 --jobs=1 \
  --out="$OUT/e2"
Rscript benchmarks/e3_calibration.R --profile=validation --n_rep=1 --n_draws=2 \
  --chains=2 --iter=40 --warmup=20 --rc_S=60 --rc_S0=20 --rc_chains=2 \
  --init=random --jobs="${NSLOTS:-3}" --out="$OUT/e3"

OUT=$OUT Rscript - <<'EOF'
out <- Sys.getenv("OUT")
ok <- TRUE
for (e in c("e1", "e2", "e3")) {
  tab <- read.csv(file.path(out, e, paste0(e, "_results.csv")))
  bad <- tab[tab$status != "OK", , drop = FALSE]
  cat(sprintf("%s: %d rows, %d not OK, peak RSS %s MB\n", toupper(e), nrow(tab),
    nrow(bad), paste(range(round(tab$peak_rss_mb)), collapse = "-")))
  if (nrow(bad)) print(bad[, intersect(c("task_id", "status", "message"), names(bad))])
  if (nrow(bad) || anyNA(tab$peak_rss_mb)) ok <- FALSE
}
if (!ok) stop("Smoke test failed: see the rows above and the logs in ", out)
cat("Smoke test passed. Submit the paper jobs (see benchmarks/hoffman2/README.md).\n")
EOF
