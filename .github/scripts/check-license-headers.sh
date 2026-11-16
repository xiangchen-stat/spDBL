#!/usr/bin/env bash
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Fail if a hand-written source file lacks the SPDX license header.
# Rcpp regenerates R/RcppExports.R and src/RcppExports.cpp, so they are exempt.
# Run from the package root: bash .github/scripts/check-license-headers.sh
set -euo pipefail

missing=0
while IFS= read -r f; do
  case "$f" in
    R/RcppExports.R|src/RcppExports.cpp) continue ;;
  esac
  if ! head -n 5 "$f" | grep -q "SPDX-License-Identifier: MIT"; then
    echo "Missing license header: $f"
    missing=1
  fi
done < <(find R src inst/stan inst/examples tests data-raw benchmarks vignettes \
           -type f \( -name '*.R' -o -name '*.cpp' -o -name '*.h' -o -name '*.stan' -o -name '*.sh' -o -name '*.Rmd' -o -name '*.qmd' \) \
           2>/dev/null | sort)

if [ "$missing" -ne 0 ]; then
  cat <<'MSG'

Add these two lines (with // instead of # in C++ and Stan files) at the top:
# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT
MSG
  exit 1
fi
echo "All source files carry the MIT license header."
