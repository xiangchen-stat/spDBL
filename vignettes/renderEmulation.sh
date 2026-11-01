#!/usr/bin/env bash
set -euo pipefail

# Usage: bash vignettes/renderEmulation.sh quarto [calibration ...]
# Without vignette names, render all three. Run from any working directory.
ENGINE=${1:-}
if [[ -z "$ENGINE" ]]; then
  read -r -p "Render with [r]markdown or [q]uarto? (r/q): " ENGINE
else
  shift
fi
case "$ENGINE" in
  r|R) ENGINE=r ;;
  q|Q|quarto) ENGINE=quarto ;;
  *) echo "Unknown engine: $ENGINE. Use 'r' or 'quarto'." >&2; exit 1 ;;
esac

VIGNETTES=(emulation calibration calibration-details)
if [[ $# -gt 0 ]]; then
  VIGNETTES=("$@")
fi
for NAME in "${VIGNETTES[@]}"; do
  case "$NAME" in
    emulation|calibration|calibration-details) ;;
    *) echo "Unknown vignette: $NAME" >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
LOCAL_LIBRARY="$PACKAGE_DIR/.r-library"
mkdir -p "$LOCAL_LIBRARY"
export R_LIBS="$LOCAL_LIBRARY${R_LIBS:+:$R_LIBS}"

echo "Installing this spDBL checkout into $LOCAL_LIBRARY"
R CMD INSTALL --library="$LOCAL_LIBRARY" "$PACKAGE_DIR"
cd "$PACKAGE_DIR"
Rscript data-raw/sync-vignettes.R
cd "$SCRIPT_DIR"

for NAME in "${VIGNETTES[@]}"; do
  echo "Rendering $NAME with $ENGINE"
  Rscript -e 'args <- commandArgs(TRUE); knitr::purl(args[1], output = args[2], documentation = 1)' \
    "$NAME.Rmd" "$NAME.R"
  if [[ "$ENGINE" == r ]]; then
    Rscript -e 'args <- commandArgs(TRUE); rmarkdown::render(args[1], output_format = "html_document", output_file = args[2])' \
      "$NAME.Rmd" "$NAME.html"
  else
    quarto render "$NAME.qmd" --to html --embed-resources
  fi
done
