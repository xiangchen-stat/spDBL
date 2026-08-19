#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# 1. Determine the rendering engine
# Check if a command-line argument was passed (e.g., "r" or "quarto")
ENGINE=$1

# If no argument is provided, prompt the user interactively
if [[ -z "$ENGINE" ]]; then
  echo "No rendering engine specified."
  read -p "Do you want to render with [r]markdown or [q]uarto? (r/q): " CHOICE
  case "$CHOICE" in 
    r|R ) ENGINE="r" ;;
    q|Q ) ENGINE="quarto" ;;
    * ) echo "Invalid choice. Exiting."; exit 1 ;;
  esac
fi

# List your Rmd files here
RMD_FILES=(
  "emulation.Rmd"
)

for RMD_FILE in "${RMD_FILES[@]}"; do
  # Extract base name without .Rmd extension
  BASE_NAME="${RMD_FILE%.Rmd}"
  R_OUTPUT="${BASE_NAME}.R"
  HTML_OUTPUT="${BASE_NAME}.html"

  echo "=================================================="
  echo "Processing: $RMD_FILE"
  echo "=================================================="

  # 1. Extract pure R code from the Rmd file[cite: 9]
  echo "[1/2] Extracting pure R code via knitr::purl()..."
  Rscript -e "knitr::purl(input = '$RMD_FILE', output = '$R_OUTPUT', documentation = 1)"

  # 2. Render Rmd file into an HTML report based on user choice
  if [[ "$ENGINE" == "r" || "$ENGINE" == "R" ]]; then
    echo "[2/2] Rendering HTML document via rmarkdown::render()..."
    Rscript -e "rmarkdown::render(input = '$RMD_FILE', output_format = 'html_document', output_file = '$HTML_OUTPUT')"
  elif [[ "$ENGINE" == "quarto" || "$ENGINE" == "q" ]]; then
    echo "[2/2] Rendering HTML document via quarto render..."
    quarto render "$RMD_FILE" --to html --embed-resources
  else
    echo "Unknown engine: $ENGINE. Please pass 'r' or 'quarto'."
    exit 1
  fi

  echo "=================================================="
  echo "Successfully generated:"
  echo " - HTML Report: $HTML_OUTPUT"
  echo " - Pure R Code: $R_OUTPUT"
  echo "=================================================="
  echo ""
done
