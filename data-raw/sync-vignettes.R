# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Run from the package directory after editing the canonical R Markdown files.
# Quarto copies are generated, never edited independently.
for (path in list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)) {
  text <- paste(readLines(path), collapse = "\n")
  if (grepl("\nformat:", text, fixed = TRUE)) {
    text <- sub("output: rmarkdown::html_vignette\n", "", text, fixed = TRUE)
  } else {
    text <- sub("output: rmarkdown::html_vignette\n", "format: html\n", text, fixed = TRUE)
  }
  text <- sub("vignette: >\n(  %[^\n]*\n)+", "", text, perl = TRUE)
  local_setup <- paste(
    "<!-- Generated from the corresponding .Rmd; do not edit. -->",
    "",
    "```{r source-library, include=FALSE}",
    "# The rendering helper installs this checkout into a local R library.",
    "local({",
    "  input_dir <- dirname(normalizePath(knitr::current_input(dir = TRUE)))",
    "  lib <- file.path(input_dir, \"..\", \".r-library\")",
    "  if (file.exists(file.path(lib, \"spDBL\", \"DESCRIPTION\"))) {",
    "    .libPaths(c(lib, .libPaths()))",
    "  }",
    "})",
    "```",
    "", sep = "\n"
  )
  text <- sub("\n---\n", paste0("\n---\n\n", local_setup), text, fixed = TRUE)
  writeLines(text, sub("[.]Rmd$", ".qmd", path))
}
