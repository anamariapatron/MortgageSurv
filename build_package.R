#!/usr/bin/env Rscript
# =============================================================================
#  build_package.R
#  Run this script from the PARENT directory of MortgageSurv/ to:
#    1. Generate man/ pages with roxygen2
#    2. Check the package
#    3. Build a .tar.gz ready for installation
#
#  Usage (from terminal):
#    Rscript build_package.R
#
#  Or from the R console:
#    source("build_package.R")
# =============================================================================

# Install build tools if needed
required <- c("devtools", "roxygen2")
missing  <- required[!required %in% rownames(installed.packages())]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

pkg_path <- "MortgageSurv"

# ── 1. Generate documentation ─────────────────────────────────────────────────
message("\n── Generating documentation with roxygen2 ──")
devtools::document(pkg_path)

# ── 2. Check package ──────────────────────────────────────────────────────────
message("\n── Running R CMD check ──")
devtools::check(pkg_path, cran = FALSE)

# ── 3. Build .tar.gz ──────────────────────────────────────────────────────────
message("\n── Building package tarball ──")
built <- devtools::build(pkg_path)
message("Package built: ", built)

# ── 4. Install locally (optional) ─────────────────────────────────────────────
message("\n── Installing package ──")
devtools::install(pkg_path)
message("\nDone!  Load with:  library(MortgageSurv)")
