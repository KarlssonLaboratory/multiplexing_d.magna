#!/usr/bin/env Rscript

#---- Install CRAN packages ---------------------------------------------------
pkgs_CRAN <- c(
  "tidyverse",
  "data.table",
  "readxl",
  "rstatix",
  "stringr",
  "DescTools",
  "patchwork",
  "DHARMa",
  "car",
  "emmeans",
  "gtools",
  "scales",
  "lme4"
)

install.packages(
  pkgs_CRAN,
  repos = "https://cloud.r-project.org",
  Ncpus = parallel::detectCores()
)

#---- Verify everything installed ---------------------------------------------
missing <- pkgs_CRAN[!pkgs_CRAN %in% installed.packages()[, "Package"]]
if (length(missing)) {
  stop("Failed to install: ", paste(missing, collapse = ", "))
}

cat("All packages installed successfully.\n")