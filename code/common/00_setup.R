## Install dependencies for mdgm simulation studies

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_if_missing("remotes")

if (!requireNamespace("mdgm", quietly = TRUE)) {
  remotes::install_github("jbcart/mdgm")
}

# Gaussian simulation study dependencies
install_if_missing("mclust")      # adjustedRandIndex
install_if_missing("bayesImageS") # smcPotts, getNeighbors, getBlocks

# Disorder analysis and plotting dependencies
install_if_missing("tidyverse")   # tables_plots.R, disorder/00_analysis.R
install_if_missing("sf")          # disorder/00_analysis.R, disorder/01_cv_single_run.R
install_if_missing("igraph")      # disorder/00_analysis.R
install_if_missing("coda")        # disorder/00_analysis.R
install_if_missing("patchwork")   # disorder/00_analysis.R
