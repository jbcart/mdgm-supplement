## Install dependencies for mdgm simulation studies

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("jbcart/mdgm")

# Gaussian simulation study dependencies
install.packages("mclust")      # adjustedRandIndex
install.packages("bayesImageS") # smcPotts, getNeighbors, getBlocks
