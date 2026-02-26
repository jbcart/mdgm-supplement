# mdgm-supplement

Supplementary code for reproducing simulation studies in:

> Carter, J. B. and Calder, C. A. (2024). "Mixture of Directed Graphical Models for Discrete Spatial Random Fields."
> [arXiv:2406.15700](https://arxiv.org/abs/2406.15700)

## Requirements

- R >= 4.1
- [mdgm](https://github.com/jbcart/mdgm) R package (installed from GitHub)

## Quick Start

```r
# Install mdgm
source("code/00_setup.R")

# Run a single scenario (quick check)
source("code/run_scenario.R")
```

## Directory Structure

```
code/          R scripts for simulation studies
output/        Generated results (not tracked by git)
```

## Reproduction

1. `Rscript code/00_setup.R` — install dependencies
2. `Rscript code/run_scenario.R` — run a single scenario for quick verification
3. For full simulation studies, see `code/README.md`
