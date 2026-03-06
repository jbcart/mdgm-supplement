# mdgm-supplement

Supplementary code for reproducing simulation studies in:

> Carter, J. B. and Calder, C. A. (2024). "Mixture of Directed Graphical Models for Discrete Spatial Random Fields."
> [arXiv:2406.15700](https://arxiv.org/abs/2406.15700)

## Requirements

- R >= 4.1
- [mdgm](https://github.com/jbcart/mdgm) R package (installed from GitHub)

## Quick Start

```bash
# Install dependencies
Rscript code/binary/00_setup.R

# Run a single binary scenario (quick check)
./scripts/binary/run_scenario.sh 1

# Run a single Gaussian scenario
./scripts/gaussian/run_scenario.sh 1
```

## Directory Structure

```
code/
  common/helpers.R           # Shared utilities
  binary/                    # Binary study (16x16, Bernoulli emissions, 32 scenarios)
  gaussian/                  # Gaussian study (100x100, 12 scenarios, k=3-6)
scripts/
  binary/                    # Shell scripts for binary batch execution
  gaussian/                  # Shell scripts for Gaussian batch execution
output/                      # Generated results (not tracked by git)
```

## Reproduction

1. `Rscript code/binary/00_setup.R` -- install dependencies
2. Run individual scenarios or server batches (see `code/README.md` for details)
3. Aggregate results:
   ```bash
   Rscript code/binary/03_aggregate_all.R
   Rscript code/gaussian/03_aggregate_all.R
   ```
