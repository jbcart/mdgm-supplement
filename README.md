# mdgm-supplement

Supplementary code for reproducing simulation studies and analyses in:

> Carter, J. B. and Calder, C. A. (2024). "Mixture of Directed Graphical Models for Discrete Spatial Random Fields."
> [arXiv:2406.15700](https://arxiv.org/abs/2406.15700)

## Requirements

- R >= 4.1 (developed with R 4.5.3)
- [mdgm](https://github.com/jbcart/mdgm) R package (>= 0.1.0, installed from GitHub)
- CRAN packages: bayesImageS, sf, igraph, tidyverse, ggplot2, coda, patchwork, mclust

Install all dependencies:
```bash
Rscript code/common/00_setup.R
```

## Directory Structure

```
code/
  common/
    00_setup.R               # Install dependencies
    helpers.R                # Shared utilities
  binary/                    # Binary simulation study (16x16, Bernoulli emissions, 32 scenarios)
  gaussian/                  # Gaussian simulation study (100x100, 12 scenarios, k=3-6)
  disorder/                  # Columbus physical disorder analysis and cross-validation
    00_analysis.R            # Main analysis (MDGM-ST and aMRF)
    01_cv_single_run.R       # Single cross-validation replicate
    02_cv_aggregate.R        # Aggregate cross-validation results
    clean_data.Rbin          # Anonymized AHDC garbage ratings data (see Data section below)
  tables_plots.R             # Generate simulation study figures
scripts/
  binary/                    # Shell scripts for binary batch execution
  gaussian/                  # Shell scripts for Gaussian batch execution
  disorder/                  # Shell scripts for cross-validation
output/
  binary/                    # Binary study: per-scenario .rds + figures/
  gaussian/                  # Gaussian study: per-scenario .rds + figures/
  disorder/                  # Disorder application: cv_summary.csv + figures/
```

## Figure and Table Reproduction

The paper is submitted as two documents: the main paper and an
online appendix with its own figure/table numbering.
References below use "Main Paper" or "Appendix" to disambiguate.

| Paper element | Script |
|---------------|--------|
| Main Paper Figure 1 (binary sim, 16x16, missing data) | `code/tables_plots.R` → `sim_o1g16_lambda_bs.jpeg` |
| Main Paper Figure 2 (Columbus posterior maps and comparison) | `code/disorder/00_analysis.R` → `posterior_mdgm_st.jpeg`, `posterior_difference.jpeg`, `posterior_comparison.jpeg` |
| Main Paper Table 1 (MCMC algorithm summary) | Hand-written; not code output |
| Appendix Figure 6 (binary sim, 16x16, complete data) | `code/tables_plots.R` → `sim_o1g16_error_bs.jpeg` |
| Appendix Figure 7 (binary sim, 32x32, complete data) | `code/tables_plots.R` → `sim_o1g32_error_bs.jpeg` |
| Appendix Figure 8 (binary sim, 32x32, missing data) | `code/tables_plots.R` → `sim_o1g32_lambda_bs.jpeg` |
| Appendix Table 1 (spanning tree generation timings) | Hand-written; not code output |
| Appendix Table 2 (Gaussian simulation scenarios) | Parameters defined in `code/gaussian/01_run_single_rep.R` |
| Appendix Figure 9 (Gaussian simulation results) | `code/tables_plots.R` → `sim_gaussian_bs.jpeg` |
| Appendix Figure 10 (edge inclusion probabilities) | `code/disorder/00_analysis.R` → `edge_inclusion.jpeg` |
| Appendix J (cross-validation results) | `code/disorder/01_cv_single_run.R` + `02_cv_aggregate.R` |

## Quick Start

To regenerate every figure and table above from the pre-computed simulation
output already bundled in `output/` (fast, ~10-20 minutes; does not re-run
the simulation studies themselves):

```bash
Rscript code/run_all.R
```

This wrapper runs each analysis step in the order the results appear in the
paper, printing which figure/table file each step produces. See
`code/run_all.R` for details, including how to re-run the full simulation
pipeline from scratch (`Rscript code/run_all.R --full`; requires GNU
parallel and takes hours to days).

## Reproduction

### Simulation studies

1. Install dependencies:
   ```bash
   Rscript code/common/00_setup.R
   ```

2. Run simulation batches (default 4 parallel jobs; pass a number to change):
   ```bash
   # Binary 16x16: complete data scenarios (1-16)
   ./scripts/binary/batch_complete.sh [n_jobs]

   # Binary 16x16: missing data scenarios (17-32)
   ./scripts/binary/batch_missing.sh [n_jobs]

   # Binary 32x32: all scenarios
   ./scripts/binary/batch_g32.sh [n_jobs]

   # Gaussian: all 12 scenarios
   ./scripts/gaussian/batch_all.sh [n_jobs]
   ```

   Or run a single scenario:
   ```bash
   # Binary 16x16: single scenario (index 1-32)
   ./scripts/binary/run_scenario.sh <scenario_index> [n_replicates] [n_jobs]

   # Binary 32x32: single scenario (index 1-16)
   ./scripts/binary/run_scenario_g32.sh <scenario_index> [n_replicates] [n_jobs]

   # Gaussian: single scenario (index 1-12)
   ./scripts/gaussian/run_scenario.sh <scenario_index> [n_replicates] [n_jobs]
   ```

3. Aggregate results:
   ```bash
   Rscript code/binary/03_aggregate_all.R
   Rscript code/gaussian/03_aggregate_all.R
   ```

4. Generate figures:
   ```bash
   Rscript code/tables_plots.R
   ```

See `code/README.md` for scenario index mappings.

### Columbus physical disorder analysis

```bash
Rscript code/disorder/00_analysis.R
```

### Cross-validation study

```bash
./scripts/disorder/run_cv.sh [n_runs]
```

## Data

The file `code/disorder/clean_data.Rbin` contains anonymized data from the Adolescent Health and Development in Context (AHDC) Study. Ratings are aggregated at the block group level with no cross-block-group respondent tracking. The file contains three R objects:

- **`bg270`**: An `sf` data frame (615 rows × 15 columns) containing census block group geometries and identifiers for the 615 block groups within the I-270 belt loop of Columbus, Ohio. See the data dictionary below.
- **`nug_bg270`**: A list of length 615 defining the first-order neighborhood structure (NUG). Each element `nug_bg270[[i]]` is an integer vector of the indices of block group `i`'s neighbors (rook contiguity).
- **`y_bg270`**: A list of length 615 containing the binary garbage ratings for each block group. Each element `y_bg270[[i]]` is a numeric vector of 0/1 ratings (1 = garbage is "a big problem" or "somewhat of a problem", 0 = "not a problem"). Block groups with no ratings have `NA`.

### Data dictionary for `bg270`

| Column | Type | Description |
|--------|------|-------------|
| `STATEFP` | character | FIPS state code (39 = Ohio) |
| `COUNTYFP` | character | FIPS county code |
| `TRACTCE` | character | Census tract code |
| `BLKGRPCE` | character | Block group code within tract |
| `GEOID` | character | Full block group FIPS identifier (state + county + tract + block group) |
| `NAMELSAD` | character | Census name description (e.g., "Block Group 1") |
| `MTFCC` | character | MAF/TIGER feature class code |
| `FUNCSTAT` | character | Functional status code |
| `ALAND` | numeric | Land area (square meters) |
| `AWATER` | numeric | Water area (square meters) |
| `INTPTLAT` | character | Latitude of the internal point |
| `INTPTLON` | character | Longitude of the internal point |
| `Area` | numeric | Area (square kilometers) |
| `in270` | integer | Indicator for block group within I-270 belt loop (all 1) |
| `geometry` | sfc_MULTIPOLYGON | Block group boundary polygon |
