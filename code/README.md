# Code

## Scripts

| File | Description |
|------|-------------|
| `00_setup.R` | Install dependencies (mdgm from GitHub) |
| `helpers.R` | Utility functions: data generation, metrics, diagnostics |
| `run_scenario.R` | Quick single-scenario runner for verification |
| `01_simulate_data.R` | Generate MRF data via `sample_mrf()` |
| `02_fit_models.R` | Fit all models for a scenario |
| `03_aggregate.R` | Collect results across replicates |

## Quick Verification

```r
source("code/run_scenario.R")
```

This runs a single 16x16 grid scenario with 4 models (mrf_exact, mrf_pl, mdgm_st, mdgm_ao)
and prints a comparison table.
