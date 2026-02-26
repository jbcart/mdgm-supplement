# Code

## Scripts

| File | Description |
|------|-------------|
| `00_setup.R` | Install dependencies (mdgm from GitHub) |
| `helpers.R` | Utility functions: data generation, metrics, TV distance, diagnostics |
| `run_scenario.R` | Quick single-scenario runner (100 reps, psi=0.5, theta=0.2/0.8) |
| `run_all_g16.R` | Full simulation study: all 32 g16 scenarios from the paper |
| `01_simulate_data.R` | Generate MRF data via `sample_mrf()` |
| `02_fit_models.R` | Fit all models for a scenario |
| `03_aggregate.R` | Collect results across replicates |

## Full Simulation Study

All 32 scenarios on a 16x16 grid, 100 replicates each, 4 models
(mrf_exact, mrf_pl, mdgm_st, mdgm_ao):

```bash
# Run all 32 scenarios sequentially
Rscript code/run_all_g16.R

# Run a single scenario by index (1-32) for parallelization
Rscript code/run_all_g16.R 1
Rscript code/run_all_g16.R 2
# ...
Rscript code/run_all_g16.R 32
```

### Scenario index mapping

| Index | Scenario | psi | theta | Obs type |
|-------|----------|-----|-------|----------|
| 1-8 | `p0_{1..8}e0_2r2o1g16` | 0.1-0.8 | (0.20, 0.80) | 2 reps/site |
| 9-16 | `p0_{1..8}e0_05r2o1g16` | 0.1-0.8 | (0.05, 0.95) | 2 reps/site |
| 17-24 | `p0_{1..8}e0_1l1_39o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(1.39), ~25% missing |
| 25-32 | `p0_{1..8}e0_1l2_3o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(2.3), ~10% missing |
