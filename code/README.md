# Code

## Directory structure

```
code/
  common/
    helpers.R              # Shared utilities: data generation, metrics, TV distance, diagnostics
  binary/
    00_setup.R             # Install dependencies (mdgm from GitHub, mclust, bayesImageS)
    01_run_single_rep.R    # Fit 4 models for one replicate of one scenario
    02_aggregate_scenario.R  # Aggregate per-rep files for a single scenario
    03_aggregate_all.R     # Aggregate all scenarios into summary CSV
    run_all_g16.R          # Full sequential runner: all 32 g16 scenarios
    run_scenario.R         # Quick single-scenario runner (100 reps, psi=0.5, theta=0.2/0.8)
  gaussian/
    helpers.R              # Gaussian-specific: generate_gaussian_obs, metrics, bayesImageS PFAB
    01_run_single_rep.R    # Fit 4 methods for one replicate of one Gaussian scenario
    02_aggregate_scenario.R  # Aggregate per-rep files for a single Gaussian scenario
    03_aggregate_all.R     # Aggregate all Gaussian scenarios into summary CSV
```

## Binary study (16x16 grid, Bernoulli emissions)

32 scenarios, 100 replicates each, 4 models (mrf_exact, mrf_pl, mdgm_st, mdgm_ao):

```bash
# Run all 32 scenarios sequentially
Rscript code/binary/run_all_g16.R

# Run a single scenario by index (1-32)
Rscript code/binary/run_all_g16.R 1

# Run via GNU parallel (recommended for servers)
./scripts/binary/run_scenario.sh 1        # scenario 1, 100 reps
./scripts/binary/batch_server1.sh         # scenarios 5-8, 13-16
```

### Scenario index mapping

| Index | Scenario | psi | theta | Obs type |
|-------|----------|-----|-------|----------|
| 1-8 | `p0_{1..8}e0_2r2o1g16` | 0.1-0.8 | (0.20, 0.80) | 2 reps/site |
| 9-16 | `p0_{1..8}e0_05r2o1g16` | 0.1-0.8 | (0.05, 0.95) | 2 reps/site |
| 17-24 | `p0_{1..8}e0_1l1_39o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(1.39), ~25% missing |
| 25-32 | `p0_{1..8}e0_1l2_3o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(2.3), ~10% missing |

## Gaussian study (100x100 grid, Gaussian emissions)

18 scenarios (k=2,3,6 x 6 beta values), 50 replicates each, 4 methods (mdgm_st, mdgm_ao, mrf_pl, bis_pfab):

```bash
# Run via GNU parallel (recommended for servers)
./scripts/gaussian/run_scenario.sh 1      # scenario 1, 50 reps
./scripts/gaussian/batch_server1.sh       # k=6 scenarios (13-18)
```

### Scenario index mapping

| Index | k | beta values |
|-------|---|-------------|
| 1-6 | 2 | 0.15, 0.30, 0.45, 0.60, 0.75, 0.85 |
| 7-12 | 3 | 0.17, 0.33, 0.50, 0.67, 0.83, 1.00 |
| 13-18 | 6 | 0.20, 0.40, 0.60, 0.80, 1.00, 1.20 |
