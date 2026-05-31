# Code

## Directory structure

```
code/
  common/
    00_setup.R             # Install dependencies (mdgm from GitHub, mclust, bayesImageS)
    helpers.R              # Shared utilities: data generation, metrics, TV distance, diagnostics
  binary/
    01_run_single_rep.R    # Fit 4 models for one replicate of one scenario
    02_aggregate_scenario.R  # Aggregate per-rep files for a single scenario
    03_aggregate_all.R     # Aggregate all scenarios into summary CSV
  gaussian/
    helpers.R              # Gaussian-specific: generate_gaussian_obs, metrics, bayesImageS PFAB
    01_run_single_rep.R    # Fit 4 methods for one replicate of one Gaussian scenario
    02_aggregate_scenario.R  # Aggregate per-rep files for a single Gaussian scenario
    03_aggregate_all.R     # Aggregate all Gaussian scenarios into summary CSV
  disorder/
    00_analysis.R          # Columbus physical disorder analysis (MDGM-ST and aMRF)
    01_cv_single_run.R     # Single cross-validation replicate
    02_cv_aggregate.R      # Aggregate cross-validation results
    clean_data.Rbin        # Anonymized AHDC garbage ratings data
  tables_plots.R           # Generate all figures from simulation and analysis output
```

## Binary study (16x16 grid, Bernoulli emissions)

32 scenarios, 100 replicates each, 4 models (mrf_exact, mrf_pl, mdgm_st, mdgm_ao):

```bash
# Run all complete data scenarios (1-16)
./scripts/binary/batch_complete.sh [n_jobs]

# Run all missing data scenarios (17-32)
./scripts/binary/batch_missing.sh [n_jobs]

# Or run a single scenario
./scripts/binary/run_scenario.sh <index> [n_replicates] [n_jobs]
```

### Scenario index mapping

| Index | Scenario | psi | theta | Obs type | Batch |
|-------|----------|-----|-------|----------|-------|
| 1-8 | `p0_{1..8}e0_2r2o1g16` | 0.1-0.8 | (0.20, 0.80) | 2 reps/site | `batch_complete.sh` |
| 9-16 | `p0_{1..8}e0_05r2o1g16` | 0.1-0.8 | (0.05, 0.95) | 2 reps/site | `batch_complete.sh` |
| 17-24 | `p0_{1..8}e0_1l1_39o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(1.39), ~25% missing | `batch_missing.sh` |
| 25-32 | `p0_{1..8}e0_1l2_3o1g16` | 0.1-0.8 | (0.10, 0.90) | Pois(2.3), ~10% missing | `batch_missing.sh` |

## Binary study (32x32 grid, Bernoulli emissions)

16 scenarios (high uncertainty only), 100 replicates each, 2 models (mdgm_st, mrf_pl):

```bash
# Run all 32x32 scenarios
./scripts/binary/batch_g32.sh [n_jobs]

# Or run a single scenario
./scripts/binary/run_scenario_g32.sh <index> [n_replicates] [n_jobs]
```

### Scenario index mapping

| Index | Scenario | psi | theta | Obs type |
|-------|----------|-----|-------|----------|
| 1-8 | `p0_{1..8}e0_2r2o1g32` | 0.1-0.8 | (0.20, 0.80) | 2 reps/site |
| 9-16 | `p0_{1..8}e0_1l1_39o1g32` | 0.1-0.8 | (0.10, 0.90) | Pois(1.39), ~25% missing |

## Gaussian study (Gaussian emissions)

12 scenarios, 20 replicates each, 4 methods (mdgm_st, mdgm_ao, mrf_pl, bis_pfab):

```bash
# Run all 12 Gaussian scenarios
./scripts/gaussian/batch_all.sh [n_jobs]

# Or run a single scenario
./scripts/gaussian/run_scenario.sh <index> [n_replicates] [n_jobs]
```

### Scenario index mapping

| Idx | k | grid | sigma | beta |
|-----|---|------|-------|------|
| 1 | 3 | 100 | 0.50 | 0.35 |
| 2 | 3 | 100 | 0.50 | 0.65 |
| 3 | 3 | 100 | 0.50 | 0.95 |
| 4 | 4 | 100 | 0.50 | 0.38 |
| 5 | 4 | 100 | 0.50 | 0.71 |
| 6 | 4 | 100 | 0.50 | 1.04 |
| 7 | 5 | 100 | 0.50 | 0.41 |
| 8 | 5 | 100 | 0.50 | 0.76 |
| 9 | 5 | 100 | 0.50 | 1.12 |
| 10 | 6 | 100 | 0.50 | 0.43 |
| 11 | 6 | 100 | 0.50 | 0.80 |
| 12 | 6 | 100 | 0.50 | 1.18 |

Tag format: `gauss_k<K>_b<BETA>_s<SIGMA>_g<GRID>` (e.g., `gauss_k3_b0_35_s0_50_g100`)
