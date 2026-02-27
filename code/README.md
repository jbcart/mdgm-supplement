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
  gaussian/
    helpers.R              # Gaussian-specific: generate_gaussian_obs, metrics, bayesImageS PFAB
    01_run_single_rep.R    # Fit 4 methods for one replicate of one Gaussian scenario
    02_aggregate_scenario.R  # Aggregate per-rep files for a single Gaussian scenario
    03_aggregate_all.R     # Aggregate all Gaussian scenarios into summary CSV
```

## Binary study (16x16 grid, Bernoulli emissions)

32 scenarios, 100 replicates each, 4 models (mrf_exact, mrf_pl, mdgm_st, mdgm_ao):

```bash
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

## Gaussian study (Gaussian emissions)

15 scenarios, 20 replicates each, 4 methods (mdgm_st, mdgm_ao, mrf_pl, bis_pfab):

```bash
# Run via GNU parallel (recommended for servers)
./scripts/gaussian/run_scenario.sh 1      # scenario 1, 20 reps
./scripts/gaussian/batch_server1.sh       # k=3 + k=4, 100x100 (1-6)
./scripts/gaussian/batch_server2.sh       # k=5 low + high noise, 100x100 (7-12)
./scripts/gaussian/batch_server3.sh       # k=5, 1000x1000 (13-15)
```

### Scenario index mapping

| Idx | k | grid | sigma | beta |
|-----|---|------|-------|------|
| 1 | 3 | 100 | 0.20 | 0.35 |
| 2 | 3 | 100 | 0.20 | 0.65 |
| 3 | 3 | 100 | 0.20 | 0.95 |
| 4 | 4 | 100 | 0.20 | 0.38 |
| 5 | 4 | 100 | 0.20 | 0.71 |
| 6 | 4 | 100 | 0.20 | 1.04 |
| 7 | 5 | 100 | 0.20 | 0.41 |
| 8 | 5 | 100 | 0.20 | 0.76 |
| 9 | 5 | 100 | 0.20 | 1.12 |
| 10 | 5 | 100 | 0.50 | 0.41 |
| 11 | 5 | 100 | 0.50 | 0.76 |
| 12 | 5 | 100 | 0.50 | 1.12 |
| 13 | 5 | 1000 | 0.20 | 0.41 |
| 14 | 5 | 1000 | 0.20 | 0.76 |
| 15 | 5 | 1000 | 0.20 | 1.12 |

Tag format: `gauss_k<K>_b<BETA>_s<SIGMA>_g<GRID>` (e.g., `gauss_k3_b0_35_s0_20_g100`)
