#!/usr/bin/env bash
## Run all replicates for a single Gaussian scenario using GNU parallel
##
## Usage: ./scripts/gaussian/run_scenario.sh <scenario_index> [n_replicates]
##
## Runs `parallel --jobs 25%` over replicates, then aggregates.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <scenario_index> [n_replicates]"
  exit 1
fi

SCENARIO_IDX=$1
N_REPS=${2:-50}

# Get scenario tag by running R to look it up
SCENARIO_TAG=$(Rscript -e '
beta_k2 <- c(0.15, 0.30, 0.45, 0.60, 0.75, 0.85)
beta_k3 <- c(0.17, 0.33, 0.50, 0.67, 0.83, 1.00)
beta_k6 <- c(0.20, 0.40, 0.60, 0.80, 1.00, 1.20)
scenarios <- character(0)
for (b in beta_k2) scenarios <- c(scenarios, sprintf("gauss_k2_b%s_g100", gsub("\\\\.", "_", sprintf("%.2f", b))))
for (b in beta_k3) scenarios <- c(scenarios, sprintf("gauss_k3_b%s_g100", gsub("\\\\.", "_", sprintf("%.2f", b))))
for (b in beta_k6) scenarios <- c(scenarios, sprintf("gauss_k6_b%s_g100", gsub("\\\\.", "_", sprintf("%.2f", b))))
cat(scenarios['"$SCENARIO_IDX"'])
')

echo "=== Scenario $SCENARIO_IDX: $SCENARIO_TAG ($N_REPS replicates) ==="

# Create output directory
mkdir -p "output/$SCENARIO_TAG"

# Run replicates in parallel (25% of available cores)
seq 1 "$N_REPS" | parallel --jobs 25% --progress \
  "Rscript code/gaussian/01_run_single_rep.R $SCENARIO_IDX {}"

echo ""
echo "All $N_REPS replicates complete for $SCENARIO_TAG"

# Aggregate per-rep files into combined result
Rscript code/gaussian/02_aggregate_scenario.R "$SCENARIO_TAG"
