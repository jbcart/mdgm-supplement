#!/usr/bin/env bash
## Run all replicates for a single scenario using GNU parallel
##
## Usage: ./scripts/binary/run_scenario.sh <scenario_index> [n_replicates]
##
## Runs `parallel --jobs 25%` over replicates, then aggregates.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <scenario_index> [n_replicates]"
  exit 1
fi

SCENARIO_IDX=$1
N_REPS=${2:-100}

# Get scenario tag by running R to look it up
SCENARIO_TAG=$(Rscript -e '
psi_vec <- seq(0.1, 0.8, by = 0.1)
scenarios <- character(0)
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_2r2o1g16", round(psi * 10)))
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_05r2o1g16", round(psi * 10)))
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_1l1_39o1g16", round(psi * 10)))
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_1l2_3o1g16", round(psi * 10)))
cat(scenarios['"$SCENARIO_IDX"'])
')

echo "=== Scenario $SCENARIO_IDX: $SCENARIO_TAG ($N_REPS replicates) ==="

# Create output directory
mkdir -p "output/$SCENARIO_TAG"

# Run replicates in parallel (25% of available cores)
seq 1 "$N_REPS" | parallel --jobs 25% --progress \
  "Rscript code/binary/01_run_single_rep.R $SCENARIO_IDX {}"

echo ""
echo "All $N_REPS replicates complete for $SCENARIO_TAG"

# Aggregate per-rep files into combined result
Rscript code/binary/02_aggregate_scenario.R "$SCENARIO_TAG"
