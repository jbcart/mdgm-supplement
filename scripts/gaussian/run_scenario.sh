#!/usr/bin/env bash
## Run all replicates for a single Gaussian scenario using GNU parallel
##
## Usage: ./scripts/gaussian/run_scenario.sh <scenario_index> [n_replicates] [n_jobs]
##
## Runs replicates in parallel, then aggregates.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <scenario_index> [n_replicates] [n_jobs]"
  exit 1
fi

SCENARIO_IDX=$1
N_REPS=${2:-20}
N_JOBS=${3:-6}

# Get scenario tag by running R to look it up
SCENARIO_TAG=$(Rscript -e '
make_tag <- function(k, beta, sigma, grid) {
  sprintf("gauss_k%d_b%s_s%s_g%d",
          k,
          gsub("\\.", "_", sprintf("%.2f", beta)),
          gsub("\\.", "_", sprintf("%.2f", sigma)),
          grid)
}
scenarios <- character(0)
for (b in c(0.35, 0.65, 0.95)) scenarios <- c(scenarios, make_tag(3, b, 0.20, 100))
for (b in c(0.38, 0.71, 1.04)) scenarios <- c(scenarios, make_tag(4, b, 0.20, 100))
for (b in c(0.41, 0.76, 1.12)) scenarios <- c(scenarios, make_tag(5, b, 0.20, 100))
for (b in c(0.41, 0.76, 1.12)) scenarios <- c(scenarios, make_tag(5, b, 0.50, 100))
for (b in c(0.41, 0.76, 1.12)) scenarios <- c(scenarios, make_tag(5, b, 0.20, 1000))
cat(scenarios['"$SCENARIO_IDX"'])
')

echo "=== Scenario $SCENARIO_IDX: $SCENARIO_TAG ($N_REPS replicates) ==="

# Create output directory
mkdir -p "output/$SCENARIO_TAG"

# Run replicates in parallel
seq 1 "$N_REPS" | parallel --jobs "$N_JOBS" --progress \
  "Rscript code/gaussian/01_run_single_rep.R $SCENARIO_IDX {}"

echo ""
echo "All $N_REPS replicates complete for $SCENARIO_TAG"

# Aggregate per-rep files into combined result
Rscript code/gaussian/02_aggregate_scenario.R "$SCENARIO_TAG"
