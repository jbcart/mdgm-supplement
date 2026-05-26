#!/usr/bin/env bash
## Run all replicates for a single 32x32 scenario using GNU parallel
##
## Usage: ./scripts/binary/run_scenario_g32.sh <scenario_index> [n_replicates] [n_jobs]

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <scenario_index> [n_replicates] [n_jobs]"
  exit 1
fi

SCENARIO_IDX=$1
N_REPS=${2:-100}
N_JOBS=${3:-4}

# Get scenario tag
SCENARIO_TAG=$(Rscript -e '
psi_vec <- seq(0.1, 0.8, by = 0.1)
scenarios <- character(0)
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_2r2o1g32", round(psi * 10)))
for (psi in psi_vec) scenarios <- c(scenarios, sprintf("p0_%de0_1l1_39o1g32", round(psi * 10)))
cat(scenarios['"$SCENARIO_IDX"'])
')

echo "=== Scenario $SCENARIO_IDX: $SCENARIO_TAG ($N_REPS replicates, $N_JOBS jobs) ==="

mkdir -p "output/$SCENARIO_TAG"

seq 1 "$N_REPS" | parallel --jobs "$N_JOBS" --progress \
  "Rscript code/binary/01_run_single_rep_g32.R $SCENARIO_IDX {}"

echo ""
echo "All $N_REPS replicates complete for $SCENARIO_TAG"

Rscript code/binary/02_aggregate_scenario_g32.R "$SCENARIO_TAG"
