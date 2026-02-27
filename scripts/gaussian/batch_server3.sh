#!/usr/bin/env bash
## Server 3 (112 cores): k=5, 1000x1000 grid (heaviest compute)
## Scenarios 13-15
##
## Memory note: each 1000x1000 rep uses ~60-80 GB peak (z + DAG matrices).
## Limit to 2 concurrent jobs to stay within ~120 GB total.
set -euo pipefail
cd "$(dirname "$0")/../.."

N_REPS=${1:-20}

for i in 13 14 15; do
  echo "=== Scenario $i ($N_REPS replicates, 2 jobs max) ==="
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
cat(scenarios['"$i"'])
')
  mkdir -p "output/$SCENARIO_TAG"
  seq 1 "$N_REPS" | parallel --jobs 2 --progress \
    "Rscript code/gaussian/01_run_single_rep.R $i {}"
  Rscript code/gaussian/02_aggregate_scenario.R "$SCENARIO_TAG"
done

echo "Server 3 Gaussian batch complete."
