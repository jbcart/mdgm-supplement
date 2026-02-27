#!/usr/bin/env bash
## Server 1 (72 cores): k=3 + k=4, 100x100 grid
## Scenarios 1-6
set -euo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=6

for i in 1 2 3 4 5 6; do
  ./scripts/gaussian/run_scenario.sh $i 20 $N_JOBS
done

echo "Server 1 Gaussian batch complete."
