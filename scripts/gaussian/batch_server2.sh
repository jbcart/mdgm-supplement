#!/usr/bin/env bash
## Server 2: k=5 + k=6, 100x100 grid, sigma=0.50
## Scenarios 7-12
set -euo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=6

for i in 7 8 9 10 11 12; do
  ./scripts/gaussian/run_scenario.sh $i 20 $N_JOBS
done

echo "Server 2 Gaussian batch complete."
