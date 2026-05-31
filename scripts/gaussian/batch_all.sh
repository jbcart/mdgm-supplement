#!/usr/bin/env bash
## Gaussian study: all 12 scenarios (k=3,4,5,6 x 3 beta values)
set -uo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=${1:-8}

for i in $(seq 1 12); do
  ./scripts/gaussian/run_scenario.sh $i 20 $N_JOBS || echo "WARNING: scenario $i failed, continuing..."
done

echo "Gaussian batch done."
