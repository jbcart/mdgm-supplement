#!/usr/bin/env bash
## Binary 16x16: complete data scenarios (eta=0.2 and eta=0.05)
## Scenarios 1-16
set -uo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=${1:-8}

for i in $(seq 1 16); do
  ./scripts/binary/run_scenario.sh $i 100 $N_JOBS || echo "WARNING: scenario $i failed, continuing..."
done

echo "Binary 16x16 complete data batch done."
