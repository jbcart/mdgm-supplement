#!/usr/bin/env bash
## Server 3 (112 cores): missing-data scenarios (generally faster)
## Scenarios 17-32 (lambda=1.39 and lambda=2.3)
set -uo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=8

for i in $(seq 17 32); do
  ./scripts/binary/run_scenario.sh $i 100 $N_JOBS || echo "WARNING: scenario $i failed, continuing..."
done

echo "Server 3 batch complete."
