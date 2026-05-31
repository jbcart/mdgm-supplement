#!/usr/bin/env bash
## Binary 16x16: missing data scenarios (lambda=1.39 and lambda=2.3)
## Scenarios 17-32
set -uo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=${1:-8}

for i in $(seq 17 32); do
  ./scripts/binary/run_scenario.sh $i 100 $N_JOBS || echo "WARNING: scenario $i failed, continuing..."
done

echo "Binary 16x16 missing data batch done."
