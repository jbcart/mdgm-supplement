#!/usr/bin/env bash
## Server 1 (72 cores): high-psi scenarios (both emission configs)
## Scenarios 5-8 (e0_2, psi=0.5-0.8) and 13-16 (e0_05, psi=0.5-0.8)
set -uo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=8

for i in 5 6 7 8 13 14 15 16; do
  ./scripts/binary/run_scenario.sh $i 100 $N_JOBS || echo "WARNING: scenario $i failed, continuing..."
done

echo "Server 1 batch complete."
