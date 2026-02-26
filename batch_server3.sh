#!/usr/bin/env bash
## Server 3 (72 cores): missing-data scenarios (generally faster)
## Scenarios 17-32 (lambda=1.39 and lambda=2.3)
set -euo pipefail
cd "$(dirname "$0")"

for i in $(seq 17 32); do
  ./run_scenario.sh $i
done

echo "Server 3 batch complete."
