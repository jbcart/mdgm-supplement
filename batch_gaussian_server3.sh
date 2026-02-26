#!/usr/bin/env bash
## Server 3 (72 cores): k=2 scenarios
## Scenarios 1-6
set -euo pipefail
cd "$(dirname "$0")"

for i in 1 2 3 4 5 6; do
  ./run_scenario_gaussian.sh $i
done

echo "Server 3 Gaussian batch complete."
