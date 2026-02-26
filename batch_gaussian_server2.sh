#!/usr/bin/env bash
## Server 2 (72 cores): k=3 scenarios
## Scenarios 7-12
set -euo pipefail
cd "$(dirname "$0")"

for i in 7 8 9 10 11 12; do
  ./run_scenario_gaussian.sh $i
done

echo "Server 2 Gaussian batch complete."
