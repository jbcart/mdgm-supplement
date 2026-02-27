#!/usr/bin/env bash
## Server 3 (112 cores): k=5, 1000x1000 grid (heaviest compute)
## Scenarios 13-15
set -euo pipefail
cd "$(dirname "$0")/../.."

for i in 13 14 15; do
  ./scripts/gaussian/run_scenario.sh $i
done

echo "Server 3 Gaussian batch complete."
