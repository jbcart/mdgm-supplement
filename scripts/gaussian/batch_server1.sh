#!/usr/bin/env bash
## Server 1 (112 cores): k=6 scenarios (hardest, most compute)
## Scenarios 13-18
set -euo pipefail
cd "$(dirname "$0")/../.."

for i in 13 14 15 16 17 18; do
  ./scripts/gaussian/run_scenario.sh $i
done

echo "Server 1 Gaussian batch complete."
