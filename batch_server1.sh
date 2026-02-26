#!/usr/bin/env bash
## Server 1 (112 cores): high-psi scenarios (both emission configs)
## Scenarios 5-8 (e0_2, psi=0.5-0.8) and 13-16 (e0_05, psi=0.5-0.8)
set -euo pipefail
cd "$(dirname "$0")"

for i in 5 6 7 8 13 14 15 16; do
  ./run_scenario.sh $i
done

echo "Server 1 batch complete."
