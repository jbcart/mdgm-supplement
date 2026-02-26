#!/usr/bin/env bash
## Server 2 (72 cores): low-psi scenarios (both emission configs)
## Scenarios 1-4 (e0_2, psi=0.1-0.4) and 9-12 (e0_05, psi=0.1-0.4)
set -euo pipefail
cd "$(dirname "$0")"

for i in 1 2 3 4 9 10 11 12; do
  ./run_scenario.sh $i
done

echo "Server 2 batch complete."
