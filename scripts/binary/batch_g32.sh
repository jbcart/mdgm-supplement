#!/usr/bin/env bash
## 32x32 binary study: all 16 scenarios
## 1-8:  complete data (eta=0.2, m=2)
## 9-16: missing data (eta=0.1, lambda=1.39)
set -euo pipefail
cd "$(dirname "$0")/../.."

N_JOBS=${1:-4}

for i in $(seq 1 16); do
  ./scripts/binary/run_scenario_g32.sh $i 100 $N_JOBS
done

echo "32x32 batch complete."
