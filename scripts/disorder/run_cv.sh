#!/bin/bash
## Run cross-validation study for disorder analysis
##
## Usage: bash code/disorder/run_cv.sh [N_RUNS]
## (Run from the mdgm-supplement root directory)

N_RUNS=${1:-100}

mkdir -p output/disorder

echo "Running ${N_RUNS} CV iterations..."
seq 1 "$N_RUNS" | parallel --jobs 25% Rscript code/disorder/01_cv_single_run.R {}

echo "Aggregating results..."
Rscript code/disorder/02_cv_aggregate.R
