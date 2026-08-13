## Master wrapper: reproduce every figure and table in the paper and appendix.
##
## Usage: Rscript code/run_all.R
## (Run from the mdgm-supplement root directory)
##
## This script runs in two stages:
##
##   Stage 1 (default, ~10-20 minutes): regenerates every figure/table from the
##   bundled, pre-computed simulation output already committed under output/.
##   This is what most readers want -- it reproduces the exact figures in the
##   paper without re-running the underlying simulation studies.
##
##   Stage 2 (opt-in, hours to days): re-runs the simulation studies themselves
##   from scratch via the shell scripts in scripts/, regenerating the .rds
##   files that Stage 1 reads. This requires GNU parallel and is only needed
##   if you want to verify the simulation pipeline itself, not just the
##   figures. Set RUN_FULL_SIMULATION <- TRUE below (or pass --full) to enable.
##
## Each step below is labeled with the exact figure/table it reproduces, using
## the numbering in the split main paper (main.pdf) and online appendix
## (supplement.pdf).

args <- commandArgs(trailingOnly = TRUE)
RUN_FULL_SIMULATION <- "--full" %in% args

stopifnot(
  "Run this script from the mdgm-supplement root directory (Rscript code/run_all.R)" =
    dir.exists("code") && dir.exists("output") && dir.exists("scripts")
)

section <- function(msg) cat(sprintf("\n=== %s ===\n", msg))

# -----------------------------------------------------------------------
# Stage 2 (optional): full re-simulation from scratch
# -----------------------------------------------------------------------
if (RUN_FULL_SIMULATION) {
  section("STAGE 2: Full re-simulation from scratch (this will take a long time)")
  cat("Running binary 16x16 complete-data batch (scenarios 1-16)...\n")
  system2("./scripts/binary/batch_complete.sh")
  cat("Running binary 16x16 missing-data batch (scenarios 17-32)...\n")
  system2("./scripts/binary/batch_missing.sh")
  cat("Running binary 32x32 batch (all scenarios)...\n")
  system2("./scripts/binary/batch_g32.sh")
  cat("Running Gaussian batch (all 12 scenarios)...\n")
  system2("./scripts/gaussian/batch_all.sh")
  cat("Aggregating binary and Gaussian results...\n")
  system2("Rscript", "code/binary/03_aggregate_all.R")
  system2("Rscript", "code/gaussian/03_aggregate_all.R")
  cat("Running cross-validation study (100 replicates)...\n")
  system2("./scripts/disorder/run_cv.sh", "100")
} else {
  section("STAGE 2 skipped: using pre-computed simulation output in output/")
  cat("Pass --full to Rscript code/run_all.R to re-run the full simulation\n")
  cat("pipeline from scratch (requires GNU parallel; hours to days).\n")
}

# -----------------------------------------------------------------------
# Stage 1: regenerate figures/tables from (possibly just-refreshed) output
# -----------------------------------------------------------------------
section("STAGE 1: Regenerating figures from simulation output")

section("Main Paper Figure 1 (binary sim, 16x16, missing data) and")
cat("Appendix Figures 6, 7, 8 (binary sim, complete/32x32) and\n")
cat("Appendix Figure 9 (Gaussian sim)\n")
source("code/tables_plots.R")

section("Main Paper Figure 2 (Columbus posterior maps) and")
cat("Appendix Figure 10 (posterior edge inclusion probabilities)\n")
source("code/disorder/00_analysis.R")

section("Appendix J (cross-validation results)")
cv_summary_path <- file.path("output", "disorder", "cv_summary.csv")
cv_raw_files <- list.files(file.path("output", "disorder"),
                           pattern = "^run_cv_.*\\.rds$")
if (length(cv_raw_files) > 0) {
  # Raw per-replicate CV files are present (e.g., after a --full run): re-aggregate.
  source("code/disorder/02_cv_aggregate.R")
} else if (file.exists(cv_summary_path)) {
  # Raw per-replicate files (100 replicates) are not bundled in this repo, only
  # the final aggregated summary. Nothing to regenerate; report the cached result.
  cat("Raw per-replicate CV files not bundled in this repo (only the final\n")
  cat("aggregated summary is committed). Using the existing cached summary.\n")
  print(read.csv(cv_summary_path), digits = 4)
} else {
  cat("No cross-validation results found. Run ./scripts/disorder/run_cv.sh\n")
  cat("[n_runs] to generate them (see Stage 2 / --full), then re-run this step.\n")
}

section("Not produced by a script")
cat("Main Paper Table 1 (MCMC algorithm summary): hand-written, not code output.\n")
cat("Appendix Table 1 (spanning tree generation timings): hand-written, not code output.\n")
cat("Appendix Table 2 (Gaussian simulation scenarios): parameters listed in\n")
cat("  code/gaussian/01_run_single_rep.R; not written to a separate file.\n")

section("Done")
cat("All figures/tables reproduced. See README.md for the full mapping from\n")
cat("paper elements to scripts.\n")
