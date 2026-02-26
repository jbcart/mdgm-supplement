## Aggregate per-rep files for a single scenario
##
## Usage: Rscript code/04_aggregate_scenario.R <scenario_tag>
##
## Reads all output/<tag>/rep_*.rds files, computes TV distances
## (mrf_exact as reference), and saves combined output/<tag>_100rep.rds
## in the same format as run_all_g16.R output.

source(file.path("code", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/04_aggregate_scenario.R <scenario_tag>")
}
sc_tag <- args[1]

# --- Read per-rep files ---
rep_dir <- file.path("output", sc_tag)
rep_files <- sort(list.files(rep_dir, pattern = "^rep_.*\\.rds$",
                             full.names = TRUE))

if (length(rep_files) == 0) {
  stop(sprintf("No rep files found in %s", rep_dir))
}

n_reps <- length(rep_files)
cat(sprintf("Aggregating %d replicates for %s\n", n_reps, sc_tag))

model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")
all_metrics <- vector("list", n_reps)

for (i in seq_along(rep_files)) {
  rep_metrics <- readRDS(rep_files[i])

  # Compute TV distances vs mrf_exact
  ref_eq <- rep_metrics[["mrf_exact"]]$eq_vec
  for (name in model_names) {
    rep_metrics[[name]]$tv_vs_mrf <- tv_distance(
      ref_eq, rep_metrics[[name]]$eq_vec
    )
  }

  all_metrics[[i]] <- rep_metrics
}

# --- Print summary table ---
cat(sprintf("\n--- %s: aggregated across %d replicates ---\n", sc_tag, n_reps))
cat(sprintf("%-12s %10s %10s %10s %10s %10s %10s\n",
            "Model", "psi_pm", "th1_pm", "th2_pm",
            "acc_pred", "TV_dist", "time"))
cat(paste(rep("-", 74), collapse = ""), "\n")

for (name in model_names) {
  vals <- lapply(all_metrics, function(r) r[[name]])
  cat(sprintf("%-12s %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.1f(%4.1f)\n",
              name,
              mean(vapply(vals, function(v) v$psi_pm, double(1))),
              1000 * sd(vapply(vals, function(v) v$psi_pm, double(1))),
              mean(vapply(vals, function(v) v$theta_pm[1], double(1))),
              1000 * sd(vapply(vals, function(v) v$theta_pm[1], double(1))),
              mean(vapply(vals, function(v) v$theta_pm[2], double(1))),
              1000 * sd(vapply(vals, function(v) v$theta_pm[2], double(1))),
              mean(vapply(vals, function(v) v$acc_pred, double(1))),
              1000 * sd(vapply(vals, function(v) v$acc_pred, double(1))),
              mean(vapply(vals, function(v) v$tv_vs_mrf, double(1))),
              1000 * sd(vapply(vals, function(v) v$tv_vs_mrf, double(1))),
              mean(vapply(vals, function(v) v$elapsed, double(1))),
              sd(vapply(vals, function(v) v$elapsed, double(1)))))
}
cat("(parentheses = SD x 1000)\n\n")

# --- Save combined result (same format as run_all_g16.R) ---
outfile <- sprintf("output/%s_%drep.rds", sc_tag, n_reps)
saveRDS(all_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
