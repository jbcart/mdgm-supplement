## Aggregate per-rep files for a single Gaussian scenario
##
## Usage: Rscript code/gaussian/02_aggregate_scenario.R <scenario_tag>
##
## Reads all output/<tag>/rep_*.rds files and saves combined
## output/<tag>_50rep.rds.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/gaussian/02_aggregate_scenario.R <scenario_tag>")
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

model_names <- c("mdgm_st", "mdgm_ao", "mrf_pl", "bis_pfab")
all_metrics <- vector("list", n_reps)

for (i in seq_along(rep_files)) {
  all_metrics[[i]] <- readRDS(rep_files[i])
}

# --- Print summary table ---
cat(sprintf("\n--- %s: aggregated across %d replicates ---\n", sc_tag, n_reps))
cat(sprintf("%-12s %10s %10s %10s\n",
            "Model", "ARI", "Misclass", "Time"))
cat(paste(rep("-", 46), collapse = ""), "\n")

for (name in model_names) {
  vals <- lapply(all_metrics, function(r) r[[name]])
  ari_vals <- vapply(vals, function(v) v$ari, double(1))
  mc_vals <- vapply(vals, function(v) v$misclass, double(1))
  t_vals <- vapply(vals, function(v) v$elapsed, double(1))

  # Handle NAs (e.g., bayesImageS failures)
  ari_mean <- mean(ari_vals, na.rm = TRUE)
  ari_sd <- sd(ari_vals, na.rm = TRUE)
  mc_mean <- mean(mc_vals, na.rm = TRUE)
  mc_sd <- sd(mc_vals, na.rm = TRUE)
  t_mean <- mean(t_vals, na.rm = TRUE)
  t_sd <- sd(t_vals, na.rm = TRUE)

  cat(sprintf("%-12s %5.3f(%3.0f) %5.3f(%3.0f) %5.1f(%4.1f)\n",
              name,
              ari_mean, 1000 * ari_sd,
              mc_mean, 1000 * mc_sd,
              t_mean, t_sd))
}
cat("(parentheses = SD x 1000 for ARI/misclass, SD for time)\n\n")

# --- Save combined result ---
outfile <- sprintf("output/%s_%drep.rds", sc_tag, n_reps)
saveRDS(all_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
