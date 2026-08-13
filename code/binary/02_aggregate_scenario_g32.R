## Aggregate per-rep files for a single 32x32 scenario
##
## Usage: Rscript code/binary/02_aggregate_scenario_g32.R <scenario_tag>
##
## Reads all output/<tag>/rep_*.rds files and saves combined
## output/<tag>_<N>rep.rds. No TV distance (no mrf_exact reference).

source(file.path("code", "common", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/binary/02_aggregate_scenario_g32.R <scenario_tag>")
}
sc_tag <- args[1]

# --- Read per-rep files ---
rep_dir <- file.path("output", "binary", sc_tag)
rep_files <- sort(list.files(rep_dir, pattern = "^rep_.*\\.rds$",
                             full.names = TRUE))

if (length(rep_files) == 0) {
  stop(sprintf("No rep files found in %s", rep_dir))
}

n_reps <- length(rep_files)
cat(sprintf("Aggregating %d replicates for %s\n", n_reps, sc_tag))

model_names <- c("mrf_pl", "mdgm_st")
all_metrics <- vector("list", n_reps)

for (i in seq_along(rep_files)) {
  all_metrics[[i]] <- readRDS(rep_files[i])
}

# --- Print summary table ---
cat(sprintf("\n--- %s: aggregated across %d replicates ---\n", sc_tag, n_reps))
cat(sprintf("%-12s %10s %10s %10s %10s %10s\n",
            "Model", "psi_pm", "th1_pm", "th2_pm", "acc_pred", "time"))
cat(paste(rep("-", 62), collapse = ""), "\n")

for (name in model_names) {
  vals <- lapply(all_metrics, function(r) r[[name]])
  cat(sprintf("%-12s %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.1f(%4.1f)\n",
              name,
              mean(vapply(vals, function(v) v$psi_pm, double(1))),
              1000 * sd(vapply(vals, function(v) v$psi_pm, double(1))),
              mean(vapply(vals, function(v) v$theta_pm[1], double(1))),
              1000 * sd(vapply(vals, function(v) v$theta_pm[1], double(1))),
              mean(vapply(vals, function(v) v$theta_pm[2], double(1))),
              1000 * sd(vapply(vals, function(v) v$theta_pm[2], double(1))),
              mean(vapply(vals, function(v) v$acc_pred, double(1))),
              1000 * sd(vapply(vals, function(v) v$acc_pred, double(1))),
              mean(vapply(vals, function(v) v$elapsed, double(1))),
              sd(vapply(vals, function(v) v$elapsed, double(1)))))
}
cat("(parentheses = SD x 1000 except time = SD)\n\n")

# --- Save combined result ---
outfile <- sprintf("output/binary/%s_%drep.rds", sc_tag, n_reps)
saveRDS(all_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
