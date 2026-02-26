## Aggregate results across replicates
##
## Usage: Rscript code/03_aggregate.R [results_dir]
##
## Reads all output/results_*.rds files and produces a summary table.

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "output"

files <- list.files(results_dir, pattern = "^results_.*\\.rds$",
                    full.names = TRUE)

if (length(files) == 0) {
  stop("No result files found in ", results_dir)
}

cat(sprintf("Found %d result files\n", length(files)))

all_results <- lapply(files, readRDS)

# Collect model names from first file
model_names <- names(all_results[[1]])

# Build summary data frame
rows <- list()
for (mname in model_names) {
  psi_pms <- vapply(all_results, function(r) r[[mname]]$psi_pm, double(1))
  psi_pmses <- vapply(all_results, function(r) r[[mname]]$psi_pmse, double(1))
  acc_preds <- vapply(all_results, function(r) r[[mname]]$acc_pred, double(1))
  acc_pms <- vapply(all_results, function(r) r[[mname]]$acc_pm, double(1))
  times <- vapply(all_results, function(r) r[[mname]]$elapsed, double(1))

  rows[[mname]] <- data.frame(
    model = mname,
    psi_pm_mean = mean(psi_pms),
    psi_pm_sd = sd(psi_pms),
    psi_pmse_mean = mean(psi_pmses),
    acc_pred_mean = mean(acc_preds),
    acc_pm_mean = mean(acc_pms),
    time_mean = mean(times, na.rm = TRUE),
    n_reps = length(all_results),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, rows)
print(summary_df, row.names = FALSE)

outfile <- file.path(results_dir, "summary.csv")
write.csv(summary_df, outfile, row.names = FALSE)
cat("\nSaved to", outfile, "\n")
