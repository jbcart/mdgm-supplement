## Aggregate all scenario results into a final summary CSV
##
## Usage: Rscript code/binary/03_aggregate_all.R [output_dir]
##
## Reads all output/binary/*_100rep.rds (or *_Nrep.rds) files and produces
## output/binary/summary_all.csv with one row per (scenario, model).

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "output/binary"

files <- sort(list.files(output_dir, pattern = "_\\d+rep\\.rds$",
                         full.names = TRUE))

if (length(files) == 0) {
  stop(sprintf("No *_Nrep.rds files found in %s", output_dir))
}

cat(sprintf("Found %d scenario files\n", length(files)))

model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")
rows <- list()

for (f in files) {
  sc_tag <- sub("_\\d+rep\\.rds$", "", basename(f))
  all_metrics <- readRDS(f)
  n_reps <- length(all_metrics)

  for (name in model_names) {
    vals <- lapply(all_metrics, function(r) r[[name]])

    psi_pms   <- vapply(vals, function(v) v$psi_pm, double(1))
    th1_pms   <- vapply(vals, function(v) v$theta_pm[1], double(1))
    th2_pms   <- vapply(vals, function(v) v$theta_pm[2], double(1))
    acc_preds <- vapply(vals, function(v) v$acc_pred, double(1))
    tv_dists  <- vapply(vals, function(v) v$tv_vs_mrf, double(1))
    psi_pmses <- vapply(vals, function(v) v$psi_pmse, double(1))
    times     <- vapply(vals, function(v) v$elapsed, double(1))

    rows[[length(rows) + 1L]] <- data.frame(
      scenario = sc_tag,
      model = name,
      n_reps = n_reps,
      psi_pm_mean = mean(psi_pms),
      psi_pm_sd = sd(psi_pms),
      psi_pmse_mean = mean(psi_pmses),
      theta1_pm_mean = mean(th1_pms),
      theta1_pm_sd = sd(th1_pms),
      theta2_pm_mean = mean(th2_pms),
      theta2_pm_sd = sd(th2_pms),
      acc_pred_mean = mean(acc_preds),
      acc_pred_sd = sd(acc_preds),
      tv_dist_mean = mean(tv_dists),
      tv_dist_sd = sd(tv_dists),
      time_mean = mean(times),
      time_sd = sd(times),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, rows)

cat(sprintf("\n%d rows (%d scenarios x %d models)\n",
            nrow(summary_df), length(files), length(model_names)))
print(summary_df[, c("scenario", "model", "psi_pm_mean", "acc_pred_mean",
                      "tv_dist_mean", "time_mean")],
      row.names = FALSE)

outfile <- file.path(output_dir, "summary_all.csv")
write.csv(summary_df, outfile, row.names = FALSE)
cat(sprintf("\nSaved %s\n", outfile))
