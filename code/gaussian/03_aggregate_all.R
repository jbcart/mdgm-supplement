## Aggregate all Gaussian scenario results into a final summary CSV
##
## Usage: Rscript code/gaussian/03_aggregate_all.R [output_dir]
##
## Reads all output/gauss_*_Nrep.rds files and produces
## output/summary_gaussian.csv with one row per (scenario, model).

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "output"

files <- sort(list.files(output_dir, pattern = "^gauss_.*_\\d+rep\\.rds$",
                         full.names = TRUE))

if (length(files) == 0) {
  stop(sprintf("No gauss_*_Nrep.rds files found in %s", output_dir))
}

cat(sprintf("Found %d Gaussian scenario files\n", length(files)))

model_names <- c("mdgm_st", "mdgm_ao", "mrf_pl", "bis_pfab")
rows <- list()

for (f in files) {
  sc_tag <- sub("_\\d+rep\\.rds$", "", basename(f))
  all_metrics <- readRDS(f)
  n_reps <- length(all_metrics)

  # Parse scenario tag to extract k and beta
  # Format: gauss_k<K>_b<BETA>_g100
  parts <- regmatches(sc_tag, regexec("gauss_k(\\d+)_b(.+)_g100", sc_tag))[[1]]
  k_val <- as.integer(parts[2])
  beta_val <- as.numeric(gsub("_", ".", parts[3]))

  for (name in model_names) {
    vals <- lapply(all_metrics, function(r) r[[name]])

    ari_vals <- vapply(vals, function(v) v$ari, double(1))
    mc_vals <- vapply(vals, function(v) v$misclass, double(1))
    t_vals <- vapply(vals, function(v) v$elapsed, double(1))

    rows[[length(rows) + 1L]] <- data.frame(
      scenario = sc_tag,
      k = k_val,
      beta = beta_val,
      model = name,
      n_reps = n_reps,
      ARI_mean = mean(ari_vals, na.rm = TRUE),
      ARI_sd = sd(ari_vals, na.rm = TRUE),
      misclass_mean = mean(mc_vals, na.rm = TRUE),
      misclass_sd = sd(mc_vals, na.rm = TRUE),
      time_mean = mean(t_vals, na.rm = TRUE),
      time_sd = sd(t_vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, rows)

cat(sprintf("\n%d rows (%d scenarios x %d models)\n",
            nrow(summary_df), length(files), length(model_names)))
print(summary_df[, c("scenario", "k", "beta", "model",
                      "ARI_mean", "misclass_mean", "time_mean")],
      row.names = FALSE)

outfile <- file.path(output_dir, "summary_gaussian.csv")
write.csv(summary_df, outfile, row.names = FALSE)
cat(sprintf("\nSaved %s\n", outfile))
