## Aggregate all Gaussian scenario results into a final summary CSV
##
## Usage: Rscript code/gaussian/03_aggregate_all.R [output_dir]
##
## Reads all output/gaussian/gauss_*_Nrep.rds files and produces
## output/gaussian/summary_gaussian.csv with one row per (scenario, model).

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "output/gaussian"

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

  # Parse scenario tag to extract k, beta, sigma, grid
  # Format: gauss_k<K>_b<BETA>_s<SIGMA>_g<GRID>
  parts <- regmatches(sc_tag, regexec("gauss_k(\\d+)_b(.+)_s(.+)_g(\\d+)", sc_tag))[[1]]
  k_val <- as.integer(parts[2])
  beta_val <- as.numeric(gsub("_", ".", parts[3]))
  sigma_val <- as.numeric(gsub("_", ".", parts[4]))
  grid_val <- as.integer(parts[5])

  for (name in model_names) {
    vals <- lapply(all_metrics, function(r) r[[name]])

    safe_extract <- function(v, field) {
      x <- v[[field]]
      if (is.null(x)) NA_real_ else x
    }

    ari_vals <- vapply(vals, function(v) v$ari, double(1))
    mc_vals <- vapply(vals, function(v) v$misclass, double(1))
    psi_pm_vals <- vapply(vals, function(v) v$psi_pm, double(1))
    rhat_vals <- vapply(vals, function(v) safe_extract(v, "rhat_psi"), double(1))
    eq_pm_vals <- vapply(vals, function(v) v$eq_pm, double(1))
    eq_true_vals <- vapply(vals, function(v) v$eq_true, double(1))
    eq_pmse_vals <- vapply(vals, function(v) v$eq_pmse, double(1))
    brier_vals <- vapply(vals, function(v) v$brier, double(1))
    accept_psi_vals <- vapply(vals, function(v) safe_extract(v, "accept_psi"), double(1))
    accept_graph_vals <- vapply(vals, function(v) safe_extract(v, "accept_graph"), double(1))
    t_vals <- vapply(vals, function(v) v$elapsed, double(1))
    t_pre_vals <- vapply(vals, function(v) safe_extract(v, "elapsed_pre"), double(1))

    rows[[length(rows) + 1L]] <- data.frame(
      scenario = sc_tag,
      k = k_val,
      beta = beta_val,
      sigma = sigma_val,
      grid = grid_val,
      model = name,
      n_reps = n_reps,
      ARI_mean = mean(ari_vals, na.rm = TRUE),
      ARI_sd = sd(ari_vals, na.rm = TRUE),
      misclass_mean = mean(mc_vals, na.rm = TRUE),
      misclass_sd = sd(mc_vals, na.rm = TRUE),
      psi_pm_mean = mean(psi_pm_vals, na.rm = TRUE),
      psi_pm_sd = sd(psi_pm_vals, na.rm = TRUE),
      rhat_psi_mean = mean(rhat_vals, na.rm = TRUE),
      rhat_psi_max = max(rhat_vals, na.rm = TRUE),
      eq_pm_mean = mean(eq_pm_vals, na.rm = TRUE),
      eq_true_mean = mean(eq_true_vals, na.rm = TRUE),
      eq_pmse_mean = mean(eq_pmse_vals, na.rm = TRUE),
      eq_pmse_sd = sd(eq_pmse_vals, na.rm = TRUE),
      brier_mean = mean(brier_vals, na.rm = TRUE),
      brier_sd = sd(brier_vals, na.rm = TRUE),
      accept_psi_mean = mean(accept_psi_vals, na.rm = TRUE),
      accept_graph_mean = mean(accept_graph_vals, na.rm = TRUE),
      time_mean = mean(t_vals, na.rm = TRUE),
      time_sd = sd(t_vals, na.rm = TRUE),
      time_pre_mean = mean(t_pre_vals, na.rm = TRUE),
      time_pre_sd = sd(t_pre_vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, rows)

cat(sprintf("\n%d rows (%d scenarios x %d models)\n",
            nrow(summary_df), length(files), length(model_names)))
print(summary_df[, c("scenario", "k", "beta", "sigma", "grid", "model",
                      "ARI_mean", "misclass_mean", "eq_pmse_mean",
                      "brier_mean", "time_mean")],
      row.names = FALSE)

outfile <- file.path(output_dir, "summary_gaussian.csv")
write.csv(summary_df, outfile, row.names = FALSE)
cat(sprintf("\nSaved %s\n", outfile))
