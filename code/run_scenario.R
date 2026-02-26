## Simulation scenario: 100 replicates on a 16x16 grid with Bernoulli emissions
## Each replicate generates fresh exact MRF data and fits 4 models.

source(file.path("code", "helpers.R"))

# --- Scenario parameters ---
grid_rows <- 16L
grid_cols <- 16L
psi_true <- 0.5
theta_true <- c(0.2, 0.8)
n_reps <- 2L
n_iter <- 5000L
burnin <- 1000L
n_replicates <- 100L

model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")

cat("=== mdgm simulation scenario ===\n")
cat(sprintf("Grid: %dx%d (%d sites)\n", grid_rows, grid_cols,
            grid_rows * grid_cols))
cat(sprintf("psi = %.2f, theta = (%.2f, %.2f), reps = %d\n",
            psi_true, theta_true[1], theta_true[2], n_reps))
cat(sprintf("MCMC: %d iterations, %d burnin, %d replicates\n\n",
            n_iter, burnin, n_replicates))

# --- Storage ---
# Each element: list of metrics per model per replicate
all_metrics <- vector("list", n_replicates)

# --- Replication loop ---
for (rep in seq_len(n_replicates)) {
  cat(sprintf("--- Replicate %d/%d ---\n", rep, n_replicates))

  # Fresh data per replicate (different seed each time)
  data_seed <- rep
  nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)
  z_true <- sample_mrf(nug, psi = psi_true, seed = data_seed)
  y <- generate_observations(z_true, theta_true, n_reps)

  # Model configurations
  model_configs <- list(
    mrf_exact = list(
      spatial = mrf(method = "exchange"),
      psi_tune = 0.15, psi_init = psi_true
    ),
    mrf_pl = list(
      spatial = mrf(method = "pseudo_likelihood"),
      psi_tune = 0.15, psi_init = psi_true
    ),
    mdgm_st = list(
      spatial = mdgm(dag_type = "spanning_tree"),
      psi_tune = 0.4, psi_init = psi_true * 2
    ),
    mdgm_ao = list(
      spatial = mdgm(dag_type = "acyclic_orientation"),
      psi_tune = 0.3, psi_init = psi_true * 1.25
    )
  )

  # Random z_init (matching original sim_study.r)
  n <- grid_rows * grid_cols
  z_init <- sample(0:1, n, replace = TRUE)

  rep_metrics <- list()
  for (name in model_names) {
    cfg <- model_configs[[name]]
    cat(sprintf("  %s...", name))

    model <- srf_model(nug, spatial = cfg$spatial, emission = "bernoulli")

    t0 <- proc.time()
    fit <- mcmc(model, y = y, z_init = z_init,
                psi_init = cfg$psi_init, theta_init = theta_true,
                n_iter = n_iter, psi_tune = cfg$psi_tune,
                seed = rep * 1000L + match(name, model_names),
                nug = nug)
    elapsed <- (proc.time() - t0)[["elapsed"]]

    rep_metrics[[name]] <- compute_metrics(
      fit, z_true, psi_true, theta_true, nug, burnin, elapsed
    )
    cat(sprintf(" %.1fs", elapsed))
  }

  # Compute TV distances vs mrf_exact for this replicate
  ref_eq <- rep_metrics[["mrf_exact"]]$eq_vec
  for (name in model_names) {
    rep_metrics[[name]]$tv_vs_mrf <- tv_distance(
      ref_eq, rep_metrics[[name]]$eq_vec
    )
  }

  all_metrics[[rep]] <- rep_metrics
  cat("\n")
}

# --- Aggregate across replicates ---
cat("\n=== Aggregated results (mean +/- sd across", n_replicates, "replicates) ===\n\n")

header <- sprintf("%-12s %10s %10s %10s %10s %10s %10s %10s",
                  "Model", "psi_pm", "th1_pm", "th2_pm",
                  "acc_pred", "TV_dist", "psi_ar", "time")
cat(header, "\n")
cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

for (name in model_names) {
  vals <- lapply(all_metrics, function(r) r[[name]])

  psi_pms   <- vapply(vals, function(v) v$psi_pm, double(1))
  th1_pms   <- vapply(vals, function(v) v$theta_pm[1], double(1))
  th2_pms   <- vapply(vals, function(v) v$theta_pm[2], double(1))
  acc_preds <- vapply(vals, function(v) v$acc_pred, double(1))
  tvds      <- vapply(vals, function(v) v$tv_vs_mrf, double(1))
  psi_ars   <- vapply(vals, function(v) v$accept_psi, double(1))
  times     <- vapply(vals, function(v) v$elapsed, double(1))

  cat(sprintf("%-12s %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.3f(%3.0f) %5.1f(%4.1f)\n",
              name,
              mean(psi_pms), 1000 * sd(psi_pms),
              mean(th1_pms), 1000 * sd(th1_pms),
              mean(th2_pms), 1000 * sd(th2_pms),
              mean(acc_preds), 1000 * sd(acc_preds),
              mean(tvds), 1000 * sd(tvds),
              mean(psi_ars), 1000 * sd(psi_ars),
              mean(times), sd(times)))
}

cat("\nTrue values: psi =", psi_true,
    ", theta = (", paste(theta_true, collapse = ", "), ")\n")
cat("(values in parentheses are SD x 1000)\n")

# --- Save full results ---
outfile <- sprintf("output/scenario_%dx%d_%drep.rds",
                   grid_rows, grid_cols, n_replicates)
saveRDS(all_metrics, outfile)
cat("\nFull results saved to", outfile, "\n")
