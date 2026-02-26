## Quick single-scenario runner for verification
## Fits 4 models on a 16x16 grid with Bernoulli emissions

source(file.path("code", "helpers.R"))

# --- Scenario parameters ---
grid_rows <- 16L
grid_cols <- 16L
psi_true <- 0.5
theta_true <- c(0.2, 0.8)
n_reps <- 2L
n_iter <- 5000L
burnin <- 1000L
data_seed <- 42L

cat("=== mdgm simulation scenario ===\n")
cat(sprintf("Grid: %dx%d (%d sites)\n", grid_rows, grid_cols,
            grid_rows * grid_cols))
cat(sprintf("psi = %.2f, theta = (%.2f, %.2f), reps = %d\n",
            psi_true, theta_true[1], theta_true[2], n_reps))
cat(sprintf("MCMC: %d iterations, %d burnin\n\n", n_iter, burnin))

# --- Generate data ---
cat("Generating data...\n")
nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)
z_true <- sample_mrf(nug, psi = psi_true, seed = data_seed)
y <- generate_observations(z_true, theta_true, n_reps)

cat(sprintf("  Sites with z=0: %d, z=1: %d\n",
            sum(z_true == 0), sum(z_true == 1)))
cat(sprintf("  Sufficient stat (same-color edges): %d\n",
            sufficient_stat(z_true, nug)))

# --- Model configurations ---
models <- list(
  mrf_exact = list(
    spatial = mrf(method = "exchange"),
    psi_tune = 0.15,
    psi_init = psi_true,
    seed = 1L
  ),
  mrf_pl = list(
    spatial = mrf(method = "pseudo_likelihood"),
    psi_tune = 0.15,
    psi_init = psi_true,
    seed = 2L
  ),
  mdgm_st = list(
    spatial = mdgm(dag_type = "spanning_tree"),
    psi_tune = 0.4,
    psi_init = psi_true * 2,
    seed = 3L
  ),
  mdgm_ao = list(
    spatial = mdgm(dag_type = "acyclic_orientation"),
    psi_tune = 0.3,
    psi_init = psi_true * 1.25,
    seed = 4L
  )
)

# --- Fit models ---
results <- list()
metrics <- list()

for (name in names(models)) {
  cfg <- models[[name]]
  cat(sprintf("Fitting %s...", name))

  model <- srf_model(nug, spatial = cfg$spatial, emission = "bernoulli")

  t0 <- proc.time()
  fit <- mcmc(model,
              y = y,
              z_init = z_true,
              psi_init = cfg$psi_init,
              theta_init = theta_true,
              n_iter = n_iter,
              psi_tune = cfg$psi_tune,
              seed = cfg$seed,
              nug = nug)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  results[[name]] <- fit
  metrics[[name]] <- compute_metrics(
    fit, z_true, psi_true, theta_true, nug, burnin, elapsed
  )

  cat(sprintf(" done (%.1fs)\n", elapsed))
}

# --- Print comparison table ---
cat("\n=== Results ===\n\n")

header <- sprintf("%-12s %8s %8s %8s %8s %8s %8s %8s %8s %8s",
                  "Model", "psi_pm", "psi_psd", "th1_pm", "th2_pm",
                  "acc_pred", "acc_pm", "eq_pm", "psi_ar", "time")
cat(header, "\n")
cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

for (name in names(metrics)) {
  m <- metrics[[name]]
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.1f %8.3f %8.1f\n",
              name,
              m$psi_pm, m$psi_psd,
              m$theta_pm[1], m$theta_pm[2],
              m$acc_pred, m$acc_pm,
              m$eq_pm,
              m$accept_psi,
              m$elapsed))
}

cat("\nTrue values: psi =", psi_true,
    ", theta = (", paste(theta_true, collapse = ", "), ")\n")

# --- R-hat diagnostics ---
cat("\n=== Diagnostics (split R-hat) ===\n\n")
cat(sprintf("%-12s %8s %8s %8s\n", "Model", "psi", "theta_1", "theta_2"))
for (name in names(metrics)) {
  m <- metrics[[name]]
  cat(sprintf("%-12s %8.3f %8.3f %8.3f\n",
              name, m$rhat_psi, m$rhat_theta[1], m$rhat_theta[2]))
}
