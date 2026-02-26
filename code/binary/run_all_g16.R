## Full simulation study: all 32 g16 scenarios from batch_master
## Each scenario: 100 replicates, 4 models (mrf_exact, mrf_pl, mdgm_st, mdgm_ao)
##
## Usage: Rscript code/binary/run_all_g16.R
##   or:  Rscript code/binary/run_all_g16.R <scenario_index>   # run a single scenario (1-32)

source(file.path("code", "common", "helpers.R"))

# --- Fixed parameters ---
grid_rows <- 16L
grid_cols <- 16L
n_iter <- 5000L
burnin <- 1000L
n_replicates <- 100L

model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")

# n_aux_sweeps for exchange algorithm, scales with psi (from simulation_setup.zsh)
get_n_aux_sweeps <- function(psi) {
  if (psi <= 0.1) return(200L)
  if (psi <= 0.3) return(300L)
  if (psi <= 0.6) return(400L)
  return(500L)
}

# --- Define all 32 scenarios ---
psi_vec <- seq(0.1, 0.8, by = 0.1)

scenarios <- list()

# Complete data: e0_2 (theta = 0.2, 0.8), 2 reps/site
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_2r2o1g16", round(psi * 10))
  scenarios[[tag]] <- list(
    psi = psi, theta = c(0.2, 0.8), n_reps = 2L, lambda = 0
  )
}

# Complete data: e0_05 (theta = 0.05, 0.95), 2 reps/site
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_05r2o1g16", round(psi * 10))
  scenarios[[tag]] <- list(
    psi = psi, theta = c(0.05, 0.95), n_reps = 2L, lambda = 0
  )
}

# Missing data: lambda = 1.39 (~25% missing), theta = 0.1, 0.9
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_1l1_39o1g16", round(psi * 10))
  scenarios[[tag]] <- list(
    psi = psi, theta = c(0.1, 0.9), n_reps = 0L, lambda = 1.39
  )
}

# Missing data: lambda = 2.3 (~10% missing), theta = 0.1, 0.9
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_1l2_3o1g16", round(psi * 10))
  scenarios[[tag]] <- list(
    psi = psi, theta = c(0.1, 0.9), n_reps = 0L, lambda = 2.3
  )
}

# --- Optionally run a single scenario by index ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  idx <- as.integer(args[1])
  if (idx < 1 || idx > length(scenarios)) {
    stop(sprintf("Scenario index must be between 1 and %d", length(scenarios)))
  }
  scenarios <- scenarios[idx]
  cat(sprintf("Running single scenario: %s (index %d)\n\n",
              names(scenarios), idx))
}

# --- Model config factory ---
make_model_configs <- function(psi_true, n_aux_sweeps) {
  list(
    mrf_exact = list(
      spatial = mrf(method = "exchange", n_aux_sweeps = n_aux_sweeps),
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
}

# --- Run all scenarios ---
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  psi_true <- sc$psi
  theta_true <- sc$theta
  n_aux_sweeps <- get_n_aux_sweeps(psi_true)

  cat(sprintf("====== Scenario: %s ======\n", sc_name))
  cat(sprintf("psi=%.1f, theta=(%.2f, %.2f), lambda=%.2f, reps=%d, M=%d\n",
              psi_true, theta_true[1], theta_true[2],
              sc$lambda, sc$n_reps, n_aux_sweeps))
  cat(sprintf("%d replicates x %d models\n\n", n_replicates, length(model_names)))

  model_configs <- make_model_configs(psi_true, n_aux_sweeps)
  all_metrics <- vector("list", n_replicates)

  for (rep in seq_len(n_replicates)) {
    cat(sprintf("  Rep %3d/%d:", rep, n_replicates))

    data_seed <- rep
    nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)
    z_true <- sample_mrf(nug, psi = psi_true, seed = data_seed)
    set.seed(data_seed)
    y <- generate_observations(z_true, theta_true, sc$n_reps,
                               lambda = sc$lambda)
    n <- grid_rows * grid_cols
    missing_sites <- vapply(y, function(yi) length(yi) == 0, logical(1))

    # Random z_init (matching original sim_study.r)
    z_init <- sample(0:1, n, replace = TRUE)

    rep_metrics <- list()
    for (name in model_names) {
      cfg <- model_configs[[name]]
      model <- srf_model(nug, spatial = cfg$spatial, emission = "bernoulli")

      t0 <- proc.time()
      fit <- mcmc(model, y = y, z_init = z_init,
                  psi_init = cfg$psi_init, theta_init = theta_true,
                  n_iter = n_iter, psi_tune = cfg$psi_tune,
                  seed = rep * 1000L + match(name, model_names),
                  nug = nug)
      elapsed <- (proc.time() - t0)[["elapsed"]]

      rep_metrics[[name]] <- compute_metrics(
        fit, z_true, psi_true, theta_true, nug, burnin, elapsed,
        missing_sites = missing_sites
      )
      cat(sprintf(" %s=%.0fs", name, elapsed))
    }

    # TV distances vs mrf_exact
    ref_eq <- rep_metrics[["mrf_exact"]]$eq_vec
    for (name in model_names) {
      rep_metrics[[name]]$tv_vs_mrf <- tv_distance(
        ref_eq, rep_metrics[[name]]$eq_vec
      )
    }

    all_metrics[[rep]] <- rep_metrics
    cat("\n")
  }

  # --- Print aggregated results for this scenario ---
  cat(sprintf("\n--- %s: aggregated across %d replicates ---\n",
              sc_name, n_replicates))

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
  cat(sprintf("True: psi=%.1f, theta=(%.2f, %.2f)\n",
              psi_true, theta_true[1], theta_true[2]))
  cat("(parentheses = SD x 1000)\n\n")

  # Save per-scenario results
  outfile <- sprintf("output/%s_%drep.rds", sc_name, n_replicates)
  saveRDS(all_metrics, outfile)
  cat(sprintf("Saved to %s\n\n", outfile))
}

cat("All scenarios complete.\n")
