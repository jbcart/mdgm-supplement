## Fit all 4 models for a single replicate of a single scenario
##
## Usage: Rscript code/binary/01_run_single_rep.R <scenario_index> <rep>
##
## Saves per-rep results (without TV distance) to:
##   output/<scenario_tag>/rep_<NNN>.rds

source(file.path("code", "common", "helpers.R"))

# --- Parse args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript code/binary/01_run_single_rep.R <scenario_index> <rep>")
}
sc_idx <- as.integer(args[1])
rep <- as.integer(args[2])

# --- Fixed parameters (must match run_all_g16.R) ---
grid_rows <- 16L
grid_cols <- 16L
n_iter <- 5000L
burnin <- 1000L

model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")

get_n_aux_sweeps <- function(psi) {
  if (psi <= 0.1) return(200L)
  if (psi <= 0.3) return(300L)
  if (psi <= 0.6) return(400L)
  return(500L)
}

# --- Define all 32 scenarios (same order as run_all_g16.R) ---
psi_vec <- seq(0.1, 0.8, by = 0.1)
scenarios <- list()

for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_2r2o1g16", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.2, 0.8), n_reps = 2L, lambda = 0)
}
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_05r2o1g16", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.05, 0.95), n_reps = 2L, lambda = 0)
}
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_1l1_39o1g16", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.1, 0.9), n_reps = 0L, lambda = 1.39)
}
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_1l2_3o1g16", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.1, 0.9), n_reps = 0L, lambda = 2.3)
}

if (sc_idx < 1 || sc_idx > length(scenarios)) {
  stop(sprintf("Scenario index must be between 1 and %d", length(scenarios)))
}

sc_name <- names(scenarios)[sc_idx]
sc <- scenarios[[sc_idx]]
psi_true <- sc$psi
theta_true <- sc$theta
n_aux_sweeps <- get_n_aux_sweeps(psi_true)

cat(sprintf("[%s] rep %d: psi=%.1f, theta=(%.2f,%.2f), lambda=%.2f, M=%d\n",
            sc_name, rep, psi_true, theta_true[1], theta_true[2],
            sc$lambda, n_aux_sweeps))

# --- Model config factory (same as run_all_g16.R) ---
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

model_configs <- make_model_configs(psi_true, n_aux_sweeps)

# --- Generate data (reproducible, non-sequential seed) ---
# Hash scenario index and rep to get well-separated seeds
data_seed <- as.integer((sc_idx * 100003L + rep * 1000033L) %% .Machine$integer.max)
nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)
z_true <- sample_mrf(nug, psi = psi_true, seed = data_seed)
set.seed(data_seed)
y <- generate_observations(z_true, theta_true, sc$n_reps, lambda = sc$lambda)
n <- grid_rows * grid_cols
missing_sites <- vapply(y, function(yi) length(yi) == 0, logical(1))

z_init <- sample(0:1, n, replace = TRUE)

# --- Fit all 4 models ---
rep_metrics <- list()
for (name in model_names) {
  cfg <- model_configs[[name]]
  model <- srf_model(nug, spatial = cfg$spatial, emission = "bernoulli")

  t0 <- proc.time()
  fit <- mcmc(model, y = y, z_init = z_init,
              psi_init = cfg$psi_init, theta_init = theta_true,
              n_iter = n_iter, psi_tune = cfg$psi_tune,
              store_z = TRUE,
              seed = as.integer((data_seed + match(name, model_names) * 999983L) %% .Machine$integer.max),
              nug = nug)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  rep_metrics[[name]] <- compute_metrics(
    fit, z_true, psi_true, theta_true, nug, burnin, elapsed,
    missing_sites = missing_sites
  )
  cat(sprintf("  %s=%.0fs", name, elapsed))
}
cat("\n")

# --- Save per-rep result (TV distance computed at aggregation) ---
outdir <- file.path("output", "binary", sc_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("rep_%03d.rds", rep))
saveRDS(rep_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
