## Fit 2 models (mdgm_st, mrf_pl) for a single replicate of a 32x32 scenario
##
## Usage: Rscript code/binary/01_run_single_rep_g32.R <scenario_index> <rep>
##
## 24 scenarios (see below); 100 reps each.
## Saves per-rep results to:
##   output/<scenario_tag>/rep_<NNN>.rds

source(file.path("code", "common", "helpers.R"))

# --- Parse args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript code/binary/01_run_single_rep_g32.R <scenario_index> <rep>")
}
sc_idx <- as.integer(args[1])
rep <- as.integer(args[2])

# --- Fixed parameters ---
grid_rows <- 32L
grid_cols <- 32L
n_iter <- 5000L
burnin <- 1000L

model_names <- c("mrf_pl", "mdgm_st")

# --- Define 16 scenarios ---
# Complete data: eta=0.2, m=2 (indices 1-8)
# Missing data:  eta=0.1, lambda=1.39 (indices 9-16)
psi_vec <- seq(0.1, 0.8, by = 0.1)
scenarios <- list()

for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_2r2o1g32", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.2, 0.8), n_reps = 2L, lambda = 0)
}
for (psi in psi_vec) {
  tag <- sprintf("p0_%de0_1l1_39o1g32", round(psi * 10))
  scenarios[[tag]] <- list(psi = psi, theta = c(0.1, 0.9), n_reps = 0L, lambda = 1.39)
}

if (sc_idx < 1L || sc_idx > length(scenarios)) {
  stop(sprintf("Scenario index must be between 1 and %d", length(scenarios)))
}

sc_name <- names(scenarios)[sc_idx]
sc <- scenarios[[sc_idx]]
psi_true <- sc$psi
theta_true <- sc$theta

cat(sprintf("[%s] rep %d: psi=%.1f, theta=(%.2f,%.2f), lambda=%.2f\n",
            sc_name, rep, psi_true, theta_true[1], theta_true[2], sc$lambda))

# --- Model config factory ---
make_model_configs <- function(psi_true) {
  list(
    mrf_pl = list(
      spatial = mrf(method = "pseudo_likelihood"),
      psi_tune = 0.15, psi_init = psi_true
    ),
    mdgm_st = list(
      spatial = mdgm(dag_type = "spanning_tree"),
      psi_tune = 0.4, psi_init = psi_true * 2
    )
  )
}

model_configs <- make_model_configs(psi_true)

# --- Generate data (reproducible, non-sequential seed) ---
data_seed <- as.integer((sc_idx * 100003L + rep * 1000033L) %% .Machine$integer.max)
nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)
z_true <- sample_mrf(nug, psi = psi_true, seed = data_seed)
set.seed(data_seed)
y <- generate_observations(z_true, theta_true, sc$n_reps, lambda = sc$lambda)
n <- grid_rows * grid_cols
missing_sites <- vapply(y, function(yi) length(yi) == 0, logical(1))

z_init <- sample(0:1, n, replace = TRUE)

# --- Fit models ---
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

  rm(fit, model); gc()
}
cat("\n")

# --- Save per-rep result ---
outdir <- file.path("output", sc_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("rep_%03d.rds", rep))
saveRDS(rep_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
