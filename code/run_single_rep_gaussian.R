## Fit all 4 methods for a single replicate of one Gaussian scenario
##
## Usage: Rscript code/run_single_rep_gaussian.R <scenario_index> <rep>
##
## Saves per-rep results to:
##   output/<scenario_tag>/rep_<NNN>.rds

source(file.path("code", "helpers.R"))
source(file.path("code", "helpers_gaussian.R"))

# --- Parse args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript code/run_single_rep_gaussian.R <scenario_index> <rep>")
}
sc_idx <- as.integer(args[1])
rep <- as.integer(args[2])

# --- Fixed parameters ---
grid_rows <- 100L
grid_cols <- 100L
n_iter <- 5000L
burnin <- 1000L

model_names <- c("mdgm_st", "mdgm_ao", "mrf_pl", "bis_smc")

# --- Define 18 scenarios ---
# k=2 (Ising): 6 scenarios
# k=3 (Potts): 6 scenarios
# k=6 (Potts): 6 scenarios
scenarios <- list()

# k=2: beta = 0.22, 0.44, 0.66, 0.88, 1.10, 1.32
beta_k2 <- c(0.22, 0.44, 0.66, 0.88, 1.10, 1.32)
for (beta in beta_k2) {
  tag <- sprintf("gauss_k2_b%s_g100", gsub("\\.", "_", sprintf("%.2f", beta)))
  scenarios[[tag]] <- list(
    k = 2L, beta = beta,
    mu = c(-1, 1), sigma2 = c(1, 1)
  )
}

# k=3: beta = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2
beta_k3 <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (beta in beta_k3) {
  tag <- sprintf("gauss_k3_b%s_g100", gsub("\\.", "_", sprintf("%.1f", beta)))
  scenarios[[tag]] <- list(
    k = 3L, beta = beta,
    mu = c(-1, 0, 1), sigma2 = c(1, 1, 1)
  )
}

# k=6: beta = 0.3, 0.5, 0.7, 0.9, 1.1, 1.3
beta_k6 <- c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3)
for (beta in beta_k6) {
  tag <- sprintf("gauss_k6_b%s_g100", gsub("\\.", "_", sprintf("%.1f", beta)))
  scenarios[[tag]] <- list(
    k = 6L, beta = beta,
    mu = c(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5), sigma2 = rep(1, 6)
  )
}

if (sc_idx < 1L || sc_idx > length(scenarios)) {
  stop(sprintf("Scenario index must be between 1 and %d", length(scenarios)))
}

sc_name <- names(scenarios)[sc_idx]
sc <- scenarios[[sc_idx]]
k <- sc$k
beta_true <- sc$beta
mu_true <- sc$mu
sigma2_true <- sc$sigma2

cat(sprintf("[%s] rep %d: k=%d, beta=%.2f, mu=(%s), sigma2=(%s)\n",
            sc_name, rep, k, beta_true,
            paste(mu_true, collapse = ","),
            paste(sigma2_true, collapse = ",")))

# --- Auxiliary sweeps for data generation ---
get_n_sw_sweeps <- function(beta) {
  if (beta <= 0.4) return(200L)
  if (beta <= 0.8) return(300L)
  return(500L)
}

# --- Model config factory ---
make_model_configs <- function(beta_true, k) {
  n_aux <- get_n_sw_sweeps(beta_true)
  list(
    mdgm_st = list(
      spatial = mdgm(dag_type = "spanning_tree"),
      psi_tune = 0.4, psi_init = beta_true * 2
    ),
    mdgm_ao = list(
      spatial = mdgm(dag_type = "acyclic_orientation"),
      psi_tune = 0.3, psi_init = beta_true * 1.25
    ),
    mrf_pl = list(
      spatial = mrf(method = "pseudo_likelihood"),
      psi_tune = 0.15, psi_init = beta_true
    )
  )
}

model_configs <- make_model_configs(beta_true, k)

# --- Generate data (reproducible, non-sequential seed) ---
data_seed <- as.integer((sc_idx * 100003L + rep * 1000033L) %% .Machine$integer.max)
nug <- nug_from_grid(grid_rows, grid_cols, seed = data_seed)

# Generate true latent field via Swendsen-Wang
n_sw_sweeps <- get_n_sw_sweeps(beta_true)
z_true <- sample_mrf(nug, psi = beta_true, n_colors = k,
                      method = "swendsen_wang", n_sweeps = n_sw_sweeps,
                      seed = data_seed)

# Generate Gaussian observations
set.seed(data_seed)
y <- generate_gaussian_obs(z_true, mu_true, sigma2_true)
y_vec <- vapply(y, `[`, double(1), 1L)  # flat vector for bayesImageS
n <- grid_rows * grid_cols

z_init <- sample(0L:(k - 1L), n, replace = TRUE)

# theta_init: concatenate means then variances (Gaussian format)
theta_init <- c(mu_true, sigma2_true)

# --- Fit mdgm/mrf methods ---
rep_metrics <- list()

for (name in c("mdgm_st", "mdgm_ao", "mrf_pl")) {
  cfg <- model_configs[[name]]
  model <- srf_model(nug, spatial = cfg$spatial, emission = "gaussian",
                     n_colors = k)

  t0 <- proc.time()
  fit <- mcmc(model, y = y, z_init = z_init,
              psi_init = cfg$psi_init, theta_init = theta_init,
              n_iter = n_iter, psi_tune = cfg$psi_tune,
              seed = as.integer((data_seed + match(name, model_names) * 999983L) %% .Machine$integer.max),
              nug = nug)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  rep_metrics[[name]] <- compute_metrics_gaussian(
    fit, z_true, nug, burnin, elapsed, k
  )
  cat(sprintf("  %s: ARI=%.3f misclass=%.3f t=%.0fs",
              name, rep_metrics[[name]]$ari,
              rep_metrics[[name]]$misclass, elapsed))
}

# --- Fit bayesImageS smcPotts ---
t0 <- proc.time()
smc_result <- tryCatch(
  fit_bayesImageS_smc(y_vec, grid_rows, grid_cols, k),
  error = function(e) {
    cat(sprintf("\n  bis_smc ERROR: %s", e$message))
    NULL
  }
)
elapsed <- (proc.time() - t0)[["elapsed"]]

if (!is.null(smc_result)) {
  rep_metrics[["bis_smc"]] <- compute_metrics_bayesImageS(
    smc_result, z_true, elapsed, k
  )
  cat(sprintf("  bis_smc: ARI=%.3f misclass=%.3f t=%.0fs",
              rep_metrics[["bis_smc"]]$ari,
              rep_metrics[["bis_smc"]]$misclass, elapsed))
} else {
  rep_metrics[["bis_smc"]] <- list(
    ari = NA_real_, misclass = NA_real_,
    psi_pm = NA_real_, rhat_psi = NA_real_,
    eq_pm = NA_real_, eq_true = NA_real_,
    accept_psi = NA_real_, accept_graph = NA_real_,
    elapsed = elapsed
  )
}
cat("\n")

# --- Save per-rep result ---
outdir <- file.path("output", sc_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("rep_%03d.rds", rep))
saveRDS(rep_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
