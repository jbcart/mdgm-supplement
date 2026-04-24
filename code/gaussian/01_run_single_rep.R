## Fit all 4 methods for a single replicate of one Gaussian scenario
##
## Usage: Rscript code/gaussian/01_run_single_rep.R <scenario_index> <rep>
##
## 12 scenarios (see scenario table below); 20 reps each.
## Saves per-rep results to:
##   output/<scenario_tag>/rep_<NNN>.rds

source(file.path("code", "common", "helpers.R"))
source(file.path("code", "gaussian", "helpers.R"))

# --- Parse args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript code/gaussian/01_run_single_rep.R <scenario_index> <rep>")
}
sc_idx <- as.integer(args[1])
rep <- as.integer(args[2])

# --- Fixed parameters ---
n_iter <- 10000L
burnin <- 2000L

model_names <- c("mdgm_st", "mdgm_ao", "mrf_pl", "bis_pfab")

# --- Define 12 scenarios ---
# Tag format: gauss_k<K>_b<BETA>_s<SIGMA>_g<GRID>
# All scenarios use 100x100 grid, sigma=0.50 (sigma2=0.25).
# Beta values at ~0.35/0.65/0.95 × beta_c where beta_c = log(1 + sqrt(k)):
#   k=3: beta_c ~= 1.005 → 0.35, 0.65, 0.95
#   k=4: beta_c ~= 1.099 → 0.38, 0.71, 1.04
#   k=5: beta_c ~= 1.179 → 0.41, 0.76, 1.12
#   k=6: beta_c ~= 1.238 → 0.43, 0.80, 1.18
make_tag <- function(k, beta, sigma, grid) {
  sprintf("gauss_k%d_b%s_s%s_g%d",
          k,
          gsub("\\.", "_", sprintf("%.2f", beta)),
          gsub("\\.", "_", sprintf("%.2f", sigma)),
          grid)
}

scenarios <- list()

# k=3, 100x100, sigma=0.50 (indices 1-3)
for (beta in c(0.35, 0.65, 0.95)) {
  tag <- make_tag(3, beta, 0.50, 100)
  scenarios[[tag]] <- list(
    k = 3L, beta = beta, grid = 100L,
    mu = c(-1, 0, 1), sigma2 = rep(0.25, 3)
  )
}

# k=4, 100x100, sigma=0.50 (indices 4-6)
for (beta in c(0.38, 0.71, 1.04)) {
  tag <- make_tag(4, beta, 0.50, 100)
  scenarios[[tag]] <- list(
    k = 4L, beta = beta, grid = 100L,
    mu = c(-1, 0, 1, 2), sigma2 = rep(0.25, 4)
  )
}

# k=5, 100x100, sigma=0.50 (indices 7-9)
for (beta in c(0.41, 0.76, 1.12)) {
  tag <- make_tag(5, beta, 0.50, 100)
  scenarios[[tag]] <- list(
    k = 5L, beta = beta, grid = 100L,
    mu = c(-2, -1, 0, 1, 2), sigma2 = rep(0.25, 5)
  )
}

# k=6, 100x100, sigma=0.50 (indices 10-12)
for (beta in c(0.43, 0.80, 1.18)) {
  tag <- make_tag(6, beta, 0.50, 100)
  scenarios[[tag]] <- list(
    k = 6L, beta = beta, grid = 100L,
    mu = c(-2, -1, 0, 1, 2, 3), sigma2 = rep(0.25, 6)
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
grid_rows <- sc$grid
grid_cols <- sc$grid

cat(sprintf("[%s] rep %d: k=%d, beta=%.2f, mu=(%s), sigma2=(%s), grid=%dx%d\n",
            sc_name, rep, k, beta_true,
            paste(mu_true, collapse = ","),
            paste(sigma2_true, collapse = ","),
            grid_rows, grid_cols))

# --- Auxiliary sweeps for data generation ---
get_n_sw_sweeps <- function(beta) {
  if (beta <= 0.4) return(200L)
  if (beta <= 0.8) return(300L)
  return(500L)
}

# --- Model config factory ---
make_model_configs <- function(beta_true, k) {
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

# --- Priors for Gaussian emission parameters ---
# Match bayesImageS: mu_k ~ N(0, 100^2), sigma_k^2 ~ InvGamma(1.5, 1.5)
# bayesImageS parameterizes as Inv-Chi^2(nu=3, s^2=1) = InvGamma(3/2, 3/2)
gaussian_prior <- c(0, 10000, 1.5, 1.5)

# --- Fit mdgm/mrf methods ---
# Process one model at a time, extracting metrics and freeing the fit object
# to keep memory usage manageable.
rep_metrics <- list()

for (name in c("mdgm_st", "mdgm_ao", "mrf_pl")) {
  cfg <- model_configs[[name]]
  model <- srf_model(nug, spatial = cfg$spatial, emission = "gaussian",
                     n_colors = k)

  t0 <- proc.time()
  fit <- mcmc(model, y = y, z_init = z_init,
              psi_init = cfg$psi_init, theta_init = theta_init,
              n_iter = n_iter, psi_tune = cfg$psi_tune,
              emission_prior_params = gaussian_prior,
              store_z = TRUE,
              seed = as.integer((data_seed + match(name, model_names) * 999983L) %% .Machine$integer.max),
              nug = nug)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  rep_metrics[[name]] <- compute_metrics_gaussian(
    fit, z_true, nug, burnin, elapsed, k
  )
  cat(sprintf("  %s: ARI=%.3f misclass=%.3f brier=%.3f t=%.0fs",
              name, rep_metrics[[name]]$ari,
              rep_metrics[[name]]$misclass,
              rep_metrics[[name]]$brier, elapsed))

  rm(fit, model); gc()
}

# --- Fit bayesImageS PFAB (mcmcPotts with algorithm="aux") ---
# Load or compute PFAB surrogate parameters (cached per grid-size x k combo)
pfab_cache_dir <- file.path("output", "pfab_cache")
dir.create(pfab_cache_dir, showWarnings = FALSE, recursive = TRUE)
pfab_cache_file <- file.path(pfab_cache_dir,
                              sprintf("k%d_g%d.rds", k, grid_rows))

t0_pre <- proc.time()
if (file.exists(pfab_cache_file)) {
  mh_params <- readRDS(pfab_cache_file)
  cat(sprintf("  PFAB: loaded cached params from %s\n", pfab_cache_file))
} else {
  cat("  PFAB: running precomputation (swNoData)...\n")
  mh_params <- pfab_precompute(grid_rows, grid_cols, k)
  saveRDS(mh_params, pfab_cache_file)
  cat(sprintf("  PFAB: cached params to %s\n", pfab_cache_file))
}
elapsed_pre <- (proc.time() - t0_pre)[["elapsed"]]
cat(sprintf("  PFAB precompute: %.0fs\n", elapsed_pre))

t0 <- proc.time()
pfab_result <- tryCatch(
  fit_bayesImageS_pfab(y_vec, grid_rows, grid_cols, k, mh_params),
  error = function(e) {
    cat(sprintf("\n  bis_pfab ERROR: %s", e$message))
    NULL
  }
)
elapsed <- (proc.time() - t0)[["elapsed"]]

if (!is.null(pfab_result)) {
  rep_metrics[["bis_pfab"]] <- compute_metrics_bayesImageS(
    pfab_result, z_true, elapsed, k, nug
  )
  rep_metrics[["bis_pfab"]]$elapsed_pre <- elapsed_pre
  cat(sprintf("  bis_pfab: ARI=%.3f misclass=%.3f brier=%.3f t=%.0fs (pre=%.0fs)",
              rep_metrics[["bis_pfab"]]$ari,
              rep_metrics[["bis_pfab"]]$misclass,
              rep_metrics[["bis_pfab"]]$brier, elapsed, elapsed_pre))
} else {
  rep_metrics[["bis_pfab"]] <- list(
    ari = NA_real_, misclass = NA_real_,
    psi_pm = NA_real_, rhat_psi = NA_real_,
    eq_pm = NA_real_, eq_true = NA_real_,
    eq_pmse = NA_real_, brier = NA_real_,
    accept_psi = NA_real_, accept_graph = NA_real_,
    elapsed = elapsed, elapsed_pre = elapsed_pre
  )
}
cat("\n")

# --- Save per-rep result ---
outdir <- file.path("output", sc_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("rep_%03d.rds", rep))
saveRDS(rep_metrics, outfile)
cat(sprintf("Saved %s\n", outfile))
