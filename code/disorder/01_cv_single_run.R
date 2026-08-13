## Cross-validation: single iteration for disorder analysis
##
## Usage: Rscript code/disorder/01_cv_single_run.R <iter>
## (Run from the mdgm-supplement root directory)
##
## Holds out 60 block groups, fits MDGM-ST and MRF-PL, computes MAE metrics.
## Saves result to output/disorder/run_cv_<NNN>.rds

library(mdgm)
library(sf)

sf_use_s2(FALSE)

# --- Parse args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript code/disorder/01_cv_single_run.R <iter>")
}
iter <- as.integer(args[1])
cat(sprintf("CV iteration %d\n", iter))

# --- Settings ---
n_iter <- 10000L
burnin <- 2000L

# --- Data ---
load(file.path("code", "disorder", "clean_data.Rbin"))
n <- length(y_bg270)

# Rook adjacency
adj_mat <- 1L * as.matrix(st_relate(bg270, pattern = "F***1****"))
diag(adj_mat) <- 0L
nug <- nug_from_adj_mat(adj_mat, seed = 42L)

# Convert observations
y_full <- lapply(y_bg270, function(x) {
  if (all(is.na(x))) numeric(0) else as.numeric(x)
})

# --- Hold out 60 sites with observations ---
which_obs <- which(vapply(y_full, function(x) length(x) > 0, logical(1)))
set.seed(iter)
test_set <- sample(which_obs, size = 60)

# Ground truth for test sites: mean of held-out observations
y_test <- sapply(y_full[test_set], mean)
y_w    <- sapply(y_full[test_set], length)

# Training data: set held-out sites to missing
y_train <- y_full
y_train[test_set] <- lapply(seq_along(test_set), function(i) numeric(0))

# --- Shared initial values ---
seed_i <- as.integer((iter * 100003L) %% .Machine$integer.max)
set.seed(seed_i)
z_init     <- rbinom(n, 1, 0.5)
theta_init <- c(0.25, 0.75)

# --- Fit MDGM-ST ---
model_st <- srf_model(nug, spatial = mdgm(dag_type = "spanning_tree"),
                       emission = "bernoulli")
fit_st <- mcmc(model_st, y = y_train, z_init = z_init,
               psi_init = 0.5, theta_init = theta_init,
               n_iter = n_iter, psi_tune = 0.5,
               store_z = TRUE, seed = seed_i, nug = nug)

# --- Fit MRF-PL ---
model_pl <- srf_model(nug, spatial = mrf(method = "pseudo_likelihood"),
                       emission = "bernoulli")
fit_pl <- mcmc(model_pl, y = y_train, z_init = z_init,
               psi_init = 0.1, theta_init = theta_init,
               n_iter = n_iter, psi_tune = 0.1,
               store_z = TRUE, seed = seed_i + 1L, nug = nug)

# --- Posterior predictive and z posterior mean ---
post_idx <- (burnin + 1L):n_iter

compute_cv_metrics <- function(result, test_set, y_test, y_w) {
  z_mat <- result$z()[, post_idx, drop = FALSE]
  p_mat <- result$emission_params()$p  # 2 x n_iter

  # Posterior predictive mean: E[p_{z_i+1}] averaged over post-burnin iterations
  pp1 <- sapply(test_set, function(i) {
    mean(vapply(seq_along(post_idx), function(j) {
      p_mat[z_mat[i, j] + 1L, post_idx[j]]
    }, double(1)))
  })

  # Posterior mean of z
  z_pm <- rowMeans(z_mat)
  z_pm_test <- z_pm[test_set]

  # MAE metrics (8 total: unweighted/weighted x pp1/zpm)
  mae_pp1     <- mean(abs(pp1 - y_test))
  mae_zpm     <- mean(abs(z_pm_test - y_test))
  mae_pp1_w   <- sum(abs(pp1 - y_test) * y_w) / sum(y_w)
  mae_zpm_w   <- sum(abs(z_pm_test - y_test) * y_w) / sum(y_w)

  c(mae_pp1 = mae_pp1, mae_zpm = mae_zpm,
    mae_pp1_w = mae_pp1_w, mae_zpm_w = mae_zpm_w)
}

metrics_st <- compute_cv_metrics(fit_st, test_set, y_test, y_w)
metrics_pl <- compute_cv_metrics(fit_pl, test_set, y_test, y_w)

result <- c(
  setNames(metrics_st, paste0(names(metrics_st), "_st")),
  setNames(metrics_pl, paste0(names(metrics_pl), "_pl"))
)

# --- Save ---
outdir <- file.path("output", "disorder")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("run_cv_%03d.rds", iter))
saveRDS(result, outfile)
cat(sprintf("Saved %s\n", outfile))
