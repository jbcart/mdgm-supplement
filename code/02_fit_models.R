## Fit all models for a given data file
##
## Usage: Rscript code/02_fit_models.R <data_file> [n_iter] [burnin]
##
## Outputs: output/results_<basename>.rds

source(file.path("code", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript code/02_fit_models.R <data_file>")

data_file <- args[1]
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 5000L
burnin <- if (length(args) >= 3) as.integer(args[3]) else 1000L

dat <- readRDS(data_file)
nug <- dat$nug
z_true <- dat$z_true
y <- dat$y
psi_true <- dat$psi_true
theta_true <- dat$theta_true

model_configs <- list(
  mrf_exact = list(spatial = mrf(method = "exchange"),
                   psi_tune = 0.15, psi_init = psi_true),
  mrf_pl    = list(spatial = mrf(method = "pseudo_likelihood"),
                   psi_tune = 0.15, psi_init = psi_true),
  mdgm_st   = list(spatial = mdgm(dag_type = "spanning_tree"),
                   psi_tune = 0.4, psi_init = psi_true * 2),
  mdgm_ao   = list(spatial = mdgm(dag_type = "acyclic_orientation"),
                   psi_tune = 0.3, psi_init = psi_true * 1.25)
)

all_metrics <- list()

for (name in names(model_configs)) {
  cfg <- model_configs[[name]]
  cat(sprintf("Fitting %s... ", name))

  model <- srf_model(nug, spatial = cfg$spatial, emission = "bernoulli")

  t0 <- proc.time()
  fit <- mcmc(model, y = y, z_init = z_true,
              psi_init = cfg$psi_init, theta_init = theta_true,
              n_iter = n_iter, psi_tune = cfg$psi_tune,
              nug = nug)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  all_metrics[[name]] <- compute_metrics(
    fit, z_true, psi_true, theta_true, nug, burnin, elapsed
  )
  cat(sprintf("done (%.1fs)\n", elapsed))
}

basename_noext <- tools::file_path_sans_ext(basename(data_file))
outfile <- sprintf("output/results_%s.rds", basename_noext)
saveRDS(all_metrics, outfile)
cat("Saved to", outfile, "\n")
