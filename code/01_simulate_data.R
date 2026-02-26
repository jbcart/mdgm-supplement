## Generate MRF data for simulation studies
##
## Usage: Rscript code/01_simulate_data.R [grid_size] [psi] [seed]
##
## Outputs: output/data_<grid_size>_psi<psi>_seed<seed>.rds

source(file.path("code", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
grid_size <- if (length(args) >= 1) as.integer(args[1]) else 16L
psi_true <- if (length(args) >= 2) as.numeric(args[2]) else 0.5
seed <- if (length(args) >= 3) as.integer(args[3]) else 42L
theta_true <- c(0.2, 0.8)
n_reps <- 2L

cat(sprintf("Generating data: %dx%d grid, psi=%.2f, seed=%d\n",
            grid_size, grid_size, psi_true, seed))

nug <- nug_from_grid(grid_size, grid_size, seed = seed)
z_true <- sample_mrf(nug, psi = psi_true, seed = seed)
y <- generate_observations(z_true, theta_true, n_reps)

out <- list(
  nug = nug, z_true = z_true, y = y,
  psi_true = psi_true, theta_true = theta_true,
  grid_size = grid_size, seed = seed, n_reps = n_reps
)

outfile <- sprintf("output/data_%d_psi%.2f_seed%d.rds",
                   grid_size, psi_true, seed)
saveRDS(out, outfile)
cat("Saved to", outfile, "\n")
