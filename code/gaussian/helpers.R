## Utility functions for Gaussian emission simulation study

library(mdgm)

#' Generate Gaussian observations from a latent field
#'
#' Produces a single observation per site: y[i] ~ N(mu[z[i]+1], sigma2[z[i]+1]).
#'
#' @param z Integer vector of latent labels (0-indexed).
#' @param mu Numeric vector of class means (length k).
#' @param sigma2 Numeric vector of class variances (length k).
#' @return List of length-1 numeric vectors (one per site), matching mdgm's y
#'   format.
generate_gaussian_obs <- function(z, mu, sigma2) {
  n <- length(z)
  lapply(seq_len(n), function(i) {
    k <- z[i] + 1L
    rnorm(1, mean = mu[k], sd = sqrt(sigma2[k]))
  })
}

#' Greedy label alignment for misclassification rate
#'
#' Finds label permutation that maximizes overlap between estimated and true
#' labels. Uses greedy matching on the confusion matrix.
#'
#' @param z_est Integer vector of estimated labels (0-indexed).
#' @param z_true Integer vector of true labels (0-indexed).
#' @param k Number of classes.
#' @return Integer vector of relabeled z_est (0-indexed).
label_align <- function(z_est, z_true, k) {
  # Build confusion matrix: overlap[i+1, j+1] = count of sites where
  # z_est==i and z_true==j
  overlap <- matrix(0L, nrow = k, ncol = k)
  for (s in seq_along(z_est)) {
    overlap[z_est[s] + 1L, z_true[s] + 1L] <- overlap[z_est[s] + 1L, z_true[s] + 1L] + 1L
  }

  # Greedy matching: iteratively pick the (est, true) pair with highest overlap
  perm <- rep(NA_integer_, k)  # perm[est_label + 1] = true_label
  used_true <- logical(k)
  used_est  <- logical(k)

  for (step in seq_len(k)) {
    best_val <- -1L
    best_i <- 0L
    best_j <- 0L
    for (i in seq_len(k)) {
      if (used_est[i]) next
      for (j in seq_len(k)) {
        if (used_true[j]) next
        if (overlap[i, j] > best_val) {
          best_val <- overlap[i, j]
          best_i <- i
          best_j <- j
        }
      }
    }
    perm[best_i] <- best_j - 1L  # back to 0-indexed
    used_est[best_i] <- TRUE
    used_true[best_j] <- TRUE
  }

  # Apply permutation
  perm[z_est + 1L]
}

#' Compute metrics from an mdgm Gaussian MCMC result
#'
#' @param result An MdgmResult object from mcmc().
#' @param z_true True latent field (0-indexed integer vector).
#' @param nug NaturalUndirectedGraph object.
#' @param burnin Number of burn-in iterations to discard.
#' @param elapsed Elapsed time in seconds.
#' @param k Number of classes.
#' @return Named list of metrics.
compute_metrics_gaussian <- function(result, z_true, nug, burnin, elapsed, k) {
  n_iter <- length(result$psi())
  post_idx <- (burnin + 1L):n_iter
  n <- length(z_true)

  # Posterior mode per site
  z_mat <- result$z()[, post_idx, drop = FALSE]
  z_mode <- apply(z_mat, 1, function(row) {
    tab <- tabulate(row + 1L, nbins = k)
    which.max(tab) - 1L
  })

  # ARI (no label alignment needed for ARI)
  ari <- mclust::adjustedRandIndex(z_mode, z_true)

  # Misclassification rate (with label alignment)
  z_aligned <- label_align(z_mode, z_true, k)
  misclass <- mean(z_aligned != z_true)

  # Sufficient statistic
  eq_true <- sufficient_stat(z_true, nug)
  eq_vec <- apply(z_mat, 2, function(z_j) sufficient_stat(z_j, nug))
  eq_pm <- mean(eq_vec)
  eq_pmse <- mean((eq_vec - eq_true)^2)

  # Brier score
  n_post <- ncol(z_mat)
  alloc_prop <- t(apply(z_mat, 1, function(row) tabulate(row + 1L, nbins = k))) / n_post
  o_mat <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(n)) o_mat[i, z_true[i] + 1L] <- 1
  brier <- mean(rowSums((alloc_prop - o_mat)^2))

  # Psi diagnostics
  psi_chain <- result$psi()[post_idx]
  psi_pm <- mean(psi_chain)
  rhat_psi <- split_rhat(psi_chain)

  # Acceptance rates
  ar <- result$acceptance_rates()

  list(
    ari = ari,
    misclass = misclass,
    psi_pm = psi_pm,
    rhat_psi = rhat_psi,
    eq_pm = eq_pm,
    eq_true = eq_true,
    eq_pmse = eq_pmse,
    brier = brier,
    accept_psi = ar[["psi"]],
    accept_graph = ar[["graph"]],
    elapsed = elapsed
  )
}

#' Compute metrics from a bayesImageS mcmcPotts (PFAB) result
#'
#' @param pfab_result Result from mcmcPotts() with algorithm="aux".
#' @param z_true True latent field (0-indexed integer vector).
#' @param elapsed Elapsed time in seconds.
#' @param k Number of classes.
#' @return Named list of metrics.
compute_metrics_bayesImageS <- function(pfab_result, z_true, elapsed, k,
                                         nug, nburn = 2000L) {
  # mcmcPotts returns $alloc: n_sites x k matrix of post-burn-in allocation
  # counts. MAP estimate = column with highest count per row.
  z_est <- max.col(pfab_result$alloc) - 1L  # 0-indexed
  n <- length(z_true)

  ari <- mclust::adjustedRandIndex(z_est, z_true)

  z_aligned <- label_align(z_est, z_true, k)
  misclass <- mean(z_aligned != z_true)

  # Extract beta posterior mean from MCMC chain
  psi_pm <- mean(pfab_result$beta)

  # Sufficient statistic: eq_pmse from full sum chain (trim to post-burnin)
  eq_true <- sufficient_stat(z_true, nug)
  sum_chain <- as.numeric(pfab_result$sum)
  sum_post <- sum_chain[(nburn + 1L):length(sum_chain)]
  eq_pm <- mean(sum_post)
  eq_pmse <- mean((sum_post - eq_true)^2)

  # Brier score from allocation proportions
  alloc_prop <- pfab_result$alloc / rowSums(pfab_result$alloc)
  o_mat <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(n)) o_mat[i, z_true[i] + 1L] <- 1
  brier <- mean(rowSums((alloc_prop - o_mat)^2))

  list(
    ari = ari,
    misclass = misclass,
    psi_pm = psi_pm,
    rhat_psi = NA_real_,
    eq_pm = eq_pm,
    eq_true = eq_true,
    eq_pmse = eq_pmse,
    brier = brier,
    accept_psi = NA_real_,
    accept_graph = NA_real_,
    elapsed = elapsed
  )
}

#' PFAB precomputation: run swNoData and fit parametric surrogate
#'
#' Runs Swendsen-Wang simulations without data at a grid of beta values, then
#' fits a piecewise parametric model to the sufficient statistic curve to obtain
#' the surrogate parameters needed by mcmcPotts(algorithm="aux").
#'
#' @param grid_rows Number of grid rows.
#' @param grid_cols Number of grid columns.
#' @param k Number of classes.
#' @param n_sw Number of Swendsen-Wang iterations per beta (default 800).
#' @param burn Burn-in iterations to discard (default 201).
#' @return Named list of mh parameters for mcmcPotts().
pfab_precompute <- function(grid_rows, grid_cols, k, n_sw = 800L, burn = 201L) {
  mask <- matrix(1L, nrow = grid_rows, ncol = grid_cols)
  neigh <- bayesImageS::getNeighbors(mask, c(2, 2, 0, 0))
  block <- bayesImageS::getBlocks(mask, 2)

  # Analytical constants
  n_edges <- 2L * grid_rows * grid_cols - grid_rows - grid_cols
  bcrit <- log(1 + sqrt(k))
  E0 <- n_edges / k
  V0 <- n_edges * (1 / k) * (1 - 1 / k)
  maxS <- n_edges

  # Beta grid with density near bcrit (following PFAB vignette)
  beta_grid <- sort(unique(c(
    seq(0, 1, by = 0.1),
    seq(1.05, 1.15, by = 0.05),
    bcrit - 0.05, bcrit - 0.02, bcrit + 0.02,
    seq(1.3, 1.4, by = 0.05),
    seq(1.5, 2, by = 0.1),
    2.5, 3
  )))

  cat(sprintf("  PFAB precompute: k=%d, grid=%dx%d, %d beta values\n",
              k, grid_rows, grid_cols, length(beta_grid)))

  # Run swNoData at each beta
  emp_mean <- numeric(length(beta_grid))
  emp_var <- numeric(length(beta_grid))

  for (i in seq_along(beta_grid)) {
    res <- bayesImageS::swNoData(beta_grid[i], k, neigh, block, niter = n_sw)
    s <- res$sum[burn:n_sw, 1]
    emp_mean[i] <- mean(s)
    emp_var[i] <- var(s)
  }

  # Fit piecewise mean curve to extract Ecrit, phi1, phi2
  # For beta <= bcrit: E(beta) = E0 + (Ecrit - E0) * (beta/bcrit)^phi1
  # For beta >  bcrit: E(beta) = maxS - (maxS - Ecrit) * ((bmax - beta)/(bmax - bcrit))^phi2
  bmax <- max(beta_grid)

  fit_mean <- function(par) {
    Ecrit <- par[1]
    phi1 <- exp(par[2])  # ensure positive
    phi2 <- exp(par[3])

    pred <- numeric(length(beta_grid))
    for (i in seq_along(beta_grid)) {
      b <- beta_grid[i]
      if (b <= bcrit) {
        pred[i] <- E0 + (Ecrit - E0) * (b / bcrit)^phi1
      } else {
        pred[i] <- maxS - (maxS - Ecrit) * ((bmax - b) / (bmax - bcrit))^phi2
      }
    }
    sum((pred - emp_mean)^2)
  }

  # Initial values
  Ecrit_init <- emp_mean[which.min(abs(beta_grid - bcrit))]
  opt <- optim(c(Ecrit_init, log(2), log(2)), fit_mean, method = "Nelder-Mead",
               control = list(maxit = 10000))
  Ecrit <- opt$par[1]
  phi1 <- exp(opt$par[2])
  phi2 <- exp(opt$par[3])

  # Variance peaks below/above critical temperature
  below <- beta_grid <= bcrit
  above <- beta_grid > bcrit
  Vmax1 <- max(emp_var[below])
  Vmax2 <- if (any(above)) max(emp_var[above]) else Vmax1

  cat(sprintf("  PFAB params: Ecrit=%.1f, phi1=%.3f, phi2=%.3f, Vmax1=%.0f, Vmax2=%.0f\n",
              Ecrit, phi1, phi2, Vmax1, Vmax2))

  list(
    algorithm = "aux",
    bandwidth = 0.02,
    Vmax1 = Vmax1,
    Vmax2 = Vmax2,
    E0 = E0,
    Ecrit = Ecrit,
    phi1 = phi1,
    phi2 = phi2,
    factor = 1,
    bcrit = bcrit,
    V0 = V0
  )
}

#' Fit bayesImageS mcmcPotts with PFAB algorithm
#'
#' Wrapper around bayesImageS::mcmcPotts with algorithm="aux" using
#' precomputed surrogate parameters.
#'
#' @param y_vec Numeric vector of observations (one per site).
#' @param grid_rows Number of grid rows.
#' @param grid_cols Number of grid columns.
#' @param k Number of classes.
#' @param mh_params Precomputed PFAB parameters from pfab_precompute().
#' @param niter Total MCMC iterations (default 10000).
#' @param nburn Burn-in iterations (default 2000).
#' @return Result from mcmcPotts().
fit_bayesImageS_pfab <- function(y_vec, grid_rows, grid_cols, k, mh_params,
                                 niter = 10000L, nburn = 2000L) {
  mask <- matrix(1L, nrow = grid_rows, ncol = grid_cols)
  neigh <- bayesImageS::getNeighbors(mask, c(2, 2, 0, 0))
  block <- bayesImageS::getBlocks(mask, 2)

  priors <- list()
  priors$k <- k
  priors$mu <- rep(0, k)
  priors$mu.sd <- rep(100, k)
  priors$sigma <- rep(1, k)
  priors$sigma.nu <- rep(3, k)
  priors$beta <- c(0, log(1 + sqrt(k)) + 0.5)

  bayesImageS::mcmcPotts(y_vec, neigh, block, priors, mh_params,
                          niter = niter, nburn = nburn)
}
