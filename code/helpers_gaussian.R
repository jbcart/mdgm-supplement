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
    accept_psi = ar[["psi"]],
    accept_graph = ar[["graph"]],
    elapsed = elapsed
  )
}

#' Compute metrics from a bayesImageS smcPotts result
#'
#' @param smc_result Result from smcPotts().
#' @param z_true True latent field (0-indexed integer vector).
#' @param elapsed Elapsed time in seconds.
#' @param k Number of classes.
#' @return Named list of metrics.
compute_metrics_bayesImageS <- function(smc_result, z_true, elapsed, k) {
  # smcPotts returns alloc: matrix of posterior allocation probabilities
  # (n_sites x k). Get MAP estimate.
  z_est <- max.col(smc_result$alloc) - 1L  # 0-indexed

  ari <- mclust::adjustedRandIndex(z_est, z_true)

  z_aligned <- label_align(z_est, z_true, k)
  misclass <- mean(z_aligned != z_true)

  list(
    ari = ari,
    misclass = misclass,
    psi_pm = NA_real_,
    rhat_psi = NA_real_,
    eq_pm = NA_real_,
    eq_true = NA_real_,
    accept_psi = NA_real_,
    accept_graph = NA_real_,
    elapsed = elapsed
  )
}

#' Fit bayesImageS smcPotts model
#'
#' Wrapper around bayesImageS::smcPotts with standard settings.
#'
#' @param y_vec Numeric vector of observations (one per site).
#' @param grid_rows Number of grid rows.
#' @param grid_cols Number of grid columns.
#' @param k Number of classes.
#' @return Result from smcPotts().
fit_bayesImageS_smc <- function(y_vec, grid_rows, grid_cols, k) {
  mask <- matrix(1L, nrow = grid_rows, ncol = grid_cols)
  neighbors <- bayesImageS::getNeighbors(mask, c(2, 2, 0, 0))
  blocks <- bayesImageS::getBlocks(mask, 2)

  priors <- list()
  priors$k <- k
  priors$mu <- rep(0, k)
  priors$mu.sd <- rep(100, k)
  priors$sigma <- rep(1, k)
  priors$sigma.nu <- rep(3, k)
  priors$beta <- c(0, 1.5)  # uniform prior on beta in [0, 1.5]

  param <- list(npart = 10000L, nstat = 50L)

  bayesImageS::smcPotts(y_vec, neighbors, blocks, priors, param = param)
}
