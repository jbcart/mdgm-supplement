## Utility functions for mdgm simulation studies

library(mdgm)

#' Generate Bernoulli observations from a latent field
#'
#' @param z Integer vector of latent labels (0-indexed).
#' @param theta Numeric vector of success probabilities per class.
#' @param n_reps Number of replicates per site.
#' @return List of numeric vectors (one per site).
generate_observations <- function(z, theta, n_reps) {
  lapply(z, function(zi) {
    rbinom(n_reps, size = 1, prob = theta[zi + 1L])
  })
}

#' Count same-color edges (sufficient statistic for Ising/Potts)
#'
#' @param z Integer vector of labels (0-indexed).
#' @param nug A NaturalUndirectedGraph object.
#' @return Integer count of edges where both endpoints share the same color.
sufficient_stat <- function(z, nug) {
  n <- nug$nvertices()
  count <- 0L
  for (i in seq_len(n)) {
    nbrs <- nug$neighbors(i)  # 1-indexed
    for (nb in nbrs) {
      if (nb > i) {
        count <- count + (z[i] == z[nb])
      }
    }
  }
  count
}

#' Split-chain R-hat diagnostic
#'
#' @param chain Numeric vector of MCMC samples.
#' @return Split R-hat statistic.
split_rhat <- function(chain) {
  n <- length(chain)
  if (n < 4L) return(NA_real_)
  mid <- n %/% 2L
  chains <- list(chain[1:mid], chain[(mid + 1L):n])
  m <- length(chains)
  ns <- vapply(chains, length, integer(1))
  means <- vapply(chains, mean, double(1))
  vars <- vapply(chains, var, double(1))
  n_min <- min(ns)
  grand_mean <- mean(means)
  B <- n_min * sum((means - grand_mean)^2) / (m - 1L)
  W <- mean(vars)
  if (W < .Machine$double.eps) return(NA_real_)
  var_hat <- (1 - 1 / n_min) * W + B / n_min
  sqrt(var_hat / W)
}

#' Compute metrics from a fitted model result
#'
#' @param result An MdgmResult object from mcmc().
#' @param z_true True latent field (0-indexed integer vector).
#' @param psi_true True psi value.
#' @param theta_true True theta vector.
#' @param nug NaturalUndirectedGraph object.
#' @param burnin Number of burn-in iterations to discard.
#' @param elapsed Elapsed time in seconds.
#' @return Named list of metrics.
compute_metrics <- function(result, z_true, psi_true, theta_true, nug,
                            burnin, elapsed = NA_real_) {
  n_iter <- length(result$psi())
  post_idx <- (burnin + 1L):n_iter
  n <- length(z_true)
  n_colors <- length(theta_true)

  # Psi
  psi_chain <- result$psi()[post_idx]
  psi_pm <- mean(psi_chain)
  psi_psd <- sd(psi_chain)
  psi_pmse <- (psi_pm - psi_true)^2

  # Theta (Bernoulli: emission_params()$p is n_colors x J matrix)
  theta_mat <- result$emission_params()$p
  theta_pm <- numeric(n_colors)
  theta_psd <- numeric(n_colors)
  theta_pmse <- numeric(n_colors)
  for (k in seq_len(n_colors)) {
    chain_k <- theta_mat[k, post_idx]
    theta_pm[k] <- mean(chain_k)
    theta_psd[k] <- sd(chain_k)
    theta_pmse[k] <- (theta_pm[k] - theta_true[k])^2
  }

  # Classification accuracy per iteration, then summarize
  z_mat <- result$z()[, post_idx, drop = FALSE]
  acc_vec <- apply(z_mat, 2, function(z_j) mean(z_j == z_true))
  acc_pm <- mean(acc_vec)
  acc_psd <- sd(acc_vec)

  # Prediction accuracy (posterior mode per site)
  z_mode <- apply(z_mat, 1, function(row) {
    tab <- tabulate(row + 1L, nbins = n_colors)
    which.max(tab) - 1L
  })
  acc_pred <- mean(z_mode == z_true)

  # Sufficient statistic (eq = edge quality)
  eq_true <- sufficient_stat(z_true, nug)
  eq_vec <- apply(z_mat, 2, function(z_j) sufficient_stat(z_j, nug))
  eq_pm <- mean(eq_vec)
  eq_psd <- sd(eq_vec)
  eq_pmse <- (eq_pm - eq_true)^2

  # Acceptance rates
  ar <- result$acceptance_rates()
  accept_psi <- ar[["psi"]]
  accept_graph <- ar[["graph"]]

  # R-hat
  rhat_psi <- split_rhat(psi_chain)
  rhat_theta <- vapply(seq_len(n_colors), function(k) {
    split_rhat(as.numeric(theta_mat[k, post_idx]))
  }, double(1))

  list(
    psi_pm = psi_pm, psi_psd = psi_psd, psi_pmse = psi_pmse,
    theta_pm = theta_pm, theta_psd = theta_psd, theta_pmse = theta_pmse,
    acc_pred = acc_pred, acc_pm = acc_pm, acc_psd = acc_psd,
    eq_pm = eq_pm, eq_psd = eq_psd, eq_pmse = eq_pmse,
    accept_psi = accept_psi, accept_graph = accept_graph,
    rhat_psi = rhat_psi, rhat_theta = rhat_theta,
    elapsed = elapsed
  )
}
