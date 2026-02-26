## Utility functions for mdgm simulation studies

library(mdgm)

#' Generate Bernoulli observations from a latent field
#'
#' @param z Integer vector of latent labels (0-indexed).
#' @param theta Numeric vector of success probabilities per class
#'   (theta[1] = P(y=1|z=0), theta[2] = P(y=1|z=1)).
#' @param n_reps Number of replicates per site (used when lambda = 0).
#' @param lambda If > 0, the number of observations per site is drawn from
#'   Poisson(lambda) instead of using n_reps. Sites with 0 draws get
#'   numeric(0) (missing data).
#' @return List of numeric vectors (one per site).
generate_observations <- function(z, theta, n_reps, lambda = 0) {
  n <- length(z)
  if (lambda > 0) {
    site_reps <- rpois(n, lambda)
  } else {
    site_reps <- rep(as.integer(n_reps), n)
  }
  lapply(seq_len(n), function(i) {
    m <- site_reps[i]
    if (m == 0L) return(numeric(0))
    rbinom(m, size = 1, prob = theta[z[i] + 1L])
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

#' Total variation distance between two empirical discrete distributions
#'
#' Constructs empirical PMFs over the union of observed values from two
#' integer-valued MCMC sample vectors, then computes
#' TV(P, Q) = 0.5 * sum_x |P(x) - Q(x)|.
#'
#' @param x Integer vector of samples from distribution P.
#' @param y Integer vector of samples from distribution Q.
#' @return Scalar total variation distance in [0, 1].
tv_distance <- function(x, y) {
  all_vals <- union(x, y)
  px <- tabulate(match(x, all_vals), nbins = length(all_vals)) / length(x)
  py <- tabulate(match(y, all_vals), nbins = length(all_vals)) / length(y)
  0.5 * sum(abs(px - py))
}

#' Compute metrics from a fitted model result
#'
#' Matches the statistics computed in the original sim_study.r from
#' AHDC-carter/spatial_ddd/supplement.
#'
#' @param result An MdgmResult object from mcmc().
#' @param z_true True latent field (0-indexed integer vector).
#' @param psi_true True psi value.
#' @param theta_true True theta vector.
#' @param nug NaturalUndirectedGraph object.
#' @param burnin Number of burn-in iterations to discard.
#' @param elapsed Elapsed time in seconds.
#' @param missing_sites Logical vector indicating sites with no observations
#'   (length n). If NULL, no missing-data accuracy is computed.
#' @return Named list of metrics.
compute_metrics <- function(result, z_true, psi_true, theta_true, nug,
                            burnin, elapsed = NA_real_,
                            missing_sites = NULL) {
  n_iter <- length(result$psi())
  post_idx <- (burnin + 1L):n_iter
  n <- length(z_true)
  n_colors <- length(theta_true)

  # Psi — posterior mean squared error = variance + bias^2
  # (matches original: mean((psi_true - out$psi[burnin:J])^2))
  psi_chain <- result$psi()[post_idx]
  psi_pm <- mean(psi_chain)
  psi_psd <- sd(psi_chain)
  psi_pmse <- mean((psi_true - psi_chain)^2)

  # Theta (Bernoulli: emission_params()$p is n_colors x J matrix)
  theta_mat <- result$emission_params()$p
  theta_pm <- numeric(n_colors)
  theta_psd <- numeric(n_colors)
  theta_pmse <- numeric(n_colors)
  for (k in seq_len(n_colors)) {
    chain_k <- theta_mat[k, post_idx]
    theta_pm[k] <- mean(chain_k)
    theta_psd[k] <- sd(chain_k)
    theta_pmse[k] <- mean((theta_true[k] - chain_k)^2)
  }

  # Classification accuracy per iteration, then summarize
  z_mat <- result$z()[, post_idx, drop = FALSE]
  acc_vec <- apply(z_mat, 2, function(z_j) mean(z_j == z_true))
  acc_pm <- mean(acc_vec)
  acc_psd <- sd(acc_vec)

  # Missing-data site accuracy (original: mcmc_acc_NA)
  if (!is.null(missing_sites) && any(missing_sites)) {
    acc_na_vec <- apply(z_mat[missing_sites, , drop = FALSE], 2,
                        function(z_j) mean(z_j == z_true[missing_sites]))
    acc_na_pm <- mean(acc_na_vec)
    acc_na_psd <- sd(acc_na_vec)
  } else {
    acc_na_pm <- NA_real_
    acc_na_psd <- NA_real_
  }

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
  eq_pmse <- mean((eq_vec - eq_true)^2)

  # Acceptance rates
  ar <- result$acceptance_rates()
  accept_psi <- ar[["psi"]]
  accept_graph <- ar[["graph"]]

  # R-hat (including sufficient stat chain, matching original)
  rhat_psi <- split_rhat(psi_chain)
  rhat_theta <- vapply(seq_len(n_colors), function(k) {
    split_rhat(as.numeric(theta_mat[k, post_idx]))
  }, double(1))
  rhat_eq <- split_rhat(as.numeric(eq_vec))

  list(
    psi_pm = psi_pm, psi_psd = psi_psd, psi_pmse = psi_pmse,
    theta_pm = theta_pm, theta_psd = theta_psd, theta_pmse = theta_pmse,
    acc_pred = acc_pred, acc_pm = acc_pm, acc_psd = acc_psd,
    acc_na_pm = acc_na_pm, acc_na_psd = acc_na_psd,
    eq_vec = eq_vec, eq_pm = eq_pm, eq_psd = eq_psd, eq_pmse = eq_pmse,
    accept_psi = accept_psi, accept_graph = accept_graph,
    rhat_psi = rhat_psi, rhat_theta = rhat_theta, rhat_eq = rhat_eq,
    elapsed = elapsed
  )
}
