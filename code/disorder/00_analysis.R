## Disorder analysis: MDGM-ST and aMRF models on Columbus block group data
##
## Usage: Rscript code/disorder/00_analysis.R
## (Run from the mdgm-supplement root directory)

library(mdgm)
library(sf)
library(ggplot2)
library(igraph)
library(coda)
library(tidyverse)

sf_use_s2(FALSE)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
load(file.path("code", "disorder", "clean_data.Rbin"))

n <- length(y_bg270)

# Rook adjacency matrix
adj_mat <- 1L * as.matrix(st_relate(bg270, pattern = "F***1****"))
diag(adj_mat) <- 0L

nug <- nug_from_adj_mat(adj_mat, seed = 42L)

bg_layout <- st_coordinates(st_centroid(bg270))

# Convert observations: all-NA entries -> numeric(0); keep non-NA as-is
y <- lapply(y_bg270, function(x) {
  if (all(is.na(x))) numeric(0) else as.numeric(x)
})

which_missing <- vapply(y, function(x) length(x) == 0, logical(1))
cat(sprintf("Sites: %d, Missing: %d, Edges: %d\n",
            n, sum(which_missing), nug$nedges()))

# ---------------------------------------------------------------------------
# MCMC settings
# ---------------------------------------------------------------------------
n_iter  <- 10000L
burnin  <- 2000L
n_chains <- 4L

# ---------------------------------------------------------------------------
# Fit MDGM-ST (4 chains)
# ---------------------------------------------------------------------------
cat("Fitting MDGM-ST ...\n")

model_st <- srf_model(nug, spatial = mdgm(dag_type = "spanning_tree"),
                       emission = "bernoulli")

fit_st <- parallel::mclapply(seq_len(n_chains), function(i) {
  seed_i <- as.integer((42L * 100003L + i * 1000033L) %% .Machine$integer.max)
  set.seed(seed_i)
  z_init     <- rbinom(n, 1, 0.5)
  theta_init <- c(runif(1, 0, 0.5), runif(1, 0.5, 1))
  psi_init   <- runif(1, 0, 5)

  mdgm::mcmc(model_st, y = y, z_init = z_init, psi_init = psi_init,
             theta_init = theta_init, n_iter = n_iter, psi_tune = 0.5,
             store_z = TRUE, seed = seed_i, nug = nug)
}, mc.cores = n_chains)

cat("  Acceptance rates (psi):",
    paste(round(sapply(fit_st, function(r) r$acceptance_rates()[["psi"]]), 3),
          collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Fit MRF-PL (4 chains)
# ---------------------------------------------------------------------------
cat("Fitting MRF-PL ...\n")

model_pl <- srf_model(nug, spatial = mrf(method = "pseudo_likelihood"),
                       emission = "bernoulli")

fit_pl <- parallel::mclapply(seq_len(n_chains), function(i) {
  seed_i <- as.integer((43L * 100003L + i * 1000033L) %% .Machine$integer.max)
  set.seed(seed_i)
  z_init     <- rbinom(n, 1, 0.5)
  theta_init <- c(runif(1, 0, 0.5), runif(1, 0.5, 1))
  psi_init   <- runif(1, 0, 0.9)

  mdgm::mcmc(model_pl, y = y, z_init = z_init, psi_init = psi_init,
             theta_init = theta_init, n_iter = n_iter, psi_tune = 0.1,
             store_z = TRUE, seed = seed_i, nug = nug)
}, mc.cores = n_chains)

cat("  Acceptance rates (psi):",
    paste(round(sapply(fit_pl, function(r) r$acceptance_rates()[["psi"]]), 3),
          collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Diagnostics (coda)
# ---------------------------------------------------------------------------
source(file.path("code", "common", "helpers.R"))

make_coda_chain <- function(result, burnin, nug) {
  post_idx <- (burnin + 1L):length(result$psi())
  p_mat    <- result$emission_params()$p   # n_colors x n_iter
  z_mat    <- result$z()[, post_idx, drop = FALSE]
  M        <- nug$nedges()
  tz       <- apply(z_mat, 2, function(z_j) sufficient_stat(z_j, nug))
  coda::mcmc(data.frame(
    p_0      = p_mat[1, post_idx],
    p_1      = p_mat[2, post_idx],
    psi      = result$psi()[post_idx],
    tz_over_M = tz / M
  ))
}

coda_st <- coda::mcmc.list(lapply(fit_st, make_coda_chain, burnin = burnin, nug = nug))
coda_pl <- coda::mcmc.list(lapply(fit_pl, make_coda_chain, burnin = burnin, nug = nug))

cat("\nGelman-Rubin diagnostics (MDGM-ST):\n")
print(coda::gelman.diag(coda_st, autoburnin = FALSE))

cat("\nGelman-Rubin diagnostics (MRF-PL):\n")
print(coda::gelman.diag(coda_pl, autoburnin = FALSE))

# ---------------------------------------------------------------------------
# Posterior inference
# ---------------------------------------------------------------------------
post_idx <- (burnin + 1L):n_iter

z_pm_mdgm_st <- Reduce(`+`, lapply(fit_st, function(r) {
  rowMeans(r$z()[, post_idx, drop = FALSE])
})) / n_chains

z_pm_mrf_pl <- Reduce(`+`, lapply(fit_pl, function(r) {
  rowMeans(r$z()[, post_idx, drop = FALSE])
})) / n_chains

# ---------------------------------------------------------------------------
# Posterior comparison data frame
# ---------------------------------------------------------------------------
make_draws_df <- function(results, model_label, burnin, nug) {
  M <- nug$nedges()
  dfs <- lapply(results, function(r) {
    post_idx <- (burnin + 1L):length(r$psi())
    p_mat <- r$emission_params()$p
    z_mat <- r$z()[, post_idx, drop = FALSE]
    tz    <- apply(z_mat, 2, function(z_j) sufficient_stat(z_j, nug))
    data.frame(
      eta_0     = p_mat[1, post_idx],
      eta_1     = p_mat[2, post_idx],
      tz_over_M = tz / M
    )
  })
  df <- bind_rows(dfs)
  df %>%
    pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
    mutate(model = model_label)
}

df_st <- make_draws_df(fit_st, "MDGM-ST", burnin, nug)
df_pl <- make_draws_df(fit_pl, "aMRF", burnin, nug)
df_draws <- bind_rows(df_st, df_pl)

df_draws$stat <- factor(df_draws$stat,
                        levels = c("eta_0", "eta_1", "tz_over_M"),
                        labels = c("eta[0]", "eta[1]", "T(z)/M"))
df_draws$model <- factor(df_draws$model, levels = c("MDGM-ST", "aMRF"))

# ---------------------------------------------------------------------------
# Edge inclusion probabilities (MDGM-ST)
# ---------------------------------------------------------------------------
eip_list <- lapply(fit_st, function(r) {
  r$edge_inclusion_probs(nug = nug, burnin = burnin)
})

eip <- eip_list[[1]]
eip$prob <- Reduce(`+`, lapply(eip_list, `[[`, "prob")) / n_chains

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
figdir <- file.path("output", "disorder", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

h <- 600 / 1.25
w <- 600 / 1.25
q <- 90

thm <- theme_bw() +
  theme(text = element_text(size = 20),
        axis.text  = element_blank(),
        axis.ticks = element_blank())
scale_fill_bw  <- scale_fill_continuous(low = "white", high = "grey20")
scale_fill_div <- scale_fill_gradient2(low = "#882255", mid = "white",
                                       high = "#117733", midpoint = 0)

# 1. Posterior MDGM-ST
jpeg(file.path(figdir, "posterior_mdgm_st.jpeg"), height = h, width = w, quality = q)
print(ggplot(bg270) + geom_sf(aes(fill = z_pm_mdgm_st)) + thm +
        scale_fill_bw + labs(fill = "E(z|y)", title = "MDGM-ST"))
dev.off()

# 2. Posterior MRF-PL
jpeg(file.path(figdir, "posterior_mrf_pl.jpeg"), height = h, width = w, quality = q)
print(ggplot(bg270) + geom_sf(aes(fill = z_pm_mrf_pl)) + thm +
        scale_fill_bw + labs(fill = "E(z|y)", title = "aMRF"))
dev.off()

# 3. Difference map
jpeg(file.path(figdir, "posterior_difference.jpeg"), height = h, width = w, quality = q)
print(ggplot(bg270) + geom_sf(aes(fill = z_pm_mdgm_st - z_pm_mrf_pl)) + thm +
        scale_fill_div + labs(fill = "", title = "MDGM-ST - aMRF") +
        geom_sf(data = bg270[which_missing, ], linewidth = 0.6, alpha = 0,
                color = "black"))
dev.off()

# 4. Posterior comparison boxplots
jpeg(file.path(figdir, "posterior_comparison.jpeg"),
     height = h, width = 1000 / 1.25, quality = q)
print(ggplot(df_draws) +
        geom_boxplot(aes(x = model, y = value, fill = model), outliers = FALSE) +
        facet_wrap(~stat, scales = "free", labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x    = element_blank(),
              axis.title.x   = element_blank(),
              axis.title.y   = element_blank(),
              text = element_text(size = 20)) +
        scale_fill_manual(values = c("#117733", "#AA4499")) +
        ggtitle(" "))
dev.off()

# 5. Edge inclusion probability plot
jpeg(file.path(figdir, "edge_inclusion.jpeg"), height = 800, width = 950, quality = q)
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
par(mar = c(1, 1, 4, 0))
el <- as.matrix(eip[, c("vertex1", "vertex2")])
g  <- igraph::graph_from_edgelist(el, directed = FALSE)
igraph::E(g)$width <- eip$prob * 3
igraph::E(g)$color <- gray(1 - eip$prob)
plot(g, layout = bg_layout, vertex.size = 1, vertex.label = NA,
     vertex.color = "grey30", vertex.frame.color = "grey30",
     rescale = TRUE, asp = 1, main = "")
title(main = "Edge Inclusion Probability", cex.main = 2)
# Line thickness legend
par(mar = c(20, 0, 20, 2))
ref_probs <- c(0.2, 0.4, 0.6, 0.8, 1.0)
n_ref <- length(ref_probs)
plot(NULL, xlim = c(0, 1), ylim = c(0.5, n_ref + 0.7),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
for (i in seq_along(ref_probs)) {
  p <- ref_probs[i]
  segments(0.05, i, 0.55, i, lwd = p * 3, col = gray(1 - p))
  text(0.75, i, labels = sprintf("%.1f", p), cex = 1.2, adj = 0)
}
text(0.4, n_ref + 0.5, "Prob.", cex = 1.5, font = 2)
dev.off()

# 6. Edge inclusion cutoff panels (2x2: quantile graphs + histogram)
qs <- quantile(eip$prob, probs = c(0.25, 0.50, 0.75))
all_verts <- sort(unique(c(eip$vertex1, eip$vertex2)))
n_verts <- max(all_verts)

jpeg(file.path(figdir, "edge_inclusion_cutoffs.jpeg"),
     height = 1200, width = 1200, quality = q)
par(mfrow = c(2, 2))
for (i in seq_along(qs)) {
  par(mar = c(2, 2, 4, 1))
  cutoff <- qs[i]
  keep <- eip$prob >= cutoff
  g_cut <- igraph::make_empty_graph(n = n_verts, directed = FALSE)
  if (sum(keep) > 0) {
    el_cut <- as.matrix(eip[keep, c("vertex1", "vertex2")])
    g_cut <- igraph::add_edges(g_cut, as.vector(t(el_cut)))
    igraph::E(g_cut)$width <- 1
    igraph::E(g_cut)$color <- "grey30"
  }
  pct_label <- c("25th", "50th", "75th")[i]
  plot(g_cut, layout = bg_layout, vertex.size = 1, vertex.label = NA,
       vertex.color = "grey30", vertex.frame.color = "grey30",
       rescale = TRUE, asp = 0, main = "")
  title(main = sprintf("%s percentile (>= %.2f)", pct_label, cutoff), cex.main = 2)
}
par(mar = c(5, 5, 4, 1))
hist(eip$prob, breaks = 30, col = "grey70", border = "white",
     main = "", xlab = "Edge Inclusion Probability",
     ylab = "Frequency", cex.lab = 1.4, cex.axis = 1.2)
title(main = "EIP Distribution", cex.main = 2)
abline(v = qs, lty = 2, lwd = 2, col = c("#332288", "#44AA99", "#882255"))
legend("topright",
       legend = sprintf("%s: %.2f", c("Q25", "Q50", "Q75"), qs),
       col = c("#332288", "#44AA99", "#882255"),
       lty = 2, lwd = 2, cex = 1.3)
dev.off()

# 7. EIP over block group outlines (continuous scale)
seg_df <- data.frame(
  x    = bg_layout[eip$vertex1, 1],
  y    = bg_layout[eip$vertex1, 2],
  xend = bg_layout[eip$vertex2, 1],
  yend = bg_layout[eip$vertex2, 2],
  prob = eip$prob
)

jpeg(file.path(figdir, "eip_over_bg.jpeg"), height = h, width = w, quality = q)
print(
  ggplot() +
    geom_sf(data = bg270, fill = NA, color = "grey80", linewidth = 0.3) +
    geom_segment(data = seg_df[seg_df$prob >= 0.1, ],
                 aes(x = x, y = y, xend = xend, yend = yend,
                     color = prob),
                 linewidth = 0.6, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno", direction = 1,
                          name = "EIP", limits = c(0, 1)) +
    theme_void() +
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "Edge Inclusion Probability")
)
dev.off()

# 8. EIP cutoffs over block group outlines (2x2)
library(patchwork)

cutoff_panels <- lapply(seq_along(qs), function(i) {
  cutoff <- qs[i]
  pct_label <- c("25th", "50th", "75th")[i]
  seg_cut <- seg_df[seg_df$prob >= cutoff, ]
  p <- ggplot() +
    geom_sf(data = bg270, fill = NA, color = "grey80", linewidth = 0.3)
  if (nrow(seg_cut) > 0) {
    p <- p + geom_segment(data = seg_cut,
                          aes(x = x, y = y, xend = xend, yend = yend),
                          color = "#882255", linewidth = 0.6, alpha = 0.7)
  }
  p + theme_void() +
    theme(text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(sprintf("%s pctl (>= %.2f)", pct_label, cutoff))
})

hist_panel <- ggplot(eip, aes(x = prob)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = qs, linetype = 2, linewidth = 0.8,
             color = c("#332288", "#44AA99", "#882255")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(x = "Edge Inclusion Probability", y = "Frequency",
       title = "EIP Distribution")

jpeg(file.path(figdir, "eip_cutoffs_over_bg.jpeg"),
     height = h * 2, width = w * 2, quality = q)
print((cutoff_panels[[1]] | cutoff_panels[[2]]) /
      (cutoff_panels[[3]] | hist_panel))
dev.off()

cat("\nPlots saved to", figdir, "\n")
cat("Done.\n")
