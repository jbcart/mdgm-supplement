## Simulation study tables and plots
##
## Usage: Rscript code/tables_plots.R
## (Run from the mdgm-supplement root directory)
##
## Produces faceted comparison plots (bootstrap CIs) for both the
## binary and Gaussian simulation studies.

library(tidyverse)

figdir <- file.path("output", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

n_boot <- 1000

compute_boot <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(data.frame(mean = NA, bsp05 = NA, bsp95 = NA))
  bs <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  data.frame(
    mean  = mean(x),
    bsp05 = unname(quantile(bs, 0.05)),
    bsp95 = unname(quantile(bs, 0.95))
  )
}

# ===========================================================================
# Binary study
# ===========================================================================
bin_files <- sort(list.files("output",
  pattern = "^p0_.*_\\d+rep\\.rds$", full.names = TRUE))

if (length(bin_files) > 0) {
  cat(sprintf("\n=== Binary study: %d scenario files ===\n", length(bin_files)))

  bin_model_names <- c("mrf_exact", "mrf_pl", "mdgm_st", "mdgm_ao")

  rows <- list()
  for (f in bin_files) {
    sc_tag <- sub("_\\d+rep\\.rds$", "", basename(f))

    # Parse tag: p0_<PSI>e0_<E>r<R>o1g<G> or p0_<PSI>e0_<E>l<L>o1g<G>
    # Extract psi
    psi_match <- regmatches(sc_tag, regexec("^p0_(\\d+)", sc_tag))[[1]]
    psi_val <- as.numeric(paste0("0.", psi_match[2]))

    # Detect lambda vs fixed-rep scenario
    is_lambda <- grepl("l\\d", sc_tag)
    if (is_lambda) {
      lam_match <- regmatches(sc_tag, regexec("l(\\d+_\\d+)", sc_tag))[[1]]
      lambda_val <- as.numeric(gsub("_", ".", lam_match[2]))
      eta_val <- NA_real_
    } else {
      eta_match <- regmatches(sc_tag, regexec("e0_(\\d+)", sc_tag))[[1]]
      eta_val <- as.numeric(paste0("0.", eta_match[2]))
      lambda_val <- NA_real_
    }

    all_reps <- readRDS(f)

    for (r in seq_along(all_reps)) {
      for (m in names(all_reps[[r]])) {
        v <- all_reps[[r]][[m]]
        if (is.null(v)) next
        rows[[length(rows) + 1L]] <- data.frame(
          psi    = psi_val,
          eta    = eta_val,
          lambda = lambda_val,
          model  = m,
          rep    = r,
          acc_pm   = v$acc_pm,
          eq_pmse  = v$eq_pmse,
          tv       = if (!is.null(v$tv_vs_mrf)) v$tv_vs_mrf else NA_real_,
          acc_na_pm = if (!is.null(v$acc_na_pm)) v$acc_na_pm else NA_real_,
          rhat_psi = v$rhat_psi,
          rhat_eq  = v$rhat_eq,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df_bin_raw <- do.call(rbind, rows)
  cat(sprintf("  %d rows\n", nrow(df_bin_raw)))

  # Filter to converged chains (rhat < 1.1), matching original
  df_bin_raw <- df_bin_raw %>%
    filter(is.na(rhat_psi) | rhat_psi < 1.1,
           is.na(rhat_eq)  | rhat_eq < 1.1)

  # Pivot to long
  df_bin_long <- df_bin_raw %>%
    mutate(rmse_eq = sqrt(eq_pmse)) %>%
    pivot_longer(cols = c(acc_pm, rmse_eq, tv),
                 names_to = "stat", values_to = "value") %>%
    filter(!is.na(value))

  # Determine faceting variable
  has_lambda <- any(!is.na(df_bin_raw$lambda))
  has_eta    <- any(!is.na(df_bin_raw$eta))

  if (has_lambda) {
    df_bin_summary <- df_bin_long %>%
      group_by(psi, lambda, model, stat) %>%
      summarise(compute_boot(value), .groups = "drop")

    df_bin_summary$facet_var <- factor(
      df_bin_summary$lambda,
      levels = sort(unique(df_bin_summary$lambda), decreasing = TRUE),
      labels = paste0("Rate Missing\u2248",
        ifelse(sort(unique(df_bin_summary$lambda), decreasing = TRUE) < 2,
               "25%", "10%")))
  } else {
    df_bin_summary <- df_bin_long %>%
      group_by(psi, eta, model, stat) %>%
      summarise(compute_boot(value), .groups = "drop")

    df_bin_summary$facet_var <- factor(
      df_bin_summary$eta,
      levels = sort(unique(df_bin_summary$eta)),
      labels = paste0("\u03B7=", sort(unique(df_bin_summary$eta))))
  }

  df_bin_summary$model <- factor(df_bin_summary$model,
    levels = c("mdgm_ao", "mdgm_st", "mrf_exact", "mrf_pl"),
    labels = c("MDGM-AO", "MDGM-ST", "MRF", "aMRF"))

  df_bin_summary$stat <- factor(df_bin_summary$stat,
    levels = c("acc_pm", "rmse_eq", "tv"),
    labels = c("Accuracy", "RMSE", "TVD"))

  colors_bin <- c("#88CCEE", "#117733", "#882255", "#AA4499")

  h_bin <- 800
  w_bin <- 1200
  q <- 90
  ebw <- 0.6
  ps  <- 0.7
  dw  <- 0.7
  text_size <- 25

  pic_name <- file.path(figdir, "sim_binary_bs.jpeg")
  jpeg(pic_name, height = h_bin, width = w_bin, quality = q)
  print(
    ggplot(df_bin_summary) +
      geom_errorbar(aes(ymin = bsp05, ymax = bsp95,
                         x = as.factor(psi), color = model),
                    position = position_dodge(width = dw),
                    width = ebw, linewidth = ebw) +
      geom_point(aes(y = mean, x = as.factor(psi), color = model),
                 size = ps, position = position_dodge(width = dw)) +
      facet_wrap(facet_var ~ stat, scales = "free") +
      theme_bw() +
      scale_color_manual(values = colors_bin) +
      theme(text = element_text(size = text_size)) +
      labs(color = "Model", y = "", x = "\u03B2")
  )
  dev.off()
  cat(sprintf("Saved %s\n", pic_name))
} else {
  cat("\nNo binary study output found, skipping.\n")
}

# ===========================================================================
# Gaussian study
# ===========================================================================
gauss_files <- sort(list.files("output",
  pattern = "^gauss_k\\d+_b.+_s0_50_g\\d+_\\d+rep\\.rds$",
  full.names = TRUE))

if (length(gauss_files) > 0) {
  cat(sprintf("\n=== Gaussian study: %d scenario files ===\n", length(gauss_files)))

  gauss_model_names <- c("mdgm_st", "mdgm_ao", "mrf_pl", "bis_pfab")

  rows <- list()
  for (f in gauss_files) {
    sc_tag <- sub("_\\d+rep\\.rds$", "", basename(f))

    parts <- regmatches(sc_tag,
      regexec("gauss_k(\\d+)_b(.+)_s(.+)_g(\\d+)", sc_tag))[[1]]
    k_val    <- as.integer(parts[2])
    beta_val <- as.numeric(gsub("_", ".", parts[3]))

    all_reps <- readRDS(f)

    for (r in seq_along(all_reps)) {
      for (m in gauss_model_names) {
        v <- all_reps[[r]][[m]]
        if (is.null(v)) next
        rows[[length(rows) + 1L]] <- data.frame(
          k     = k_val,
          beta  = beta_val,
          model = m,
          rep   = r,
          ari      = v$ari,
          misclass = v$misclass,
          brier    = v$brier,
          eq_pmse  = v$eq_pmse,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df_gauss_raw <- do.call(rbind, rows)
  cat(sprintf("  %d rows\n", nrow(df_gauss_raw)))

  df_gauss_long <- df_gauss_raw %>%
    mutate(rmse_eq = sqrt(eq_pmse)) %>%
    pivot_longer(cols = c(brier, rmse_eq),
                 names_to = "stat", values_to = "value")

  df_gauss_summary <- df_gauss_long %>%
    group_by(k, beta, model, stat) %>%
    summarise(compute_boot(value), .groups = "drop")

  df_gauss_summary$model <- factor(df_gauss_summary$model,
    levels = c("mdgm_ao", "mdgm_st", "mrf_pl", "bis_pfab"),
    labels = c("MDGM-AO", "MDGM-ST", "aMRF", "bayesImageS"))

  df_gauss_summary$stat <- factor(df_gauss_summary$stat,
    levels = c("brier", "rmse_eq"),
    labels = c("Brier", "RMSE"))

  df_gauss_summary$k_label <- factor(paste0("k = ", df_gauss_summary$k),
    levels = paste0("k = ", sort(unique(df_gauss_summary$k))))

  colors_gauss <- c("#88CCEE", "#117733", "#AA4499", "#882255")

  h_gauss <- 700
  w_gauss <- 1200
  q <- 90
  ebw <- 0.6
  ps  <- 1.5
  dw  <- 0.7
  text_size <- 20

  pic_name <- file.path(figdir, "sim_gaussian_bs.jpeg")
  jpeg(pic_name, height = h_gauss, width = w_gauss, quality = q)
  print(
    ggplot(df_gauss_summary) +
      geom_errorbar(aes(ymin = bsp05, ymax = bsp95,
                         x = as.factor(beta), color = model),
                    position = position_dodge(width = dw),
                    width = ebw, linewidth = ebw) +
      geom_point(aes(y = mean, x = as.factor(beta), color = model),
                 size = ps, position = position_dodge(width = dw)) +
      facet_wrap(stat ~ k_label, scales = "free", nrow = 2) +
      theme_bw() +
      scale_color_manual(values = colors_gauss) +
      theme(text = element_text(size = text_size)) +
      labs(color = "Model", y = "", x = "\u03B2")
  )
  dev.off()
  cat(sprintf("Saved %s\n", pic_name))
} else {
  cat("\nNo Gaussian study output found, skipping.\n")
}
