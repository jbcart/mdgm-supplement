## Simulation study tables and plots
##
## Usage: Rscript code/tables_plots.R
## (Run from the mdgm-supplement root directory)
##
## Produces faceted comparison plots (bootstrap CIs) for the
## binary (16x16 and 32x32) and Gaussian simulation studies.
## Output filenames match those referenced in the paper.

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

# ---------------------------------------------------------------------------
# Helper: plot a binary simulation figure
# ---------------------------------------------------------------------------
plot_binary <- function(df_summary, pic_name, h = 800, w = 1200) {
  colors_bin <- c("MDGM-AO" = "#88CCEE", "MDGM-ST" = "#117733",
                   "MRF" = "#882255", "aMRF" = "#AA4499")
  q <- 90; ebw <- 0.6; ps <- 0.7; dw <- 0.7; text_size <- 25

  jpeg(pic_name, height = h, width = w, quality = q)
  print(
    ggplot(df_summary) +
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
}

# ===========================================================================
# Binary study
# ===========================================================================
bin_files <- sort(list.files("output",
  pattern = "^p0_.*_\\d+rep\\.rds$", full.names = TRUE))

if (length(bin_files) > 0) {
  cat(sprintf("\n=== Binary study: %d scenario files ===\n", length(bin_files)))

  rows <- list()
  for (f in bin_files) {
    sc_tag <- sub("_\\d+rep\\.rds$", "", basename(f))

    # Extract psi
    psi_match <- regmatches(sc_tag, regexec("^p0_(\\d+)", sc_tag))[[1]]
    psi_val <- as.numeric(paste0("0.", psi_match[2]))

    # Extract grid size
    grid_match <- regmatches(sc_tag, regexec("g(\\d+)$", sc_tag))[[1]]
    grid_val <- as.integer(grid_match[2])

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
          grid   = grid_val,
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

  # Model factor levels
  model_levels <- c("mdgm_ao", "mdgm_st", "mrf_exact", "mrf_pl")
  model_labels <- c("MDGM-AO", "MDGM-ST", "MRF", "aMRF")

  stat_levels <- c("acc_pm", "rmse_eq", "tv")
  stat_labels <- c("Accuracy", "RMSE", "TVD")

  # --- Generate figures per grid size and scenario type ---
  grid_sizes <- sort(unique(df_bin_long$grid))

  for (g in grid_sizes) {
    df_g <- df_bin_long %>% filter(grid == g)

    # Lambda scenarios (missing data)
    df_lambda <- df_g %>% filter(!is.na(lambda))
    if (nrow(df_lambda) > 0) {
      df_summary <- df_lambda %>%
        group_by(psi, lambda, model, stat) %>%
        summarise(compute_boot(value), .groups = "drop")

      df_summary$facet_var <- factor(
        df_summary$lambda,
        levels = sort(unique(df_summary$lambda), decreasing = TRUE),
        labels = paste0("Rate Missing\u2248",
          ifelse(sort(unique(df_summary$lambda), decreasing = TRUE) < 2,
                 "25%", "10%")))

      df_summary$model <- factor(df_summary$model,
        levels = model_levels, labels = model_labels)
      df_summary$stat <- factor(df_summary$stat,
        levels = stat_levels, labels = stat_labels)

      pic_name <- file.path(figdir, sprintf("sim_o1g%d_lambda_bs.jpeg", g))
      plot_binary(df_summary, pic_name)
    }

    # Eta scenarios (complete data)
    df_eta <- df_g %>% filter(!is.na(eta))
    if (nrow(df_eta) > 0) {
      df_summary <- df_eta %>%
        group_by(psi, eta, model, stat) %>%
        summarise(compute_boot(value), .groups = "drop")

      df_summary$facet_var <- factor(
        df_summary$eta,
        levels = sort(unique(df_summary$eta)),
        labels = paste0("\u03B7=", sort(unique(df_summary$eta))))

      df_summary$model <- factor(df_summary$model,
        levels = model_levels, labels = model_labels)
      df_summary$stat <- factor(df_summary$stat,
        levels = stat_levels, labels = stat_labels)

      pic_name <- file.path(figdir, sprintf("sim_o1g%d_error_bs.jpeg", g))
      plot_binary(df_summary, pic_name)
    }
  }
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
