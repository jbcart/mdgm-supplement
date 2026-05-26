## Aggregate cross-validation results for disorder analysis
##
## Usage: Rscript code/disorder/02_cv_aggregate.R
## (Run from the mdgm-supplement root directory)

cvdir <- file.path("output", "disorder_cv")
files <- list.files(cvdir, pattern = "^run_cv_.*\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No CV result files found in ", cvdir)
}

cat(sprintf("Aggregating %d CV runs from %s\n", length(files), cvdir))

results <- do.call(rbind, lapply(files, function(f) {
  as.data.frame(t(readRDS(f)))
}))

# Summary: mean and SD of each metric
summary_df <- data.frame(
  metric = names(results),
  mean   = colMeans(results),
  sd     = apply(results, 2, sd),
  row.names = NULL
)

cat("\nCross-Validation Summary:\n")
print(summary_df, digits = 4)

# Fraction where MDGM-ST has lower MAE than MRF-PL
cat(sprintf("\nFraction pp1 MAE: ST < PL: %.3f\n",
            mean(results$mae_pp1_st < results$mae_pp1_pl)))
cat(sprintf("Fraction zpm MAE: ST < PL: %.3f\n",
            mean(results$mae_zpm_st < results$mae_zpm_pl)))

# Save
outfile <- file.path(cvdir, "cv_summary.csv")
write.csv(summary_df, outfile, row.names = FALSE)
cat(sprintf("\nSaved %s\n", outfile))
