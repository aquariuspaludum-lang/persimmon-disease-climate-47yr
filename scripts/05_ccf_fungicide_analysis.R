# =============================================================================
# 05_ccf_fungicide_analysis.R
# Cross-correlation function (CCF) analysis of fungicide application
# frequency (by FRAC group) vs. disease incidence, 1995-2025 (n = 31 years)
# (Methods 2.5 "Cross-correlation function analysis"; Supplementary Table S1)
# =============================================================================

library(dplyr)

disease_yearly <- readRDS("prepared_data/disease_yearly.rds")
fungicide      <- readRDS("prepared_data/fungicide.rds")

merged <- disease_yearly %>%
  filter(year >= 1995, year <= 2025) %>%
  inner_join(fungicide, by = "year") %>%
  arrange(year)

diseases <- c("Powdery_mildew_mean", "Angular_leaf_spot_mean",
              "Circular_leaf_spot_mean", "Anthracnose_mean")
frac_groups <- c("FRAC_3", "FRAC_M1", "FRAC_M2", "FRAC_M3", "FRAC_1M",
                  "FRAC_11", "FRAC_P3", "FRAC_7", "FRAC_U13", "total_fungicide")

# Compute correlation at a specific lag (positive lag: fungicide leads disease
# by `lag` years; negative lag: disease leads fungicide)
lagged_cor <- function(disease, fungicide_var, lag) {
  n <- length(disease)
  if (lag >= 0) {
    x <- fungicide_var[1:(n - lag)]
    y <- disease[(1 + lag):n]
  } else {
    x <- fungicide_var[(1 - lag):n]
    y <- disease[1:(n + lag)]
  }
  if (length(x) < 10) return(c(r = NA, p = NA))
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  c(r = unname(ct$estimate), p = ct$p.value)
}

ccf_table <- list()
for (dis in diseases) {
  for (fg in frac_groups) {
    for (lag in -3:3) {
      res <- lagged_cor(merged[[dis]], merged[[fg]], lag)
      ccf_table[[length(ccf_table) + 1]] <- data.frame(
        disease = dis, frac_group = fg, lag = lag,
        r = res["r"], p_value = res["p"]
      )
    }
  }
}
ccf_results <- do.call(rbind, ccf_table)
rownames(ccf_results) <- NULL

# Multiple-comparison correction: 4 diseases x 10 FRAC groups x 7 lags = 280
# tests were run without correction, which risks a substantial number of
# false positives (~14 expected by chance alone at alpha = 0.05). Apply both
# Bonferroni (conservative, family-wise error rate) and Benjamini-Hochberg
# (FDR, less conservative) corrections across all 280 tests together.
ccf_results$p_bonferroni <- p.adjust(ccf_results$p_value, method = "bonferroni")
ccf_results$p_fdr <- p.adjust(ccf_results$p_value, method = "BH")
ccf_results$significant_uncorrected <- ccf_results$p_value < 0.05
ccf_results$significant_fdr <- ccf_results$p_fdr < 0.05
ccf_results$significant_bonferroni <- ccf_results$p_bonferroni < 0.05

cat("=== CCF associations: uncorrected vs. corrected significance ===\n")
cat("Total tests:", nrow(ccf_results), " (NA r/p, e.g. from zero-variance FRAC series, excluded from counts below)\n")
cat("Significant at uncorrected p < 0.05:", sum(ccf_results$significant_uncorrected, na.rm = TRUE), "\n")
cat("Significant after FDR (Benjamini-Hochberg) correction:", sum(ccf_results$significant_fdr, na.rm = TRUE), "\n")
cat("Significant after Bonferroni correction:", sum(ccf_results$significant_bonferroni, na.rm = TRUE), "\n\n")

cat("=== Associations surviving FDR correction (q < 0.05) ===\n")
print(ccf_results[which(ccf_results$significant_fdr), ], digits = 2)

cat("\n=== Associations surviving Bonferroni correction (p < 0.05) ===\n")
print(ccf_results[which(ccf_results$significant_bonferroni), ], digits = 2)

dir.create("results", showWarnings = FALSE)
write.csv(ccf_results, "results/tableS1_ccf_fungicide_disease.csv", row.names = FALSE)

# ---- Total fungicide application trend over study period ----
total_trend <- lm(total_fungicide ~ year, data = merged)
cat("\n=== Total fungicide application trend, 1995-2025 ===\n")
cat("Mean applications/year:", round(mean(merged$total_fungicide), 1),
    " SD:", round(sd(merged$total_fungicide), 1), "\n")
cat("Linear trend slope:", round(coef(total_trend)[2], 4),
    " p =", round(summary(total_trend)$coefficients[2, 4], 3), "\n")
