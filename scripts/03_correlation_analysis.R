# =============================================================================
# 03_correlation_analysis.R
# Pearson / Spearman correlations between disease incidence and
# environmental variables (Methods 2.5 "Correlation analyses"; Table 2)
# Includes Cook's distance influential-point screening (used for
# circular leaf spot; Supplementary Figure S1).
# =============================================================================

library(dplyr)

disease_yearly <- readRDS("prepared_data/disease_yearly.rds")

diseases <- c("Powdery_mildew_mean", "Angular_leaf_spot_mean",
              "Circular_leaf_spot_mean", "Anthracnose_mean")

run_correlation <- function(data, disease_col, env_col) {
  d <- data %>% filter(!is.na(.data[[disease_col]]), !is.na(.data[[env_col]]))
  pear <- cor.test(d[[disease_col]], d[[env_col]], method = "pearson")
  spear <- cor.test(d[[disease_col]], d[[env_col]], method = "spearman")
  data.frame(
    disease = disease_col, env_var = env_col, n = nrow(d),
    pearson_r = unname(pear$estimate), pearson_p = pear$p.value,
    spearman_rho = unname(spear$estimate), spearman_p = spear$p.value
  )
}

cat("=== Temperature correlations (Table 2A) ===\n")
temp_results <- do.call(rbind, lapply(diseases, run_correlation,
                                       data = disease_yearly, env_col = "annual_mean_temp"))
print(temp_results, digits = 3)

cat("\n=== Early summer precipitation correlations (Table 2B) ===\n")
precip_results <- do.call(rbind, lapply(diseases, run_correlation,
                                         data = disease_yearly, env_col = "early_summer_precip"))
print(precip_results, digits = 3)

# ---- Cook's distance screening: circular leaf spot vs temperature ----
d_cls <- disease_yearly %>%
  filter(!is.na(Circular_leaf_spot_mean), !is.na(annual_mean_temp))

lm_cls <- lm(Circular_leaf_spot_mean ~ annual_mean_temp, data = d_cls)
cooksd <- cooks.distance(lm_cls)
threshold <- 4 / nrow(d_cls)
influential <- d_cls$year[cooksd > threshold]

cat("\n=== Cook's distance screening: circular leaf spot ===\n")
cat("Threshold (4/n):", round(threshold, 4), "\n")
cat("Influential years:", paste(influential, collapse = ", "), "\n")

# Re-run correlation excluding influential points
d_cls_clean <- d_cls %>% filter(!year %in% influential)
pear_clean <- cor.test(d_cls_clean$Circular_leaf_spot_mean, d_cls_clean$annual_mean_temp,
                        method = "pearson")
spear_clean <- cor.test(d_cls_clean$Circular_leaf_spot_mean, d_cls_clean$annual_mean_temp,
                         method = "spearman")
cat("After removing influential points:\n")
cat("  Pearson r =", round(pear_clean$estimate, 3), " p =", format.pval(pear_clean$p.value), "\n")
cat("  Spearman rho =", round(spear_clean$estimate, 3), " p =", format.pval(spear_clean$p.value), "\n")

dir.create("results", showWarnings = FALSE)
write.csv(temp_results, "results/table2A_temperature_correlations.csv", row.names = FALSE)
write.csv(precip_results, "results/table2B_precipitation_correlations.csv", row.names = FALSE)
