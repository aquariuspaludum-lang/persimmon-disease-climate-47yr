# =============================================================================
# 01_temperature_trend_analysis.R
# Regional warming trend analysis (Methods 2.3; Results 3.2, Figure 2A)
# =============================================================================

library(dplyr)
library(Kendall)   # Mann-Kendall test
library(trend)     # Sen's slope estimator
library(lmtest)    # Durbin-Watson test

annual_temp <- readRDS("prepared_data/annual_temp.rds")

# ---- Linear regression: temperature ~ year ----
lm_fit <- lm(annual_mean_temp ~ year, data = annual_temp)
lm_summary <- summary(lm_fit)

cat("=== Linear regression: annual mean temperature ~ year ===\n")
cat("Slope:", round(coef(lm_fit)[2], 4), "C/year\n")
cat("R-squared:", round(lm_summary$r.squared, 3), "\n")
cat("p-value:", format.pval(lm_summary$coefficients[2, 4]), "\n")
cat("Total increase (1979-2025):",
    round(coef(lm_fit)[2] * (max(annual_temp$year) - min(annual_temp$year)), 2), "C\n\n")

# ---- Mann-Kendall trend test + Sen's slope (robustness check) ----
mk <- MannKendall(annual_temp$annual_mean_temp)
sens <- sens.slope(annual_temp$annual_mean_temp)

cat("=== Mann-Kendall test ===\n")
cat("tau:", round(mk$tau, 3), " p-value:", format.pval(mk$sl), "\n")
cat("Sen's slope:", round(sens$estimates, 4), "C/year\n\n")

# ---- Durbin-Watson test for autocorrelation of residuals ----
dw <- dwtest(lm_fit)
cat("=== Durbin-Watson test ===\n")
cat("DW =", round(dw$statistic, 2), " p =", round(dw$p.value, 2), "\n")
