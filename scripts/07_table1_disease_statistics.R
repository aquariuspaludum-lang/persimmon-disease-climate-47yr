# =============================================================================
# 07_table1_disease_statistics.R
# Table 1: summary statistics of annual mean incidence (non-astringent
# cultivars, 1979-2025 excluding 1985) and Period 1 vs Period 2 comparison.
# Also outputs the powdery mildew share of total incidence (Results 3.1),
# decadal means (1980s vs 2020s), and the annual mean temperature comparison
# between Period 1 and Period 2 (Results 3.5; Welch two-sample t-test).
# Requires: prepared_data/disease_yearly.rds and annual_temp.rds
# (created by 00_data_preparation.R)
# =============================================================================

library(dplyr)

disease_yearly <- readRDS("prepared_data/disease_yearly.rds")

diseases <- c(Powdery_mildew = "Powdery_mildew_mean",
              Angular_leaf_spot = "Angular_leaf_spot_mean",
              Circular_leaf_spot = "Circular_leaf_spot_mean",
              Anthracnose = "Anthracnose_mean")

summarise_disease <- function(col) {
  x  <- disease_yearly[[col]]
  yr <- disease_yearly$year
  ok <- !is.na(x)
  p1 <- mean(x[yr >= 1995 & yr <= 2010], na.rm = TRUE)
  p2 <- mean(x[yr >= 2011 & yr <= 2025], na.rm = TRUE)
  data.frame(n_years = sum(ok), mean = mean(x[ok]), sd = sd(x[ok]),
             min = min(x[ok]), max = max(x[ok]), cv = 100 * sd(x[ok]) / mean(x[ok]),
             period1_mean = p1, period2_mean = p2,
             change_pct = 100 * (p2 - p1) / p1)
}

table1 <- do.call(rbind, lapply(diseases, summarise_disease))
table1 <- cbind(disease = names(diseases), table1)
rownames(table1) <- NULL

cat("=== Table 1: annual mean incidence (%), non-astringent cultivars ===\n")
print(table1, digits = 4)

cat("\nPowdery mildew share of the sum of the four disease means (%):",
    round(100 * table1$mean[1] / sum(table1$mean), 1), "\n")

pm <- disease_yearly[, c("year", "Powdery_mildew_mean")]
cat("Powdery mildew mean, 1980s (1980-1989):",
    round(mean(pm$Powdery_mildew_mean[pm$year >= 1980 & pm$year <= 1989], na.rm = TRUE), 1),
    "| 2020s (2020-2025):",
    round(mean(pm$Powdery_mildew_mean[pm$year >= 2020 & pm$year <= 2025], na.rm = TRUE), 1), "\n")

# ---- Annual mean temperature: decadal means and Period 1 vs Period 2 ----
annual_temp <- readRDS("prepared_data/annual_temp.rds")
tt <- annual_temp$annual_mean_temp
yy <- annual_temp$year
cat("\nAnnual mean temperature, 1980s (1980-1989):", round(mean(tt[yy >= 1980 & yy <= 1989]), 2),
    "| 2020s (2020-2025):", round(mean(tt[yy >= 2020 & yy <= 2025]), 2), "\n")
t1 <- tt[yy >= 1995 & yy <= 2010]   # Period 1
t2 <- tt[yy >= 2011 & yy <= 2025]   # Period 2
tp <- t.test(t2, t1)                # Welch two-sample t-test
cat("Period 1 mean temperature:", round(mean(t1), 2), "C | Period 2:", round(mean(t2), 2),
    "C | difference:", round(mean(t2) - mean(t1), 2), "C | Welch t =", round(tp$statistic, 2),
    ", df =", round(tp$parameter, 1), ", p =", round(tp$p.value, 3), "\n")

# ---- Early summer (May-July) precipitation trend (Table 2C; Results 3.2) ----
ep <- disease_yearly[, c("year", "early_summer_precip")]
pm_fit <- lm(early_summer_precip ~ year, data = ep)
cat("\nEarly summer precipitation: mean", round(mean(ep$early_summer_precip), 1),
    "mm, range", round(range(ep$early_summer_precip), 1),
    "| slope", round(coef(pm_fit)[2], 3), "mm/year, p =", round(summary(pm_fit)$coefficients[2, 4], 3), "\n")

dir.create("results", showWarnings = FALSE)
write.csv(table1, "results/table1_disease_statistics.csv", row.names = FALSE)
