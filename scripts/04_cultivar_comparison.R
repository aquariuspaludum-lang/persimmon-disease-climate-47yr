# =============================================================================
# 04_cultivar_comparison.R
# Cultivar-level resistance analysis, 1993-2025 (Methods 2.5 "Cultivar
# resistance analysis"; Table 3, Figure 3)
#
# Analysis is performed on ANNUAL MEANS per cultivar (year-mean basis,
# n = 33 years) rather than raw orchard-level observations, to avoid
# pseudoreplication from repeated sampling of the same orchards across
# years (Mann-Whitney U / Wilcoxon rank-sum test on annual means).
# =============================================================================

library(dplyr)

foliar_orchard <- readRDS("prepared_data/foliar_orchard.rds")
fruit_orchard  <- readRDS("prepared_data/fruit_orchard.rds")

YEAR_MIN <- 1993
YEAR_MAX <- 2025

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp
}

summarise_cultivar_years <- function(data, disease_col, sample_col = "sample_size") {
  data %>%
    filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
    mutate(incidence_pct = 100 * .data[[disease_col]] / .data[[sample_col]]) %>%
    group_by(year, cultivar_type) %>%
    summarise(mean_incidence = mean(incidence_pct, na.rm = TRUE), .groups = "drop")
}

compare_cultivars <- function(annual_means, disease_name) {
  na_vals <- annual_means$mean_incidence[annual_means$cultivar_type == "non_astringent"]
  a_vals  <- annual_means$mean_incidence[annual_means$cultivar_type == "astringent"]

  wt <- wilcox.test(na_vals, a_vals)
  d  <- cohens_d(na_vals, a_vals)

  data.frame(
    disease = disease_name,
    n_years = length(na_vals),
    non_astringent_mean = mean(na_vals), non_astringent_se = sd(na_vals) / sqrt(length(na_vals)),
    astringent_mean = mean(a_vals), astringent_se = sd(a_vals) / sqrt(length(a_vals)),
    pct_reduction = 100 * (mean(na_vals) - mean(a_vals)) / mean(na_vals),
    relative_risk = mean(na_vals) / mean(a_vals),
    cohens_d = d,
    p_value = wt$p.value
  )
}

pm_annual  <- summarise_cultivar_years(foliar_orchard, "Powdery_mildew")
als_annual <- summarise_cultivar_years(foliar_orchard, "Angular_leaf_spot")
cls_annual <- summarise_cultivar_years(foliar_orchard, "Circular_leaf_spot")
anth_annual <- summarise_cultivar_years(fruit_orchard, "Anthracnose")

table3 <- bind_rows(
  compare_cultivars(pm_annual, "Powdery mildew"),
  compare_cultivars(als_annual, "Angular leaf spot"),
  compare_cultivars(cls_annual, "Circular leaf spot"),
  compare_cultivars(anth_annual, "Anthracnose")
)

cat("=== Table 3: Cultivar-specific disease incidence and resistance analysis ===\n")
print(table3, digits = 3)

dir.create("results", showWarnings = FALSE)
write.csv(table3, "results/table3_cultivar_comparison.csv", row.names = FALSE)

# Save the annual-mean datasets too, for figure generation
saveRDS(list(powdery_mildew = pm_annual, angular_leaf_spot = als_annual,
             circular_leaf_spot = cls_annual, anthracnose = anth_annual),
        "results/cultivar_annual_means.rds")
