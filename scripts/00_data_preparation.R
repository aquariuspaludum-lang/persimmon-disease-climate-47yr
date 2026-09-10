# =============================================================================
# 00_data_preparation.R
# Kishi et al. -- Climate-driven disease suppression and broad-spectrum
# cultivar resistance in Japanese persimmon: A 47-year monitoring study
#
# Purpose: Load raw data files, recompute annual mean temperature directly
#          from raw pentad weather records, and prepare orchard-level and
#          year-level datasets used by all subsequent analysis scripts.
#
# IMPORTANT DATA NOTES:
#   1. The 'annual_mean_temp' column in disease_yearly_summary_1979_2025.csv
#      was found to contain a systematic error (values ~1.0-1.5C too high in
#      most years, and >2.5C too high in 2025), likely from an earlier,
#      incorrect aggregation script. This script RECOMPUTES annual mean
#      temperature directly from the raw pentad records in
#      weather_data_1970_2025.csv, which matches the statistics reported in
#      the manuscript (slope = 0.0237 C/year, R2 = 0.344, range 13.7-16.3C).
#   2. disease_yearly_summary_1979_2025.csv's disease incidence columns
#      (e.g. Powdery_mildew_mean) are POOLED across all cultivars, not
#      non-astringent (Fuyu)-only as the manuscript's Figure 1/2, Table 2,
#      and Supplementary Figures S1-S2 require. This script does not use
#      that file at all; disease incidence series are rebuilt directly from
#      the orchard-level raw data, filtered to non-astringent orchards only
#      (see Section 2 below).
#   Do NOT use disease_yearly_summary_1979_2025.csv for any purpose in this
#   pipeline -- it is retained in data/ only as a historical/reference file.
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)

data_dir <- "data"

# ---- 1. Raw weather data (pentad-level) -> correct annual mean temperature ----

weather_raw <- read_csv(file.path(data_dir, "weather_data_1970_2025.csv"),
                         locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
names(weather_raw) <- gsub("^\uFEFF", "", names(weather_raw))  # strip BOM if present

weather_raw <- weather_raw %>%
  mutate(
    tem_mean_num      = as.numeric(tem_mean),
    precipitation_num = as.numeric(precipitation),
    month = sub("_Pentad.*$", "", period)
  )

annual_temp <- weather_raw %>%
  filter(year >= 1979, year <= 2025) %>%
  group_by(year) %>%
  summarise(annual_mean_temp = mean(tem_mean_num, na.rm = TRUE), .groups = "drop")

# Early-summer (May-Jul) cumulative precipitation, per manuscript section 2.3
early_summer_precip <- weather_raw %>%
  filter(year >= 1979, year <= 2025, month %in% c("May", "June", "July")) %>%
  group_by(year) %>%
  summarise(early_summer_precip = sum(precipitation_num, na.rm = TRUE), .groups = "drop")

# ---- 2. Year-level disease summary for non-astringent (Fuyu) cultivars ----
#
# IMPORTANT: Figure 1, Figure 2, Table 2, Supplementary Figures S1-S2, and the
# CCF analysis (script 05) all report statistics for NON-ASTRINGENT (Fuyu)
# cultivars specifically (see manuscript Table 1 caption and Methods 2.2).
# The pre-computed 'disease_yearly_summary_1979_2025.csv' file's *_mean
# columns were found to be POOLED across all cultivars (non-astringent +
# astringent combined), NOT non-astringent-only -- this pools in the much
# lower astringent-cultivar incidence and dilutes the true non-astringent
# values (e.g. powdery mildew max 67.4% pooled vs 79.2% non-astringent-only,
# which matches the manuscript's Table 1 exactly). This script instead
# rebuilds the year-level disease series directly from the orchard-level
# (site-level) raw data, filtered to non-astringent (Fuyu) orchards only,
# matching the manuscript's reported statistics.

foliar_orchard_raw <- read_csv(
  file.path(data_dir, "fruit_damage_rate_powdery_mildew_and_leaf_spot_1979_2025.csv"),
  show_col_types = FALSE)
names(foliar_orchard_raw) <- gsub("^\uFEFF", "", names(foliar_orchard_raw))

fruit_orchard_raw <- read_csv(
  file.path(data_dir, "fruit_damage_rate_damage_1979_2025.csv"),
  show_col_types = FALSE)
names(fruit_orchard_raw) <- gsub("^\uFEFF", "", names(fruit_orchard_raw))

foliar_na <- foliar_orchard_raw %>%
  filter(year != 1985, cultivar_type == "Non_astringent") %>%
  filter(Powdery_mildew <= sample_size, Angular_leaf_spot <= sample_size,
         Circular_leaf_spot <= sample_size) %>%  # exclude data-entry errors (diseased > sample_size); see 02_glmm script for detail
  mutate(
    Powdery_mildew_pct    = 100 * Powdery_mildew    / sample_size,
    Angular_leaf_spot_pct = 100 * Angular_leaf_spot / sample_size,
    Circular_leaf_spot_pct = 100 * Circular_leaf_spot / sample_size
  ) %>%
  group_by(year) %>%
  summarise(
    Powdery_mildew_mean    = mean(Powdery_mildew_pct, na.rm = TRUE),
    Angular_leaf_spot_mean = mean(Angular_leaf_spot_pct, na.rm = TRUE),
    Circular_leaf_spot_mean = mean(Circular_leaf_spot_pct, na.rm = TRUE),
    .groups = "drop"
  )

fruit_na <- fruit_orchard_raw %>%
  filter(year != 1985, cultivar_type == "Non_astringent") %>%
  mutate(Anthracnose_pct = 100 * Anthracnose / sample_size) %>%
  group_by(year) %>%
  summarise(Anthracnose_mean = mean(Anthracnose_pct, na.rm = TRUE), .groups = "drop")

disease_yearly <- foliar_na %>%
  full_join(fruit_na, by = "year") %>%
  left_join(annual_temp, by = "year") %>%
  left_join(early_summer_precip, by = "year") %>%
  arrange(year)

# Sanity check against manuscript Table 1 (non-astringent, 1979-2025):
# Powdery mildew mean ~26.3%, max 79.22%
stopifnot(abs(max(disease_yearly$Powdery_mildew_mean, na.rm = TRUE) - 79.22) < 0.5)

# ---- 3. Orchard-level data (for GLMM and cultivar comparison) ----
# Reuses the raw orchard-level tables already loaded above (all cultivars).

foliar_orchard <- foliar_orchard_raw %>%
  filter(year != 1985) %>%
  mutate(cultivar_type = recode(cultivar_type,
                                 "Non_astringent" = "non_astringent",
                                 "Astringent"     = "astringent")) %>%
  left_join(annual_temp, by = "year")

fruit_orchard <- fruit_orchard_raw %>%
  filter(year != 1985) %>%
  mutate(cultivar_type = recode(cultivar_type,
                                 "Non_astringent" = "non_astringent",
                                 "Astringent"     = "astringent")) %>%
  left_join(annual_temp, by = "year")

# ---- 4. Fungicide application data (1995-2025) ----

fungicide <- read_csv(file.path(data_dir, "fungicide_frac_summary_1995_2025.csv"),
                       show_col_types = FALSE)

# ---- 5. Save prepared objects for downstream scripts ----

dir.create("prepared_data", showWarnings = FALSE)
saveRDS(disease_yearly,  "prepared_data/disease_yearly.rds")
saveRDS(foliar_orchard,  "prepared_data/foliar_orchard.rds")
saveRDS(fruit_orchard,   "prepared_data/fruit_orchard.rds")
saveRDS(fungicide,       "prepared_data/fungicide.rds")
saveRDS(annual_temp,     "prepared_data/annual_temp.rds")

cat("Data preparation complete.\n")
cat("Years covered (disease_yearly):", range(disease_yearly$year), "\n")
cat("Recomputed annual_mean_temp range:", round(range(annual_temp$annual_mean_temp), 2), "\n")
