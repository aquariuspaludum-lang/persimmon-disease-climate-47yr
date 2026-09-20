# Data and code for: Climate-driven disease suppression and broad-spectrum cultivar resistance in Japanese persimmon: A 47-year monitoring study

Manabu Kishi, Yoko Otani, Takuto Hirooka
Wakayama Prefectural Plant Protection Center, Laboratory of Persimmon and Peach, Wakayama, Japan

Submitted to *Plant Pathology*.

## Contents

```
data/       Raw input data (CSV)
scripts/    R analysis scripts, run in numeric order
run_all.R   Master script that runs the full pipeline
prepared_data/  Intermediate datasets (created by scripts/00_data_preparation.R)
results/    Statistical outputs (created when scripts are run)
figures/    Figures, colour and black-and-white versions (created when scripts are run)
```

## Data files (`data/`)

| File | Description |
|---|---|
| `disease_yearly_summary_1979_2025.csv` | Year-level summary of disease incidence and climate variables, 1979-2025. **Retained for reference only; it is not used by the analysis pipeline.** As originally compiled, (i) its `annual_mean_temp` column contained a systematic error (~1.0-1.5°C too high in most years, >2.5°C in 2025), and (ii) its disease incidence columns are pooled across all cultivars rather than restricted to non-astringent cultivars. `scripts/00_data_preparation.R` therefore recomputes annual mean temperature directly from the raw pentad records in `weather_data_1970_2025.csv` (n = 47 years; slope = 0.0230°C/year, R² = 0.336, range 13.7-16.3°C, as reported in the manuscript) and rebuilds the disease incidence series from the orchard-level files below, restricted to non-astringent cultivars. |
| `fruit_damage_rate_powdery_mildew_and_leaf_spot_1979_2025.csv` | Orchard-level (site-level) counts of diseased leaves for powdery mildew, angular leaf spot, and circular leaf spot, with sample size, cultivar type, district, and year. Incidence (%) is calculated as 100 × diseased leaves / sample size, using the sample size recorded for each survey (90 leaves in 1979-1988; 100 leaves from 1989). No records are available for 1985, which is excluded. Used for GLMM analysis (script 02) and cultivar comparison (script 04). |
| `fruit_damage_rate_damage_1979_2025.csv` | Orchard-level counts for anthracnose (and other pests not analyzed in this manuscript), with sample size, cultivar type, district, and year. The sample size recorded for anthracnose is 90 fruits in 1979-1988, 100 fruits in 1989-1993, and 50 fruits from 1994; incidence is calculated with the recorded sample size. |
| `weather_data_1970_2025.csv` | Pentad (5-day) meteorological records (mean/max/min temperature, precipitation) from the Katsuragi weather station, 1970-2025. |
| `fungicide_frac_summary_1995_2025.csv` | Annual counts of recommended fungicide applications by FRAC (Fungicide Resistance Action Committee) mode-of-action group, 1995-2025. FRAC 7 and FRAC U13 were not applied in any year (all zeros). |

## Scripts (`scripts/`), run in order via `run_all.R`

| Script | Purpose | Manuscript section |
|---|---|---|
| `00_data_preparation.R` | Loads raw data, recomputes annual mean temperature from pentad records, rebuilds the non-astringent orchard-level and year-level disease series, merges datasets | Methods 2.1-2.4 |
| `01_temperature_trend_analysis.R` | Linear regression, Mann-Kendall test, Sen's slope, Durbin-Watson test | Methods 2.3; Results 3.2 |
| `02_glmm_temperature_disease.R` | Beta-binomial GLMM (foliar diseases) and negative binomial GLMM (anthracnose), by cultivar type, with DHARMa diagnostics | Methods 2.5; Results 3.2 |
| `03_correlation_analysis.R` | Pearson/Spearman correlations (temperature, precipitation) for all four diseases; Cook's distance screening | Methods 2.5; Results 3.1-3.3; Table 2 |
| `04_cultivar_comparison.R` | Wilcoxon rank-sum tests, relative risk, percent reduction, Cohen's d on annual means (1993-2025) | Methods 2.5; Results 3.4; Table 3 |
| `05_ccf_fungicide_analysis.R` | Cross-correlation function analysis, fungicide use vs. disease incidence, lags -3 to +3 years; Benjamini-Hochberg (FDR) and Bonferroni corrections across the 212 valid tests (of 280 combinations; combinations in which one series is constant are excluded); trend in total fungicide applications | Methods 2.5; Results 3.6; Table S1 |
| `06_figures.R` | Reproduces Figures 1-3 and Supplementary Figures S1-S2 (colour and black-and-white versions) | All figures |
| `07_table1_disease_statistics.R` | Table 1 summary statistics and Period 1 vs. Period 2 comparison; powdery mildew share of total incidence and decadal means; annual mean temperature by decade and Period 1 vs. Period 2 (Welch's t-test); early summer precipitation trend | Results 3.1, 3.2, 3.5; Table 1; Table 2C |

## Requirements

R version 4.3.0 or later (tested with R 4.5.1), with packages: `dplyr`, `tidyr`, `readr`, `ggplot2`, `scales`, `cowplot`, `glmmTMB`, `DHARMa`, `Kendall`, `trend`, `lmtest`.

## Notes

- This code was written to reproduce the analyses exactly as described in the manuscript's Methods section. It was prepared with the assistance of Claude Sonnet 5 (Anthropic) and was executed end-to-end in R 4.5.1 (Windows), reproducing the values reported in the manuscript (e.g., temperature trend: slope = 0.023°C/year, R² = 0.336; powdery mildew GLMM in non-astringent cultivars: β = -1.102 per °C; Table 3 reductions of 71.8-83.6%). Please report any errors encountered when running it. If a version-mismatch warning between `glmmTMB` and `TMB` appears, reinstall `glmmTMB` from source; it did not affect the results.
- One orchard-level record in the powdery mildew file (Hashimotoshi_N_01, 1984), in which the number of diseased leaves (100) exceeded the number examined (90), was corrected to 90 after checking the original field records.
- Figure 3's black-and-white version, as submitted with the manuscript, additionally applies a diagonal hatch pattern to the astringent-cultivar bars for extra visual redundancy; this was generated in Python/matplotlib (not included here) because ggplot2 has no native hatch-fill support. The R version in this repository conveys the same underlying data using grayscale fill only.
- Fungicide application data represent recommended spray schedules from the local agricultural cooperative (JA), not verified actual applications by individual growers.

## License

Code: MIT License. Data: CC BY 4.0. See `LICENSE`.

Zenodo archive: DOI to be added after acceptance.
