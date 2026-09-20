# =============================================================================
# run_all.R
# Master script: runs the full analysis pipeline in order.
# Run this from the repository root (working directory containing
# scripts/, data/, prepared_data/, results/, figures/).
# =============================================================================

required_packages <- c("dplyr", "tidyr", "readr", "ggplot2", "scales", "cowplot",
                        "glmmTMB", "DHARMa", "Kendall", "trend", "lmtest")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))")
}

source("scripts/00_data_preparation.R")
source("scripts/01_temperature_trend_analysis.R")
source("scripts/02_glmm_temperature_disease.R")
source("scripts/03_correlation_analysis.R")
source("scripts/04_cultivar_comparison.R")
source("scripts/05_ccf_fungicide_analysis.R")
source("scripts/06_figures.R")
source("scripts/07_table1_disease_statistics.R")

cat("\n=============================================\n")
cat("Full analysis pipeline complete.\n")
cat("See results/ for statistical outputs (CSV) and figures/ for plots.\n")
cat("=============================================\n")
