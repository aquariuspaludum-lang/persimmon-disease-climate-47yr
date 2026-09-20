# =============================================================================
# 02_glmm_temperature_disease.R
# Beta-binomial / negative binomial GLMMs of disease incidence vs. annual
# mean temperature (Methods 2.5 "Generalized linear mixed models...")
#
# Model structure (foliar diseases, beta-binomial):
#   logit(p_ij) = b0 + b1 * Temperature_j + u_j + v_i
#     u_j ~ N(0, sigma2_year), v_i ~ N(0, sigma2_orchard)
#
# Model structure (anthracnose, negative binomial, with offset):
#   log(mu_ij) = b0 + b1 * Temperature_j + u_j + v_i + log(N_ij)
# =============================================================================

library(dplyr)
library(glmmTMB)
library(DHARMa)

foliar_orchard <- readRDS("prepared_data/foliar_orchard.rds")
fruit_orchard  <- readRDS("prepared_data/fruit_orchard.rds")

foliar_diseases <- c("Powdery_mildew", "Angular_leaf_spot", "Circular_leaf_spot")
cultivars <- c("non_astringent", "astringent")

fit_foliar_glmm <- function(data, disease_col, cultivar) {
  d <- data %>%
    filter(cultivar_type == cultivar) %>%
    mutate(
      diseased = .data[[disease_col]],
      healthy  = sample_size - diseased,
      temp_scaled = as.numeric(scale(annual_mean_temp))
    ) %>%
    filter(!is.na(diseased), !is.na(annual_mean_temp), sample_size > 0)

  if (nrow(d) < 20) return(NULL)

  ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))

  # Attempt 1: full model (year + orchard random effects), scaled temperature
  model <- tryCatch(
    glmmTMB(
      cbind(diseased, healthy) ~ temp_scaled + (1 | year) + (1 | site_id),
      family = betabinomial(link = "logit"),
      data = d, control = ctrl
    ),
    error = function(e) { attr(e, "d") <- d; e }
  )
  if (inherits(model, "error")) {
    message("  [", disease_col, "/", cultivar, "] full model failed: ", conditionMessage(model))
    # Attempt 2: drop the orchard random effect (often singular when many
    # orchards each contribute few observations per year)
    model <- tryCatch(
      glmmTMB(
        cbind(diseased, healthy) ~ temp_scaled + (1 | year),
        family = betabinomial(link = "logit"),
        data = d, control = ctrl
      ),
      error = function(e) { attr(e, "d") <- d; e }
    )
    if (inherits(model, "error")) {
      message("  [", disease_col, "/", cultivar, "] year-only model also failed: ", conditionMessage(model))
      return(NULL)
    }
    message("  [", disease_col, "/", cultivar, "] converged after dropping the orchard random effect (site_id variance was likely near-zero/singular).")
  }
  model
}

fit_anthracnose_glmm <- function(data, cultivar) {
  d <- data %>%
    filter(cultivar_type == cultivar) %>%
    mutate(temp_scaled = as.numeric(scale(annual_mean_temp))) %>%
    filter(!is.na(Anthracnose), !is.na(annual_mean_temp), sample_size > 0)

  if (nrow(d) < 20) return(NULL)

  ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
  model <- tryCatch(
    glmmTMB(
      Anthracnose ~ temp_scaled + (1 | year) + (1 | site_id) + offset(log(sample_size)),
      family = nbinom2(link = "log"),
      data = d, control = ctrl
    ),
    error = function(e) { message("  [Anthracnose/", cultivar, "] failed: ", conditionMessage(e)); NULL }
  )
  model
}

results <- list()

for (disease in foliar_diseases) {
  for (cv in cultivars) {
    key <- paste(disease, cv, sep = "_")
    cat("Fitting:", key, "...\n")
    m <- fit_foliar_glmm(foliar_orchard, disease, cv)
    results[[key]] <- m
    if (!is.null(m)) {
      s <- summary(m)
      coef_name <- if ("temp_scaled" %in% rownames(s$coefficients$cond)) "temp_scaled" else "annual_mean_temp"
      beta_scaled <- s$coefficients$cond[coef_name, "Estimate"]
      temp_sd <- sd(foliar_orchard$annual_mean_temp[foliar_orchard$cultivar_type == cv], na.rm = TRUE)
      cat("  beta (temperature, per SD) =", round(beta_scaled, 3),
          " | beta (per 1C, for comparison with manuscript) =", round(beta_scaled / temp_sd, 3),
          " p =", format.pval(s$coefficients$cond[coef_name, "Pr(>|z|)"]), "\n")
      cat("  Random effects:\n"); print(VarCorr(m))
    } else {
      cat("  Model failed to converge even after fallback attempts. See message() output above for details.\n")
    }
  }
}

for (cv in cultivars) {
  key <- paste("Anthracnose", cv, sep = "_")
  cat("Fitting:", key, "...\n")
  m <- fit_anthracnose_glmm(fruit_orchard, cv)
  results[[key]] <- m
  if (!is.null(m)) {
    s <- summary(m)
    cat("  beta (temperature, scaled) =", round(s$coefficients$cond["temp_scaled", "Estimate"], 3),
        " p =", format.pval(s$coefficients$cond["temp_scaled", "Pr(>|z|)"]), "\n")
  }
}

# ---- Model diagnostics (DHARMa) for the primary reported model: ----
# Powdery mildew, non-astringent cultivars (main text: beta = -1.102, p < 0.001)
if (!is.null(results[["Powdery_mildew_non_astringent"]])) {
  sim_res <- simulateResiduals(results[["Powdery_mildew_non_astringent"]])
  cat("\n=== DHARMa diagnostics: Powdery mildew, non-astringent ===\n")
  print(testUniformity(sim_res))
  print(testDispersion(sim_res))
}

dir.create("results", showWarnings = FALSE)
saveRDS(results, "results/glmm_models.rds")
cat("\nGLMM fitting complete. Models saved to results/glmm_models.rds\n")
