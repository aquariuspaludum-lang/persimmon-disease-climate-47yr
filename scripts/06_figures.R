# =============================================================================
# 06_figures.R
# Reproduces Figures 1-3 and Supplementary Figures S1-S2.
# Both color and black-and-white (print-safe) versions are generated for
# each figure, consistent with Plant Pathology's colour-figure policy.
# =============================================================================

library(dplyr)
library(ggplot2)
library(scales)

disease_yearly <- readRDS("prepared_data/disease_yearly.rds")
cultivar_means <- readRDS("results/cultivar_annual_means.rds")

dir.create("figures", showWarnings = FALSE)
TEMP_RANGE <- range(disease_yearly$annual_mean_temp, na.rm = TRUE)

# ---------------------------------------------------------------------------
# Figure 1. Long-term disease dynamics (1979-2025), non-astringent cultivars,
# with annual mean temperature on secondary axis.
# ---------------------------------------------------------------------------
scale_factor <- diff(range(disease_yearly$Powdery_mildew_mean, na.rm = TRUE)) /
  diff(TEMP_RANGE)

fig1_data <- disease_yearly %>%
  select(year, Powdery_mildew_mean, Angular_leaf_spot_mean,
         Circular_leaf_spot_mean, Anthracnose_mean, annual_mean_temp) %>%
  tidyr::pivot_longer(cols = c(Powdery_mildew_mean, Angular_leaf_spot_mean,
                                Circular_leaf_spot_mean, Anthracnose_mean),
                       names_to = "disease", values_to = "incidence") %>%
  mutate(disease = recode(disease,
    Powdery_mildew_mean = "Powdery mildew (P. kakicola)",
    Angular_leaf_spot_mean = "Angular leaf spot (C. kaki)",
    Circular_leaf_spot_mean = "Circular leaf spot (M. nawae)",
    Anthracnose_mean = "Anthracnose (C. horii)"))

make_figure1 <- function(temp_color, temp_linetype = "dashed") {
  ggplot(fig1_data, aes(x = year, y = incidence, color = disease, linetype = disease, shape = disease)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    geom_line(data = disease_yearly,
              aes(x = year, y = (annual_mean_temp - TEMP_RANGE[1]) * scale_factor),
              inherit.aes = FALSE, color = temp_color, linetype = temp_linetype, linewidth = 1) +
    scale_y_continuous(
      name = "Disease incidence (%)",
      sec.axis = sec_axis(~ . / scale_factor + TEMP_RANGE[1], name = "Annual mean temperature (\u00b0C)")
    ) +
    labs(x = "Year", color = NULL, linetype = NULL, shape = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# Color version
fig1_color <- make_figure1(temp_color = "#E69F00")
ggsave("figures/Figure1_Disease_Dynamics_COLOR.pdf", fig1_color, width = 9.5, height = 6.5)

# Black-and-white version (grayscale disease lines via manual scale + gray temp line)
fig1_bw <- make_figure1(temp_color = "black", temp_linetype = "longdash") +
  scale_color_manual(values = c("black", "gray45", "gray65", "gray80")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(16, 17, 15, 18))
ggsave("figures/Figure1_Disease_Dynamics_BW.pdf", fig1_bw, width = 9.5, height = 6.5)

# ---------------------------------------------------------------------------
# Figure 2. Temperature effects on powdery mildew.
# ---------------------------------------------------------------------------
fig2a_data <- disease_yearly %>% filter(!is.na(annual_mean_temp))

make_fig2a <- function(point_col, line_col) {
  ggplot(fig2a_data, aes(x = year, y = annual_mean_temp)) +
    geom_point(color = point_col, size = 2) +
    geom_smooth(method = "lm", color = line_col, se = TRUE) +
    labs(x = "Year", y = "Annual mean temperature (\u00b0C)") +
    theme_classic(base_size = 12)
}
make_fig2b <- function(point_col, line_col) {
  ggplot(fig2a_data, aes(x = annual_mean_temp, y = Powdery_mildew_mean)) +
    geom_point(color = point_col, size = 2) +
    geom_smooth(method = "lm", color = line_col, se = TRUE) +
    labs(x = "Annual mean temperature (\u00b0C)", y = "Powdery mildew incidence (%)") +
    theme_classic(base_size = 12)
}

fig2_color <- cowplot::plot_grid(
  make_fig2a("#E69F00", "#8B0000"), make_fig2b("#D55E00", "#08306B"),
  labels = c("A", "B"))
ggsave("figures/Figure2_Temperature_Correlation_COLOR.pdf", fig2_color, width = 10, height = 4.5)

fig2_bw <- cowplot::plot_grid(
  make_fig2a("gray30", "black"), make_fig2b("gray30", "black"),
  labels = c("A", "B"))
ggsave("figures/Figure2_Temperature_Correlation_BW.pdf", fig2_bw, width = 10, height = 4.5)

# ---------------------------------------------------------------------------
# Figure 3. Cultivar comparison bar charts (color and BW).
# ---------------------------------------------------------------------------
table3 <- read.csv("results/table3_cultivar_comparison.csv")

fig3_data <- table3 %>%
  select(disease, non_astringent_mean, non_astringent_se,
         astringent_mean, astringent_se, p_value) %>%
  tidyr::pivot_longer(cols = c(non_astringent_mean, astringent_mean),
                       names_to = "cultivar", values_to = "mean_incidence") %>%
  mutate(
    se = ifelse(cultivar == "non_astringent_mean", non_astringent_se, astringent_se),
    cultivar = recode(cultivar,
                       non_astringent_mean = "Non-astringent (Fuyu)",
                       astringent_mean = "Astringent (Tonewase/Hiratanenashi)"),
    sig_label = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.05, "*", "ns"))
  )

make_figure3 <- function(fill_values) {
  ggplot(fig3_data, aes(x = cultivar, y = mean_incidence, fill = cultivar)) +
    geom_col(color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = mean_incidence - se, ymax = mean_incidence + se), width = 0.15) +
    facet_wrap(~ disease, scales = "free_y") +
    scale_fill_manual(values = fill_values) +
    labs(x = NULL, y = "Disease incidence (%)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(size = 8))
}

fig3_color <- make_figure3(c("#E69F00", "#009E73"))
ggsave("figures/Figure3_Cultivar_Comparison_COLOR.pdf", fig3_color, width = 8, height = 7)

fig3_bw <- make_figure3(c("gray35", "gray80"))
ggsave("figures/Figure3_Cultivar_Comparison_BW.pdf", fig3_bw, width = 8, height = 7)
# NOTE: the version submitted with the manuscript additionally applies a
# diagonal hatch pattern to the astringent bars for extra redundancy; this
# was produced in Python/matplotlib (see figures/make_figure3_bw.py) because
# ggplot2 does not natively support hatched fills. Both versions convey the
# same data.

# ---------------------------------------------------------------------------
# Supplementary Figure S1. Circular leaf spot vs temperature, with
# Cook's-distance-influential points flagged.
# ---------------------------------------------------------------------------
d_cls <- disease_yearly %>% filter(!is.na(Circular_leaf_spot_mean), !is.na(annual_mean_temp))
lm_cls <- lm(Circular_leaf_spot_mean ~ annual_mean_temp, data = d_cls)
d_cls$cooksd <- cooks.distance(lm_cls)
d_cls$influential <- d_cls$cooksd > 4 / nrow(d_cls)

figS1 <- ggplot(d_cls, aes(x = annual_mean_temp, y = Circular_leaf_spot_mean)) +
  geom_smooth(method = "lm", color = "navy", se = TRUE) +
  geom_point(aes(color = influential, shape = influential), size = 3) +
  geom_text(data = filter(d_cls, influential), aes(label = year), vjust = -1, size = 3) +
  scale_color_manual(values = c(`FALSE` = "forestgreen", `TRUE` = "red")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
  labs(x = "Annual mean temperature (\u00b0C)", y = "Circular leaf spot incidence (%)") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
ggsave("figures/FigureS1_Circular_Leaf_Spot_Temperature_COLOR.pdf", figS1, width = 7, height = 5.5)
# Grayscale version: point shape (open vs filled circle) already distinguishes
# influential points independent of color, so this converts cleanly.
figS1_bw <- figS1 + scale_color_manual(values = c(`FALSE` = "gray40", `TRUE` = "black"))
ggsave("figures/FigureS1_Circular_Leaf_Spot_Temperature_BW.pdf", figS1_bw, width = 7, height = 5.5)

# ---------------------------------------------------------------------------
# Supplementary Figure S2. Period comparison (Period 1 vs Period 2).
# ---------------------------------------------------------------------------
period_data <- disease_yearly %>%
  mutate(period = ifelse(year <= 2010, "Period 1\n(1995-2010)", "Period 2\n(2011-2025)")) %>%
  filter(year >= 1995) %>%
  tidyr::pivot_longer(cols = c(Powdery_mildew_mean, Angular_leaf_spot_mean,
                                Circular_leaf_spot_mean, Anthracnose_mean),
                       names_to = "disease", values_to = "incidence") %>%
  group_by(disease, period) %>%
  summarise(mean_incidence = mean(incidence, na.rm = TRUE), .groups = "drop") %>%
  mutate(disease = recode(disease,
    Powdery_mildew_mean = "Powdery mildew", Angular_leaf_spot_mean = "Angular leaf spot",
    Circular_leaf_spot_mean = "Circular leaf spot", Anthracnose_mean = "Anthracnose"))

figS2 <- ggplot(period_data, aes(x = period, y = mean_incidence, fill = period)) +
  geom_col(color = "black", width = 0.6) +
  facet_wrap(~ disease, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("#4C72B0", "#DD8452")) +
  labs(x = NULL, y = "Disease incidence (%)") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
ggsave("figures/FigureS2_Period_Comparison_COLOR.pdf", figS2, width = 10, height = 4)

figS2_bw <- figS2 + scale_fill_manual(values = c("gray35", "gray75"))
ggsave("figures/FigureS2_Period_Comparison_BW.pdf", figS2_bw, width = 10, height = 4)

cat("All figures written to figures/\n")
