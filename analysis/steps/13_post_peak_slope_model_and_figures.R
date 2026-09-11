## ---- post-peak-slope-model-and-figures ----
slope_control <- nlme::lmeControl(
  opt = "optim", maxIter = 200, msMaxIter = 200,
  niterEM = 50, returnObject = TRUE
)
slope_model <- nlme::lme(
  eGFR ~ time_from_T1_years * cluster_slope,
  random = ~ time_from_T1_years | RKD.ID,
  data = slope_data, method = "REML", na.action = na.exclude,
  control = slope_control
)

# Sensitivity checks for the random-slope and linear-time specifications.
slope_random_intercept_model <- nlme::lme(
  eGFR ~ time_from_T1_years * cluster_slope,
  random = ~ 1 | RKD.ID,
  data = slope_data, method = "REML", na.action = na.exclude,
  control = slope_control
)
slope_random_effect_comparison <- anova(
  slope_random_intercept_model, slope_model
)
slope_linear_ml_model <- update(slope_model, method = "ML")
slope_quadratic_ml_model <- nlme::lme(
  eGFR ~ (time_from_T1_years + I(time_from_T1_years^2)) * cluster_slope,
  random = ~ time_from_T1_years | RKD.ID,
  data = slope_data, method = "ML", na.action = na.exclude,
  control = slope_control
)
slope_linearity_comparison <- anova(
  slope_linear_ml_model, slope_quadratic_ml_model
)

slope_beta <- nlme::fixef(slope_model)
slope_vcov <- vcov(slope_model)
slope_time_term <- "time_from_T1_years"
slope_interaction_term <- grep(
  "time_from_T1_years:cluster_slope|cluster_slope.*:time_from_T1_years",
  names(slope_beta), value = TRUE
)
stopifnot(length(slope_interaction_term) == 1L)
slope_df <- min(summary(slope_model)$tTable[
  c(slope_time_term, slope_interaction_term), "DF"
])

slope_contrast <- function(weights, label) {
  L <- setNames(rep(0, length(slope_beta)), names(slope_beta))
  L[names(weights)] <- weights
  estimate <- sum(L * slope_beta)
  SE <- sqrt(as.numeric(t(L) %*% slope_vcov %*% L))
  critical <- qt(0.975, slope_df)
  tibble(
    estimand = label, estimate, SE, df = slope_df,
    lower_95_CI = estimate - critical * SE,
    upper_95_CI = estimate + critical * SE,
    p_value = 2 * pt(abs(estimate / SE), slope_df, lower.tail = FALSE)
  )
}

slope_results <- bind_rows(
  slope_contrast(
    setNames(1, slope_time_term),
    "Static group annual post-peak eGFR slope"
  ),
  slope_contrast(
    setNames(c(1, 1), c(slope_time_term, slope_interaction_term)),
    "Recovered group annual post-peak eGFR slope"
  ),
  slope_contrast(
    setNames(1, slope_interaction_term),
    "Recovered minus Static slope difference"
  )
) %>%
  mutate(
    unit = "mL/min/1.73 m^2/year",
    estimate_95_CI = sprintf(
      "%.2f (%.2f to %.2f)", estimate, lower_95_CI, upper_95_CI
    ),
    p_value_formatted = vapply(p_value, format_p, character(1))
  )

fixed_effect_table <- as.data.frame(summary(slope_model)$tTable) %>%
  tibble::rownames_to_column("term")

prediction_grid <- expand_grid(
  time_from_T1_years = seq(0, max(slope_data$time_from_T1_years), length.out = 120),
  cluster_slope = factor(
    c("Static group", "Recovered group"),
    levels = c("Static group", "Recovered group")
  )
)
slope_design <- model.matrix(
  ~ time_from_T1_years * cluster_slope, prediction_grid
)
prediction_grid <- prediction_grid %>%
  mutate(
    estimate = as.numeric(slope_design %*% slope_beta),
    SE = sqrt(rowSums((slope_design %*% slope_vcov) * slope_design)),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  )

slope_legend <- slope_data %>%
  distinct(RKD.ID, cluster_slope) %>%
  count(cluster_slope, name = "n") %>%
  transmute(
    cluster_slope,
    label = paste0(cluster_slope, " (n = ", n, ")")
  ) %>%
  { setNames(.$label, .$cluster_slope) }

p_slope_trajectory <- ggplot(
  slope_data,
  aes(time_from_T1_years, eGFR, group = RKD.ID, colour = cluster_slope)
) +
  geom_line(alpha = 0.075, linewidth = 0.28) +
  geom_point(alpha = 0.12, size = 0.62) +
  geom_ribbon(
    data = prediction_grid,
    aes(
      x = time_from_T1_years, ymin = lower, ymax = upper,
      fill = cluster_slope, group = cluster_slope
    ),
    inherit.aes = FALSE, alpha = 0.16, colour = NA
  ) +
  geom_line(
    data = prediction_grid,
    aes(
      time_from_T1_years, estimate,
      colour = cluster_slope, group = cluster_slope
    ),
    inherit.aes = FALSE, linewidth = 1.18
  ) +
  scale_colour_manual(
    values = group_cols,
    breaks = c("Recovered group", "Static group"),
    labels = slope_legend[c("Recovered group", "Static group")]
  ) +
  scale_fill_manual(values = group_cols) +
  scale_x_continuous(
    breaks = 0:5, limits = c(0, max(slope_data$time_from_T1_years)),
    expand = expansion(mult = c(0.01, 0.025))
  ) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 90)) +
  coord_cartesian(ylim = c(0, 90), expand = FALSE) +
  labs(
    title = "Post-peak eGFR trajectories by recovery group",
    x = "Years since T1",
    y = expression(eGFR~"(mL/min/1.73"~m^2*")"),
    colour = NULL, fill = NULL
  ) +
  theme_classic(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "bottom",
    legend.justification = "left",
    legend.key.width = grid::unit(1.6, "lines")
  ) +
  guides(
    fill = "none",
    colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))
  )

forest_data <- slope_results %>%
  mutate(
    row = c(3, 2, 1),
    label = c(
      "Static group slope", "Recovered group slope",
      "Recovered - Static difference"
    ),
    series = c("Static group", "Recovered group", "Difference"),
    result_label = c(
      sprintf(
        "%.2f (%.2f to %.2f); P = %.3f",
        estimate[1], lower_95_CI[1], upper_95_CI[1], p_value[1]
      ),
      sprintf(
        "%.2f (%.2f to %.2f); P = %.3f",
        estimate[2], lower_95_CI[2], upper_95_CI[2], p_value[2]
      ),
      sprintf(
        "%.2f (%.2f to %.2f); interaction P = %.3f",
        estimate[3], lower_95_CI[3], upper_95_CI[3], p_value[3]
      )
    )
  )

p_slope_forest <- ggplot(forest_data, aes(y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_segment(
    aes(
      x = lower_95_CI, xend = upper_95_CI,
      yend = row, colour = series
    ),
    linewidth = 1.05, lineend = "round"
  ) +
  geom_point(aes(x = estimate, colour = series, shape = series), size = 3.4) +
  geom_text(
    aes(x = 5.35, label = result_label),
    hjust = 0, size = 3.35, colour = "black"
  ) +
  annotate(
    "text", x = 5.35, y = 3.72,
    label = "Estimate (95% CI); P value",
    hjust = 0, fontface = "bold", size = 3.45
  ) +
  scale_colour_manual(
    values = c(group_cols, "Difference" = "#222222"), guide = "none"
  ) +
  scale_shape_manual(
    values = c("Static group" = 16, "Recovered group" = 16, "Difference" = 18),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = forest_data$row, labels = forest_data$label,
    limits = c(0.55, 3.85), expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(breaks = seq(-4, 4, 2), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(-4.3, 5.15), clip = "off") +
  labs(
    title = "Annual post-peak eGFR slopes",
    x = expression(
      "Annual eGFR change or difference (mL/min/1.73"~m^2*"/year)"
    ),
    y = NULL
  ) +
  theme_classic(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(8, 245, 10, 12)
  )

p_slope_combined <- (
  p_slope_trajectory /
    patchwork::free(
      p_slope_forest,
      type = "space",
      side = "r"
    )
) +
  patchwork::plot_layout(
    heights = c(1.35, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

p_slope_combined <- p_slope_combined &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 15
    )
  )

print(slope_results %>% select(estimand, estimate_95_CI, p_value_formatted))
print(slope_random_effect_comparison)
print(slope_linearity_comparison)
print(nlme::VarCorr(slope_model))
print(p_slope_combined)

save_table(slope_flow, "post_peak_egfr_slope/RETRACE_post_peak_slope_cohort_flow.csv")
save_table(t1_disposition, "post_peak_egfr_slope/RETRACE_post_peak_T1_disposition.csv")
save_table(t1_by_patient, "post_peak_egfr_slope/RETRACE_post_peak_T1_by_patient.csv")
save_table(dialysis3_audit, "post_peak_egfr_slope/RETRACE_post_peak_dialysis3_audit.csv")
save_table(
  slope_patient_eligibility,
  "post_peak_egfr_slope/RETRACE_post_peak_slope_patient_eligibility.csv"
)
save_table(same_day_audit, "post_peak_egfr_slope/RETRACE_post_peak_same_day_audit.csv")
save_table(t1_tie_audit, "post_peak_egfr_slope/RETRACE_post_peak_T1_tie_audit.csv")
save_table(slope_results, "post_peak_egfr_slope/RETRACE_post_peak_slope_results.csv")
save_table(
  fixed_effect_table,
  "post_peak_egfr_slope/RETRACE_post_peak_slope_fixed_effects.csv"
)
save_table(
  prediction_grid,
  "post_peak_egfr_slope/RETRACE_post_peak_slope_fitted_curves.csv"
)
save_table(
  as.data.frame(slope_random_effect_comparison),
  "post_peak_egfr_slope/RETRACE_post_peak_random_effect_comparison.csv"
)
save_table(
  as.data.frame(slope_linearity_comparison),
  "post_peak_egfr_slope/RETRACE_post_peak_linearity_comparison.csv"
)
saveRDS(slope_model, file.path(slope_output_dir, "RETRACE_post_peak_slope_model.rds"))
save_plot(
  p_slope_trajectory,
  "post_peak_egfr_slope/Figure_post_peak_eGFR_trajectories.png", 10.2, 5.35
)
save_plot(
  p_slope_forest,
  "post_peak_egfr_slope/Figure_post_peak_eGFR_slope_forest.png", 10.2, 4.5
)
save_plot(
  p_slope_combined,
  "post_peak_egfr_slope/Figure_post_peak_eGFR_combined.png", 10.8, 10.3
)
save_plot(
  p_slope_combined,
  "post_peak_egfr_slope/Figure_post_peak_eGFR_combined.pdf", 10.8, 10.3
)
