## ---- ect-figures ----
time_var <- "Interval.from.diagnosis..months."

overall_grid <- setNames(
  data.frame(seq(0, 6, length.out = 121)),
  time_var
)
eff_overall_df <- predict_fixed_lme(lm2, overall_grid)

# Extract the three spline slopes and P values directly from lm2.
lm2_table <- summary(lm2)$tTable
slope_rows <- grep("^lspline", rownames(lm2_table))
stopifnot(length(slope_rows) == 3L)
overall_stats <- tibble(
  x = c(1, 3, 5),
  interval = c("0-2 months", "2-4 months", "4-6 months"),
  slope = lm2_table[slope_rows, "Value"],
  p = lm2_table[slope_rows, "p-value"],
  label = paste0(
    interval, "\nSlope = ", sprintf("%.2f", slope), "/month\n",
    vapply(p, format_p, character(1), prefix = TRUE)
  )
)

p_overall <- ggplot(eff_overall_df, aes(x = .data[[time_var]], y = fit)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_vline(xintercept = c(2, 4), linetype = "dotted", linewidth = 0.35, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.55) +
  geom_line(linewidth = 1.35, colour = "#0072B2", lineend = "round") +
  geom_label(
    data = overall_stats,
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.15,
    size = 2.8,
    lineheight = 0.95,
    label.padding = grid::unit(0.14, "lines"),
    label.r = grid::unit(0.08, "lines"),
    linewidth = 0.2,
    fill = scales::alpha("white", 0.92)
  ) +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = scales::pretty_breaks(6), expand = expansion(mult = c(0.08, 0.32))) +
  labs(x = "Months since diagnosis", y = "Change in eGFR from baseline (mL/min/1.73 m²)") +
  theme_manuscript()

induction_grid <- expand_grid(
  time_value = seq(0, 6, length.out = 121),
  inductionFinal = factor(
    c("Cyclo", "Rituximab", "Other"),
    levels = c("Cyclo", "Rituximab", "Other")
  )
) %>%
  rename(!!time_var := time_value)

eff_ind_df <- predict_fixed_lme(lm7, induction_grid) %>%
  mutate(
    induction_label = recode(
      as.character(inductionFinal),
      Cyclo = "Cyclophosphamide",
      Rituximab = "Rituximab",
      Other = "Other"
    ),
    induction_label = factor(
      induction_label,
      levels = c("Cyclophosphamide", "Rituximab", "Other")
    )
  )

p_induction <- ggplot(
  eff_ind_df,
  aes(
    x = .data[[time_var]], y = fit,
    colour = induction_label, fill = induction_label, group = induction_label
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_vline(xintercept = c(2, 4), linetype = "dotted", linewidth = 0.35, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.25, lineend = "round") +
  scale_colour_manual(values = induction_cols, name = NULL) +
  scale_fill_manual(values = induction_cols, guide = "none") +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = scales::pretty_breaks(6), expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Months since diagnosis", y = expression(Delta~eGFR~"(mL/min/1.73"~m^2*")")) +
  theme_manuscript()

Figure2_combined <- (p_overall / p_induction) +
  plot_layout(heights = c(1, 1.05), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 14))

print(Figure2_combined)
save_plot(p_overall, "Figure_2_overall_eGFR_recovery.png", 6.8, 4.6)
save_plot(p_induction, "Figure_induction_treatment_eGFR_lm7.png", 6.8, 4.8)
save_plot(Figure2_combined, "Figure_2_combined.png", 7.2, 8.5)
ggsave(
  file.path(output_dir, "Figure_2_combined.pdf"),
  Figure2_combined, width = 7.2, height = 8.5, bg = "white"
)
