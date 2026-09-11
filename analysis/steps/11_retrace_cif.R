## ---- retrace-cif ----
if (!any(RESULT_crr$status == 1L)) stop("No ESKD events available for the CIF.")

cif_fit <- cmprsk::cuminc(
  ftime = RESULT_crr$time,
  fstatus = RESULT_crr$status,
  group = RESULT_crr$cluster_label,
  cencode = 0
)

cif_df_eskd <- tidy_cuminc(cif_fit, event_code = 1L) %>%
  mutate(group = factor(group, levels = c("Recovered group", "Static group")))

fg_univariable_result <- fine_gray_stats(fit_crr)
fg_annotation <- with(
  fg_univariable_result,
  paste0(
    "Univariable Fine–Gray model\n",
    "Static group vs Recovered group: SHR ", sprintf("%.2f", SHR),
    " (95% CI ", sprintf("%.2f", lower), "–", sprintf("%.2f", upper), "); ",
    format_p(p, prefix = TRUE)
  )
)

p_cif <- ggplot(
  cif_df_eskd,
  aes(time, estimate, colour = group, fill = group)
) +
  geom_step(linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.14, colour = NA, show.legend = FALSE
  ) +
  annotate(
    "text", x = 0.20, y = 0.385, label = fg_annotation,
    hjust = 0, vjust = 1, size = 3.8, lineheight = 1.05
  ) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 0.4), expand = FALSE) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(
    breaks = seq(0, 0.4, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_colour_manual(values = group_cols, drop = FALSE) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(
    title = "Cumulative incidence of ESKD",
    x = "Time in years since diagnosis",
    y = "Cumulative incidence",
    colour = NULL,
    fill = NULL
  ) +
  theme_manuscript() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

print(p_cif)
save_plot(p_cif, "Figure_5_ESKD_CIF.png", 7.2, 5.35)
save_plot(p_cif, "Figure_5_ESKD_CIF.pdf", 7.2, 5.35)
