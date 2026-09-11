## ---- trajectory-figure ----
group_info <- Data_wide %>%
  distinct(RKD.ID, cluster_label) %>%
  count(cluster_label, name = "n") %>%
  mutate(
    pct = round(100 * n / sum(n), 1),
    legend_lab = paste0(cluster_label, " (n = ", n, ", ", pct, "%)")
  )
legend_labs <- setNames(group_info$legend_lab, group_info$cluster_label)

plot_c <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey55") +
  geom_line(
    data = Data_wide,
    aes(time_bin, mean, group = RKD.ID, colour = cluster_label),
    alpha = 0.15, linewidth = 0.35, show.legend = FALSE
  ) +
  geom_line(
    data = cluster_means,
    aes(time_bin, mean_value, colour = cluster_label, group = cluster_label),
    linewidth = 1.4
  ) +
  geom_point(
    data = cluster_means,
    aes(time_bin, mean_value, colour = cluster_label),
    size = 2.4
  ) +
  scale_colour_manual(values = group_cols, labels = legend_labs, drop = FALSE) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c("0-1", ">1-2", ">2-3", ">3-4", ">4-5", ">5-6")
  ) +
  scale_y_continuous(
    breaks = seq(-100, 100, by = 10),
    labels = function(x) ifelse(x > 0, paste0("+", x), x)
  ) +
  coord_cartesian(ylim = c(-50, 50)) +
  labs(
    title = "Early eGFR recovery groups",
    x = "Months since diagnosis",
    y = "Change in eGFR from baseline (mL/min/1.73 m²)",
    colour = NULL
  ) +
  theme_manuscript() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

print(plot_c)
save_plot(plot_c, "Figure_1_egfr_recovery_groups.png", 7.2, 5.2)
save_plot(plot_c, "Figure_1_egfr_recovery_groups.pdf", 7.2, 5.2)
