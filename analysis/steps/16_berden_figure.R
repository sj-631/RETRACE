## ---- berden-figure ----
berden_data <- RESULT_crr %>%
  transmute(
    cluster_label,
    Berden = factor(
      standardize_berden(Berden.score.on.renal.biopsy),
      levels = berden_levels
    )
  )

berden_tab <- with(berden_data, table(cluster_label, Berden))
berden_available_tab <- berden_tab[, setdiff(colnames(berden_tab), "Missing"), drop = FALSE]
if (nrow(berden_available_tab) >= 2 && ncol(berden_available_tab) >= 2) {
  set.seed(123)
  print(fisher.test(berden_available_tab, simulate.p.value = TRUE, B = 10000))
}

# Missingness is tested from patient-level rows, not from already aggregated
# category rows.
berden_missing_tab <- with(
  berden_data,
  table(cluster_label, Berden_missing = Berden == "Missing" | is.na(Berden))
)
print(berden_missing_tab)
if (all(dim(berden_missing_tab) >= 2)) print(fisher.test(berden_missing_tab))

berden_plot_data <- berden_data %>%
  count(cluster_label, Berden, .drop = FALSE, name = "n") %>%
  group_by(cluster_label) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

p_berden <- ggplot(
  berden_plot_data,
  aes(cluster_label, proportion, fill = Berden)
) +
  geom_col() +
  geom_text(
    aes(label = if_else(proportion > 0, scales::percent(proportion, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 4.2
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = berden_cols, drop = FALSE) +
  labs(
    title = "Berden biopsy class distribution",
    x = NULL,
    y = "Proportion",
    fill = "Berden biopsy class"
  ) +
  theme_manuscript(legend = "right") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.title = element_text(face = "bold")
  )

print(p_berden)
save_plot(p_berden, "Figure_3_Berden_distribution.png", 7.2, 5.2)
save_plot(p_berden, "Figure_3_Berden_distribution.pdf", 7.2, 5.2)
