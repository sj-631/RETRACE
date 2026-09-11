## ---- cluster-centroids ----
cluster_centroids_long <- Data_wide %>%
  transmute(
    RKD.ID,
    cluster_label = as.character(cluster_label),
    time_bin,
    eGFR_change = num_clean(mean)
  ) %>%
  group_by(RKD.ID, cluster_label, time_bin) %>%
  summarise(
    eGFR_change = if (all(is.na(eGFR_change))) NA_real_ else mean(eGFR_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cluster_label, time_bin) %>%
  summarise(
    centroid = mean(eGFR_change, na.rm = TRUE),
    n_available = sum(!is.na(eGFR_change)),
    .groups = "drop"
  ) %>%
  mutate(centroid = if_else(is.nan(centroid), NA_real_, centroid)) %>%
  arrange(cluster_label, time_bin)

cluster_centroids_wide <- cluster_centroids_long %>%
  select(cluster_label, time_bin, centroid) %>%
  mutate(time_bin = paste0("month_bin_", time_bin)) %>%
  pivot_wider(names_from = time_bin, values_from = centroid)

save_table(cluster_centroids_long, "eGFR_cluster_centroids_long.csv")
save_table(cluster_centroids_wide, "eGFR_cluster_centroids_wide.csv")
save_table(kml_wide, "kml_input.csv")
