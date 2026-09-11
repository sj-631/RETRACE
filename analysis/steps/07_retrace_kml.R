## ---- retrace-kml ----
RETRACE_base <- crtn %>%
  mutate(
    dateOfCensored_days = as.numeric(lastRecordedContact - dateOfDiagnosis),
    daysfrmESKD_days = as.numeric(dateOfESKD - dateOfDiagnosis)
  ) %>%
  filter(!is.na(dateOfCensored_days), dateOfCensored_days > 180)

ESKD_EARLY <- RETRACE_base %>%
  filter(!is.na(daysfrmESKD_days), daysfrmESKD_days <= 180) %>%
  distinct(RKD.ID) %>%
  pull(RKD.ID)

kmlDataMean_select <- RETRACE_base %>%
  filter(!RKD.ID %in% ESKD_EARLY) %>%
  mutate(
    time = pmax(Interval.from.diagnosis..months., 0),
    time_bin = pmin(pmax(ceiling(time), 1), 6)
  ) %>%
  group_by(RKD.ID, time_bin) %>%
  summarise(mean = mean(eGFR_delta, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(mean)) %>%
  group_by(RKD.ID) %>%
  filter(n_distinct(time_bin) > 1) %>%
  ungroup()

if (n_distinct(kmlDataMean_select$RKD.ID) < 2) {
  stop("KML cohort has fewer than two patients after filtering.")
}

month_cols <- paste0("month_", 1:6)
kml_wide <- kmlDataMean_select %>%
  mutate(time_bin = factor(time_bin, levels = 1:6, labels = month_cols)) %>%
  pivot_wider(id_cols = RKD.ID, names_from = time_bin, values_from = mean) %>%
  arrange(RKD.ID) %>%
  as.data.frame()

for (column in setdiff(month_cols, names(kml_wide))) kml_wide[[column]] <- NA_real_
kml_wide <- kml_wide[, c("RKD.ID", month_cols)]
id_order <- kml_wide$RKD.ID

kml_obj <- kml::cld(kml_wide, timeInData = 2:ncol(kml_wide))
set.seed(123)
kml::kml(
  kml_obj, nbClusters = 2:5, nbRedrawing = 20,
  parAlgo = kml::parALGO(saveFreq = Inf)
)
cl2 <- getClusters(kml_obj, 2)

cluster_df_raw <- tibble(RKD.ID = id_order, cluster_raw = as.integer(cl2))
latest_by_id <- kmlDataMean_select %>%
  group_by(RKD.ID) %>%
  slice_max(time_bin, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(cluster_df_raw, by = "RKD.ID")

# Name clusters by the mean change at each patient's final available month.
cluster_order <- latest_by_id %>%
  group_by(cluster_raw) %>%
  summarise(latest_mean_eGFR_delta = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(latest_mean_eGFR_delta)) %>%
  mutate(
    cluster = row_number(),
    cluster_name = c("Recovered group", "Static group")[cluster]
  )

cluster_df <- cluster_df_raw %>%
  left_join(cluster_order, by = "cluster_raw") %>%
  select(RKD.ID, cluster, cluster_raw, cluster_name)

cluster_count <- cluster_df %>%
  count(cluster, cluster_name, name = "n") %>%
  arrange(cluster) %>%
  mutate(
    prop = round(100 * n / sum(n), 1),
    cluster_label = cluster_name
  )

cluster_df <- cluster_df %>%
  left_join(
    cluster_count %>% select(cluster, cluster_name, n, prop, cluster_label),
    by = c("cluster", "cluster_name")
  )

Data_wide <- kmlDataMean_select %>%
  left_join(cluster_df, by = "RKD.ID") %>%
  mutate(cluster_label = factor(
    cluster_label,
    levels = c("Recovered group", "Static group")
  ))

cluster_means <- Data_wide %>%
  group_by(cluster, cluster_label, time_bin) %>%
  summarise(mean_value = mean(mean, na.rm = TRUE), .groups = "drop")

print(cluster_count)
stopifnot(
  n_distinct(RETRACE_base$RKD.ID) == 285L,
  length(ESKD_EARLY) == 34L,
  n_distinct(kmlDataMean_select$RKD.ID) == 234L,
  setequal(cluster_count$n, c(97L, 137L))
)
