## ---- supplementary-full-and-excluded-cif ----

cif_followup_years <- 5
prediagnosis_eskd_policy <- "time_zero"
stopifnot(prediagnosis_eskd_policy %in% c("time_zero", "exclude", "stop"))

cif_output_dir <- file.path(output_dir, "cif_full_and_excluded")
dir.create(cif_output_dir, recursive = TRUE, showWarnings = FALSE)

full_cif_fields <- c(
  "dateOfDiagnosis", "dateOfESKD", "dateOfDeath",
  "lastRecordedContact", "ESKDStandard", "death"
)
assert_columns(Full_cohort, full_cif_fields)

# Confirm patient-level dates/statuses do not conflict across visit rows.
full_cif_field_audit <- Full_cohort %>%
  filter(!is.na(RKD.ID)) %>%
  group_by(RKD.ID) %>%
  summarise(
    across(
      all_of(full_cif_fields),
      ~ n_distinct(.x[!is.na(.x) & trimws(as.character(.x)) != ""]),
      .names = "n_values_{.col}"
    ),
    .groups = "drop"
  )

full_cif_field_conflicts <- full_cif_field_audit %>%
  filter(if_any(starts_with("n_values_"), ~ .x > 1))
write.csv(
  full_cif_field_conflicts,
  file.path(cif_output_dir, "CIF_patient_field_conflicts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (nrow(full_cif_field_conflicts)) {
  stop("Conflicting patient-level CIF fields; review CIF_patient_field_conflicts.csv.")
}

full_cif_patient <- Full_cohort %>%
  filter(!is.na(RKD.ID)) %>%
  arrange(RKD.ID, Interval.from.diagnosis..months., Date.Of.Visit) %>%
  group_by(RKD.ID) %>%
  summarise(across(all_of(full_cif_fields), first_nonblank), .groups = "drop")

retrace_cif_ids <- RESULT_crr$RKD.ID
stopifnot(
  !anyDuplicated(retrace_cif_ids),
  all(retrace_cif_ids %in% full_cif_patient$RKD.ID)
)

full_cif_patient <- full_cif_patient %>%
  mutate(
    trajectory_eligibility = factor(
      if_else(RKD.ID %in% retrace_cif_ids, "RETRACE included", "RETRACE excluded"),
      levels = c("RETRACE included", "RETRACE excluded")
    ),
    eskd_flag = coalesce(toupper(trimws(ESKDStandard)) == "YES", FALSE),
    death_flag = coalesce(toupper(trimws(death)) == "DEAD", FALSE),
    days_to_eskd_raw = as.numeric(dateOfESKD - dateOfDiagnosis),
    days_to_death_raw = as.numeric(dateOfDeath - dateOfDiagnosis),
    days_to_censor_raw = as.numeric(lastRecordedContact - dateOfDiagnosis),
    prediagnosis_eskd = eskd_flag & !is.na(days_to_eskd_raw) & days_to_eskd_raw < 0,
    eskd_after_last_contact = eskd_flag &
      !is.na(dateOfESKD) & !is.na(lastRecordedContact) &
      dateOfESKD > lastRecordedContact
  )

prediagnosis_eskd_records <- full_cif_patient %>%
  filter(prediagnosis_eskd) %>%
  select(
    RKD.ID, trajectory_eligibility,
    dateOfDiagnosis, dateOfESKD, days_to_eskd_raw
  )
write.csv(
  prediagnosis_eskd_records,
  file.path(cif_output_dir, "CIF_prediagnosis_ESKD_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (prediagnosis_eskd_policy == "stop" && nrow(prediagnosis_eskd_records)) {
  stop("ESKD dates before diagnosis require adjudication.")
}

full_cif_patient <- full_cif_patient %>%
  mutate(
    days_to_eskd_used = case_when(
      prediagnosis_eskd & prediagnosis_eskd_policy == "time_zero" ~ 0,
      prediagnosis_eskd ~ NA_real_,
      TRUE ~ days_to_eskd_raw
    ),
    time_to_eskd_days = if_else(
      eskd_flag & !is.na(days_to_eskd_used) & days_to_eskd_used >= 0,
      days_to_eskd_used, Inf
    ),
    time_to_death_days = if_else(
      death_flag & !is.na(days_to_death_raw) & days_to_death_raw >= 0,
      days_to_death_raw, Inf
    ),
    time_to_censor_days = if_else(
      !is.na(days_to_censor_raw) & days_to_censor_raw >= 0,
      days_to_censor_raw, Inf
    ),
    first_time_days = pmin(
      time_to_eskd_days, time_to_death_days, time_to_censor_days
    ),
    exclusion_reason = case_when(
      is.na(dateOfDiagnosis) ~ "Missing diagnosis date",
      eskd_flag & is.na(dateOfESKD) ~ "ESKD status without ESKD date",
      death_flag & is.na(dateOfDeath) ~ "Death status without death date",
      !is.na(days_to_death_raw) & days_to_death_raw < 0 ~ "Death date before diagnosis",
      !is.na(days_to_censor_raw) & days_to_censor_raw < 0 ~ "Last contact before diagnosis",
      prediagnosis_eskd & prediagnosis_eskd_policy == "exclude" ~ "ESKD before diagnosis",
      !is.finite(first_time_days) ~ "No valid outcome follow-up",
      TRUE ~ NA_character_
    ),
    status_raw = case_when(
      !is.finite(first_time_days) ~ NA_integer_,
      time_to_eskd_days <= time_to_death_days &
        time_to_eskd_days <= time_to_censor_days ~ 1L,
      time_to_death_days < time_to_eskd_days &
        time_to_death_days <= time_to_censor_days ~ 2L,
      TRUE ~ 0L
    ),
    status = if_else(
      first_time_days > cif_followup_years * 365.25, 0L, status_raw
    ),
    time = pmin(first_time_days / 365.25, cif_followup_years)
  )

full_cif_quality_audit <- full_cif_patient %>%
  filter(
    is.na(dateOfDiagnosis) | prediagnosis_eskd | eskd_after_last_contact |
      (!is.na(days_to_death_raw) & days_to_death_raw < 0) |
      (!is.na(days_to_censor_raw) & days_to_censor_raw < 0) |
      !is.na(exclusion_reason)
  ) %>%
  select(
    RKD.ID, trajectory_eligibility,
    dateOfDiagnosis, dateOfESKD, dateOfDeath, lastRecordedContact,
    days_to_eskd_raw, days_to_death_raw, days_to_censor_raw,
    prediagnosis_eskd, eskd_after_last_contact, exclusion_reason
  )

full_cif_disposition <- full_cif_patient %>%
  mutate(CIF_disposition = coalesce(exclusion_reason, "Included in CIF")) %>%
  count(trajectory_eligibility, CIF_disposition, name = "n")

write.csv(
  full_cif_quality_audit,
  file.path(cif_output_dir, "CIF_patient_level_data_quality_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  full_cif_disposition,
  file.path(cif_output_dir, "CIF_analysis_disposition.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

full_cif_analysis <- full_cif_patient %>%
  filter(is.na(exclusion_reason), is.finite(time), !is.na(status))
excluded_cif_analysis <- full_cif_analysis %>%
  filter(trajectory_eligibility == "RETRACE excluded")

fit_single_cif <- function(data) {
  if (!any(data$status == 1L)) stop("A plotted cohort has no ESKD events.")
  cmprsk::cuminc(data$time, data$status, cencode = 0)
}

full_cif_curves <- bind_rows(
  tidy_cuminc(fit_single_cif(full_cif_analysis), 1L, "Full cohort"),
  tidy_cuminc(fit_single_cif(excluded_cif_analysis), 1L, "RETRACE excluded")
)

full_cif_summary <- extract_cif_times(full_cif_curves) %>%
  mutate(across(c(estimate, lower, upper), ~ 100 * .x, .names = "{.col}_percent"))

full_cif_outcome_counts <- bind_rows(
  full_cif_analysis %>%
    summarise(
      population = "Full cohort", N = n(),
      ESKD_events = sum(status == 1L),
      competing_deaths = sum(status == 2L),
      censored = sum(status == 0L)
    ),
  full_cif_analysis %>%
    group_by(trajectory_eligibility) %>%
    summarise(
      population = as.character(first(trajectory_eligibility)), N = n(),
      ESKD_events = sum(status == 1L),
      competing_deaths = sum(status == 2L),
      censored = sum(status == 0L),
      .groups = "drop"
    ) %>%
    select(-trajectory_eligibility)
)

write.csv(
  full_cif_curves %>% select(population, time, estimate, lower, upper),
  file.path(cif_output_dir, "CIF_curve_coordinates.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  full_cif_summary,
  file.path(cif_output_dir, "CIF_summary_1_3_5_years.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  full_cif_outcome_counts,
  file.path(cif_output_dir, "CIF_outcome_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

full_source_n <- nrow(full_cif_patient)
full_analysis_n <- nrow(full_cif_analysis)
excluded_source_n <- sum(full_cif_patient$trajectory_eligibility == "RETRACE excluded")
excluded_analysis_n <- nrow(excluded_cif_analysis)

cif_y_upper <- min(
  1,
  max(0.25, ceiling(max(full_cif_curves$upper, na.rm = TRUE) / 0.05) * 0.05)
)

make_cif_panel <- function(population_value, colour, title) {
  ggplot(
    full_cif_curves %>% filter(.data$population == .env$population_value),
    aes(time, estimate)
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = colour, alpha = 0.14, colour = NA
    ) +
    geom_step(linewidth = 1.05, colour = colour) +
    coord_cartesian(
      xlim = c(0, cif_followup_years),
      ylim = c(0, cif_y_upper),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = 0:cif_followup_years) +
    scale_y_continuous(
      breaks = seq(0, cif_y_upper, by = 0.05),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      title = title,
      x = "Time in years since diagnosis",
      y = "Cumulative incidence of ESKD"
    ) +
    theme_manuscript(base_size = 12, legend = "none") +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))
}

p_full_cohort_cif <- make_cif_panel(
  "Full cohort", "grey20",
  sprintf(
    "A. Full cohort (n = %d/%d analyzable)",
    full_analysis_n, full_source_n
  )
)
p_excluded_cohort_cif <- make_cif_panel(
  "RETRACE excluded", "#D55E00",
  sprintf(
    "B. RETRACE-excluded cohort (n = %d/%d analyzable)",
    excluded_analysis_n, excluded_source_n
  )
)

combined_full_excluded_cif <-
  (p_full_cohort_cif | p_excluded_cohort_cif) +
  plot_annotation(
    title = "Cumulative incidence of ESKD",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

print(combined_full_excluded_cif)
ggsave(
  file.path(cif_output_dir, "Supplementary_CIF_full_and_excluded.png"),
  combined_full_excluded_cif, width = 12, height = 5.5, dpi = 600, bg = "white"
)
ggsave(
  file.path(cif_output_dir, "Supplementary_CIF_full_and_excluded.pdf"),
  combined_full_excluded_cif, width = 12, height = 5.5, bg = "white"
)

stopifnot(
  full_source_n == 742L,
  full_analysis_n == 739L,
  excluded_source_n == 508L,
  excluded_analysis_n == 505L,
  nrow(prediagnosis_eskd_records) == 11L,
  nrow(full_cif_field_conflicts) == 0L
)

print(full_cif_outcome_counts)
print(full_cif_summary)
print(full_cif_disposition)
