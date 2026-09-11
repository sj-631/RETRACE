## ---- post-peak-slope-data ----
slope_output_dir <- file.path(output_dir, "post_peak_egfr_slope")
dir.create(slope_output_dir, recursive = TRUE, showWarnings = FALSE)

slope_ids <- cluster_df %>%
  transmute(RKD.ID, cluster_label = as.character(cluster_label)) %>%
  distinct()

patient_endpoints <- renal_raw %>%
  semi_join(slope_ids, by = "RKD.ID") %>%
  group_by(RKD.ID) %>%
  summarise(
    dateOfESKD = first_valid_date(dateOfESKD),
    dateOfDeath = first_valid_date(dateOfDeath),
    lastRecordedContact = first_valid_date(lastRecordedContact),
    .groups = "drop"
  )

transplant_dates <- if ("Date.of.transplant." %in% names(raw_data)) {
  raw_data %>%
    semi_join(slope_ids, by = "RKD.ID") %>%
    transmute(RKD.ID, dateOfTransplant = date_clean(Date.of.transplant.)) %>%
    group_by(RKD.ID) %>%
    summarise(
      dateOfTransplant = first_valid_date(sort(dateOfTransplant)),
      .groups = "drop"
    )
} else {
  slope_ids %>% transmute(RKD.ID, dateOfTransplant = as.Date(NA_character_))
}

slope_source <- renal_raw %>%
  select(
    RKD.ID, Date.Of.Visit, Interval.from.diagnosis..months.,
    egfr_calculated, Dialysis.dependent
  ) %>%
  inner_join(slope_ids, by = "RKD.ID") %>%
  left_join(patient_endpoints, by = "RKD.ID") %>%
  left_join(transplant_dates, by = "RKD.ID") %>%
  mutate(
    interval_months = num_clean(Interval.from.diagnosis..months.),
    dialysis_dependent = Dialysis.dependent == "Yes",
    eGFR = if_else(dialysis_dependent, 3, num_clean(egfr_calculated))
  ) %>%
  filter(
    !is.na(Date.Of.Visit), is.finite(interval_months), is.finite(eGFR),
    between(interval_months, 4, 60),
    is.na(dateOfTransplant) | Date.Of.Visit < dateOfTransplant,
    is.na(dateOfDeath) | Date.Of.Visit <= dateOfDeath,
    is.na(lastRecordedContact) | Date.Of.Visit <= lastRecordedContact
  )

# Collapse same-day records; any dialysis-dependent row fixes that date at 3.
same_day_audit <- slope_source %>%
  count(RKD.ID, Date.Of.Visit, name = "n_rows") %>%
  filter(n_rows > 1)

slope_visits <- slope_source %>%
  group_by(RKD.ID, Date.Of.Visit) %>%
  summarise(
    cluster_label = first(cluster_label),
    dateOfESKD = first_valid_date(dateOfESKD),
    interval_months = mean(interval_months),
    dialysis_dependent = any(dialysis_dependent),
    eGFR = if_else(dialysis_dependent, 3, mean(eGFR)),
    .groups = "drop"
  )

t1_candidates <- slope_visits %>% filter(between(interval_months, 4, 12))
t1_tie_audit <- t1_candidates %>%
  group_by(RKD.ID) %>%
  filter(eGFR == max(eGFR)) %>%
  summarise(n_tied_dates = n_distinct(Date.Of.Visit), .groups = "drop") %>%
  filter(n_tied_dates > 1)

t1_by_patient <- t1_candidates %>%
  arrange(RKD.ID, desc(eGFR), Date.Of.Visit) %>%
  group_by(RKD.ID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    RKD.ID, T1_date = Date.Of.Visit,
    T1_month = interval_months, T1_eGFR = eGFR
  )

t1_source_audit <- renal_raw %>%
  semi_join(slope_ids, by = "RKD.ID") %>%
  transmute(
    RKD.ID, Date.Of.Visit,
    interval_months = num_clean(Interval.from.diagnosis..months.),
    eGFR = if_else(Dialysis.dependent == "Yes", 3, num_clean(egfr_calculated))
  ) %>%
  group_by(RKD.ID) %>%
  summarise(
    any_record_4_12 = any(
      !is.na(interval_months) & between(interval_months, 4, 12)
    ),
    any_usable_eGFR_4_12 = any(
      !is.na(Date.Of.Visit) & !is.na(interval_months) &
        between(interval_months, 4, 12) & is.finite(eGFR)
    ),
    .groups = "drop"
  )

t1_disposition <- slope_ids %>%
  left_join(t1_source_audit, by = "RKD.ID") %>%
  mutate(
    across(c(any_record_4_12, any_usable_eGFR_4_12), ~ coalesce(.x, FALSE)),
    T1_status = case_when(
      RKD.ID %in% t1_by_patient$RKD.ID ~ "T1 available",
      !any_record_4_12 ~ "No record at 4-12 months",
      !any_usable_eGFR_4_12 ~ "Record present, but no usable eGFR",
      TRUE ~ "Excluded by endpoint/date rule"
    )
  ) %>%
  count(T1_status, name = "n")

slope_after_t1 <- slope_visits %>%
  inner_join(t1_by_patient, by = "RKD.ID") %>%
  filter(Date.Of.Visit >= T1_date, interval_months <= 60) %>%
  mutate(time_from_T1_years = as.numeric(Date.Of.Visit - T1_date) / 365.25)

slope_patient_eligibility <- slope_ids %>%
  left_join(t1_by_patient, by = "RKD.ID") %>%
  left_join(
    slope_after_t1 %>%
      group_by(RKD.ID) %>%
      summarise(
        n_slope_visits = n_distinct(Date.Of.Visit),
        has_post_T1_visit = any(time_from_T1_years > 0),
        maximum_followup_from_T1_years = max(time_from_T1_years),
        .groups = "drop"
      ),
    by = "RKD.ID"
  ) %>%
  mutate(
    has_T1 = !is.na(T1_date),
    n_slope_visits = coalesce(n_slope_visits, 0L),
    has_post_T1_visit = coalesce(has_post_T1_visit, FALSE),
    model_eligible = has_T1 & n_slope_visits >= 2 & has_post_T1_visit
  )

slope_flow <- bind_rows(
  slope_patient_eligibility %>%
    group_by(cluster_label) %>%
    summarise(
      RETRACE_N = n(), T1_available = sum(has_T1),
      model_N = sum(model_eligible), .groups = "drop"
    ),
  slope_patient_eligibility %>%
    summarise(
      cluster_label = "Overall", RETRACE_N = n(),
      T1_available = sum(has_T1), model_N = sum(model_eligible)
    )
)

slope_data <- slope_after_t1 %>%
  semi_join(
    slope_patient_eligibility %>% filter(model_eligible),
    by = "RKD.ID"
  ) %>%
  mutate(
    cluster_slope = factor(
      cluster_label, levels = c("Static group", "Recovered group")
    )
  ) %>%
  arrange(RKD.ID, time_from_T1_years)

dialysis3_audit <- bind_rows(
  slope_source %>% summarise(
    analysis_stage = "All eligible 4-60 month source rows",
    dialysis_rows = sum(dialysis_dependent),
    dialysis_patients = n_distinct(RKD.ID[dialysis_dependent])
  ),
  slope_data %>% summarise(
    analysis_stage = "Primary model data from T1 onward",
    dialysis_rows = sum(dialysis_dependent),
    dialysis_patients = n_distinct(RKD.ID[dialysis_dependent])
  )
)

stopifnot(
  nrow(slope_ids) == 234L,
  slope_flow$T1_available[slope_flow$cluster_label == "Overall"] == 166L,
  slope_flow$model_N[slope_flow$cluster_label == "Overall"] == 148L,
  sum(slope_patient_eligibility$model_eligible &
        slope_patient_eligibility$cluster_label == "Recovered group") == 58L,
  sum(slope_patient_eligibility$model_eligible &
        slope_patient_eligibility$cluster_label == "Static group") == 90L,
  nrow(slope_data) == 875L
)

print(slope_flow)
print(t1_disposition)
print(dialysis3_audit)
