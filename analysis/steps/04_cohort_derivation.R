## ---- cohort-derivation ----
Full_cohort <- renal_raw

# Full-cohort baseline rule: earliest value in [-0.5, 1] months; if absent,
# use the earliest dated value anywhere.
baseline_values <- Full_cohort %>%
  group_by(RKD.ID) %>%
  group_modify(~ {
    time <- .x$Interval.from.diagnosis..months.
    preferred <- !is.na(time) & time >= -0.5 & time <= 1
    if (any(preferred)) {
      .x[preferred, ] %>% arrange(Interval.from.diagnosis..months.) %>% slice(1)
    } else if (any(!is.na(time))) {
      .x %>%
        filter(!is.na(Interval.from.diagnosis..months.)) %>%
        arrange(Interval.from.diagnosis..months.) %>%
        slice(1)
    } else {
      slice(.x, 1)
    }
  }) %>%
  ungroup() %>%
  select(
    RKD.ID,
    creatininetest1 = creatininetest,
    egfr_calculated1 = egfr_calculated
  )

Full_cohort <- Full_cohort %>% left_join(baseline_values, by = "RKD.ID")

# ECT source: renal visits in the first six months. Dialysis-dependent visits
# are assigned eGFR=3 and creatinine=500 before baseline and change calculation.
crtn <- renal_raw %>%
  filter(
    affectedOrgan_Renal == "Checked",
    !is.na(Date.Of.Visit),
    !is.na(Interval.from.diagnosis..months.),
    Interval.from.diagnosis..months. > -0.5,
    Interval.from.diagnosis..months. <= 6
  ) %>%
  mutate(
    egfr_calculated = if_else(Dialysis.dependent == "Yes", 3, egfr_calculated),
    creatininetest = if_else(Dialysis.dependent == "Yes", 500, creatininetest)
  ) %>%
  filter(is.finite(egfr_calculated)) %>%
  arrange(RKD.ID, Date.Of.Visit, Interval.from.diagnosis..months.) %>%
  distinct(RKD.ID, Date.Of.Visit, .keep_all = TRUE) %>%
  arrange(RKD.ID, Interval.from.diagnosis..months.)

# ECT baseline rule: earliest value at 0-1 month; otherwise the closest value
# before diagnosis in (-0.5, 0).
baseline_post <- crtn %>%
  filter(between(Interval.from.diagnosis..months., 0, 1)) %>%
  arrange(RKD.ID, Interval.from.diagnosis..months.) %>%
  distinct(RKD.ID, .keep_all = TRUE)

baseline_pre <- crtn %>%
  filter(
    !RKD.ID %in% baseline_post$RKD.ID,
    Interval.from.diagnosis..months. > -0.5,
    Interval.from.diagnosis..months. < 0
  ) %>%
  arrange(RKD.ID, desc(Interval.from.diagnosis..months.)) %>%
  distinct(RKD.ID, .keep_all = TRUE)

crtn_start <- bind_rows(baseline_post, baseline_pre) %>%
  distinct(RKD.ID, .keep_all = TRUE) %>%
  transmute(
    RKD.ID,
    creatininetest1 = creatininetest,
    eGFR1 = egfr_calculated
  )

baseline_eGFR_status <- Full_cohort %>%
  distinct(RKD.ID) %>%
  mutate(has_baseline_eGFR = RKD.ID %in% crtn_start$RKD.ID)

crtn <- crtn %>%
  inner_join(crtn_start, by = "RKD.ID") %>%
  mutate(
    creatininetest_delta = creatininetest - creatininetest1,
    creatininetest_percent = if_else(
      !is.na(creatininetest1) & creatininetest1 > 0,
      100 * creatininetest_delta / creatininetest1,
      NA_real_
    ),
    eGFR_delta = egfr_calculated - eGFR1,
    eGFR_percent = if_else(
      !is.na(eGFR1) & eGFR1 > 0, 100 * eGFR_delta / eGFR1, NA_real_
    )
  ) %>%
  filter(is.finite(eGFR_delta)) %>%
  group_by(RKD.ID) %>%
  mutate(
    n_of_eGFR = n(),
    has_post_diagnosis = any(Interval.from.diagnosis..months. >= 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(n_of_eGFR > 1, has_post_diagnosis)

ECT_cohort <- crtn
Ect_data <- ECT_cohort %>%
  mutate(
    gender = factor(gender),
    ancaSpec = factor(ancaSpec, levels = c("MPO", "PR3", "Other")),
    inductionFinal = factor(inductionFinal, levels = c("Cyclo", "Rituximab", "Other")),
    Diagnosis = factor(Diagnosis, levels = c("MPA", "GPA"))
  )

cat("Full cohort IDs:", n_distinct(Full_cohort$RKD.ID), "\n")
cat("Patients without baseline eGFR:", sum(!baseline_eGFR_status$has_baseline_eGFR), "\n")
cat("ECT cohort IDs/rows:", n_distinct(Ect_data$RKD.ID), "/", nrow(Ect_data), "\n")

stopifnot(
  n_distinct(Full_cohort$RKD.ID) == 742L,
  n_distinct(Ect_data$RKD.ID) == 311L
)
