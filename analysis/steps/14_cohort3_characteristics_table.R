## ---- cohort3-characteristics-table ----
cohort3_ids <- slope_patient_eligibility %>%
  filter(model_eligible) %>%
  select(
    RKD.ID, cluster_label, T1_date, T1_month, T1_eGFR,
    n_slope_visits, maximum_followup_from_T1_years
  )

# Reuse the Cohort 2 baseline rule and its dialysis substitutions.
cohort3_baseline <- bind_rows(baseline_post, baseline_pre) %>%
  distinct(RKD.ID, .keep_all = TRUE) %>%
  transmute(
    RKD.ID,
    baseline_visit_date = Date.Of.Visit,
    baseline_month = Interval.from.diagnosis..months.,
    baseline_creatinine = creatininetest,
    baseline_eGFR = egfr_calculated,
    baseline_dialysis_dependent = Dialysis.dependent == "Yes"
  )

any_flagged <- function(x, yes = "Checked") {
  x <- trimws(as.character(x))
  if (any(x %in% yes, na.rm = TRUE)) TRUE else
    if (all(is.na(x) | x == "")) NA else FALSE
}

cohort3_features <- raw_data %>%
  semi_join(cohort3_ids, by = "RKD.ID") %>%
  group_by(RKD.ID) %>%
  summarise(
    age_at_diagnosis = first_finite(ageAtDiagnosis),
    sex = first_nonblank(gender),
    anca_specificity = {
      x <- first_nonblank(ancaSpec)
      if (is.na(x)) NA_character_ else if (x %in% c("MPO", "PR3")) x else "Other"
    },
    diagnosis = case_when(
      any(num_clean(DefiniteMPA) == 1, na.rm = TRUE) ~ "MPA",
      any(num_clean(DefiniteGPA) == 1, na.rm = TRUE) ~ "GPA",
      TRUE ~ NA_character_
    ),
    berden_class = {
      x <- trimws(as.character(Berden.score.on.renal.biopsy))
      x <- x[x %in% c("Crescentic", "Focal", "Mixed", "Sclerotic")]
      if (length(x)) x[1] else "Berden class unavailable"
    },
    rituximab_used = any_flagged(inductiontreatmentType_Rituximab),
    cyclophosphamide_used = any(
      trimws(as.character(inductiontreatmentType_Cyclophosphamide)) == "Checked" |
        trimws(as.character(
          inductiontreatmentType_.Daily.Oral.Cyclophosphamide
        )) == "Checked",
      na.rm = TRUE
    ),
    oral_corticosteroids = any_flagged(inductiontreatmentType_Oral.corticosteroids),
    pulsed_iv_corticosteroids = any_flagged(
      inductiontreatmentType_.Pulsed.IV.corticosteroids
    ),
    plasma_exchange = any_flagged(inductiontreatmentType_Plasma.exchange),
    across(all_of(organ_cols), any_flagged),
    .groups = "drop"
  ) %>%
  mutate(
    induction_treatment = case_when(
      cyclophosphamide_used & !rituximab_used ~ "Cyclophosphamide",
      rituximab_used & !cyclophosphamide_used ~ "Rituximab",
      TRUE ~ "Other"
    )
  )

cohort3 <- cohort3_ids %>%
  left_join(cohort3_features, by = "RKD.ID") %>%
  left_join(cohort3_baseline, by = "RKD.ID") %>%
  mutate(
    cluster_label = factor(
      cluster_label, levels = c("Recovered group", "Static group")
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    anca_specificity = factor(
      anca_specificity, levels = c("MPO", "PR3", "Other")
    ),
    diagnosis = factor(diagnosis, levels = c("MPA", "GPA")),
    berden_class = factor(
      berden_class,
      levels = c(
        "Crescentic", "Focal", "Mixed", "Sclerotic",
        "Berden class unavailable"
      )
    ),
    induction_treatment = factor(
      induction_treatment,
      levels = c("Cyclophosphamide", "Rituximab", "Other")
    )
  )

cohort3_groups <- c("Overall", "Recovered group", "Static group")
cohort3_n <- c(
  Overall = nrow(cohort3),
  `Recovered group` = sum(cohort3$cluster_label == "Recovered group"),
  `Static group` = sum(cohort3$cluster_label == "Static group")
)
cohort3_group_columns <- setNames(
  sprintf("%s (N = %d)", cohort3_groups, cohort3_n[cohort3_groups]),
  cohort3_groups
)

c3_subset <- function(group) {
  if (group == "Overall") cohort3 else filter(cohort3, cluster_label == group)
}
c3_median_iqr <- function(data, variable, digits = 1) {
  x <- data[[variable]]
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    q[2], q[1], q[3]
  )
}
c3_n_pct <- function(data, variable, level = TRUE, denominator = NULL) {
  x <- data[[variable]]
  if (is.null(denominator)) denominator <- sum(!is.na(x))
  n <- sum(x == level, na.rm = TRUE)
  if (!denominator) return("")
  sprintf("%d/%d (%.1f%%)", n, denominator, 100 * n / denominator)
}
c3_p <- function(p) if (is.na(p)) "" else format_p(p)
c3_values <- function(formatter) {
  formatter <- rlang::as_function(formatter)
  setNames(
    vapply(cohort3_groups, function(g) formatter(c3_subset(g)), character(1)),
    cohort3_group_columns
  )
}
c3_row <- function(section, variable, values, p = "") {
  bind_cols(
    tibble(Section = section, Variable = variable),
    tibble::as_tibble_row(values),
    tibble(`P value` = p)
  )
}
c3_continuous <- function(section, label, variable, digits = 1) {
  c3_row(
    section, label,
    c3_values(~ c3_median_iqr(.x, variable, digits)),
    c3_p(wilcox_p(cohort3[[variable]], cohort3$cluster_label))
  )
}
c3_binary <- function(section, label, variable, level = TRUE) {
  c3_row(
    section, label,
    c3_values(~ c3_n_pct(.x, variable, level)),
    c3_p(fisher_p(cohort3[[variable]], cohort3$cluster_label))
  )
}
c3_factor_header <- function(section, label, variable, data = cohort3) {
  c3_row(
    section, label, setNames(rep("", 3), cohort3_group_columns),
    c3_p(fisher_p(data[[variable]], data$cluster_label))
  )
}
c3_factor_level <- function(section, label, variable, level) {
  c3_row(
    section, paste0("  ", label),
    c3_values(~ c3_n_pct(.x, variable, level))
  )
}
c3_berden <- function(data, level) {
  denominator <- if (level == "Berden class unavailable") {
    nrow(data)
  } else {
    sum(!is.na(data$berden_class) &
          data$berden_class != "Berden class unavailable")
  }
  c3_n_pct(data, "berden_class", level, denominator)
}

c3_berden_available <- cohort3 %>%
  filter(berden_class != "Berden class unavailable") %>%
  droplevels()
c3_has_berden <- cohort3$berden_class != "Berden class unavailable"

cohort3_table <- bind_rows(
  c3_continuous("Demographic characteristics", "Age at diagnosis, years", "age_at_diagnosis"),
  c3_binary("Demographic characteristics", "Female sex", "sex", "Female"),
  c3_factor_header("Disease characteristics", "Clinical diagnosis", "diagnosis"),
  c3_factor_level("Disease characteristics", "MPA", "diagnosis", "MPA"),
  c3_factor_level("Disease characteristics", "GPA", "diagnosis", "GPA"),
  c3_factor_header("Disease characteristics", "ANCA specificity", "anca_specificity"),
  c3_factor_level("Disease characteristics", "MPO-ANCA", "anca_specificity", "MPO"),
  c3_factor_level("Disease characteristics", "PR3-ANCA", "anca_specificity", "PR3"),
  c3_factor_level("Disease characteristics", "Other ANCA result*", "anca_specificity", "Other"),
  c3_binary("Organ/system involvement (ever recorded)", "Constitutional", "affectedOrgan_General"),
  c3_binary("Organ/system involvement (ever recorded)", "Cutaneous", "affectedOrgan_Cutaneous"),
  c3_binary("Organ/system involvement (ever recorded)", "Mucous membrane/ocular", "affectedOrgan_Mucous_membranes_eyes"),
  c3_binary("Organ/system involvement (ever recorded)", "Ear, nose, and throat", "affectedOrgan_ENT"),
  c3_binary("Organ/system involvement (ever recorded)", "Pulmonary", "affectedOrgan_Chest"),
  c3_binary("Organ/system involvement (ever recorded)", "Cardiovascular", "affectedOrgan_Cardiovascular"),
  c3_binary("Organ/system involvement (ever recorded)", "Gastrointestinal", "affectedOrgan_Abdominal"),
  c3_binary("Organ/system involvement (ever recorded)", "Neurological", "affectedOrgan_Nervous_system"),
  c3_binary("Organ/system involvement (ever recorded)", "Renal", "affectedOrgan_Renal"),
  c3_binary("Organ/system involvement (ever recorded)", "Other recorded system involvement**", "affectedOrgan_Other"),
  c3_continuous("Kidney characteristics", "Baseline serum creatinine, µmol/L", "baseline_creatinine", 0),
  c3_continuous("Kidney characteristics", "Baseline eGFR, mL/min/1.73 m²", "baseline_eGFR", 0),
  c3_binary("Kidney characteristics", "Dialysis-dependent at baseline", "baseline_dialysis_dependent"),
  c3_factor_header(
    "Kidney histology", "Berden histological class", "berden_class",
    c3_berden_available
  ),
  c3_row("Kidney histology", "  Crescentic", c3_values(~ c3_berden(.x, "Crescentic"))),
  c3_row("Kidney histology", "  Focal", c3_values(~ c3_berden(.x, "Focal"))),
  c3_row("Kidney histology", "  Mixed", c3_values(~ c3_berden(.x, "Mixed"))),
  c3_row("Kidney histology", "  Sclerotic", c3_values(~ c3_berden(.x, "Sclerotic"))),
  c3_row(
    "Kidney histology", "Berden class unavailable",
    c3_values(~ c3_berden(.x, "Berden class unavailable")),
    c3_p(fisher_p(c3_has_berden, cohort3$cluster_label))
  ),
  c3_factor_header("Induction treatment", "Recorded induction regimen", "induction_treatment"),
  c3_factor_level("Induction treatment", "Cyclophosphamide only", "induction_treatment", "Cyclophosphamide"),
  c3_factor_level("Induction treatment", "Rituximab only", "induction_treatment", "Rituximab"),
  c3_factor_level("Induction treatment", "Combined CYC + RTX or other induction***", "induction_treatment", "Other"),
  c3_binary("Induction treatment", "Oral corticosteroids", "oral_corticosteroids"),
  c3_binary("Induction treatment", "Pulsed intravenous corticosteroids", "pulsed_iv_corticosteroids"),
  c3_binary("Induction treatment", "Plasma exchange", "plasma_exchange"),
  c3_continuous("Post-peak analysis", "T1, months after diagnosis", "T1_month"),
  c3_continuous("Post-peak analysis", "eGFR at T1, mL/min/1.73 m²", "T1_eGFR", 0),
  c3_continuous("Post-peak analysis", "Included eGFR measurement dates from T1, n", "n_slope_visits", 0),
  c3_continuous("Post-peak analysis", "Span from T1 to last included eGFR measurement, years", "maximum_followup_from_T1_years")
)

cohort3_notes <- c(
  "Characteristics of Cohort 3 (the post-peak eGFR slope analysis subset of Cohort 2), overall and by early eGFR recovery group.",
  "Values are median (interquartile range) or n/N (%), where N is the number with available data unless otherwise stated.",
  "For dialysis-dependent visits, eGFR was assigned 3 mL/min/1.73 m² before baseline and T1 selection; serum creatinine was assigned 500 µmol/L before baseline selection.",
  "The four Berden-class percentages use patients with an available class as the denominator (overall N=110; Recovered N=41; Static N=69); the unavailable row uses the full group denominator.",
  "* Other ANCA result comprises two ELISA-negative and two dual MPO/PR3 results.",
  "** Other recorded system involvement comprises peripheral vascular disease with gangrenous toes/amputation, a constitutional entry, and a mediastinal mass.",
  "*** This category comprises 13 patients who received both cyclophosphamide and rituximab and four who received another induction agent (three azathioprine and one mycophenolate mofetil).",
  "T1 is the earliest occurrence of the highest available eGFR between 4 and 12 months after diagnosis, after endpoint/date filtering and same-day record collapse.",
  "Slope-model eligibility required T1 and at least one later eGFR measurement on a distinct date; measurements from T1 through 60 months after diagnosis were included.",
  "P values compare the Recovered and Static groups using Wilcoxon rank-sum or Fisher's exact tests; they are descriptive, unadjusted, and were not corrected for multiple comparisons."
)

save_table(cohort3, "post_peak_egfr_slope/cohort3_post_peak_patient_level.csv")
save_table(
  cohort3_table,
  "post_peak_egfr_slope/cohort3_post_peak_characteristics_table.csv"
)
save_table(
  cohort3_table %>% filter(Section != "Post-peak analysis"),
  "post_peak_egfr_slope/cohort3_post_peak_baseline_characteristics.csv"
)
save_table(
  cohort3_table %>% filter(Section == "Post-peak analysis"),
  "post_peak_egfr_slope/cohort3_post_peak_analysis_characteristics.csv"
)
writeLines(
  cohort3_notes,
  file.path(slope_output_dir, "cohort3_post_peak_table_notes.txt"),
  useBytes = TRUE
)

stopifnot(
  nrow(cohort3) == 148L,
  sum(cohort3$cluster_label == "Recovered group") == 58L,
  sum(cohort3$cluster_label == "Static group") == 90L,
  nrow(cohort3_table) == 39L,
  !any(is.na(cohort3$baseline_eGFR))
)
print(cohort3_table, n = Inf, width = Inf)
