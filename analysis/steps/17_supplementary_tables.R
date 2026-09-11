## ---- supplementary-tables ----
summarise_characteristics <- function(data, groups = character()) {
  if (length(groups)) data <- data %>% group_by(across(all_of(groups)))
  data %>%
    summarise(
      N = n(),
      mean_age = round(median(ageAtDiagnosis, na.rm = TRUE), 1),
      age_q1 = round(quantile(ageAtDiagnosis, 0.25, na.rm = TRUE), 1),
      age_q3 = round(quantile(ageAtDiagnosis, 0.75, na.rm = TRUE), 1),
      female_n = sum(gender == "Female", na.rm = TRUE),
      mean_creatinine = round(median(creatininetest1, na.rm = TRUE), 0),
      creatinine_q1 = round(quantile(creatininetest1, 0.25, na.rm = TRUE), 0),
      creatinine_q3 = round(quantile(creatininetest1, 0.75, na.rm = TRUE), 0),
      mean_egfr = round(median(eGFR1, na.rm = TRUE), 0),
      egfr_q1 = round(quantile(eGFR1, 0.25, na.rm = TRUE), 0),
      egfr_q3 = round(quantile(eGFR1, 0.75, na.rm = TRUE), 0),
      MPA_n = sum(DefiniteMPA == 1, na.rm = TRUE),
      GPA_n = sum(DefiniteGPA == 1, na.rm = TRUE),
      Berden_Crescentic = sum(Berden.score.on.renal.biopsy == "Crescentic", na.rm = TRUE),
      Berden_Focal = sum(Berden.score.on.renal.biopsy == "Focal", na.rm = TRUE),
      Berden_Mixed = sum(Berden.score.on.renal.biopsy == "Mixed", na.rm = TRUE),
      Berden_Sclerotic = sum(Berden.score.on.renal.biopsy == "Sclerotic", na.rm = TRUE),
      Berden_Missing = sum(Berden.score.on.renal.biopsy == "Missing", na.rm = TRUE),
      MPO_n = sum(ancaSpec == "MPO", na.rm = TRUE),
      PR3_n = sum(ancaSpec == "PR3", na.rm = TRUE),
      Other_n = sum(ancaSpec == "Other", na.rm = TRUE),
      Ind_Cyclo = sum(inductionFinal == "Cyclo", na.rm = TRUE),
      Ind_Rituximab = sum(inductionFinal == "Rituximab", na.rm = TRUE),
      Ind_Others = sum(inductionFinal == "Other", na.rm = TRUE),
      Dialysis_dependent_n = sum(Dialysis.dependent == "Yes", na.rm = TRUE),
      Corti_Oral = sum(
        inductiontreatmentType_Oral.corticosteroids == "Oral.Steroids",
        na.rm = TRUE
      ),
      Corti_Intra = sum(
        inductiontreatmentType_.Pulsed.IV.corticosteroids == "IV.Steroids",
        na.rm = TRUE
      ),
      Plasma = sum(inductiontreatmentType_Plasma.exchange == "Checked", na.rm = TRUE),
      .groups = "drop"
    )
}

Full_cohort_unique <- Full_cohort %>%
  distinct(RKD.ID, .keep_all = TRUE) %>%
  rename(eGFR1 = egfr_calculated1)

Full_TABLE <- summarise_characteristics(Full_cohort_unique) %>%
  rename(Ind_other = Ind_Others)

Ect_data_unique <- Ect_data %>%
  distinct(RKD.ID, .keep_all = TRUE) %>%
  mutate(inductionFinal = as.character(inductionFinal))

ECT_TABLE <- bind_rows(
  Ect_data_unique,
  Ect_data_unique %>% mutate(inductionFinal = "Overall")
) %>%
  summarise_characteristics("inductionFinal") %>%
  select(-Ind_Cyclo, -Ind_Rituximab, -Ind_Others) %>%
  arrange(inductionFinal == "Overall")

RETRACE_summary_source <- bind_rows(
  RESULT_crr %>% mutate(cluster_label = as.character(cluster_label)),
  RESULT_crr %>% mutate(cluster_label = "Overall")
)

RETRACE_TABLE <- summarise_characteristics(RETRACE_summary_source, "cluster_label") %>%
  left_join(
    RETRACE_summary_source %>%
      group_by(cluster_label) %>%
      summarise(
        ESKD_n = sum(status == 1L),
        Death_competing_n = sum(status == 2L),
        .groups = "drop"
      ),
    by = "cluster_label"
  ) %>%
  arrange(cluster_label == "Overall")

RETRACE_TABLE_WIDE <- RETRACE_TABLE %>%
  pivot_longer(-cluster_label, names_to = "Characteristic", values_to = "Value") %>%
  pivot_wider(names_from = cluster_label, values_from = Value) %>%
  relocate(Overall, .after = last_col())

# Table 3 P-value calculations.
test_data <- RESULT_crr %>%
  mutate(
    diagnosis_group = case_when(
      DefiniteMPA == 1 ~ "MPA",
      DefiniteGPA == 1 ~ "GPA",
      TRUE ~ NA_character_
    ),
    oral_steroid = inductiontreatmentType_Oral.corticosteroids == "Oral.Steroids",
    iv_steroid = inductiontreatmentType_.Pulsed.IV.corticosteroids == "IV.Steroids",
    plasma_exchange = inductiontreatmentType_Plasma.exchange == "Checked"
  )
g <- test_data$cluster_label

p_values <- tibble(
  Characteristic = c(
    "mean_age", "female_n", "mean_creatinine", "mean_egfr", "MPA_n",
    "Berden_Crescentic", "MPO_n", "Ind_Cyclo", "Dialysis_dependent_n",
    "Corti_Oral", "Corti_Intra", "Plasma", "ESKD_n"
  ),
  p_raw = c(
    wilcox_p(test_data$ageAtDiagnosis, g),
    fisher_p(test_data$gender, g),
    wilcox_p(test_data$creatininetest1, g),
    wilcox_p(test_data$eGFR1, g),
    fisher_p(test_data$diagnosis_group, g),
    fisher_p(test_data$Berden.score.on.renal.biopsy, g),
    fisher_p(test_data$ancaSpec, g),
    fisher_p(test_data$inductionFinal, g),
    fisher_p(test_data$Dialysis.dependent, g),
    fisher_p(test_data$oral_steroid, g),
    fisher_p(test_data$iv_steroid, g),
    fisher_p(test_data$plasma_exchange, g),
    fisher_p(test_data$status, g)
  )
) %>%
  mutate(`P value` = if_else(is.na(p_raw), "", vapply(p_raw, format_p, character(1))))

RETRACE_TABLE_WIDE <- RETRACE_TABLE_WIDE %>%
  left_join(p_values %>% select(Characteristic, `P value`), by = "Characteristic") %>%
  relocate(`P value`, .after = last_col())

save_table(Full_TABLE, "supplementary_table_1_full_cohort_egfr.csv")
save_table(ECT_TABLE, "supplementary_table_2_ECT_cohort_egfr.csv")
save_table(RETRACE_TABLE, "supplementary_table_3_RETRACE_cohort_egfr.csv")
save_table(
  RETRACE_TABLE_WIDE,
  "supplementary_table_3_RETRACE_cohort_egfr_wide_with_p_values.csv"
)

print(Full_TABLE)
print(ECT_TABLE)
print(RETRACE_TABLE_WIDE)
