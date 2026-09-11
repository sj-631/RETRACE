## ---- load-filter-clean ----
raw_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = TRUE)

disease_labels <- c(
  "Granulomatosis with polyangiitis (Wegener) - Orpha:900",
  "Microscopic polyangiitis (including renal limited vasculitis) - ORPHA:727"
)

exRKD <- raw_data %>%
  filter(
    diagnosis_confidence == "Definite",
    small_vessel_vas_anca %in% disease_labels,
    Secondary.vasculitis != "Yes",
    Other != "Yes",
    is.na(Medium.vessel.vasculitis) | Medium.vessel.vasculitis == "",
    is.na(Large.vessel.vasculitis) | Large.vessel.vasculitis == "",
    is.na(Variable.vessel.vasculitis) | Variable.vessel.vasculitis == "",
    is.na(small_vessel_vas_immune) |
      small_vessel_vas_immune %in% c("", "Anti-GBM disease - ORPHA:375")
  ) %>%
  mutate(
    inductiontreatmentType_Rituximab = checked_label(
      inductiontreatmentType_Rituximab, "Rituximab"
    ),
    inductiontreatmentType_Cyclophosphamide = checked_label(
      inductiontreatmentType_Cyclophosphamide, "Cyclo.P"
    ),
    inductiontreatmentType_.Daily.Oral.Cyclophosphamide = checked_label(
      inductiontreatmentType_.Daily.Oral.Cyclophosphamide, "Cyclo.O"
    ),
    inductiontreatmentType_.Pulsed.IV.corticosteroids = checked_label(
      inductiontreatmentType_.Pulsed.IV.corticosteroids, "IV.Steroids"
    ),
    inductiontreatmentType_Oral.corticosteroids = checked_label(
      inductiontreatmentType_Oral.corticosteroids, "Oral.Steroids"
    ),
    cyclo = if_else(
      inductiontreatmentType_.Daily.Oral.Cyclophosphamide == "Cyclo.O" |
        inductiontreatmentType_Cyclophosphamide == "Cyclo.P",
      "Cyclo", " "
    ),
    corti = if_else(
      inductiontreatmentType_Oral.corticosteroids == "Oral.Steroids" |
        inductiontreatmentType_.Pulsed.IV.corticosteroids == "IV.Steroids",
      "Corti", " "
    ),
    corti_seperate = str_c(
      inductiontreatmentType_Oral.corticosteroids, ";",
      inductiontreatmentType_.Pulsed.IV.corticosteroids
    ),
    induction_trtmnt = str_c(inductiontreatmentType_Rituximab, " ", cyclo),
    inductionFinal = trimws(induction_trtmnt),
    inductionFinal = if_else(
      inductionFinal %in% c("Cyclo", "Rituximab"), inductionFinal, "Other"
    ),
    Berden.score.on.renal.biopsy = standardize_berden(
      Berden.score.on.renal.biopsy
    ),
    ancaSpec = if_else(ancaSpec %in% c("MPO", "PR3"), ancaSpec, "Other")
  )

organ_cols <- c(
  "affectedOrgan_Renal", "affectedOrgan_Cutaneous", "affectedOrgan_General",
  "affectedOrgan_Abdominal", "affectedOrgan_Mucous_membranes_eyes",
  "affectedOrgan_Cardiovascular", "affectedOrgan_Other", "affectedOrgan_ENT",
  "affectedOrgan_Chest", "affectedOrgan_Nervous_system"
)
core_cols <- c(
  "RKD.ID", "creatininetest", "egfr_calculated", "dateOfDiagnosis",
  "Interval.from.diagnosis..months.", "Date.Of.Visit", "Dialysis.dependent",
  "Berden.score.on.renal.biopsy", "dateOfDeath", "dateOfESKD", "ancaSpec",
  "gender", "DefiniteGPA", "DefiniteMPA", "ESKDStandard",
  "lastRecordedContact", "death", "ageAtDiagnosis", "inductionFinal",
  "inductiontreatmentType_.Pulsed.IV.corticosteroids",
  "inductiontreatmentType_Oral.corticosteroids",
  "inductiontreatmentType_Plasma.exchange",
  "inductiontreatmentType_Methotrexate",
  "inductiontreatmentType_Rituximab", "cyclo"
)
needed_cols <- unique(c(core_cols, organ_cols))
assert_columns(exRKD, needed_cols)

renal_raw <- exRKD %>%
  select(all_of(needed_cols)) %>%
  mutate(
    across(c(creatininetest, egfr_calculated, Interval.from.diagnosis..months.), num_clean),
    across(c(Date.Of.Visit, dateOfDiagnosis, dateOfDeath, dateOfESKD,
             lastRecordedContact), date_clean),
    Dialysis.dependent = trimws(as.character(Dialysis.dependent)),
    Dialysis.dependent = if_else(
      is.na(Dialysis.dependent) | Dialysis.dependent == "", "No", Dialysis.dependent
    ),
    Diagnosis = case_when(
      DefiniteMPA == 1 ~ "MPA",
      DefiniteGPA == 1 ~ "GPA",
      TRUE ~ NA_character_
    )
  )

# Organ fields are summarized separately so that patient-level collapsing does
# not alter the visit-level renal filter used for the ECT cohort.
organ_summary <- renal_raw %>%
  group_by(RKD.ID) %>%
  summarise(across(all_of(organ_cols), collapse_checked), .groups = "drop") %>%
  summarise(across(all_of(organ_cols), ~ sum(.x == "Checked", na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Organ", values_to = "n_checked")

cat("Original IDs:", n_distinct(raw_data$RKD.ID), "\n")
cat("Filtered AAV IDs:", n_distinct(renal_raw$RKD.ID), "\n")
print(organ_summary)
