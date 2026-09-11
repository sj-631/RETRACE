## ---- fine-gray-models ----
RETRACE_cohort <- RETRACE_base %>% filter(RKD.ID %in% cluster_df$RKD.ID)

RESULT_crr <- RETRACE_cohort %>%
  arrange(RKD.ID, Interval.from.diagnosis..months.) %>%
  distinct(RKD.ID, .keep_all = TRUE) %>%
  left_join(cluster_df, by = "RKD.ID") %>%
  mutate(
    daysfrmDeath = as.numeric(dateOfDeath - dateOfDiagnosis) / 365.25,
    dateOfCensored = dateOfCensored_days / 365.25,
    daysfrmESKD = daysfrmESKD_days / 365.25,
    time_to_eskd = if_else(!is.na(daysfrmESKD), daysfrmESKD, Inf),
    time_to_death = if_else(!is.na(daysfrmDeath), daysfrmDeath, Inf),
    time_to_censor = if_else(!is.na(dateOfCensored), dateOfCensored, Inf),
    time = pmin(time_to_eskd, time_to_death, time_to_censor, na.rm = TRUE),
    status = case_when(
      ESKDStandard == "Yes" &
        time_to_eskd <= time_to_death & time_to_eskd <= time_to_censor ~ 1L,
      death == "Dead" &
        time_to_death < time_to_eskd & time_to_death <= time_to_censor ~ 2L,
      TRUE ~ 0L
    ),
    status = if_else(time > 5, 0L, status),
    time = pmin(time, 5),
    cluster_label = factor(
      cluster_label,
      levels = c("Recovered group", "Static group")
    ),
    eGFR1_10 = eGFR1 / 10,
    creatininetest1_10 = creatininetest1 / 10
  )

X <- model.matrix(~ factor(cluster), RESULT_crr)[, -1, drop = FALSE]
fit_crr <- cmprsk::crr(
  ftime = RESULT_crr$time, fstatus = RESULT_crr$status,
  cov1 = X, failcode = 1, cencode = 0, maxiter = 100
)

fg_data <- RESULT_crr %>%
  select(
    time, status, cluster, gender, ageAtDiagnosis, ancaSpec,
    inductionFinal, eGFR1_10, Berden.score.on.renal.biopsy
  ) %>%
  mutate(
    gender = factor(gender),
    ancaSpec = factor(ancaSpec, levels = c("MPO", "PR3", "Other")),
    inductionFinal = factor(inductionFinal, levels = c("Cyclo", "Rituximab", "Other")),
    Berden.score.on.renal.biopsy = factor(Berden.score.on.renal.biopsy)
  ) %>%
  filter(complete.cases(.)) %>%
  droplevels()

X_full <- model.matrix(
  ~ factor(cluster) + gender + ageAtDiagnosis + ancaSpec +
    inductionFinal + eGFR1_10 + Berden.score.on.renal.biopsy,
  fg_data
)[, -1, drop = FALSE]

fit_crr_full <- cmprsk::crr(
  ftime = fg_data$time, fstatus = fg_data$status,
  cov1 = X_full, failcode = 1, cencode = 0, maxiter = 100
)

capture.output(
  summary(fit_crr_full),
  file = file.path(output_dir, "fit_crr_full_summary.txt")
)

outcome_check <- with(
  RESULT_crr,
  table(cluster_label, factor(status, levels = 0:2))
)
print(outcome_check)
print(summary(fit_crr))
print(summary(fit_crr_full))

stopifnot(
  nrow(RESULT_crr) == 234L,
  !anyDuplicated(RESULT_crr$RKD.ID),
  all(outcome_check == matrix(c(89L, 114L, 1L, 15L, 7L, 8L), nrow = 2L))
)
