## ---- ect-models ----
fit_lme <- function(formula, data = Ect_data, na.action = na.fail) {
  model <- nlme::lme(
    formula, random = ~ 1 | RKD.ID, data = data,
    method = "REML", na.action = na.action
  )
  model
}

lm1 <- fit_lme(eGFR_delta ~ Interval.from.diagnosis..months.)
lm2 <- fit_lme(eGFR_delta ~ lspline(Interval.from.diagnosis..months., knot = c(2, 4)))
lm3 <- fit_lme(eGFR_delta ~ gender * lspline(Interval.from.diagnosis..months., knot = c(2, 4)))

Ect_data_2 <- Ect_data %>% filter(ancaSpec != "Other") %>% droplevels()
lm5 <- fit_lme(
  eGFR_delta ~ ancaSpec * lspline(Interval.from.diagnosis..months., knot = c(2, 4)),
  Ect_data_2
)
lm6 <- fit_lme(
  eGFR_delta ~ Diagnosis * lspline(Interval.from.diagnosis..months., knot = c(2, 4)),
  na.action = na.omit
)
lm7 <- fit_lme(
  eGFR_delta ~ inductionFinal * lspline(Interval.from.diagnosis..months., knot = c(2, 4))
)
lm8 <- fit_lme(
  eGFR_delta ~ ancaSpec + gender +
    inductionFinal * lspline(Interval.from.diagnosis..months., knot = c(2, 4)),
  na.action = na.omit
)

model_exports <- list(lm3 = lm3, lm5 = lm5, lm6 = lm6, lm7 = lm7, lm8 = lm8)
invisible(Map(
  function(model, model_name) {
    model_table <- summary(model)$tTable %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term")
    save_table(model_table, paste0(model_name, "_tTable_egfr.csv"))
  },
  model_exports,
  names(model_exports)
))

print(summary(lm2)$tTable)
print(summary(lm7)$tTable)
