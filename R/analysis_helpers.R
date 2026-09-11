num_clean <- function(x) suppressWarnings(as.numeric(as.character(x)))
date_clean <- function(x) suppressWarnings(as.Date(as.character(x)))

standardize_berden <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "N/A", "na", "n/a")] <- "Missing"
  x
}

checked_label <- function(x, label) {
  x <- as.character(x)
  x[x == "Checked"] <- label
  x[x == "Unchecked"] <- " "
  x
}

collapse_checked <- function(x) {
  x <- trimws(as.character(x))
  if (any(x == "Checked", na.rm = TRUE)) "Checked" else
    if (all(is.na(x) | x == "")) NA_character_ else "Unchecked"
}

first_nonblank <- function(x) {
  keep <- !is.na(x) & trimws(as.character(x)) != ""
  if (any(keep)) x[which(keep)[1]] else x[NA_integer_][1]
}
first_valid_date <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[1] else as.Date(NA_character_)
}

first_finite <- function(x) {
  x <- num_clean(x)
  x <- x[is.finite(x)]
  if (length(x)) x[1] else NA_real_
}

assert_columns <- function(data, columns, object_name = deparse(substitute(data))) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(object_name, " is missing: ", paste(missing, collapse = ", "))
  }
  invisible(data)
}

save_plot <- function(plot, filename, width, height, dpi = 600) {
  path <- file.path(output_dir, filename)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Figure was not saved correctly: ", path)
  }
  invisible(path)
}

save_table <- function(x, filename) {
  path <- file.path(output_dir, filename)
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

theme_manuscript <- function(base_size = 13, legend = "bottom") {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend
    )
}

format_p <- function(p, prefix = FALSE) {
  label <- ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  if (prefix) paste("P", ifelse(p < 0.001, label, paste("=", label))) else label
}

wilcox_p <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (nlevels(group) != 2 || length(unique(x)) < 2) return(NA_real_)
  wilcox.test(x ~ group, exact = FALSE)$p.value
}

fisher_p <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  tab <- table(group[keep], x[keep])
  if (nrow(tab) != 2 || ncol(tab) < 2) return(NA_real_)
  fisher.test(tab)$p.value
}

tidy_cuminc <- function(fit, event_code = 1L, population = NULL) {
  curve_names <- setdiff(names(fit), "Tests")
  curve_names <- curve_names[grepl(paste0("(^| )", event_code, "$"), curve_names)]
  if (!length(curve_names)) stop("No CIF curve found for event code ", event_code)

  bind_rows(lapply(curve_names, function(curve_name) {
    curve <- fit[[curve_name]]
    data.frame(
      name = curve_name,
      time = curve$time,
      estimate = curve$est,
      variance = curve$var
    )
  })) %>%
    mutate(
      population = population,
      group = sub(paste0(" ?", event_code, "$"), "", name),
      se = sqrt(pmax(variance, 0)),
      lower = pmax(estimate - 1.96 * se, 0),
      upper = pmin(estimate + 1.96 * se, 1)
    )
}

fine_gray_stats <- function(fit, term = 1L) {
  beta <- unname(fit$coef[term])
  se <- sqrt(fit$var[term, term])
  tibble(
    beta = beta,
    SHR = exp(beta),
    lower = exp(beta - 1.96 * se),
    upper = exp(beta + 1.96 * se),
    p = 2 * pnorm(-abs(beta / se))
  )
}

predict_fixed_lme <- function(model, newdata) {
  design <- model.matrix(delete.response(terms(model)), newdata)
  beta <- nlme::fixef(model)
  covariance <- vcov(model)
  if (!setequal(colnames(design), names(beta))) {
    stop("Prediction design matrix does not match the model coefficients.")
  }
  design <- design[, names(beta), drop = FALSE]
  estimate <- as.numeric(design %*% beta)
  se <- sqrt(rowSums((design %*% covariance) * design))
  critical <- qt(0.975, df = nobs(model) - length(beta))
  bind_cols(
    as_tibble(newdata),
    tibble(
      fit = estimate,
      lower = estimate - critical * se,
      upper = estimate + critical * se
    )
  )
}

extract_cif_times <- function(curve_data, times = c(1, 3, 5)) {
  bind_rows(lapply(unique(curve_data$population), function(population_value) {
    curve <- curve_data %>% filter(population == population_value) %>% arrange(time)
    bind_rows(lapply(times, function(target_time) {
      idx <- which(curve$time <= target_time)
      if (!length(idx)) {
        return(tibble(
          population = population_value, year = target_time,
          estimate = 0, lower = 0, upper = 0
        ))
      }
      row <- curve[max(idx), ]
      tibble(
        population = population_value, year = target_time,
        estimate = row$estimate, lower = row$lower, upper = row$upper
      )
    }))
  }))
}
