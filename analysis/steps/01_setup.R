## ---- setup ----
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, error = FALSE)
options(rgl.useNULL = TRUE)

# Keep labels and exported tables reproducible when command-line R inherits a
# non-UTF-8 locale. Try common locale names without changing numeric/date rules.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (locale_name in c("en_US.UTF-8", "C.UTF-8", "C.utf8")) {
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", locale_name), silent = TRUE))
    if (isTRUE(l10n_info()[["UTF-8"]])) break
  }
}
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  warning("No UTF-8 locale is available; non-ASCII figure labels may be escaped.")
}

# Paths supplied by run_analysis.R or the YAML parameters.
data_path <- params$data_path
output_dir <- params$output_dir
if (!file.exists(data_path)) stop("Input data file not found: ", data_path)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
  stop("Use a new or empty output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "nlme", "lspline",
  "stringr", "kml", "cmprsk", "scales", "patchwork"
)
namespace_available <- function(pkg) {
  check <- function() requireNamespace(pkg, quietly = TRUE)
  if (identical(pkg, "kml")) suppressWarnings(check()) else check()
}
missing_packages <- required_packages[!vapply(
  required_packages, namespace_available, logical(1)
)]
if (length(missing_packages)) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, function(pkg) {
  load_package <- function() library(pkg, character.only = TRUE)
  if (identical(pkg, "kml")) suppressWarnings(load_package()) else load_package()
}))

group_cols <- c("Recovered group" = "#0072B2", "Static group" = "#D55E00")
induction_cols <- c(
  "Cyclophosphamide" = "#0072B2",
  "Rituximab" = "#D55E00",
  "Other" = "#009E73"
)
berden_levels <- c("Focal", "Crescentic", "Mixed", "Sclerotic", "Missing")
berden_cols <- c(
  "Focal" = "#1a9850", "Crescentic" = "#66bd63",
  "Mixed" = "#fdae61", "Sclerotic" = "#d73027", "Missing" = "grey70"
)

sys.source(file.path(project_root, "R", "analysis_helpers.R"), envir = environment())
