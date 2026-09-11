#!/usr/bin/env Rscript
# Run the RETRACE analysis.
main <- function() {
  # kml imports rgl; use its headless device for scripted/reproducible runs.
  options(rgl.useNULL = TRUE)

  # Set a UTF-8 character locale before parsing the modular source files.
  # This matters because knitr reads those files before their setup chunk runs.
  if (!isTRUE(l10n_info()[["UTF-8"]])) {
    for (locale_name in c("en_US.UTF-8", "C.UTF-8", "C.utf8")) {
      suppressWarnings(try(Sys.setlocale("LC_CTYPE", locale_name), silent = TRUE))
      if (isTRUE(l10n_info()[["UTF-8"]])) break
    }
  }
  if (!isTRUE(l10n_info()[["UTF-8"]])) {
    warning("No UTF-8 locale is available; non-ASCII labels may be escaped.")
  }

  args <- commandArgs(trailingOnly = TRUE)
  help <- function() cat(
    "Rscript run_analysis.R --check\n",
    "Rscript run_analysis.R --input /path/to/harmonized.csv --output outputs/run_01\n",
    "Relative input/output paths are relative to the invoking working directory.\n",
    "The full analysis requires the clinical cohort data and the listed R packages.\n",
    sep = "")
  if (!length(args) || identical(args, "--help")) return(help())
  check_only <- "--check" %in% args
  options <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--check") { i <- i + 1L; next }
    if (!key %in% c("--input", "--output") || i == length(args))
      stop("Invalid arguments; use --help.", call. = FALSE)
    if (!is.null(options[[key]])) stop("Repeated argument: ", key)
    options[[key]] <- args[i + 1L]
    i <- i + 2L
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  packages <- c(
    "rmarkdown", "knitr", "dplyr", "tidyr", "ggplot2", "nlme",
    "lspline", "stringr", "kml", "cmprsk", "scales", "patchwork",
    "tibble", "rlang"
  )
  namespace_available <- function(pkg) {
    check <- function() requireNamespace(pkg, quietly = TRUE)
    if (identical(pkg, "kml")) suppressWarnings(check()) else check()
  }
  available <- vapply(packages, namespace_available, logical(1))
  versions <- vapply(seq_along(packages), function(i) {
    if (available[[i]]) as.character(utils::packageVersion(packages[[i]])) else "MISSING"
  }, character(1))
  print(data.frame(package = packages, version = versions), row.names = FALSE)
  if (any(!available)) stop("Missing packages: ", paste(packages[!available], collapse = ", "))
  if (!rmarkdown::pandoc_available()) stop("Pandoc is unavailable.")
  files <- c(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
             list.files(file.path(root, "analysis", "steps"), "[.]R$", full.names = TRUE))
  for (path in files) parse(file = path)
  # Native knitr check: verify every empty report chunk can retrieve its code.
  for (path in list.files(file.path(root, "analysis", "steps"), "[.]R$", full.names = TRUE))
    knitr::read_chunk(path)
  if (check_only) {
    cat("Environment and native R parsing passed. No clinical analysis was run.\n")
    return(invisible(NULL))
  }
  if (is.null(options[["--input"]]) || is.null(options[["--output"]]))
    stop("Supply both --input and --output; use --help.")
  input <- normalizePath(options[["--input"]], mustWork = TRUE)
  output <- options[["--output"]]
  if (file.exists(output) && !dir.exists(output)) stop("Output path is an existing file.")
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)))
    stop("Output directory is not empty; choose a new run directory.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  status <- "FAILED"
  on.exit({
    writeLines(status, file.path(output, "RUN_STATUS.txt"))
    writeLines(capture.output(utils::sessionInfo()), file.path(output, "session-info.txt"))
    utils::write.csv(data.frame(package = packages, version = versions),
                     file.path(output, "package-versions.csv"), row.names = FALSE)
  }, add = TRUE)
  env <- new.env(parent = globalenv())
  run_warning_messages <- character()
  report <- withCallingHandlers(
    rmarkdown::render(
      input = file.path(root, "Workflow.Rmd"), output_file = "RETRACE_report.html",
      output_dir = output, intermediates_dir = output, knit_root_dir = root,
      params = list(project_root = root, data_path = input,
                    output_dir = file.path(output, "artifacts")),
      envir = env, encoding = "UTF-8", quiet = FALSE
    ),
    warning = function(w) {
      run_warning_messages <<- c(run_warning_messages, conditionMessage(w))
    }
  )
  writeLines(
    if (!length(run_warning_messages)) "No warnings." else unique(run_warning_messages),
    file.path(output, "warnings.txt")
  )
  status <- "COMPLETED"
  cat("Report: ", report, "\n", sep = "")
}
main()
