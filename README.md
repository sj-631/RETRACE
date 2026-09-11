# RETRACE

R code for the analysis of early kidney-function recovery trajectories and
renal outcomes in ANCA-associated vasculitis.

## Files

| Path | Purpose |
| --- | --- |
| `Workflow.Rmd` | Main report and analysis order |
| `analysis/steps/` | Analysis code for Cohorts 1–3, figures, tables, and outcome models |
| `R/analysis_helpers.R` | Shared plotting, formatting, and export functions |
| `run_analysis.R` | Command-line runner with configurable input and output paths |
| `.gitignore` | Excludes clinical data, generated outputs, and R session files |

`Workflow.Rmd` is the only report entry point. It loads the numbered scripts in
`analysis/steps/` automatically; they do not need to be run individually.

## Run the analysis

From the repository root:

```sh
Rscript run_analysis.R --input "/path/to/harmonized.csv" --output "outputs/run_01"
```

The output directory must be new or empty. The rendered report is saved as
`RETRACE_report.html`, with figures and tables in its `artifacts/` directory.

To check the R environment and parse the analysis without using clinical data:

```sh
Rscript run_analysis.R --check
```

Required R packages are checked by the runner: `rmarkdown`, `knitr`, `dplyr`,
`tidyr`, `ggplot2`, `nlme`, `lspline`, `stringr`, `kml`, `cmprsk`, `scales`,
`patchwork`, `tibble`, and `rlang`. Pandoc is also required for rendering.

Clinical data are not included in the repository.
