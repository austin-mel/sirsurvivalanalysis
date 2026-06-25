analysis_dir <- file.path("analysis", "survival_baseline")
output_dir <- file.path(analysis_dir, "outputs")
raw_dir <- file.path(output_dir, "raw")
derived_dir <- file.path(output_dir, "derived")
figure_dir <- file.path(output_dir, "figures")
report_dir <- file.path(analysis_dir, "reports")
local_library <- file.path(analysis_dir, ".r_libs")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(local_library)) {
  .libPaths(c(normalizePath(local_library), .libPaths()))
}

pandoc_candidates <- c(
  file.path(Sys.getenv("RSTUDIO_PANDOC"), "pandoc.exe"),
  file.path("C:", "Program Files", "RStudio", "resources", "app", "bin", "quarto", "bin", "tools", "pandoc.exe"),
  file.path("C:", "Program Files", "RStudio", "resources", "app", "bin", "pandoc", "pandoc.exe")
)
pandoc_path <- pandoc_candidates[file.exists(pandoc_candidates)][1]
if (!is.na(pandoc_path) && nzchar(pandoc_path)) {
  Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_path))
}

required_packages <- c("survival", "dplyr", "tidyr", "readr", "ggplot2", "broom", "rmarkdown", "knitr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Run analysis/survival_baseline/setup/install_packages.R first.",
    call. = FALSE
  )
}

load_local_sirs <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }

  if (requireNamespace("SIRSsim", quietly = TRUE)) {
    library(SIRSsim)
    return(invisible(TRUE))
  }

  source(file.path("R", "states.R"))
  source(file.path("R", "neighbors.R"))
  source(file.path("R", "matrices.R"))
  source(file.path("R", "simulation.R"))
  source(file.path("R", "plotting.R"))
  invisible(TRUE)
}

derive_survival_endpoints <- function(full_log) {
  center_row <- (max(full_log$row) + 1) / 2
  center_col <- (max(full_log$col) + 1) / 2

  dplyr::summarise(
    dplyr::group_by(full_log, cell_id),
    row = dplyr::first(row),
    col = dplyr::first(col),
    initial_state = state[step == 0][1],
    initially_infected = state[step == 0][1] == 1,
    distance_from_center = sqrt((dplyr::first(row) - center_row)^2 + (dplyr::first(col) - center_col)^2),
    first_infected_step = if (any(was_infected)) min(step[was_infected]) else NA_integer_,
    recovered_step = if (any(was_immune)) min(step[was_immune]) else NA_integer_,
    death_step = if (any(was_deceased)) min(step[was_deceased]) else NA_integer_,
    last_observed_step = max(step),
    event_infection = !is.na(first_infected_step),
    event_recovery = !is.na(recovered_step),
    event_death = !is.na(death_step),
    time_to_infection = ifelse(event_infection, first_infected_step, last_observed_step),
    time_to_recovery = ifelse(event_recovery, recovered_step, last_observed_step),
    time_to_death = ifelse(event_death, death_step, last_observed_step),
    final_state = state[step == max(step)][1],
    .groups = "drop"
  )
}

validate_outputs <- function(result, endpoints) {
  if (!"full_log" %in% names(result)) {
    stop("Simulation result does not contain full_log.", call. = FALSE)
  }

  expected_rows <- 100 * (result$steps + 1)
  if (nrow(result$full_log) != expected_rows) {
    stop("full_log row count mismatch. Expected ", expected_rows, ", got ", nrow(result$full_log), ".", call. = FALSE)
  }

  if (nrow(endpoints) != 100 || length(unique(endpoints$cell_id)) != 100) {
    stop("survival_endpoints must contain exactly one row for each of 100 cells.", call. = FALSE)
  }

  event_columns <- c("event_infection", "event_recovery", "event_death")
  if (!all(vapply(endpoints[event_columns], is.logical, logical(1)))) {
    stop("Endpoint event indicators must be logical columns.", call. = FALSE)
  }

  time_columns <- c("time_to_infection", "time_to_recovery", "time_to_death")
  if (!all(vapply(endpoints[time_columns], is.numeric, logical(1)))) {
    stop("Endpoint time columns must be numeric.", call. = FALSE)
  }
}

write_matrix_csv <- function(matrix_value, path) {
  readr::write_csv(
    data.frame(row = seq_len(nrow(matrix_value)), matrix_value, check.names = FALSE),
    path
  )
}

save_figures <- function(history, final_matrix, endpoints) {
  history_long <- tidyr::pivot_longer(
    history,
    cols = c("susceptible", "infected", "recovered"),
    names_to = "state",
    values_to = "cells"
  )

  history_plot <- ggplot2::ggplot(history_long, ggplot2::aes(step, cells, color = state)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(
      title = "Baseline SIRS Simulation History",
      x = "Simulation step",
      y = "Cells",
      color = "State"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(figure_dir, "simulation_history.png"), history_plot, width = 7, height = 4.5, dpi = 300)

  final_grid <- data.frame(
    row = rep(seq_len(nrow(final_matrix)), times = ncol(final_matrix)),
    col = rep(seq_len(ncol(final_matrix)), each = nrow(final_matrix)),
    state = factor(
      as.vector(final_matrix),
      levels = c(0, 1, 2, 3),
      labels = c("Susceptible", "Infected", "Recovered", "Deceased")
    )
  )

  final_plot <- ggplot2::ggplot(final_grid, ggplot2::aes(col, row, fill = state)) +
    ggplot2::geom_tile(color = "grey85", linewidth = 0.2) +
    ggplot2::scale_y_reverse(breaks = seq_len(nrow(final_matrix))) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Final Matrix State", x = "Column", y = "Row", fill = "State") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(figure_dir, "final_matrix.png"), final_plot, width = 5.5, height = 5.5, dpi = 300)

  save_km_plot(endpoints, "time_to_infection", "event_infection", "Time to First Infection", "km_infection.png")
  save_km_plot(endpoints, "time_to_recovery", "event_recovery", "Time to Recovery", "km_recovery.png")
  save_km_plot(endpoints, "time_to_death", "event_death", "Time to Death", "km_death.png")
  save_stratified_km_plot(
    endpoints,
    "time_to_infection",
    "event_infection",
    "initially_infected",
    "Time to First Infection by Initial Infection Status",
    "km_infection_initial_status.png"
  )
  save_stratified_km_plot(
    endpoints,
    "time_to_death",
    "event_death",
    "distance_group",
    "Time to Death by Distance from Center",
    "km_death_distance_group.png"
  )
}

save_km_plot <- function(endpoints, time_column, event_column, title, filename) {
  fit <- survival::survfit(
    survival::Surv(endpoints[[time_column]], endpoints[[event_column]]) ~ 1
  )

  plot_data <- data.frame(
    time = fit$time,
    survival = fit$surv,
    lower = fit$lower,
    upper = fit$upper
  )

  km_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(time, survival)) +
    ggplot2::geom_step(linewidth = 1, color = "#2f6f9f") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, fill = "#2f6f9f") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title, x = "Simulation step", y = "Survival probability") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(figure_dir, filename), km_plot, width = 6.5, height = 4.25, dpi = 300)
}

save_stratified_km_plot <- function(endpoints, time_column, event_column, strata_column, title, filename) {
  plot_data_source <- endpoints
  plot_data_source[[strata_column]] <- factor(plot_data_source[[strata_column]])
  fit <- survival::survfit(
    stats::as.formula(paste0("survival::Surv(", time_column, ", ", event_column, ") ~ ", strata_column)),
    data = plot_data_source
  )

  plot_data <- data.frame(
    time = fit$time,
    survival = fit$surv,
    strata = rep(names(fit$strata), fit$strata)
  )

  km_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(time, survival, color = strata)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title, x = "Simulation step", y = "Survival probability", color = "Strata") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(figure_dir, filename), km_plot, width = 6.8, height = 4.25, dpi = 300)
}

render_report <- function(input_file) {
  rmarkdown::render(
    input = file.path(report_dir, input_file),
    output_dir = output_dir,
    knit_root_dir = normalizePath("."),
    quiet = TRUE
  )
}

load_local_sirs()

initial_matrix <- create_cntr_matrix(row = 10, col = 10)
result <- simulate_sir(
  prob_infect = 0.25,
  input_matrix = initial_matrix,
  model = "SIRS",
  imm_prob = 0.70,
  allow_death = TRUE,
  fat_prob = 0.15,
  seed = 94128,
  full_log = TRUE
)

survival_endpoints <- derive_survival_endpoints(result$full_log)
survival_endpoints <- dplyr::mutate(
  survival_endpoints,
  distance_group = ifelse(distance_from_center <= median(distance_from_center), "Near center", "Far from center")
)
validate_outputs(result, survival_endpoints)

readr::write_csv(result$full_log, file.path(raw_dir, "full_log.csv"))
readr::write_csv(result$history, file.path(raw_dir, "history.csv"))
write_matrix_csv(result$final_matrix, file.path(raw_dir, "final_matrix.csv"))
readr::write_csv(survival_endpoints, file.path(derived_dir, "survival_endpoints.csv"))

save_figures(result$history, result$final_matrix, survival_endpoints)

render_report("01_simulation_run.Rmd")
render_report("02_survival_analysis.Rmd")

writeLines(
  c(
    "Rscript analysis/survival_baseline/setup/install_packages.R",
    "Rscript analysis/survival_baseline/run_all.R",
    "",
    "R console fallback:",
    "source(\"analysis/survival_baseline/setup/install_packages.R\")",
    "source(\"analysis/survival_baseline/run_all.R\")"
  ),
  file.path(output_dir, "commands_run.txt")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

message("Survival baseline analysis completed. Outputs written to ", output_dir)
