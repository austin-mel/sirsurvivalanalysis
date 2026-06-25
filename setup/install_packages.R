required_packages <- c(
  "survival",
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "broom",
  "rmarkdown",
  "knitr"
)

analysis_dir <- file.path("analysis", "survival_baseline")
local_library <- file.path(analysis_dir, ".r_libs")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(local_library), .libPaths()))
Sys.setenv(R_LIBS_USER = normalizePath(local_library))

available_packages <- available.packages(repos = "https://cloud.r-project.org")
dependency_map <- tools::package_dependencies(
  required_packages,
  db = available_packages,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)
needed_packages <- unique(c(required_packages, unlist(dependency_map, use.names = FALSE)))
needed_packages <- setdiff(needed_packages, "R")
local_packages <- rownames(installed.packages(lib.loc = normalizePath(local_library)))
missing_packages <- setdiff(needed_packages, local_packages)

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    lib = normalizePath(local_library),
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

status <- data.frame(
  package = required_packages,
  installed = vapply(required_packages, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)

print(status)

if (!all(status$installed)) {
  stop("Some required packages are still missing after installation.", call. = FALSE)
}
