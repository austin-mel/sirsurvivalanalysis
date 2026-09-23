# Baseline SIRS Survival Analysis

## Table of Contents

- [Project Summary](#project-summary)
- [Executive Summary](#executive-summary)
- [Getting Started & Replication](#getting-started--replication)
  - [1. Install Prerequisites](#1-install-prerequisites)
  - [2. Clone the Repositories](#2-clone-the-repositories)
  - [3. Prepare the Report Folder](#3-prepare-the-report-folder)
  - [4. Install R Dependencies](#4-install-r-dependencies)
  - [5. Run the Analysis](#5-run-the-analysis)
  - [6. Verify the Results](#6-verify-the-results)
- [Further Details](#further-details)
  - [Simulation Settings](#simulation-settings)
  - [Endpoint Timing](#endpoint-timing)
  - [Report Findings](#report-findings)
  - [Survival Term Definitions](#survival-term-definitions)
  - [Limitations](#limitations)

## Project Summary

This project tracks first infection, first recovery, and death across 100
simulated cells on a `10 x 10` grid. It reports one baseline run. SIRS stands
for susceptible, infected, recovered, and susceptible again.

## Executive Summary

| Area | Summary |
| --- | --- |
| Situation | One center cell starts infected in a simulated population of 100 cells. |
| Problem | Survival analysis asks whether and when infection, recovery, or death occurs during the simulation. |
| Findings and Conclusions | The infected count falls to zero at step 9, when the baseline run ends. |
| Results | The final grid contains 83 susceptible, 0 infected, 13 recovered, and 4 deceased cells. The table below counts events throughout the run. |

The [survival report](02_survival_analysis.html) records these results:

| Event | Cells with the event | Censored cells |
| --- | ---: | ---: |
| First infection | 20 | 80 |
| First recovery | 13 | 87 |
| Death | 4 | 96 |

Censored cells have no recorded event by the final step. Each row tracks a
separate outcome, or endpoint, across all 100 cells. The groups can overlap
because one cell can experience more than one type of event.

## Getting Started & Replication

The analysis depends on the `SIRSsim` R package from
[sirmodelsimulation](https://github.com/austin-mel/sirmodelsimulation).
See its [installation instructions](https://github.com/austin-mel/sirmodelsimulation#installation)
for package setup. The steps below use a compatible simulator revision and
the folder layout that the current analysis scripts require.

### 1. Install Prerequisites

Install R, Git, and Pandoc. RStudio can supply Pandoc for report rendering.
Use a terminal for the Git commands and an R console for the R commands.
These instructions use forward slashes and work across Windows, macOS, and
Linux.

### 2. Clone the Repositories

Start in a new working folder. Clone the simulator, select the compatible
revision, and clone this analysis into its expected location:

```sh
git clone https://github.com/austin-mel/sirmodelsimulation.git
cd sirmodelsimulation
git checkout --detach 43bf4b88ac1d54fd298bf24d41c1ac56b58118d2
git clone https://github.com/austin-mel/sirsurvivalanalysis.git analysis/survival_baseline
```

Revision `43bf4b8` provides `create_cntr_matrix()` and `full_log`, which
`run_all.R` requires. Later simulator revisions rename the matrix helper to
`create_center_matrix()`. Keep the pinned revision for this workflow.

Open `SIRSsim.Rproj` from the simulator folder in RStudio, or start R from
that folder. Keep this working directory and R session for steps 3 through 6.
All paths below start at the **simulator repository root**, which contains
`DESCRIPTION` and `R/`.

### 3. Prepare the Report Folder

The runner reads report templates from `analysis/survival_baseline/reports`.
Copy the two templates there from the analysis root:

```r
stopifnot(
  file.exists("DESCRIPTION"),
  dir.exists("R"),
  file.exists("analysis/survival_baseline/run_all.R")
)

report_dir <- "analysis/survival_baseline/reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_sources <- file.path(
  "analysis/survival_baseline",
  c("01_simulation_run.Rmd", "02_survival_analysis.Rmd")
)
stopifnot(all(file.copy(report_sources, report_dir, overwrite = TRUE)))
```

Repeat this copy after editing either original `.Rmd` file. It replaces the
copies that the runner renders.

### 4. Install R Dependencies

Run the [package setup script](setup/install_packages.R). It installs
`survival`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `broom`, `rmarkdown`, `knitr`,
and their required dependencies into `analysis/survival_baseline/.r_libs`.
Then install `SIRSsim` from the pinned local clone using the simulator's
documented local installation method:

```r
source("analysis/survival_baseline/setup/install_packages.R")
install.packages(
  "remotes",
  repos = "https://cloud.r-project.org",
  lib = .libPaths()[1]
)
remotes::install_local(
  ".",
  lib = .libPaths()[1],
  dependencies = NA,
  upgrade = "never",
  build = FALSE
)

library(SIRSsim)
stopifnot(
  is.function(create_cntr_matrix),
  "full_log" %in% names(formals(simulate_sir))
)
stopifnot(rmarkdown::pandoc_available())
```

The [local installer](https://remotes.r-lib.org/reference/install_local.html)
uses `build = FALSE` to skip archiving the simulator folder, which now also
contains the analysis and its package library.

If the Pandoc check fails, open the project in RStudio or install Pandoc and
add it to `PATH`, then restart R. The
[Pandoc availability check](https://rmarkdown.rstudio.com/docs/reference/pandoc_available.html)
looks for both the RStudio copy and executables on `PATH`.

### 5. Run the Analysis

From the simulator root, run:

```r
source("analysis/survival_baseline/run_all.R")
```

The runner uses seed `94128` and the settings in [Simulation Settings](#simulation-settings).
It generates data, figures, and both HTML reports under
`analysis/survival_baseline/outputs`. It also records the R session and package
versions in `session_info.txt` in that folder. Rerunning replaces these outputs.

The simulator root matters: the runner loads local source through `pkgload`
when available, otherwise it tries installed `SIRSsim`, then local `R/` files.
Both reports also locate the simulator root through `DESCRIPTION` and `R/`.

### 6. Verify the Results

Run these checks in the same R session:

```r
output_dir <- "analysis/survival_baseline/outputs"
endpoints <- read.csv(file.path(output_dir, "derived", "survival_endpoints.csv"))
history <- read.csv(file.path(output_dir, "raw", "history.csv"))
full_log <- read.csv(file.path(output_dir, "raw", "full_log.csv"))
event_counts <- colSums(endpoints[c("event_infection", "event_recovery", "event_death")])

stopifnot(
  nrow(endpoints) == 100,
  length(unique(endpoints$cell_id)) == 100,
  nrow(full_log) == 1000,
  max(history$step) == 9,
  tail(history$infected, 1) == 0,
  all(event_counts == c(20, 13, 4)),
  all(file.exists(file.path(
    output_dir,
    c("01_simulation_run.html", "02_survival_analysis.html", "session_info.txt")
  )))
)
event_counts
```

The checks compare your run with the committed reports: 100 cells, 1,000 log
rows, 9 steps, 20 first infections, 13 recoveries, and 4 deaths. If a check
fails, confirm the simulator revision, working directory, and simulation
settings before interpreting the results.

Open `01_simulation_run.html` and `02_survival_analysis.html` from the new
`outputs` folder in a browser. CRAN packages are not version-locked, so report
formatting and model output can vary with dependency versions. Use the
numerical checks and `session_info.txt` when comparing runs.

## Further Details

### Simulation Settings

| Setting | Value |
| --- | --- |
| Grid size | `10 x 10` cells |
| Initial infection | `create_cntr_matrix(row = 10, col = 10)` seeds one center cell |
| Model | `SIRS` |
| Infection probability | `0.25` |
| Immunity probability | `0.70` |
| Mortality enabled | `TRUE` |
| Fatality probability | `0.15` |
| Random seed | `94128` |
| Full per-cell log | `TRUE` |

The final matrix codes states as `0` susceptible, `1` infected, `2` recovered,
and `3` deceased.

### Endpoint Timing

All three endpoints measure time from simulation step zero across all 100
cells. Recovery time marks the first step with `was_immune = TRUE`. It does
not measure illness duration among infected cells. For each endpoint without
an event, the code uses the cell's last observed step as its censoring time.

### Report Findings

The [simulation report](01_simulation_run.html) records 1,000 log rows across
steps 0 through 9 and 100 endpoint rows. Its checks pass for the expected log
size, one endpoint row per cell, and logical event indicators.

Fewer than half the cells experience each event, so the analysis cannot
estimate a Kaplan-Meier median. The survival report's separate `median_time = 9`
includes censoring times. It is not a median event time.

The infection log-rank test reports `Chisq = 99` and `p < 2e-16`. The comparison
reflects the starting conditions: one cell starts infected at step zero and
the other 99 start uninfected.

The death log-rank test reports rounded values of `Chisq = 0` and `p = 1`.
The test detected no difference in death timing in this run.
Each distance group records two deaths, among 52 near-center cells and 48
far-from-center cells. The code splits groups at the median distance from
the grid center.

### Survival Term Definitions

**Kaplan-Meier curves** estimate the share of cells that have not yet
experienced an event at each step. The infection curve tracks cells with
no previous infection. The recovery curve tracks cells that have not yet
recovered. The death curve tracks cells that remain alive.

**Censoring** keeps cells without a recorded event in the analysis through
their last observed step. It does not tell us whether an event occurs later.

**Log-rank tests** compare event timing between groups using their survival
curves. A test that detects no difference does not establish that groups
behave identically.

**Cox models** relate event rates to row, column, distance from center, and
initial infection status. A hazard ratio above `1` indicates a higher event
rate among cells still awaiting the event. A value below `1` indicates a lower
rate, holding the other model variables constant.

### Limitations

- One baseline run does not establish how results change with seed, grid size,
  or model probabilities.
- The project uses simulated data. It provides no epidemiological evidence
  about a real disease or population.
- The Cox models describe associations within this run and do not support
  causal claims. Four deaths provide limited information for the death model
  and group comparison.
- The Cox output includes missing estimates for initial infection status in
  the infection and recovery models. The death model reports a confidence
  interval from zero to infinity for that variable, which provides no useful
  bound on the association.
