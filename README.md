# Baseline SIRS Survival Analysis

## Project Title And Short Summary

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

**Replication blockers:** The current checkout cannot complete the full
workflow without changes to its path and simulator setup.

- The scripts use `analysis/survival_baseline`. The runner expects report
  sources in `analysis/survival_baseline/reports`, but both `.Rmd` files sit at
  the repository root.
- Both `.Rmd` files search for a simulator package root containing
  `DESCRIPTION` and `R/`. This checkout contains neither.
- The runner selects `pkgload::load_all(".")` whenever `pkgload` is available.
  If that call fails, execution stops. Only when `pkgload` is unavailable does
  the runner try installed `SIRSsim`, then local `R/` sources.

Installing `SIRSsim` alone does not resolve these blockers. The runner needs
its `create_cntr_matrix()` and `simulate_sir()` functions. This checkout does
not supply the simulator source or a package installer.

**Requirements:** R, the `SIRSsim` dependency, and Pandoc for HTML rendering.
The [package setup script](setup/install_packages.R) installs `survival`,
`dplyr`, `tidyr`, `readr`, `ggplot2`, `broom`, `rmarkdown`, `knitr`, and their
required CRAN dependencies into `analysis/survival_baseline/.r_libs`. It does
not install `SIRSsim` or Pandoc.

From the repository root, install the R packages with:

```sh
Rscript setup/install_packages.R
```

After resolving the replication blockers, run the analysis from the same
directory:

```sh
Rscript run_all.R
```

These commands work in PowerShell and other shells when `Rscript` is on
`PATH`. Alternatively, use an R console with the repository root as its
working directory. The same replication blockers apply:

```r
source("setup/install_packages.R")
source("run_all.R")
```

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
