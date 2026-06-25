# Reproducible Baseline SIRS Survival Analysis

This directory contains a reproducible baseline SIRS simulation with mortality
and a complete cell-level survival analysis.

## Baseline Scenario

- Matrix: `10 x 10`
- Initial infection: center seeded with `create_cntr_matrix(10, 10)`
- Model: `SIRS`
- Infection probability: `0.25`
- Immunity probability: `0.70`
- Mortality enabled: `TRUE`
- Fatality probability: `0.15`
- Seed: `94128`
- Full per-cell log: `TRUE`

Each cell is treated as a subject. The survival endpoints are time to first
infection, time to first recovery, and time to death. Cells that do not
experience an endpoint before the simulation ends are censored at their last
observed step.

## Setup

From the project root, install or check the analysis packages:

```powershell
Rscript analysis/survival_baseline/setup/install_packages.R
```

The setup script installs required analysis packages into the project-local
library `analysis/survival_baseline/.r_libs` so the setup does not require admin
rights.

Then run the full simulation and render both reports:

```powershell
Rscript analysis/survival_baseline/run_all.R
```

If `Rscript` is not on your terminal `PATH`, open RStudio or an R console at the
project root and run:

```r
source("analysis/survival_baseline/setup/install_packages.R")
source("analysis/survival_baseline/run_all.R")
```

On this Windows machine, an explicit Rscript path may also work:

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" analysis/survival_baseline/run_all.R
```

The runner auto-detects the RStudio-bundled Pandoc path used for rendering
`.Rmd` reports on this machine.

## Outputs

Running `run_all.R` writes:

- `outputs/raw/full_log.csv`
- `outputs/raw/history.csv`
- `outputs/raw/final_matrix.csv`
- `outputs/derived/survival_endpoints.csv`
- `outputs/figures/simulation_history.png`
- `outputs/figures/final_matrix.png`
- `outputs/figures/km_infection.png`
- `outputs/figures/km_recovery.png`
- `outputs/figures/km_death.png`
- `outputs/figures/km_infection_initial_status.png`
- `outputs/figures/km_death_distance_group.png`
- `outputs/01_simulation_run.html`
- `outputs/02_survival_analysis.html`
- `outputs/commands_run.txt`
- `outputs/session_info.txt`

The command log and session info are written on every run so the exact
execution context is preserved alongside the simulation outputs.
