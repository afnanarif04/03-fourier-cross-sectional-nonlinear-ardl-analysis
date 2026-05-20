# Cross-Sectionally Augmented Nonlinear Distributed-Lag Estimation with Fourier Approximation

This repository provides the code, data, and results to **(1) replicate the Monte Carlo simulation** and **(2) apply the estimator to panel data**. It accompanies a working paper on cross-sectionally augmented nonlinear distributed-lag panel estimation under cross-sectional dependence, asymmetric adjustment, and smooth structural change.

The materials are self-contained and require only a standard R installation.

---

## Contents

| File | Description |
|------|-------------|
| `01_monte_carlo_simulation.R` | Monte Carlo simulation. Generates panel data with a common-factor structure, asymmetric long-run effects, and smooth structural change, then reports bias and RMSE across several panel sizes. No external data required. The actual results are documented at the bottom of the script. |
| `02_empirical_estimation.R` | Applies the three estimators of the distributed-lag family (CS-DL, CS-NDL, FCS-NDL) to a panel dataset in an Excel file. Reports mean-group long-run elasticities, standard errors, and a Wald test of long-run symmetry. |
| `mc_results.csv` | Output of the Monte Carlo script (the values it produces). |
| `empirical_results.csv` | Output of the empirical script for the example dataset. |
| `data.xlsx` | Example panel dataset in the required four-column format (`id`, `time`, `y`, `x`). |
| `DATA_SOURCE.md` | The public source of the data used in the paper, so results can be independently verified. |
| `LICENSE` | Terms of use. |

---

## Requirements

- **R** version 4.0 or later — <https://cran.r-project.org/>
- Optional: **RStudio** — <https://posit.co/download/rstudio-desktop/>

The scripts install any missing packages automatically the first time they run. The only packages used are `MASS` (simulation) and `readxl` (reading Excel). No figures are produced; the scripts report numerical results only.

---

## How to reproduce

### 1. Monte Carlo simulation

```r
source("01_monte_carlo_simulation.R")
```

Runtime is a few minutes. Results print to the console and are written to `mc_results.csv`. The same numbers are recorded in the `REFERENCE OUTPUT` block at the bottom of the script.

### 2. Empirical estimation

Place your data in `data.xlsx` (four columns: `id`, `time`, `y`, `x`), then run:

```r
source("02_empirical_estimation.R")
```

Results print to the console and are written to `empirical_results.csv`.

---

## Data format

The Excel file must contain four columns:

| Column | Meaning |
|--------|---------|
| `id`   | Identifier for each cross-sectional unit (e.g., country) |
| `time` | Time period (e.g., year) |
| `y`    | Dependent variable |
| `x`    | Explanatory variable of interest |

To apply the method to a different dataset, replace `data.xlsx` (keeping the same column structure) and re-run `02_empirical_estimation.R`. If your columns have different names, edit the `USER SETTINGS` block at the top of the script.

> The `data.xlsx` shipped here is a small illustrative panel so the scripts run instantly. Replace it with the full panel built from the source in `DATA_SOURCE.md` to reproduce the paper's reported values.

---

## License

Released for academic, non-commercial use. See `LICENSE`.
