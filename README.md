# Analysis code: discharge motor FIM in patients aged 90 years and over

R scripts for the analyses reported in

> Motor Functional Independence Measure scores at discharge in patients aged 90 years and over in a Japanese convalescent rehabilitation ward: a retrospective cohort study.

The study compared the discharge motor Functional Independence Measure (FIM) scores of patients aged ≥90 years with those of patients aged 65–74 and 75–89 years admitted to a convalescent rehabilitation ward in Japan (April 2017 to March 2024). The statistical analysis plan (version 1, 20 September 2026; version 2, 22 September 2026) is available on Zenodo: [DOI to be added].

## Data availability

Individual patient data are **not** included in this repository and cannot be shared, because of the scope of the ethical approval and the protection of personal information. The scripts can be read to check every analysis step, and they can be run on a dataset with the same structure (see below).

## Repository contents

| Path | Description |
|---|---|
| `script/` | R scripts (see the table below) |
| `sessionInfo.md` | R session information for the analyses |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT License |
| `.gitignore` | Excludes `output/` and other files derived from patient-level data |

## Requirements

- R 4.6.1 (`07_1_sens_ranktest.R` was run on Ubuntu 24.04.4 LTS and `12_figures_tables.R` on Windows 11 x64; see `sessionInfo.md`)
- Packages: `ordinal` (2026.7-26), `rankFD` (0.1.1), `splines` (base R), `ggplot2` (4.0.3; ≥3.4 required), `patchwork` (1.3.2), `flextable` (0.10.1), `officer` (0.7.6)
- Optional: `parallel` (parallel bootstrap), `ragg`, `systemfonts`, `magick` (figure output); `MASS` (only for `08_po_thresholds.R`)

See `sessionInfo.md` for the full list of package versions.

## Input data

The scripts expect a data frame named `Alldata` in the R workspace (one row per admission) with the following columns:

| Column | Content |
|---|---|
| `id` | Patient identifier |
| `age` | Age at admission (years) |
| `sex` | Sex (`男` male, `女` female) |
| `day_in`, `day_out` | Dates of admission and discharge (`day_out` is optional; it is used for checks and for the length of stay in Table 1) |
| `disposition` | Discharge destination (`病院・診療所へ転院` and `医療機関`, transfer; `終了（死亡等）`, death; other values, discharge to the community) |
| `class` | Disease category (`脳血管` cerebrovascular disease, `運動器` musculoskeletal disorder, `廃用` disuse syndrome) |
| `support_in` | Pre-admission care need (`あり` certified as requiring long-term care, `なし` not certified or requiring support) |
| `mFIM_in`, `cFIM_in` | Motor and cognitive FIM scores at admission |
| `mFIM_out` | Motor FIM score at discharge |

Level values are in Japanese, as in the source database; English labels are defined in `script/01_labels.R`.

## How to run

1. Set the working directory to `script/` and load `Alldata` into the workspace.
2. Run the scripts in the order below with `source("<file>", encoding = "UTF-8")`. Each script reads the `.rda` files written by earlier scripts and writes its results, checks and a log to `output/`. The `output/` folder contains patient-level data and must not be committed (it is listed in `.gitignore`).

| Step | Script | Purpose | Main outputs used in the manuscript |
|---|---|---|---|
| — | `01_labels.R` | English labels for all outputs (sourced by the other scripts) | — |
| 1 | `02_preprocess.R` | Data cleaning, selection of participants, exclusions, sensitivity analysis dataset | Figure 1 counts |
| 2 | `03_m0_ranktest.R` | Two-sample relative effects and permutation Brunner–Munzel test | Figure 2 |
| 3 | `04_ridit_spline.R` | Missing-data rules, admission period, ridits and spline knots | Appendix 1 |
| 4 | `05_m2_clm.R` | Partial proportional odds and proportional odds models, profile-likelihood intervals | Table 2, Table S4 |
| 5 | `06_standardize_boot.R` | Marginal standardisation and stratified bootstrap | Table 2, Table S2 |
| 6 | `07_1_sens_ranktest.R` | Sensitivity analysis: relative effects | Table S2 |
| 7 | `07_2_sens_clm_boot.R` | Sensitivity analysis: models and standardisation | Table S2 |
| 8 | `09_positivity.R` | Overlap of covariate distributions | Table S3 |
| 9 | `10_age_continuous.R` | Proportional odds model with age as a continuous variable | Figure S1 |
| 10 | `11_supp_tables.R` | Missing data, excluded patients and model coefficients | Tables S1 and S4 |
| 11 | `12_figures_tables.R` | Figures and tables (Word and TIFF) | Tables 1–2, Figure 2, Tables S1–S4, Figure S1 |

`08_po_thresholds.R` implements the check of the proportional odds assumption specified in version 1 of the analysis plan (threshold-specific binary logistic models at 40, 65 and 80 points). This check was replaced when the plan was revised (version 2); the script is not needed to reproduce the final tables and figures and is kept for transparency. Figure 1 (flow diagram) was drawn separately from the counts produced by `02_preprocess.R`.

Comments in the scripts are written in Japanese, and section numbers (§) in the comments refer to the analysis plan. All labels in the output tables and figures are in English.

## Citation

Please cite the article and this repository (see `CITATION.cff`). Zenodo DOI: [DOI to be added].

## License

The code is released under the MIT License (see `LICENSE`).
