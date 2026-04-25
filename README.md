# CHIP-mutant-burden
# CHIP burden pipeline (script overview)

## Data preparation
- **1. Comparison of baseline characteristics.r**: Merges CHIP gene matrix + modifiable factors + cognitive data; summarizes baseline characteristics by CHIP status; runs univariable logistic regressions and correlation visualization; exports summary tables.
- **2. Search terms data gene.r**: Parses ICD-10 diagnosis strings and corresponding diagnosis dates; reshapes to a wide table (one column per ICD code); exports a code-by-participant date table.
- **3. CHIP data creation.r**: Builds the core CHIP cohort by merging participant data with death registry; derives `Death_event` and `Death_time`; filters to eligible ancestry/CHIP availability; outputs the main CHIP dataset.
- **4. Future disease split.r**: Splits diagnosis dates into *past* vs *future* relative to enrollment date; writes past/future datasets and a non-missing summary.
- **5. ICD10 combine.r**: Collapses duplicated ICD-10 columns (same prefix) by taking the earliest diagnosis date; outputs a combined future-disease dataset.
- **6. CHIP gene matrix data creation.r**: Extracts gene mutation tokens and constructs a binary (0/1) gene matrix per participant; outputs a CHIP gene matrix dataset.

## Survival analysis
- **7. Gene burden cox.r**: Runs Cox models for multiple future disease endpoints using CHIP burden as exposure and baseline covariates; saves Cox results.
- **8. Gene burden with gene cox.r**: Runs Cox models including CHIP burden and multiple gene indicators simultaneously; applies basic checks for variable variation and events-to-predictors; saves results.
- **9 cox survival death.R**: Performs survival analysis for all-cause death: KM curves by CHIP groups and Cox models with CHIP burden (continuous and binary).

## Result review & visualization
- **10. Result check.r**: Filters Cox results by event count threshold; applies multiple-testing correction; maps outcomes to ICD chapters using `coding19.tsv`; produces summary exports and plots.

## Gene contribution / overlap
- **11. Special contribution of genes.R**: Summarizes mutation frequencies for selected genes; analyzes significant outcomes per gene and visualizes intersections with an UpSet plot (optionally colored by ICD chapter).

## Protein analyses
- **12. Protein continuous variables.r**: Tests associations between proteins and CHIP (binary logistic and continuous partial-correlation style analysis); adjusts p-values; exports results and produces volcano/Venn plots.
- **13. protein function prediction.r**: Performs GO/KEGG/DO enrichment analysis for significant proteins; produces enrichment plots.

## Modifiable factors
- **14. Comparative adjustable factors.r**: Screens modifiable factors for association with CHIP (binary logistic + continuous partial correlation); adjusts p-values; exports results and produces volcano/Venn plots.

## Mediation / SEM-style analysis
- **15. SEM.R**: Runs mediation analysis (exposure → CHIP → disease) using a fixed follow-up window; loops over multiple exposures and outcomes; exports mediation results.
