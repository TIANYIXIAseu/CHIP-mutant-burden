library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(mediation)

file_chip <- "20251211_chip.csv"
file_future <- "20251211_future_data_combination.csv"
file_regular <- "regular indicators and modifiable behaviors.csv"

dat_chip <- fread(file_chip)
dat_future <- fread(file_future)
dat_regular <- fread(file_regular)

if ("Participant ID" %in% names(dat_chip)) setnames(dat_chip, "Participant ID", "ID")

if ("Participant_ID" %in% names(dat_future)) setnames(dat_future, "Participant_ID", "ID")

if ("Participant ID" %in% names(dat_regular)) setnames(dat_regular, "Participant ID", "ID")

cat("Merging datasets...\n")

merged_dat <- merge(dat_chip, dat_future, by = "ID", all = FALSE)

merged_dat <- merge(merged_dat, dat_regular, by = "ID", all = FALSE)

cat("Merge Complete.\n")
cat("Dimensions of merged data:", dim(merged_dat)[1], "rows,", dim(merged_dat)[2], "columns.\n")

dupes <- grep("\\.x$|\\.y$", names(merged_dat), value = TRUE)
if (length(dupes) > 0) {
  cat("Warning: Duplicate column names found (suffixed with .x/.y):\n")
  print(dupes)
}

cat("\nSelecting target variables...\n")

target_cols <- c(
  "ID", 
  "Age at recruitment", 
  "Sex", 
  "Body mass index (BMI) | Instance 0",
  
  "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants", 
  
  "enrollment_date",
  "A41", "J18", "I22", "I48", "I50", "I60", "C34", "C15", "M18", "N39", "J44", "E87", 
  "D46", "C92", "C93", "D45", "D47", "D61", "D69", "D70", "D64", "C83", 
  
  "Smoking status",
  "Cereal intake", 
  "Salt added to food", 
  "Coffee intake", 
  "Time spent watching television (TV)", 
  "Frequency of stair climbing in last 4 weeks"
)

missing_cols <- setdiff(target_cols, names(merged_dat))
if (length(missing_cols) > 0) {
  cat("Warning: The following requested columns were NOT found in the merged dataset:\n")
  print(missing_cols)
  cat("Proceeding with available columns.\n")
}

available_cols <- intersect(target_cols, names(merged_dat))
final_dat <- merged_dat[, ..available_cols]

cat("Final dataset dimensions:", dim(final_dat)[1], "rows,", dim(final_dat)[2], "columns.\n")

rm(dat_chip, dat_future, dat_regular, merged_dat)
gc()
merged_dat <- final_dat
cat("Environment cleaned. Only 'final_dat' remains.\n")

if ("Sex" %in% names(merged_dat)) {
  merged_dat[, Sex := fcase(
    Sex == "Male", "0 Male",
    Sex == "Female", "1 Female",
    default = Sex
  )]
}

char_cols <- names(merged_dat)[sapply(merged_dat, is.character)]

for (col in char_cols) {
  vals <- na.omit(unique(merged_dat[[col]]))
  vals <- vals[vals != ""]
  if (length(vals) > 0) {
    if (all(grepl("^\\d+\\s+.*|^\\d+$", vals))) {
      cat("Converting column to numeric (extracting prefix):", col, "\n")
      num_vals <- as.numeric(sub("^(\\d+).*", "\\1", merged_dat[[col]]))
      set(merged_dat, j = col, value = num_vals)
    }
  }
}

for (col in names(merged_dat)) {
  if (is.numeric(merged_dat[[col]])) {
    mean_val <- mean(merged_dat[[col]], na.rm = TRUE)
    if (!is.na(mean_val) && any(is.na(merged_dat[[col]]))) {
      set(merged_dat, i = which(is.na(merged_dat[[col]])), j = col, value = mean_val)
      cat("Imputed column", col, "with mean:", mean_val, "\n")
    }
  }
}


if ("Age at recruitment" %in% names(merged_dat)) setnames(merged_dat, "Age at recruitment", "age")
if ("Sex" %in% names(merged_dat)) setnames(merged_dat, "Sex", "sex")
if ("Body mass index (BMI) | Instance 0" %in% names(merged_dat)) setnames(merged_dat, "Body mass index (BMI) | Instance 0", "bmi")
if ("Clonal haematopoiesis of indeterminate potential (CHIP) number of variants" %in% names(merged_dat)) setnames(merged_dat, "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants", "CHIP_nvar")
if ("Smoking status" %in% names(merged_dat)) setnames(merged_dat, "Smoking status", "smoking")
if ("Cereal intake" %in% names(merged_dat)) setnames(merged_dat, "Cereal intake", "cereal")
if ("Salt added to food" %in% names(merged_dat)) setnames(merged_dat, "Salt added to food", "salt")
if ("Coffee intake" %in% names(merged_dat)) setnames(merged_dat, "Coffee intake", "coffee")
if ("Time spent watching television (TV)" %in% names(merged_dat)) setnames(merged_dat, "Time spent watching television (TV)", "tv")
if ("Frequency of stair climbing in last 4 weeks" %in% names(merged_dat)) setnames(merged_dat, "Frequency of stair climbing in last 4 weeks", "stairs")
if ("enrollment_date" %in% names(merged_dat)) setnames(merged_dat, "enrollment_date", "enroll_date")

if ("CHIP_nvar" %in% names(merged_dat)) merged_dat[, CHIP_bin := as.integer(CHIP_nvar > 0)]

if ("enroll_date" %in% names(merged_dat)) {
  merged_dat[, enroll_date := as.Date(enroll_date)]

  disease_codes <- c("A41", "J18", "I22", "I48", "I50", "I60", "C34", "C15", "M18", "N39", "J44", "E87",
                     "D46", "C92", "C93", "D45", "D47", "D61", "D69", "D70", "D64", "C83")

  cutoff_date <- as.Date("2023-03-31")

  for (code in disease_codes) {
    if (code %in% names(merged_dat)) {
      date_col <- paste0(code, "_date")
      time_col <- paste0("time_", code)
      event_col <- paste0("event_", code)
      merged_dat[, (date_col) := as.Date(get(code))]
      merged_dat[, (time_col) := as.numeric(fifelse(!is.na(get(date_col)), get(date_col), cutoff_date) - enroll_date)]
      merged_dat[, (event_col) := as.integer(!is.na(get(date_col)))]
    }
  }

  time_windows <- c(10)

  for (window in time_windows) {
    for (code in disease_codes) {
      time_col <- paste0("time_", code)
      event_col <- paste0("event_", code)
      eventW_col <- paste0("event", window, "y_", code)
      if (all(c(time_col, event_col) %in% names(merged_dat))) {
        merged_dat[, (eventW_col) := as.integer(get(event_col) == 1 & get(time_col) <= window * 365.25)]
      }
    }
  }
}

exposures <- c("smoking", "cereal", "salt", "coffee", "tv", "stairs")
mediator <- "CHIP_bin"
outcomes <- c("A41", "J18", "I22", "I48", "I50", "I60", "C34", "C15", "M18", "N39", "J44", "E87",
              "D46", "C92", "C93", "D45", "D47", "D61", "D69", "D70", "D64", "C83")
covariates <- c("age", "sex", "bmi")

run_mediation_analysis <- function(
  data,
  exposures,
  mediator,
  outcomes,
  time_windows = c(10),
  covariates = c("age", "sex", "bmi"),
  sims = 100,
  boot = FALSE
) {
  results_list <- list()

  for (window in time_windows) {
    cat("\n=== Time window:", window, "years ===\n")

    for (exp in exposures) {
      if (!exp %in% names(data)) next
      cat("Exposure:", exp, "\n")

      for (outcome in outcomes) {
        event_col <- paste0("event", window, "y_", outcome)
        if (!event_col %in% names(data)) next

        formula_M <- as.formula(
          paste(mediator, "~", paste(c(covariates, exposures), collapse = " + "))
        )

        formula_Y <- as.formula(
          paste(event_col, "~", mediator, "+", paste(c(covariates, exposures), collapse = " + "))
        )

        fit_M <- tryCatch(glm(formula_M, family = binomial, data = data),
                          error = function(e) NULL)
        fit_Y <- tryCatch(glm(formula_Y, family = binomial, data = data),
                          error = function(e) NULL)

        if (is.null(fit_M) || is.null(fit_Y)) next

        med_result <- tryCatch(
          mediate(
            model.m = fit_M,
            model.y = fit_Y,
            treat = exp,
            mediator = mediator,
            boot = boot,
            sims = sims
          ),
          error = function(e) NULL
        )

        if (is.null(med_result)) next

        df <- data.frame(
          exposure = exp,
          mediator = mediator,
          outcome = outcome,
          time_window = window,
          ACME_est = med_result$d0,
          ACME_lower = med_result$d0.ci[1],
          ACME_upper = med_result$d0.ci[2],
          ACME_p = med_result$d0.p,      
          ADE_est = med_result$z0,
          ADE_lower = med_result$z0.ci[1],
          ADE_upper = med_result$z0.ci[2],
          ADE_p = med_result$z0.p,
          total_effect = med_result$tau.coef,
          total_lower = med_result$tau.ci[1],
          total_upper = med_result$tau.ci[2],
          prop_mediated = med_result$n0
        )
        
        df <- df %>%
          mutate(
            ACME_OR = exp(ACME_est),
            ADE_OR = exp(ADE_est),
            Total_OR = exp(total_effect),
            prop_mediated_percent = prop_mediated * 100
          )

        cat(
          "  Outcome:", outcome,
          "| ACME:", round(df$ACME_est, 4), "(p=", signif(df$ACME_p, 3), ")",
          "| ADE:", round(df$ADE_est, 4),
          "| Total:", round(df$total_effect, 4),
          "| Prop_med:", round(df$prop_mediated_percent, 4), "%",
          "\n"
        )

        results_list[[paste(exp, outcome, window, sep = "_")]] <- df
      }
    }
  }

  if (length(results_list) == 0) {
    cat("No mediation results were generated.\n")
    return(NULL)
  }

  results <- bind_rows(results_list)
  return(results)
}

med_results <- run_mediation_analysis(
  data = merged_dat,
  exposures = exposures,
  mediator = mediator,
  outcomes = outcomes,
  time_windows = c(10),
  covariates = covariates,
  sims = 100,
  boot = FALSE
)

if (!is.null(med_results)) {
  if (!dir.exists("results")) dir.create("results", recursive = TRUE)
  fwrite(med_results, file.path("results", "med_results.csv"))
}

