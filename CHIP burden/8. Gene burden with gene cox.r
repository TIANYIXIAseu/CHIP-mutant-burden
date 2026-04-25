library(survival)
library(data.table)
library(dplyr)

dat1 <- fread('20251211_future_data.csv')
dat <- read.csv('20251211_chip_gene_matrix.csv')

median_value <- median(dat$Body.mass.index..BMI....Instance.0, na.rm = TRUE)
dat$Body.mass.index..BMI....Instance.0[is.na(dat$Body.mass.index..BMI....Instance.0)] <- median_value
dat <- dat %>%
  mutate(Smoking.status...Instance.0 = case_when(
    Smoking.status...Instance.0 == "Current" ~ "Current",
    Smoking.status...Instance.0 == "Never" ~ "Never",
    Smoking.status...Instance.0 == "Previous" ~ "Previous",
    TRUE ~ "Others"
  ))
dat <- dat %>%
  mutate(Alcohol.drinker.status...Instance.0 = case_when(
    Alcohol.drinker.status...Instance.0 == "Current" ~ "Current",
    Alcohol.drinker.status...Instance.0 == "Never" ~ "Never",
    Alcohol.drinker.status...Instance.0 == "Previous" ~ "Previous",
    TRUE ~ "Others"
  ))

dat2 <- dat[, c('Participant.ID', 'Age.at.recruitment', 'Sex',
                'Smoking.status...Instance.0', 'Alcohol.drinker.status...Instance.0', 'Body.mass.index..BMI....Instance.0',
                'Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants', names(dat)[14:21])]

dat2 <- dat2[complete.cases(dat2), ]
dat <- merge(dat2, dat1, by.x = "Participant.ID", by.y = "Participant_ID", all.x = TRUE)

exposure_columns <- names(dat)[7:15]
covariates <- 2:6
diag_cols <- setdiff(names(dat1), c("Participant_ID", "enrollment_date"))
outcome_columns <- match(diag_cols, names(dat))

dat_clean <- dat

results <- data.frame(Outcome = character(),
                      Exposure = character(),
                      HR = numeric(),
                      CI_Lower = numeric(),
                      CI_Upper = numeric(),
                      p_value = numeric(),
                      N = integer(),
                      Events = integer(),
                      AIC = numeric(),
                      Converged = logical(),
                      stringsAsFactors = FALSE)

for (i in outcome_columns) {
  outcome <- names(dat)[i]
  cat("Processing outcome:", outcome, " with all exposures 7:16\n")
  dat_clean <- dat
  dat_clean$Time <- as.integer(!is.na(dat_clean[[outcome]]))
  dat_clean$Time_data <- ifelse(dat_clean$Time == 1,
                                as.numeric(difftime(dat_clean[[outcome]], dat_clean$enrollment_date, units = "days")),
                                as.numeric(difftime(as.Date("2023-03-31"), dat_clean$enrollment_date, units = "days")))
  var_ok <- sapply(exposure_columns, function(col) {
    v <- dat_clean[[col]]
    u <- unique(v[!is.na(v)])
    length(u) >= 2
  })
  candidate_exposures <- exposure_columns[var_ok]
  cox_model <- NULL
  summary_model <- NULL
  hazard_ratios <- numeric(0)
  conf_int <- matrix(numeric(0), nrow = 0, ncol = 2)
  p_values <- numeric(0)
  if (length(candidate_exposures) > 0) {
    design_formula <- as.formula(paste("~ -1 +", paste(c(candidate_exposures, names(dat)[covariates]), collapse = " + ")))
    mm <- model.matrix(design_formula, data = dat_clean)
    events_count <- sum(dat_clean$Time == 1, na.rm = TRUE)
    predictors_count <- ncol(mm)
    if (events_count < 10L * predictors_count) {
      cat("Skipping outcome:", outcome, " — events=", events_count, " predictors=", predictors_count, " (< 10x).\n", sep = "")
    } else {
      formula_cox <- as.formula(paste("Surv(Time_data, Time) ~",
                                      paste(c(candidate_exposures, names(dat)[covariates]), collapse = " + ")))
      try({
        cox_model <- coxph(formula_cox, data = dat_clean, na.action = na.omit)
        summary_model <- summary(cox_model)
        hazard_ratios <- exp(coef(cox_model))
        conf_int <- exp(confint.default(cox_model, level = 0.95))
        p_values <- summary_model$coefficients[, 5]
      }, silent = TRUE)
    }
  } else {
    cat("Skipping outcome:", outcome, " — no varying exposures among 7:16.\n")
  }
  for (exposure_column in exposure_columns) {
    hr <- if (!is.null(cox_model) && exposure_column %in% names(hazard_ratios)) hazard_ratios[exposure_column] else NA_real_
    ci_l <- if (!is.null(cox_model) && exposure_column %in% rownames(conf_int)) conf_int[exposure_column, 1] else NA_real_
    ci_u <- if (!is.null(cox_model) && exposure_column %in% rownames(conf_int)) conf_int[exposure_column, 2] else NA_real_
    pv <- if (!is.null(cox_model) && exposure_column %in% names(p_values)) p_values[exposure_column] else NA_real_
    results <- rbind(results, data.frame(
      Outcome = outcome,
      Exposure = exposure_column,
      HR = hr,
      CI_Lower = ci_l,
      CI_Upper = ci_u,
      p_value = pv,
      N = if (!is.null(cox_model)) nobs(cox_model) else NA_integer_,
      Events = if (!is.null(summary_model)) summary_model$nevent else NA_integer_,
      AIC = if (!is.null(cox_model)) AIC(cox_model) else NA_real_,
      Converged = if (!is.null(cox_model)) isTRUE(cox_model$converged) else FALSE
    ))
  }
}

results$FDR_BH <- p.adjust(results$p_value, method = "BH")
results$P_bonf <- p.adjust(results$p_value, method = "bonferroni")
write.csv(results, file = "cox_results_chip_gene.csv", row.names = FALSE)
saveRDS(results, file = "cox_results_chip_gene.rds")