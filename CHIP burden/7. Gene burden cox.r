library(survival)
library(data.table)
library(dplyr)

dat1 <- fread('20251211_future_data.csv')
dat <- read.csv('20251211_chip.csv')

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
                'Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants')]

dat2 <- dat2[complete.cases(dat2), ]
dat <- merge(dat2, dat1, by.x = "Participant.ID", by.y = "Participant_ID", all.x = TRUE)

exposure_columns <- names(dat)[7]
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

for (exposure_column in exposure_columns) {
  dat_clean <- dat[!is.na(dat[[exposure_column]]), ]
  for (i in outcome_columns) {
    outcome <- names(dat)[i]
    cat("Processing outcome:", outcome, " with exposure:", exposure_column, " index:", i, "\n")
    
    dat_clean$Time <- as.integer(!is.na(dat_clean[[outcome]]))
    dat_clean$Time_data <- ifelse(dat_clean$Time == 1,
                                  as.numeric(difftime(dat_clean[[outcome]], dat_clean$enrollment_date, units = "days")),
                                  as.numeric(difftime(as.Date("2023-03-31"), dat_clean$enrollment_date, units = "days")))
    
    formula_cox <- as.formula(paste("Surv(Time_data, Time) ~", exposure_column, "+", paste(names(dat)[covariates], collapse = " + ")))
    
    cox_model <- NULL
    try({
      cox_model <- coxph(formula_cox, data = dat_clean, na.action = na.omit)
      summary_model <- summary(cox_model)
      
      hazard_ratios <- exp(coef(cox_model))
      conf_int <- exp(confint.default(cox_model, level = 0.95))
      p_values <- summary_model$coefficients[, 5]
      
      if (exposure_column %in% rownames(summary_model$coefficients)) {
        results <- rbind(results, data.frame(
          Outcome = outcome,
          Exposure = exposure_column,
          HR = hazard_ratios[exposure_column],
          CI_Lower = conf_int[exposure_column, 1],
          CI_Upper = conf_int[exposure_column, 2],
          p_value = p_values[exposure_column],
          N = nobs(cox_model),
          Events = summary_model$nevent,
          AIC = AIC(cox_model),
          Converged = isTRUE(cox_model$converged)
        ))
      }
    }, silent = TRUE)
    
    if (is.null(cox_model)) {
      results <- rbind(results, data.frame(
        Outcome = outcome,
        Exposure = exposure_column,
        HR = NA,
        CI_Lower = NA,
        CI_Upper = NA,
        p_value = NA,
        N = NA_integer_,
        Events = NA_integer_,
        AIC = NA_real_,
        Converged = FALSE
      ))
      cat("Skipping outcome:", outcome, " with exposure:", exposure_column, " due to insufficient sample size or other issues.\n")
    }
  }
}

saveRDS(results, file = "result/cox_results_data_combination.rds")