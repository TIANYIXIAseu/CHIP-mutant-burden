library(survival)
library(data.table)
library(dplyr)
library(survminer)
library(ggplot2)

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

dat <- dat[, c('Participant.ID', 'Age.at.recruitment', 'Sex',
                'Smoking.status...Instance.0', 'Alcohol.drinker.status...Instance.0', 'Body.mass.index..BMI....Instance.0',
                'Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants', "Death_event", "Death_time")]

dat <- dat[complete.cases(dat), ]

if ("Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants" %in% names(dat)) {
  dat$CHIP_Count <- dat$Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants
  dat$CHIP_Binary <- ifelse(dat$CHIP_Count > 0, "CHIP+", "CHIP-")
  
  if (all(c("Death_time", "Death_event") %in% names(dat))) {
    dat$Death_time_Years <- dat$Death_time / 365.25
    
    dat$CHIP_Group <- ifelse(dat$CHIP_Count >= 3, "3+", as.character(dat$CHIP_Count))
    dat$CHIP_Group <- factor(dat$CHIP_Group, levels = c("0", "1", "2", "3+"))
    
    if (length(unique(na.omit(dat$CHIP_Group))) > 1) {
      cat("Checking sample sizes for CHIP Groups:\n")
      print(table(CHIP_Group = dat$CHIP_Group, Death_event = dat$Death_event))
      
      fit_count <- survfit(Surv(Death_time_Years, Death_event) ~ CHIP_Group, data = dat)
      p1 <- ggsurvplot(
        fit_count,
        data = dat,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        xlab = "Time (Years)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(),
        palette = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
      )
      print(p1)
    }
    
    cat("Plotting KM Curve: CHIP- vs CHIP+ (Time in Years)...\n")
    print(table(CHIP_Binary = dat$CHIP_Binary, Death_event = dat$Death_event))
    fit_binary <- survfit(Surv(Death_time_Years, Death_event) ~ CHIP_Binary, data = dat)
    p2 <- ggsurvplot(
      fit_binary,
      data = dat,
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      xlab = "Time (Years)",
      ylab = "Survival Probability",
      ggtheme = theme_minimal(),
      linetype = "solid",
      palette = c("#2C6DB2", "#C6295C")
    )
    print(p2)
    cat("Survival analysis complete.\n")
  } else {
    cat("Survival columns (Death_time, Death_event) missing in dat, skipping KM plots.\n")
  }
} else {
  cat("CHIP Count column not found, skipping survival analysis.\n")
}

covars_ok <- all(c("Age.at.recruitment","Sex","Body.mass.index..BMI....Instance.0",
                   "CHIP_Count","Smoking.status...Instance.0",
                   "Alcohol.drinker.status...Instance.0","Death_time","Death_event") %in% names(dat))
if (covars_ok) {
  cox_cont <- coxph(Surv(Death_time, Death_event) ~ CHIP_Count +
                      Age.at.recruitment + Sex +
                      Body.mass.index..BMI....Instance.0 +
                      Alcohol.drinker.status...Instance.0 +
                      Smoking.status...Instance.0,
                    data = dat)
  cox_bin <- coxph(Surv(Death_time, Death_event) ~ CHIP_Binary +
                     Age.at.recruitment + Sex +
                     Body.mass.index..BMI....Instance.0 +
                     Alcohol.drinker.status...Instance.0 +
                     Smoking.status...Instance.0,
                   data = dat)
  s1 <- summary(cox_cont)
  s2 <- summary(cox_bin)
  cat("\nCox model with CHIP_Count (continuous):\n")
  print(s1)
  cat("\nCox model with CHIP_Binary (0/1):\n")
  print(s2)
}
