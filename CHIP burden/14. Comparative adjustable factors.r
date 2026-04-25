library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtext)

if(!dir.exists("results")) dir.create("results")

cov_file <- "20251211_chip.csv"
if (!file.exists(cov_file)) stop("20251211_chip.csv not found")

cov_dat <- fread(cov_file, select = c("Participant ID", "Age at recruitment", "Sex", "Body mass index (BMI) | Instance 0", "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants"))

setnames(cov_dat, 
         old = c("Participant ID", "Age at recruitment", "Sex", "Body mass index (BMI) | Instance 0", "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants"), 
         new = c("ID", "Age", "Sex", "BMI", "CHIP_Raw"))

pred_file <- "regular indicators and modifiable behaviors.csv"
if (!file.exists(pred_file)) stop("regular indicators and modifiable behaviors.csv not found")
pred_dat <- fread(pred_file)

if ("Participant ID" %in% names(pred_dat)) {
  setnames(pred_dat, "Participant ID", "ID")
}

merged_dat <- merge(cov_dat, pred_dat, by = "ID")

if (is.character(merged_dat$Sex)) {
  merged_dat[, Sex := ifelse(Sex == "Female", 0, 1)]
}

if (any(is.na(merged_dat$BMI))) {
  merged_dat[is.na(BMI), BMI := mean(BMI, na.rm = TRUE)]
}

merged_dat[, CHIP_Binary := ifelse(CHIP_Raw > 0, 1, 0)]

exclude_cols <- c("ID", "Age", "Sex", "BMI", "CHIP_Raw", "CHIP_Binary")
target_indices <- 84:123
valid_indices <- target_indices[target_indices <= ncol(merged_dat)]

if (length(valid_indices) < length(target_indices)) {
  warning("Some requested predictor indices (84:123) are out of bounds. Using available columns.")
}

if (length(valid_indices) == 0) {
  stop("No valid predictor columns found in range 84:123.")
}

predictor_cols <- names(merged_dat)[valid_indices]
predictor_cols <- setdiff(predictor_cols, exclude_cols)

cat("Selected", length(predictor_cols), "predictors from indices 84:123.\n")

for (col in predictor_cols) {
  if (is.character(merged_dat[[col]])) {
    vals <- na.omit(unique(merged_dat[[col]]))
    vals <- vals[vals != ""]
    
    if (length(vals) > 0 && all(grepl("^\\d", vals))) {
      merged_dat[[col]] <- as.numeric(sub("^(\\d+).*", "\\1", merged_dat[[col]]))
    }
  }
}

results_list_bin <- list()
cat("Starting Logistic Regression (Binary)...\n")

for (x_name in predictor_cols) {
  
  tmp_dat <- merged_dat[, c("CHIP_Binary", x_name, "Age", "Sex", "BMI"), with = FALSE]
  
  setnames(tmp_dat, x_name, "X_Var")
  
  tmp_dat <- na.omit(tmp_dat)
  
  if (nrow(tmp_dat) < 50) next
  if (length(unique(tmp_dat$X_Var)) < 2) next
  if (min(table(tmp_dat$CHIP_Binary)) < 5) next
  
  fit <- try(glm(CHIP_Binary ~ X_Var + Age + Sex + BMI, data = tmp_dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)
    coefs <- s$coefficients
    
    x_rows <- grep("X_Var", rownames(coefs))
    
    for (i in x_rows) {
      term_name <- rownames(coefs)[i]
      clean_name <- sub("X_Var", x_name, term_name)
      
      est <- coefs[i, 1]
      se <- coefs[i, 2]
      pval <- coefs[i, 4]
      
      results_list_bin[[length(results_list_bin) + 1]] <- data.frame(
        Variable = clean_name,
        Base_Var = x_name,
        Beta = est,
        OR = exp(est),
        CI_Lo = exp(est - 1.96 * se),
        CI_Hi = exp(est + 1.96 * se),
        P_Value = pval,
        N = nrow(tmp_dat),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(results_list_bin) > 0) {
  res_bin <- do.call(rbind, results_list_bin)
  res_bin$FDR <- p.adjust(res_bin$P_Value, method = "BH")
  res_bin <- res_bin[order(res_bin$P_Value), ]
  
  write.csv(res_bin, "results/Adjustable_Factors_Binary_Logistic.csv", row.names = FALSE)
  cat("Saved Binary Results.\n")
}

results_list_cont <- list()
cat("Starting Partial Correlation (Continuous)...\n")

for (x_name in predictor_cols) {
  
  if (!is.numeric(merged_dat[[x_name]])) next
  
  tmp_dat <- merged_dat[, c("CHIP_Raw", x_name, "Age", "Sex", "BMI"), with = FALSE]
  setnames(tmp_dat, x_name, "X_Var")
  tmp_dat <- na.omit(tmp_dat)
  
  if (nrow(tmp_dat) < 50) next
  if (var(tmp_dat$X_Var) == 0) next
  
  fit_y <- lm(CHIP_Raw ~ Age + Sex + BMI, data = tmp_dat)
  res_y <- residuals(fit_y)
  
  fit_x <- lm(X_Var ~ Age + Sex + BMI, data = tmp_dat)
  res_x <- residuals(fit_x)
  
  ct <- cor.test(res_y, res_x, method = "pearson")
  
  results_list_cont[[length(results_list_cont) + 1]] <- data.frame(
    Variable = x_name,
    Partial_Cor = ct$estimate,
    P_Value = ct$p.value,
    N = nrow(tmp_dat),
    stringsAsFactors = FALSE
  )
}

if (length(results_list_cont) > 0) {
  res_cont <- do.call(rbind, results_list_cont)
  res_cont$FDR <- p.adjust(res_cont$P_Value, method = "BH")
  res_cont <- res_cont[order(res_cont$P_Value), ]
  
  write.csv(res_cont, "results/Adjustable_Factors_Continuous_PartialCorr.csv", row.names = FALSE)
  cat("Saved Continuous Results.\n")
}

library(VennDiagram)
library(grid)

final_res <- if(exists("res_bin")) res_bin else NULL
final_res_cont <- if(exists("res_cont")) res_cont else NULL
out_dir <- "results"

cat("\nStarting visualization...\n")

plot_volcano_enhanced <- function(data, x_col, p_col, fdr_col, title, x_label, fc_cutoff = 0, fdr_cutoff = 0.05) {
  
  plot_dat <- data.table(data)
  
  plot_dat[, nlog10FDR := -log10(get(fdr_col))]
  plot_dat[, P_Val := get(p_col)]
  plot_dat[, X_Val := get(x_col)]
  
  plot_dat <- plot_dat %>% 
    mutate(trend = case_when(
      X_Val > fc_cutoff & get(fdr_col) < fdr_cutoff ~ "up",
      X_Val < -fc_cutoff & get(fdr_col) < fdr_cutoff ~ "down",
      TRUE ~ "non-sig"
    ))
  
  p <- ggplot(plot_dat, aes(x = X_Val, y = nlog10FDR, color = trend)) +
    geom_point(alpha = 0.7, aes(size = nlog10FDR)) +
    scale_color_manual(values = c("up" = "#C6295C", "down" = "#2C6DB2", "non-sig" = "grey")) +
    scale_size_continuous(range = c(1, 3)) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = plot_dat %>% 
                      filter(trend != "non-sig") %>%
                      arrange(P_Val) %>%
                      slice_head(n = 15),
                    aes(label = Variable),
                    size = 3.5, segment.size = 0.2, max.overlaps = 20,
                    box.padding = 0.5) +
    guides(size = guide_legend(order = 1)) +
    labs(title = title,
         x = x_label,
         y = "-log10(FDR)", color = "Risk Trend") +
    theme_test() +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_markdown(),
          plot.title = element_text(hjust = 0.5)) +
    xlim(min(plot_dat$X_Val, na.rm = TRUE) - 0.5, max(plot_dat$X_Val, na.rm = TRUE) + 0.5)
  
  return(p)
}

p1 <- plot_volcano_enhanced(final_res, 
                            x_col = "Beta", 
                            p_col = "P_Value", 
                            fdr_col = "FDR", 
                            title = "Binary (Logistic)", 
                            x_label = "Beta (log OR)")
print(p1)


p2 <- plot_volcano_enhanced(final_res_cont, 
                            x_col = "Partial_Cor", 
                            p_col = "P_Value", 
                            fdr_col = "FDR", 
                            title = "Continuous (Partial Corr)", 
                            x_label = "Partial Correlation (r)")
print(p2)

cat("Plotting Venn Diagram...\n")

sig_binary <- if(!is.null(final_res)) final_res$Variable[final_res$FDR < 0.05] else character(0)
sig_cont   <- if(!is.null(final_res_cont)) final_res_cont$Variable[final_res_cont$FDR < 0.05] else character(0)

venn_list <- list(
  "Binary (Logistic)" = sig_binary,
  "Continuous (Partial Corr)" = sig_cont
)

if (length(sig_binary) > 0 || length(sig_cont) > 0) {
  
  venn.plot <- venn.diagram(
    x = venn_list,
    filename = NULL,
    fill = c("#C6295C", "#2C6DB2"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.2,
    main = "Overlap of Significant Factors (FDR < 0.05)",
    margin = 0.05
  )
  
  pdf(file.path(out_dir, "Venn_Significant_Factors.pdf"))
  grid.draw(venn.plot)
  dev.off()
  
  overlap_prots <- intersect(sig_binary, sig_cont)
  if(length(overlap_prots) > 0) {
    write.csv(data.frame(Variable = overlap_prots), 
              file.path(out_dir, "Common_Significant_Factors.csv"), 
              row.names = FALSE)
    cat("Found", length(overlap_prots), "common significant factors.\n")
  } else {
    cat("No common significant factors found.\n")
  }
  
} else {
  cat("No significant factors found in either analysis, skipping Venn diagram.\n")
}

cat("Visualization complete. Check 'results' folder.\n")
