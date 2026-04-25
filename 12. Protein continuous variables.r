library(data.table)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(ggrepel)
library(ggtext)

out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir)

dat1_path <- "20251211_chip.csv"
if (!file.exists(dat1_path)) stop("20251211_chip.csv not found")
dat1 <- fread(dat1_path)

prot_path <- "Proteins.csv"
if (!file.exists(prot_path)) stop("Proteins.csv not found")
prot <- fread(prot_path, check.names = FALSE)

cols_needed <- c(
  "Participant ID",
  "Age at recruitment",
  "Sex",
  "Smoking status | Instance 0",
  "Alcohol drinker status | Instance 0",
  "Body mass index (BMI) | Instance 0",
  "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants"
)

available_cols <- intersect(cols_needed, names(dat1))
if (length(available_cols) < length(cols_needed)) {
  warning("Some specified columns are missing in dat1: ", 
          paste(setdiff(cols_needed, names(dat1)), collapse = ", "))
}
df_cov <- dat1[, ..available_cols]

setnames(df_cov, 
         old = c("Participant ID", 
                 "Age at recruitment", 
                 "Sex", 
                 "Smoking status | Instance 0", 
                 "Alcohol drinker status | Instance 0", 
                 "Body mass index (BMI) | Instance 0", 
                 "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants"),
         new = c("ID", "Age", "Sex", "Smoke", "Alcohol", "BMI", "CHIP_Raw"),
         skip_absent = TRUE)

df_cov[, CHIP_Status := ifelse(CHIP_Raw > 0, 1, 0)]

if (is.character(df_cov$Sex)) {
  df_cov[, Sex := ifelse(Sex == "Female", 0, 1)]
}

recode_four <- function(x) {
  x[is.na(x)] <- "Unknown"
  ifelse(x %in% c("Current", "Never", "Previous"), x, "Others")
}

df_cov[, Smoke := factor(recode_four(Smoke), levels = c("Never", "Previous", "Current", "Others", "Unknown"))]
df_cov[, Alcohol := factor(recode_four(Alcohol), levels = c("Never", "Previous", "Current", "Others", "Unknown"))]

if (any(is.na(df_cov$BMI))) {
  df_cov[is.na(BMI), BMI := mean(BMI, na.rm = TRUE)]
}

if ("Participant ID" %in% names(prot)) setnames(prot, "Participant ID", "ID")

prot_cols_orig <- setdiff(names(prot), "ID")
prot_names_clean <- sapply(strsplit(prot_cols_orig, ";"), `[`, 1)

setnames(prot, old = prot_cols_orig, new = prot_names_clean)

unique_prot_cols <- unique(prot_names_clean)
prot <- prot[, c("ID", unique_prot_cols), with = FALSE]

merged_dat <- merge(df_cov, prot, by = "ID")

cat("Data loaded and merged.\n")
cat("Participants: ", nrow(merged_dat), "\n")
cat("Proteins to test: ", length(unique_prot_cols), "\n")

results_list <- list()

covars_formula <- "Age + Sex + factor(Smoke) + factor(Alcohol) + BMI"

cat("Starting regression analysis...\n")

for (pname in unique_prot_cols) {
  idx <- which(unique_prot_cols == pname)
  cat(sprintf("[%d/%d] Processing Continuous Model for: %s\n", idx, length(unique_prot_cols), pname))

  idx <- which(unique_prot_cols == pname)
  cat(sprintf("[%d/%d] Processing Binary Model for: %s\n", idx, length(unique_prot_cols), pname))

  
  if (!pname %in% names(merged_dat)) next
  
  tmp_dat <- merged_dat[, c("CHIP_Status", pname, "Age", "Sex", "Smoke", "Alcohol", "BMI"), with = FALSE]
  
  tmp_dat <- na.omit(tmp_dat)
  
  if (nrow(tmp_dat) < 50 || sum(tmp_dat$CHIP_Status == 1) < 5 || sum(tmp_dat$CHIP_Status == 0) < 5) next
  
  fml <- as.formula(paste0("CHIP_Status ~ `", pname, "` + ", covars_formula))
  
  fit <- try(glm(fml, data = tmp_dat, family = binomial), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    coefs <- summary(fit)$coefficients
    
    row_idx <- grep(paste0("^`?", pname, "`?$"), rownames(coefs))
    
    if (length(row_idx) > 0) {
      est <- coefs[row_idx, 1]
      se  <- coefs[row_idx, 2]
      pval <- coefs[row_idx, 4]
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Protein = pname,
        Beta = est,
        OR = exp(est),
        SE = se,
        P_Value = pval,
        N = nrow(tmp_dat),
        N_Cases = sum(tmp_dat$CHIP_Status == 1),
        stringsAsFactors = FALSE
      )
    }
  }
}

cat("\nAnalysis complete.\n")

if (length(results_list) > 0) {
  final_res <- do.call(rbind, results_list)
  
  final_res <- final_res[order(final_res$P_Value), ]
  final_res$FDR <- p.adjust(final_res$P_Value, method = "BH")
  final_res$Bonferroni <- p.adjust(final_res$P_Value, method = "bonferroni")
  
  write.csv(final_res, file.path(out_dir, "CHIP_Binary_Protein_Association_Results.csv"), row.names = FALSE)
  
  sig_res <- final_res[final_res$FDR < 0.05, ]
  if (nrow(sig_res) > 0) {
    write.csv(sig_res, file.path(out_dir, "CHIP_Binary_Protein_Association_Significant.csv"), row.names = FALSE)
  }
  
  cat("Binary analysis results saved to:", out_dir, "\n")
  
} else {
  cat("No valid binary results generated.\n")
}

results_list_cont <- list()

cat("Starting partial correlation analysis (Continuous CHIP)...\n")

for (pname in unique_prot_cols) {
  
  if (!pname %in% names(merged_dat)) next
  
  tmp_dat <- merged_dat[, c("CHIP_Raw", pname, "Age", "Sex", "Smoke", "Alcohol", "BMI"), with = FALSE]
  
  tmp_dat <- na.omit(tmp_dat)
  
  if (nrow(tmp_dat) < 50) next
  
  fml_chip <- as.formula(paste0("CHIP_Raw ~ ", covars_formula))
  fit_chip <- lm(fml_chip, data = tmp_dat)
  res_chip <- residuals(fit_chip)
  
  fml_prot <- as.formula(paste0("`", pname, "` ~ ", covars_formula))
  fit_prot <- lm(fml_prot, data = tmp_dat)
  res_prot <- residuals(fit_prot)
  
  ct <- cor.test(res_chip, res_prot, method = "pearson")
  
  results_list_cont[[length(results_list_cont) + 1]] <- data.frame(
    Protein = pname,
    Partial_Cor = ct$estimate,
    P_Value = ct$p.value,
    N = nrow(tmp_dat),
    stringsAsFactors = FALSE
  )
}

if (length(results_list_cont) > 0) {
  final_res_cont <- do.call(rbind, results_list_cont)
  
  final_res_cont <- final_res_cont[order(final_res_cont$P_Value), ]
  final_res_cont$FDR <- p.adjust(final_res_cont$P_Value, method = "BH")
  final_res_cont$Bonferroni <- p.adjust(final_res_cont$P_Value, method = "bonferroni")
  
  write.csv(final_res_cont, file.path(out_dir, "CHIP_Continuous_Protein_Association_Results.csv"), row.names = FALSE)
  
  sig_res_cont <- final_res_cont[final_res_cont$FDR < 0.05, ]
  if (nrow(sig_res_cont) > 0) {
    write.csv(sig_res_cont, file.path(out_dir, "CHIP_Continuous_Protein_Association_Significant.csv"), row.names = FALSE)
  }
  
  cat("Continuous analysis results saved to:", out_dir, "\n")
  
} else {
  cat("No valid continuous results generated.\n")
}

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
                    aes(label = Protein), 
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
                            title = "Binary", 
                            x_label = "Beta (log OR)")
p1

p2 <- plot_volcano_enhanced(final_res_cont, 
                            x_col = "Partial_Cor", 
                            p_col = "P_Value", 
                            fdr_col = "FDR", 
                            title = "Continuous (Partial Corr)", 
                            x_label = "Partial Correlation (r)")
p2

cat("Plotting Venn Diagram...\n")

sig_binary <- if(exists("final_res")) final_res$Protein[final_res$FDR < 0.05] else character(0)
sig_cont   <- if(exists("final_res_cont")) final_res_cont$Protein[final_res_cont$FDR < 0.05] else character(0)

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
    main = "Overlap of Significant Proteins (FDR < 0.05)",
    margin = 0.05
  )
  
  pdf(file.path(out_dir, "Venn_Significant_Proteins.pdf"))
  grid.draw(venn.plot)
  dev.off()
  
  overlap_prots <- intersect(sig_binary, sig_cont)
  if(length(overlap_prots) > 0) {
    write.csv(data.frame(Protein = overlap_prots), 
              file.path(out_dir, "Common_Significant_Proteins.csv"), 
              row.names = FALSE)
    cat("Found", length(overlap_prots), "common significant proteins.\n")
  } else {
    cat("No common significant proteins found.\n")
  }
  
} else {
  cat("No significant proteins found in either analysis, skipping Venn diagram.\n")
}

cat("Visualization complete. Check 'results' folder.\n")


