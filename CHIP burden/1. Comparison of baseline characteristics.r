library(data.table)

dat1 <- fread("20251211_chip_gene_matrix.csv")
dat2 <- fread("regular indicators and modifiable behaviors.csv")
dat3 <- fread("cognitive.csv")
merged_dat <- merge(dat1, dat2, by = "Participant ID", all.x = TRUE)
merged_dat <- merge(merged_dat, dat3, by = "Participant ID", all.x = TRUE)

cols_to_keep <- c(1:3, 13, 68:145, 191:195)
cols_to_keep <- cols_to_keep[cols_to_keep <= ncol(merged_dat)]
merged_dat <- merged_dat[, ..cols_to_keep]

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

outcome_col <- names(merged_dat)[4]

max_col <- ncol(merged_dat)
pred_indices <- c(2:3, 5:max_col)
pred_indices <- pred_indices[pred_indices <= max_col]
predictor_cols <- names(merged_dat)[pred_indices]

results_list <- list()

y <- merged_dat[[outcome_col]]

for (x_name in predictor_cols) {
  x <- merged_dat[[x_name]]
  
  valid <- !is.na(x) & !is.na(y)
  x_clean <- x[valid]
  y_clean <- y[valid]

  if (length(x_clean) < 3 || length(unique(x_clean)) < 2) {
    next
  }
  
  res_item <- list(Variable = x_name, 
                   CHIP_Stat = NA_character_, 
                   Non_CHIP_Stat = NA_character_,
                   OR = NA_character_, 
                   CI_95 = NA_character_, 
                   P_Value = NA_real_,
                   Note = NA_character_)
  
  groups_list <- list(
    "CHIP" = which(y_clean >= 1),
    "Non_CHIP" = which(y_clean == 0)
  )
  
  for (grp_name in names(groups_list)) {
    idx <- groups_list[[grp_name]]
    col_name <- paste0(grp_name, "_Stat")
    
    if (length(idx) == 0) {
      res_item[[col_name]] <- "0 (0.0%)"
      next
    }
    
    x_sub <- x_clean[idx]
    
    if (is.numeric(x_clean) && !is.factor(x_clean)) {
      m_val <- mean(x_sub, na.rm = TRUE)
      sd_val <- sd(x_sub, na.rm = TRUE)
      stat_str <- sprintf("%.2f (%.2f)", m_val, sd_val)
    } else {
      tbl <- table(x_sub)
      props <- prop.table(tbl) * 100
      stat_parts <- paste0(names(tbl), ": ", as.integer(tbl), " (", sprintf("%.1f", props), "%)")
      stat_str <- paste(stat_parts, collapse = "; ")
    }
    res_item[[col_name]] <- stat_str
  }
  
  y_bin <- ifelse(y_clean >= 1, 1, 0)
  
  tryCatch({
    if (length(unique(y_bin)) < 2) {
      res_item$Note <- "Outcome has only 1 level"
    } else if (length(unique(x_clean)) < 2) {
      res_item$Note <- "Predictor has only 1 level"
    } else {
      
      if (is.numeric(x_clean) && !is.factor(x_clean) && length(unique(x_clean)) > 2) {
        m <- glm(y_bin ~ x_clean, family = binomial)
        s <- summary(m)
        
        est <- coef(m)[2]
        se <- s$coefficients[2, 2]
        p_val <- s$coefficients[2, 4]
        
        or_val <- exp(est)
        ci_lo <- exp(est - 1.96 * se)
        ci_hi <- exp(est + 1.96 * se)
        
        res_item$OR <- sprintf("%.2f", or_val)
        res_item$CI_95 <- sprintf("%.2f - %.2f", ci_lo, ci_hi)
        res_item$P_Value <- p_val
        res_item$Note <- "Continuous"
        
      } else {
        x_fact <- as.factor(x_clean)
        m <- glm(y_bin ~ x_fact, family = binomial)
        
        if (length(levels(x_fact)) == 2) {
          s <- summary(m)
          est <- coef(m)[2]
          se <- s$coefficients[2, 2]
          p_val <- s$coefficients[2, 4]
          
          or_val <- exp(est)
          ci_lo <- exp(est - 1.96 * se)
          ci_hi <- exp(est + 1.96 * se)
          
          res_item$OR <- sprintf("%.2f", or_val)
          res_item$CI_95 <- sprintf("%.2f - %.2f", ci_lo, ci_hi)
          res_item$P_Value <- p_val
          res_item$Note <- paste("Binary:", levels(x_fact)[2], "vs", levels(x_fact)[1])
          
        } else {
          m0 <- glm(y_bin ~ 1, family = binomial)
          lrt <- anova(m0, m, test = "LRT")
          p_val <- lrt[["Pr(>Chi)"]][2]
          
          res_item$OR <- "Ref (Mixed)"
          res_item$CI_95 <- "-"
          res_item$P_Value <- p_val
          res_item$Note <- "Categorical (Overall P)"
        }
      }
    }
  }, error = function(e) {
    res_item$Note <- paste("Error:", e$message)
  })
  
  results_list[[length(results_list) + 1]] <- res_item
}

results_df <- rbindlist(results_list)

if (nrow(results_df) > 0) {
  results_df[, FDR := p.adjust(P_Value, method = "BH")]
}

output_file <- "results/baseline_characteristics_logistic.csv"
if(!dir.exists("results")) dir.create("results")
fwrite(results_df, output_file)
cat("Baseline characteristics with Logistic Regression results saved to:", output_file, "\n")

if (requireNamespace("ggplot2", quietly = TRUE) && 
    requireNamespace("ggrepel", quietly = TRUE) && 
    requireNamespace("ggtext", quietly = TRUE) &&
    requireNamespace("dplyr", quietly = TRUE)) {
  
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(dplyr)
  
  plot_data <- results_df[OR != "Ref (Mixed)" & !is.na(P_Value) & !is.na(OR)]
  plot_data[, OR_num := as.numeric(OR)]
  plot_data[, log2FC := log2(OR_num)]
  plot_data[, nlog10FDR := -log10(FDR)]
  
  fc_cutoff <- 0 
  fdr_cutoff <- 0.05
  
  plot_data <- plot_data %>% 
    mutate(trend = case_when(
      log2FC > fc_cutoff & FDR < fdr_cutoff ~ "up",
      log2FC < -fc_cutoff & FDR < fdr_cutoff ~ "down",
      TRUE ~ "non-sig"
    ))
  
  p_vol <- ggplot(plot_data, aes(x = log2FC, y = nlog10FDR, color = trend)) +
    geom_point(alpha = 0.7, aes(size = nlog10FDR)) +
    scale_color_manual(values = c("up" = "#C6295C", "down" = "#2C6DB2", "non-sig" = "grey")) +
    scale_size_continuous(range = c(1, 3)) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = plot_data %>% 
                      filter(trend != "non-sig") %>%
                      arrange(P_Value) %>%
                      slice_head(n = 15),
                    aes(label = Variable),
                    size = 3.5, segment.size = 0, max.overlaps = 15,
                    box.padding = 0.5) +
    guides(size = guide_legend(order = 1)) +
    labs(x = "log2(Odds Ratio)",
         y = "-log10(FDR)", color = "Risk Trend") +
    theme_test() +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_markdown()) +
    xlim(-3, 3)
  
  print(p_vol)
  cat("Enhanced Volcano plot displayed.\n")
} else {
  cat("Please install 'ggplot2', 'ggrepel', 'ggtext', and 'dplyr' to view the Enhanced Volcano Plot.\n")
}

if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ggplot2)

chip_y <- merged_dat[[outcome_col]]

pred_num_cols <- predictor_cols[sapply(merged_dat[, ..predictor_cols], is.numeric)]
if (length(pred_num_cols) == 0) {
  stop("No numeric predictors available for correlation plot.")
}

corr_dt <- merged_dat[, c(outcome_col, pred_num_cols), with = FALSE]
cc <- stats::complete.cases(corr_dt)
corr_dt <- corr_dt[cc]

x_mat <- as.matrix(corr_dt[, ..pred_num_cols])
y_df <- data.frame(CHIP = as.numeric(corr_dt[[outcome_col]]))

pearson <- psych::corr.test(x_mat, y_df, method = "pearson", adjust = "none")

r <- data.table::as.data.table(pearson$r)
p <- data.table::as.data.table(pearson$p)

r[, env := rownames(pearson$r)]
p[, env := rownames(pearson$p)]

r_long <- data.table::melt(r, id.vars = "env", variable.name = "spe", value.name = "pearson_correlation")
p_long <- data.table::melt(p, id.vars = "env", variable.name = "spe", value.name = "p.value")

pearson_long <- merge(r_long, p_long, by = c("env", "spe"), all = TRUE)

pearson_long[, p_fdr := p.adjust(p.value, method = "BH")]
pearson_long[, sig := ""]
pearson_long[p_fdr < 0.001, sig := "***"]
pearson_long[p_fdr < 0.01 & p_fdr >= 0.001, sig := "**"]
pearson_long[p_fdr < 0.05 & p_fdr >= 0.01, sig := "*"]

output_corr_file <- "results/baseline_characteristics_correlation.csv"
if(!dir.exists("results")) dir.create("results")
data.table::fwrite(
  pearson_long[, .(
    Variable = env,
    Outcome = spe,
    Pearson_r = pearson_correlation,
    P_Value = p.value,
    FDR = p_fdr,
    Sig = sig
  )],
  output_corr_file
)
cat("Correlation results saved to:", output_corr_file, "\n")

pearson_long <- pearson_long[!is.na(p_fdr) & p_fdr < 0.05]
if (nrow(pearson_long) == 0) {
  stop("No correlations pass FDR < 0.05.")
}

pearson_long[, spe := factor(spe, levels = "CHIP")]

wrap_label <- function(x, width = 18) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

pearson_long <- pearson_long[order(-abs(pearson_correlation))]
top_n <- min(30, nrow(pearson_long))
pearson_plot <- pearson_long[seq_len(top_n)]

pearson_plot[, env_label := wrap_label(env, 18)]
pearson_plot[, label := paste0(env_label, "  ", sig)]
pearson_plot[, env_f := factor(env_label, levels = env_label)]

p_corr <- ggplot(pearson_plot, aes(x = env_f, y = spe, fill = pearson_correlation)) +
  geom_tile(aes(border = after_stat(fill)), color = "grey") +
  scale_fill_gradient2(
    low = "#2C6DB2",
    mid = "white",
    high = "#C6295C",
    limits = c(-1, 1)
  ) +
  coord_polar() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(25, 180, 25, 180)
  ) +
  geom_text(
    aes(label = label, y = 1.3),
    size = 3,
    angle = 0,
    lineheight = 0.9,
    check_overlap = TRUE
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0)))

print(p_corr)
