if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)

df <- readRDS("results/cox_results_data.rds")
row.names(df) <- NULL; df$Exposure <- "CHIP"

filtered <- df[df$Events >= 60, ]


filtered$Exposure <- as.character(filtered$Exposure)
filtered <- filtered[!is.na(filtered$Exposure) & nzchar(trimws(filtered$Exposure)), ]

filtered$fdr <- ave(filtered$p_value, filtered$Exposure, FUN = function(p) p.adjust(p, method = "BH"))
filtered$bonf <- ave(filtered$p_value, filtered$Exposure, FUN = function(p) p.adjust(p, method = "bonferroni"))

counts <- tapply(filtered$fdr < 0.05, filtered$Exposure, function(x) sum(x, na.rm = TRUE))
print(data.frame(Exposure = names(counts), fdr_lt_0_05 = as.integer(counts)))
write.csv(filtered, "results/cox_events_ge60_adjusted_subclass.csv", row.names = FALSE)

df2 <- filtered
if (ncol(df2) == 12) {
  names(df2) <- c("name", "Covariate", "effect", "CI_Lower", "CI_Upper", "pval", "N", "Events", "AIC", "Converged", "fdr_pval", "bonf_pval")
} else {
  warning(paste("Column count mismatch: Expected 12, got", ncol(df2), ". Keeping original names."))
}

index <- read.csv("coding19.tsv", sep = "\t")
if (!all(c("meaning", "source") %in% names(index))) stop("Index file missing 'meaning' or 'source' columns.")

df2 <- df2 %>%
  rowwise() %>%
  mutate(
    match_index = which(startsWith(index$meaning, name))[1],
    name = dplyr::if_else(!is.na(match_index), as.character(index$meaning[match_index]), as.character(name)),
    source = dplyr::if_else(!is.na(match_index), as.character(index$source[match_index]), NA_character_)
  ) %>%
  ungroup() %>%
  select(-match_index)

prepare_plot_data <- function(covariate) {
  df2 %>%
    filter(Covariate == covariate) %>%
    rowwise() %>%
    mutate(
      match_index = which(startsWith(index$meaning, name))[1],
      name = dplyr::if_else(!is.na(match_index), as.character(index$meaning[match_index]), as.character(name)),
      source = dplyr::if_else(!is.na(match_index), as.character(index$source[match_index]), NA_character_)
    ) %>%
    ungroup() %>%
    select(-match_index) %>%
    mutate(source = as.character(source)) %>%
    mutate(
      NegativeLogP = -log10(pval),
      level = case_when(
        fdr_pval >= 0.05 ~ "not_sig",
        effect > 1 ~ "enriched_up",
        effect < 1 ~ "depleted_down",
        TRUE ~ "not_sig"
      ),
      Group = case_when(
        source %in% c(
          "Chapter XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
          "Chapter XIX Injury, poisoning and certain other consequences of external causes",
          "Chapter XXII Codes for special purposes",
          "Chapter XX External causes of morbidity and mortality",
          "Chapter XXI Factors influencing health status and contact with health services"
        ) ~ "Others",
        TRUE ~ as.character(source)
      )
    ) %>%
    arrange(Group) %>%
    mutate(
      name = factor(name, levels = unique(name[order(source)])),
      Group = if(requireNamespace("forcats", quietly=TRUE)) forcats::fct_inorder(Group) else factor(Group, levels = unique(Group))
    )
}

getBreak <- function(x, y) {
  freq <- as.vector(table(y))
  half_freq <- freq %/% 2
  for (i in seq(2, length(freq))) {
    new_num <- freq[i] + freq[i-1]
    freq[i] <- new_num
  }
  pos <- freq - half_freq
  break_point <- as.vector(x[pos])
  return(break_point)
}

mycol <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#d62728", "#ff9896",
  "#c5b0d5", "#c49c94", "#f4a582", "#9c9ede", "#dbdb8d"
)

shape_values <- c(
  "not_sig" = 16,
  "enriched_up" = 24,
  "depleted_down" = 25
)

target_covariates <- unique(na.omit(df2$Covariate))
target_covariates <- head(target_covariates, 4)
if (length(target_covariates) == 0) stop("No valid Covariate values found in df2")

plot_list <- lapply(target_covariates, function(cov) {
  plot_data <- prepare_plot_data(cov)
  gtext <- plot_data %>% 
    filter(fdr_pval < 0.05) %>% 
    arrange(pval) %>% 
    head(10)
  
  p <- ggplot(plot_data, aes(x = name, y = NegativeLogP)) +
    geom_point(
      aes(color = Group, size = effect, shape = level),
      alpha = 0.7,
      stroke = 0.3,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = gtext,
      aes(label = name),
      color = "black",
      size = 3.0,
      box.padding = 0.15,
      point.padding = 0.05,
      segment.size = 0.2,
      segment.color = "grey60",
      max.overlaps = 30,
      min.segment.length = 0.1,
      force = 0.5,
      direction = "both",
      seed = 123,
      nudge_x = 0.05,
      nudge_y = 0.05
    ) +
    scale_shape_manual(
      values = shape_values,
      labels = c("Not Significant", "Enriched (HR > 1)", "Depleted (HR < 1)")
    ) +
    scale_size_continuous(
      range = c(0.5, 4),
      breaks = c(1, 2, 4)
    ) +
    scale_color_manual(values = mycol) +
    scale_x_discrete(
      breaks = if(nlevels(plot_data$Group) > 1) getBreak(plot_data$name, plot_data$Group) else NULL,
      labels = NULL
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "red",
      size = 0.3
    ) +
    labs(
      x = "",
      y = {
        idx <- match(cov, target_covariates)
        if (!is.na(idx) && idx %in% c(1, min(4, length(target_covariates)))) "-log10(p-value)" else ""
      }
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(color = "black", size = 0.3),
      axis.line = element_line(color = "black", size = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 9),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
  
  return(p)
})

dummy_data <- prepare_plot_data(target_covariates[1])

legend_plot <- ggplot(dummy_data, aes(x = name, y = NegativeLogP)) +
  geom_point(
    aes(color = Group, size = effect, shape = level),
    alpha = 0.7,
    stroke = 0.3
  ) +
  scale_shape_manual(
    values = shape_values,
    breaks = c("enriched_up", "depleted_down", "not_sig"),
    labels = c("Enriched", "Depleted", "Not Sig."),
    drop = FALSE,
    name = NULL
  ) +
  guides(shape = guide_legend(override.aes = list(fill = NA))) +
  scale_size_continuous(
    range = c(1, 3),
    breaks = c(1, 2, 4),
    name = "HR"
  ) +
  scale_color_manual(
    values = mycol,
    name = NULL,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme(
    legend.spacing = grid::unit(0.1, "cm"),
    legend.key.size = grid::unit(0.4, "cm"),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7),
    legend.box.margin = margin(0, 0, 0, 0)
  )

safe_get_legend <- function(p) {
  tryCatch(suppressWarnings(cowplot::get_legend(p)), error = function(e) NULL)
}

legend_chapter <- safe_get_legend(
  legend_plot + guides(shape = "none", size = "none")
)

l1 <- safe_get_legend(legend_plot + guides(color = "none", size = "none"))
l2 <- safe_get_legend(legend_plot + guides(color = "none", shape = "none"))

legend_other <- if (!is.null(l1) || !is.null(l2)) {
  cowplot::plot_grid(
    l1, l2,
    ncol = 1,
    align = "v",
    rel_heights = c(1, 0.8)
  )
} else {
  NULL
}

legend_combined <- if (!is.null(legend_chapter) || !is.null(legend_other)) {
  cowplot::plot_grid(
    legend_chapter,
    legend_other,
    nrow = 1,
    rel_widths = c(1.2, 0.8)
  )
} else {
  NULL
}

n_plots <- length(plot_list)
final_grid <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol = min(2, n_plots)
)

if (!is.null(legend_combined)) {
  final_plot <- cowplot::plot_grid(final_grid, legend_combined, ncol = 2, rel_widths = c(4, 1))
} else {
  final_plot <- final_grid
}

print(final_plot)

chapter_summary <- df2 %>%
  filter(!is.na(source) & source != "") %>%
  mutate(
    short_source = case_when(
      source == "Others" ~ "Others",
      grepl("^Chapter [IVXLCDM]+", source) ~ sub("^(Chapter [IVXLCDM]+).*", "\\1", source),
      TRUE ~ source
    )
  ) %>%
  group_by(short_source) %>%
  summarise(
    Total = n(),
    Sig = sum(fdr_pval < 0.05, na.rm = TRUE)
  ) %>%
  arrange(desc(Total))

chapter_summary$short_source <- factor(chapter_summary$short_source, levels = chapter_summary$short_source)

p_chapter <- ggplot(chapter_summary, aes(x = rev(short_source))) +
  geom_col(aes(y = Total), fill = "#2C6DB2", width = 0.7) +
  geom_col(aes(y = Sig), fill = "#C6295C", width = 0.7) +
  geom_text(
    aes(y = Total, label = paste0(Sig, "/", Total)),
    hjust = -0.1, 
    size = 3.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  coord_flip() +
  labs(
    x = NULL,
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9)
  )

print(p_chapter)

top20_df <- df2 %>%
  filter(fdr_pval < 0.05) %>%
  mutate(abs_log_hr = abs(log(effect))) %>%
  arrange(desc(abs_log_hr)) %>%
  head(20) %>%
  mutate(
    name = factor(name, levels = rev(name)),
    RiskGroup = ifelse(effect > 1, "High Risk (>1)", "Low Risk (<1)")
  )

p_forest <- ggplot(top20_df, aes(x = effect, y = name, color = RiskGroup)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("High Risk (>1)" = "#C6295C", "Low Risk (<1)" = "#2C6DB2")) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    color = "Effect Direction"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_forest)
