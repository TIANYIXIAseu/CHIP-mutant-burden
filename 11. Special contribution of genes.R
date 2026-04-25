library(data.table)
library(ggplot2)
library(dplyr)

file_path <- "20251211_chip_gene_matrix.csv" 
if (!file.exists(file_path)) {
  stop(paste("File not found:", file_path))
}

dat <- fread(file_path)

gene_cols <- c("DNMT3A","TET2","ASXL1","PPM1D","SRSF2", "TP53", "JAK2","SF3B1",
               "GNB1","NF1","PRPF40B","GNAS","IDH2","CBL","U2AF1", "ASXL2", 
               "PRPF8","CUX1","SETDB1","BRCC3","BCORL1","CREBBP","KRAS","BCOR",
               "CEBPA", "EP300", "PHIP","RAD21","STAG2","NRAS","ETV6","KDM6A",
               "MPL","PTPN11","EZH2", "PDS5B", "RUNX1","SMC3","PHF6","ZRSR2",
               "KIT","SUZ12","WT1","BRAF","ETNK1", "CTCF", "IDH1","SETD2",
               "SETBP1","CBLB","CSF3R","GATA2","IKZF1","SMC1A")

available_genes <- intersect(gene_cols, names(dat))
missing_genes <- setdiff(gene_cols, names(dat))

if (length(missing_genes) > 0) {
  cat("Warning: The following genes were not found in the dataset:\n")
  print(missing_genes)
}

if (length(available_genes) == 0) {
  stop("None of the specified genes were found in the dataset.")
}

gene_counts <- dat[, lapply(.SD, function(x) sum(x > 0, na.rm = TRUE)), .SDcols = available_genes]

plot_dat <- data.table(
  Gene = names(gene_counts),
  Count = as.numeric(gene_counts[1, ])
)

plot_dat <- plot_dat[order(-Count)]

plot_dat_filtered <- plot_dat[Count > 0]

cat("Top 10 mutated genes:\n")
print(head(plot_dat_filtered, 10))

plot_dat_filtered <- plot_dat[Count > 20]

top_n <- 5 
plot_dat_filtered[, Color_Group := ifelse(rank(-Count) <= top_n, "Top", "Other")]

anno_df <- data.frame(
  Gene = c("DNMT3A", "TET2", "ASXL1", "PPM1D", "SRSF2", "TP53", "JAK2", "SF3B1"),
  Class = c(12, 0, 1, 0, 14, 0, 24, 5),
  Subclass = c(6, 1, 1, 0, 12, 0, 24, 4)
)

if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

anno_df_t <- t(anno_df[, -1])
colnames(anno_df_t) <- anno_df$Gene

tt <- ttheme_default(core = list(fg_params = list(cex = 0.9)),
                     colhead = list(fg_params = list(cex = 0.9)),
                     rowhead = list(fg_params = list(cex = 0.9)))
tbl <- tableGrob(anno_df_t, theme = tt)

p <- ggplot(plot_dat_filtered, aes(x = reorder(Gene, -Count), y = Count, fill = Color_Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Top" = "#C6295C", "Other" = "#2C6DB2")) +
  labs(
    x = "Gene",
    y = "Number of Participants with Mutation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    plot.margin = margin(10, 10, 10, 30),
    legend.position = "none"
  ) +
  geom_text(aes(label = Count), vjust = -0.3, size = 3) +
  annotation_custom(tbl, 
                    xmin = length(unique(plot_dat_filtered$Gene)) * 0.4, 
                    xmax = length(unique(plot_dat_filtered$Gene)), 
                    ymin = max(plot_dat_filtered$Count) * 0.5, 
                    ymax = max(plot_dat_filtered$Count))

print(p)

venn_file <- "results/cox_results_chip_gene_combination.csv"
if(!file.exists(venn_file)) {
  stop(paste("Venn input file not found:", venn_file))
}

venn_dat <- fread(venn_file)

exclude_exp <- "Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants"

venn_dat_sig <- venn_dat[FDR_BH < 0.05 & Exposure != exclude_exp]

venn_list <- split(venn_dat_sig$Outcome, venn_dat_sig$Exposure)

venn_list <- lapply(venn_list, unique)

cat("Exposures found:", paste(names(venn_list), collapse = ", "), "\n")

if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
library(UpSetR)
library(grid)

all_outcomes <- unique(unlist(venn_list))
all_genes <- names(venn_list)

upset_data <- data.frame(row.names = all_outcomes)
for (gene in all_genes) {
  upset_data[[gene]] <- as.integer(rownames(upset_data) %in% venn_list[[gene]])
}

upset_data$pattern <- apply(upset_data[, all_genes], 1, function(x) {
  genes_in <- names(x)[x == 1]
  if (length(genes_in) == 0) return(NA)
  paste(sort(genes_in), collapse = "-")
})

upset_data_plot <- upset_data[!is.na(upset_data$pattern), ]

outcomes_list <- split(rownames(upset_data_plot), upset_data_plot$pattern)

counts <- sapply(outcomes_list, length)
sorted_patterns <- names(counts)[order(-counts)]

intersections_list <- lapply(sorted_patterns, function(pat) {
  unlist(strsplit(pat, "-"))
})

outcome_labels <- sapply(sorted_patterns, function(pat) {
  outcomes <- outcomes_list[[pat]]
  paste(outcomes, collapse = "\n")
})

index_file <- "coding19.tsv"
upset_queries <- list()

if (file.exists(index_file)) {
  index <- read.csv(index_file, sep = "\t")
  
  mycol <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#e6ab02", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#fdcdac", "#ff9896",
    "#c5b0d5", "#c49c94", "#f4a582", "#9c9ede", "#dbdb8d"
  )
  
  get_group_for_single_outcome <- function(outcome_name) {
    outcome_name <- trimws(outcome_name)
    match_index <- match(outcome_name, trimws(index$meaning))
    
    if (is.na(match_index)) {
      safe_outcome <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", outcome_name)
      matches <- grep(safe_outcome, index$meaning, ignore.case = TRUE)
      if (length(matches) > 0) match_index <- matches[1]
    }
    
    if (is.na(match_index)) {
      matches <- which(startsWith(index$meaning, outcome_name))
      if (length(matches) > 0) match_index <- matches[1]
    }
    
    if (is.na(match_index)) return("Unknown")
    
    source_val <- as.character(index$source[match_index])
    
    if (source_val %in% c(
      "Chapter XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
      "Chapter XIX Injury, poisoning and certain other consequences of external causes",
      "Chapter XXII Codes for special purposes",
      "Chapter XX External causes of morbidity and mortality",
      "Chapter XXI Factors influencing health status and contact with health services"
    )) {
      return("Others")
    } else {
      return(source_val)
    }
  }
  
  upset_data$Group <- sapply(rownames(upset_data), get_group_for_single_outcome)
  
  elements <- function(row, column, value) {
    return(row[[column]] == value)
  }
  
  upset_data$Group <- factor(upset_data$Group, 
                             levels = names(sort(table(upset_data$Group), decreasing = TRUE)))
  
  unique_groups <- levels(upset_data$Group)
  
  mycol <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#e6ab02", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#fdcdac", "#ff9896",
    "#c5b0d5", "#c49c94", "#f4a582", "#9c9ede", "#dbdb8d"
  )
  
  group_colors <- setNames(mycol[seq_along(unique_groups)], unique_groups)
  
  if ("Unknown" %in% unique_groups) group_colors["Unknown"] <- "black"
  
  cat("Group distribution:\n")
  print(table(upset_data$Group))
  cat("\nGroup colors assigned:\n")
  print(group_colors)
  
  upset_queries <- lapply(levels(upset_data$Group), function(grp) {
    list(
      query = elements,
      params = list("Group", grp),
      color = group_colors[[grp]],
      active = TRUE,
      query.name = grp
    )
  })
  
} else {
  warning("coding19.tsv not found. Stacked bars will not be generated.")
}

grid.newpage()

p_upset <- grid::grid.grabExpr({
  print(upset(upset_data, 
        sets = all_genes, 
        intersections = intersections_list,
        queries = upset_queries,
        query.legend = "none",
        main.bar.color = "gray90",
        sets.bar.color = "#C6295C", 
        set_size.show = TRUE,
        sets.x.label = "Set Size",
        order.by = "freq",
        text.scale = c(1.4, 1.4, 1, 1, 1.5, 1.2)
  ))
}, wrap.grobs = TRUE)

grid.draw(p_upset)
grid.force()

sorted_counts <- counts[sorted_patterns]
target_seq <- as.character(sorted_counts)
n_target <- length(target_seq)

grob_list <- grid.ls(print = FALSE)
text_indices <- grep("text", grob_list$name)
text_names <- grob_list$name[text_indices]

all_labels <- sapply(text_names, function(nm) {
  g <- grid.get(nm)
  if (length(g$label) > 0) as.character(g$label)[1] else ""
})

match_idx <- integer(0)
n_all <- length(all_labels)

if (n_all >= n_target) {
  for (i in 1:(n_all - n_target + 1)) {
    sub_seq <- all_labels[i:(i + n_target - 1)]
    if (all(sub_seq == target_seq)) {
      match_idx <- c(match_idx, i)
    }
  }
}

final_names <- NULL
if (length(match_idx) > 0) {
  best_i <- match_idx[length(match_idx)]
  
  for (i in match_idx) {
    curr_names <- text_names[i:(i + n_target - 1)]
    if (any(grepl("geom_text", curr_names))) {
      best_i <- i
      break
    }
  }
  final_names <- text_names[best_i:(best_i + n_target - 1)]
  
  cat("Successfully identified bar labels. Replacing with outcomes...\n")
  
  for (i in seq_along(final_names)) {
    t_name <- final_names[i]
    grid.edit(t_name, label = outcome_labels[i], 
              gp = gpar(fontsize = 7, col = "black"), 
              rot = 90, hjust = 0)
  }
  
} else {
  cat("Warning: Could not identify bar labels by sequence matching.\n")
  cat("Target sequence:", paste(head(target_seq), collapse=", "), "...\n")
}

if (exists("group_colors") && length(group_colors) > 0) {
  
  group_counts <- table(upset_data$Group)
  
  ordered_groups <- levels(upset_data$Group)
  
  legend_df <- data.frame(
    Group = ordered_groups,
    Color = as.character(group_colors[ordered_groups]),
    Count = as.numeric(group_counts[ordered_groups]),
    stringsAsFactors = FALSE
  )
  
  cat("Legend order and colors:\n")
  print(legend_df)
  
  pushViewport(viewport(x = 0.6, y = 0.95, width = 0.35, height = 0.5, just = c("left", "top")))
  
  grid.text("Disease chapters:", x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
  
  start_y <- 0.9
  gap_y <- 0.08
  
  for (i in 1:nrow(legend_df)) {
    cur_y <- start_y - (i-1) * gap_y
    
    if (cur_y < 0) break
    
    grid.rect(x = 0, y = cur_y, width = 0.05, height = 0.05, just = c("left", "center"),
              gp = gpar(fill = legend_df$Color[i], col = NA))
    
    grid.text(legend_df$Group[i], x = 0.07, y = cur_y, just = c("left", "center"),
              gp = gpar(fontsize = 8))
  }
  
  popViewport()
}

cat("UpSet plot displayed with outcomes on bars.\n")
  
