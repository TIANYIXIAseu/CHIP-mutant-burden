infile  <- "20251211_chip.csv"
outfile <- "20251211_chip_gene_matrix.csv"

parse_gene_tokens <- function(v) {
  if (length(v) == 0 || is.na(v) || v == "") return(character(0))
  toks <- unlist(strsplit(v, "\\|", fixed = TRUE))
  genes <- sub("^\\s*([^:]+).*", "\\1", toks)
  genes <- trimws(genes)
  genes[genes != "" & !is.na(genes)]
}

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

keep_first_13 <- df[, 1:13, drop = FALSE]
gene_cols_idx <- intersect(14:20, seq_len(ncol(df)))

row_genes_list <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  g_i <- character(0)
  for (j in gene_cols_idx) {
    g_i <- c(g_i, parse_gene_tokens(df[i, j]))
  }
  row_genes_list[[i]] <- unique(g_i)
}

gene_set <- sort(unique(unlist(row_genes_list)))

if (length(gene_set) == 0) {
  gene_mat <- matrix(integer(0), nrow = nrow(df), ncol = 0)
  colnames(gene_mat) <- character(0)
} else {
  rows_fac  <- factor(rep(seq_len(nrow(df)), lengths(row_genes_list)),
                      levels = seq_len(nrow(df)))
  genes_fac <- factor(unlist(row_genes_list), levels = gene_set)

  tab <- xtabs(~ rows_fac + genes_fac)

  gene_mat <- (tab > 0) + 0L
  colnames(gene_mat) <- gene_set
}

if (ncol(gene_mat) > 0) {
  ord_cols <- order(colSums(gene_mat), decreasing = TRUE)
  gene_mat <- gene_mat[, ord_cols, drop = FALSE]
}
out_df <- cbind(keep_first_13, as.data.frame(gene_mat, check.names = FALSE))


if (ncol(gene_mat) > 0) {
  gc <- colSums(gene_mat)
  ord <- order(gc, decreasing = TRUE)
  gc_df <- data.frame(Gene = colnames(gene_mat)[ord], Count = as.integer(gc[ord]),
                      check.names = FALSE, stringsAsFactors = FALSE)
  print(gc_df)
} else {
  print(data.frame(Gene = character(0), Count = integer(0)))
}

write.csv(out_df, outfile, row.names = FALSE)
