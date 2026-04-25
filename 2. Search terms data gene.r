library(dplyr)

dat <- read.csv("Diagnosis.csv", stringsAsFactors = FALSE)
gene <- read.csv("GENE.csv", stringsAsFactors = FALSE)

valid_id <- gene %>%
  filter(!is.na(`Clonal.haematopoiesis.of.indeterminate.potential..CHIP..number.of.variants`),
         `Genomic.ancestry` == "European ancestry (EUR)") %>%
  pull(Participant.ID) %>%          
  unique()                          
dat <- dat %>% filter(Participant.ID %in% valid_id)

participant_ids <- dat$Participant.ID
diagnoses <- dat$Diagnoses...ICD10
date_columns <- dat[, grep("Date.of.first.in.patient.diagnosis...ICD10...Array.", names(dat), fixed = TRUE)]

diagnoses <- as.character(dat$Diagnoses...ICD10)
cat("Splitting diagnoses for", length(diagnoses), "participants...\n")

diag_tokens_list <- strsplit(diagnoses, "\\|")
code_list <- lapply(diag_tokens_list, function(x) sub("^(\\S+).*", "\\1", x))
cat("Extracted", sum(lengths(code_list)), "diagnosis tokens; preparing indices...\n")

row_index <- rep(seq_len(nrow(dat)), lengths(code_list))
pos_index <- sequence(lengths(code_list))

date_matrix <- as.matrix(date_columns)
valid <- pos_index <= ncol(date_matrix)
date_values <- rep(NA_character_, length(pos_index))
date_values[valid] <- as.character(date_matrix[cbind(row_index[valid], pos_index[valid])])
cat("Mapped date values for", sum(valid), "tokens.\n")

long_df <- data.frame(
  Participant_ID = participant_ids[row_index],
  code = unlist(code_list),
  date = date_values,
  stringsAsFactors = FALSE
)

long_df$date[long_df$date == ""] <- NA
long_df$code <- trimws(long_df$code)
long_df <- long_df[!is.na(long_df$code) & long_df$code != "", ]

manual_terms <- NULL
if (!is.null(manual_terms) && length(manual_terms) > 0) {
  cat("Filtering to manual terms:", paste(manual_terms, collapse=", "), "\n")
  long_df <- long_df[long_df$code %in% manual_terms, ]
  for (term in manual_terms) {
    cnt <- sum(long_df$code == term & !is.na(long_df$date))
    cat("Matched term", term, "with", cnt, "non-NA dates.\n")
  }
} else {
  cat("No manual terms provided; using all observed codes (", length(unique(long_df$code)), ").\n", sep = "")
}

if (nrow(long_df) > 0) {
  agg <- aggregate(date ~ Participant_ID + code, data = long_df, FUN = function(x) {
    idx <- which(!is.na(x))
    if (length(idx) > 0) x[idx[1]] else NA_character_
  })
  cat("Aggregated to", nrow(agg), "Participant_ID+code pairs.\n")
  wide <- reshape(agg, v.names = "date", idvar = "Participant_ID", timevar = "code", direction = "wide")
  names(wide) <- sub("^date\\.", "", names(wide))
} else {
  cat("No matches after filtering; creating empty wide table.\n")
  wide <- data.frame(Participant_ID = participant_ids, stringsAsFactors = FALSE)
}

idx <- match(participant_ids, wide$Participant_ID)
wide <- wide[idx, , drop = FALSE]
wide$Participant_ID <- participant_ids

if (!is.null(manual_terms) && length(manual_terms) > 0) {
  existing <- setdiff(names(wide), "Participant_ID")
  missing <- setdiff(manual_terms, existing)
  if (length(missing) > 0) {
    for (m in missing) wide[[m]] <- NA_character_
  }
  wide <- wide[, c("Participant_ID", manual_terms), drop = FALSE]
}

cat("Writing output to search_terms_data_gene.csv with", max(ncol(wide)-1, 0), "code columns.\n")
write.csv(wide, "search_terms_data_gene.csv", row.names = FALSE)
