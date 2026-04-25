library(stringr)

dat <- read.csv("20251211_future_data.csv")
dat_sub <- dat[, 3:ncol(dat)]
sorted_names <- sort(colnames(dat_sub))
dat_sub <- dat_sub[, sorted_names]

parse_formats <- c("%Y-%m-%d", "%y-%m-%d", "%Y/%m/%d", "%d/%m/%Y", "%m/%d/%Y")

merge_keys <- sub("\\.[^.]+$", "", colnames(dat_sub))
unique_keys <- unique(merge_keys)
merged_cols <- vector("list", length(unique_keys))
names(merged_cols) <- unique_keys
total_keys <- length(unique_keys)

for (i in seq_along(unique_keys)) {
  key <- unique_keys[i]
  cols_to_merge <- which(merge_keys == key)
  
  cat(sprintf("Processing %s (%d of %d): merging %d column(s)...\n", 
              key, i, total_keys, length(cols_to_merge)))
  
  if (length(cols_to_merge) == 1) {
    x <- dat_sub[[cols_to_merge]]
    merged_cols[[key]] <- as.Date(as.character(x), tryFormats = parse_formats)
  } else {
    date_matrix <- as.data.frame(dat_sub[, cols_to_merge])
    date_matrix[] <- lapply(date_matrix, function(x) as.Date(as.character(x), tryFormats = parse_formats))
    nums <- lapply(date_matrix, as.numeric)
    res_num <- do.call(pmin, c(nums, na.rm = TRUE))
    merged_cols[[key]] <- as.Date(res_num, origin = "1970-01-01")
  }
}

merged_df <- as.data.frame(merged_cols)
dat_final <- cbind(dat[, 1:2], merged_df)

convert_to_date <- function(x) {
  as.Date(as.character(x), tryFormats = parse_formats)
}

date_cols <- names(dat_final)[3:ncol(dat_final)]
dat_final[date_cols] <- lapply(dat_final[date_cols], convert_to_date)

counts <- sapply(dat_final[date_cols], function(x) sum(!is.na(x)))
counts_df <- data.frame(column = names(counts), count = as.integer(counts))

write.csv(dat_final,'20251211_future_data_combination.csv', row.names = FALSE)
