library(data.table)

infile  <- "/public/home/b_tyxia/CHIP/search_terms_data_gene.csv"
chipfile <- "/public/home/b_tyxia/CHIP/20251211_chip.csv"
past_out <- "/public/home/b_tyxia/CHIP/20251211_past_data.csv"
future_out <- "/public/home/b_tyxia/CHIP/20251211_future_data.csv"
summary_out <- "/public/home/b_tyxia/CHIP/20251211_non_na_summary.csv"

header <- fread(infile, nrows = 0)
col_names <- names(header)

dat2 <- fread(chipfile)
dat2[, enrollment_date := as.IDate(`Date of attending assessment centre | Instance 0`, format = "%Y-%m-%d")]
setkey(dat2, `Participant ID`)

chunk_size <- 5000L
processed <- 0L
first_write <- TRUE

past_counts <- integer(0)
future_counts <- integer(0)

repeat {
  start <- 2L + processed
  end   <- start + chunk_size - 1L
  cmd <- sprintf("sed -n '%d,%dp' %s", start, end, infile)
  chunk <- try(fread(cmd = cmd, header = FALSE, col.names = col_names, showProgress = FALSE), silent = TRUE)
  if (inherits(chunk, "try-error") || nrow(chunk) == 0) break
  
  chunk <- merge(chunk,
                 dat2[, .(`Participant ID`, enrollment_date)],
                 by.x = "Participant_ID", by.y = "Participant ID",
                 all.x = TRUE)
  chunk <- chunk[!is.na(enrollment_date)]
  
  diag_cols <- setdiff(names(chunk), c("Participant_ID", "enrollment_date"))
  chunk[, (diag_cols) := lapply(.SD, function(x) as.IDate(x, format = "%Y/%m/%d")), .SDcols = diag_cols]
  
  past_chunk <- chunk[, c("Participant_ID", "enrollment_date", diag_cols), with = FALSE]
  future_chunk <- copy(past_chunk)
  
  past_chunk[, (diag_cols) := Map(function(x) fifelse(!is.na(x) & x <= enrollment_date, x, as.IDate(NA)), .SD), .SDcols = diag_cols]
  future_chunk[, (diag_cols) := Map(function(x) fifelse(!is.na(x) & x >  enrollment_date, x, as.IDate(NA)), .SD), .SDcols = diag_cols]
  
  pc <- vapply(diag_cols, function(col) sum(!is.na(past_chunk[[col]])), integer(1))
  fc <- vapply(diag_cols, function(col) sum(!is.na(future_chunk[[col]])), integer(1))
  
  if (length(past_counts) == 0) {
    past_counts <- pc; names(past_counts) <- diag_cols
    future_counts <- fc; names(future_counts) <- diag_cols
  } else {
    past_counts[diag_cols] <- past_counts[diag_cols] + pc
    future_counts[diag_cols] <- future_counts[diag_cols] + fc
  }
  
  fwrite(past_chunk, past_out, append = !first_write)
  fwrite(future_chunk, future_out, append = !first_write)
  first_write <- FALSE
  
  processed <- processed + nrow(chunk)
  cat("Processed rows:", processed, "\n")
}

non_na_summary <- data.table(
  Diagnosis = names(past_counts),
  Future_Non_NA_Count = as.integer(future_counts[names(past_counts)]),
  Past_Non_NA_Count   = as.integer(past_counts)
)[order(-Future_Non_NA_Count)]

fwrite(non_na_summary, summary_out)
cat("Done. Past:", past_out, "Future:", future_out, "Summary:", summary_out, "\n")
