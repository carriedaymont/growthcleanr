#' split_input
#'
#' \code{split_input} Splits input based on keepcol specified, yielding csv files each with at least the mininum
#' number of rows that are written and saved separately (except for the last split file written, which may be
#' smaller). Allows splitting input data while ensuring all records for each individual subject will stay together
#' in one file. Pads split filenames with zeros out to five digits for consistency, assuming < 100,000 file count
#' result.
#'
#' @param df data frame to split
#' @param fname new name for each of the split files to start with
#' @param fdir directory to put each of the split files (default working directory)
#' @param min_row minimum number of rows for each split file (default 10000)
#' @param keepcol the column name (default "subjid") to use to keep records with the same values together in the same single split file
#'
#' @return the count number referring to the last split file written
#'
#' @export
#' @examples
#' \donttest{
#' # Run on given data
#' df <- as.data.frame(syngrowth)
#'
#' # Run with all defaults
#' split_input(df)
#'
#' # Specifying the name, directory and minimum row size
#' split_input(df, fname = "syngrowth", fdir = tempdir(), min_nrow = 5000)
#'
#' # Specifying a different subject ID column
#' colnames(df)[colnames(df) == "subjid"] <- "sub_id"
#' split_input(df, keepcol = "sub_id")
#' }
split_input <- function(df,
                       fname = deparse(substitute(df)),
                       fdir = ".",
                       min_nrow = 10000,
                       keepcol = "subjid") {
  # first, check if the given directory exists
  fdir <- normalizePath(fdir, mustWork = TRUE)
  if (!dir.exists(fdir)) {
    stop("invalid directory")
  }

  fname_counter <- 0
  # row_count <- 0
  split_df <- data.frame()

  # split the data frame by the grouping that user specifies
  split_sample <- split(df, df[[keepcol]])

  # grab the individual grouping names
  split_sample_names <- names(split_sample)

  for (name in split_sample_names) {
    # append the rows from name, store new total row_count for current split file
    split_df <- rbind(split_df, split_sample[[name]])
    current_nrow <- nrow(split_df)

    # check if updated row count will exceed min row count,
    # if min nrow is exceeded, then write.csv the current split file and clear the split dataframe starter (start from 0)
    if (current_nrow > min_nrow) {
      fname_counter_str <- sprintf("%05d", fname_counter) # pad 0s
      fwrite(
        x = split_df,
        file = file.path(fdir, paste(fname, fname_counter_str, "csv", sep = ".")),
        row.names = FALSE
      )

      split_df <- data.frame() # reset split_df
      fname_counter <- fname_counter + 1
    } else if (name == split_sample_names[nrow(split_sample_names), ]) {
      # for last part, just write
      fname_counter_str <- sprintf("%05d", fname_counter) # pad 0s
      fwrite(
        x = split_df,
        file = file.path(fdir, paste(fname, fname_counter_str, "csv", sep = ".")),
        row.names = FALSE,
        na = ""
      )
    }
  }

  return(fname_counter)
}
