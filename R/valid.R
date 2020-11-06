#' Helper function for cleanbatch to identify subset of observations that are either "included" or a "temporary duplicate"
#'
#' @keywords internal
#' @noRd
valid <- function(df,
                  include.temporary.duplicates = FALSE,
                  include.duplicates = FALSE,
                  include.carryforward = FALSE) {
  exclude <- if (is.data.frame(df)) df$exclude else df
  return(
    exclude < "Exclude" |
      include.temporary.duplicates & exclude == "Exclude-Temporary-Duplicate" |
      include.duplicates & exclude == "Exclude-Duplicate" |
      include.carryforward & exclude == "Exclude-Carried-Forward"
  )
}