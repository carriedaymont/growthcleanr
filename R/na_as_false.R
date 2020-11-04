#' Helper function to treat NA values as FALSE
#'
#' @keywords internal
#' @noRd
na_as_false <- function(v) {
  v[is.na(v)] <- FALSE
  v
}
