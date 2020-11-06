#' Function to convert vector of ages, in days, into matrix
#'
#' @keywords internal
#' @noRd
as_matrix_delta <- function(agedays) {
  n <- length(agedays)
  abs(matrix(rep(agedays, n), n, byrow = TRUE) - agedays)
}
