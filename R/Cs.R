#' Character strings from unquoted names
#'
#' Makes a vector of character strings from a list of valid S names
#'
#' @param ...	any number of names separated by commas
#'
#' @return character string vector
#'
#' @keywords internal
#' @noRd
Cs <- function(...) as.character(sys.call())[-1]
