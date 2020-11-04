#' Function to identify switches (step 6):
#'  7.  Identify switches (weight and height each recorded as the other parameter)
#' Utility function to find the value associated with the opposite parameter
#' For example, if param.1='WEIGHTKG' and param.2='HEIGHTCM', then a vector of height values observed on the same day as the weight values is returned and vice versa.
#'
#' @import data.table
#'
#' @keywords internal
#' @noRd
swap_parameters <- function(param.1 = "WEIGHTKG",
                            param.2 = "HEIGHTCM",
                            field.name = "tbc.sd",
                            df) {
  valid.rows <- valid(df)
  # copy swap field to a new value for convenience in the code below
  df$swap <- df[, field.name, with = FALSE]
  # construct an indexed table with valid parameters for speed in swapping
  valid.other <- df[valid.rows, list(
    subjid.other = subjid,
    param.other = ifelse(
      param == param.1,
      param.2,
      ifelse(param == param.2, param.1, NA)
    ),
    agedays.other = agedays,
    swap
  )]
  setkeyv(valid.other, c("subjid.other", "param.other", "agedays.other"))
  # clear swap field (should retain original type)
  df[, swap := NA]
  # use indexed table for speed -- assign all values in one statement
  df[valid.rows, swap := valid.other[list(subjid, param, agedays), swap]]
  return(df$swap)
}
