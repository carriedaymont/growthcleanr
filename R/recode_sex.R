#' recode_sex
#'
#' \code{recode_sex} recodes a binary sex variable for a given source column in a data frame or data table.
#' Useful in transforming output from growthcleanr::cleangrowth() into a format suitable for growthcleanr::ext_bmiz().
#'
#' @param input_data a data frame or data table to be transformed. Expects a source column containing a binary sex variable.
#' @param sourcecol name of sex descriptor column. Defaults to "sex"
#' @param sourcem variable indicating "male" sex in input data. Defaults to "0"
#' @param sourcef variable indicating "female" sex in input data. Defaults to "1"
#' @param targetcol desired name of recoded sex descriptor column. Defaults to "sex_recoded"
#' @param targetm desired name of recoded sex variable indicating "male" sex in output data. Defaults to 1
#' @param targetf desired name of recoded sex variable indicating "female" sex in output data. Defaults to 2
#'
#' @return Returns a data table with recoded sex variables.
#'
#' @export
#' @examples
#' # Run on given data
#' df <- as.data.frame(syngrowth)
#'
#' # Run with all defaults
#' df_r <- recode_sex(df)
#'
#' # Specify different targets
#' df_rt <- recode_sex(df, targetcol = "sexr", targetm = "Male", targetf = "Female")
#'
#' # Specify different inputs
#' df_ri <- recode_sex(df_rt, sourcecol = "sexr", sourcem = "Male", sourcef = "Female")
recode_sex <- function(input_data,
                       sourcecol = "sex",
                       sourcem = "0",
                       sourcef = "1",
                       targetcol = "sex_recoded",
                       targetm = 1L,
                       targetf = 2L) {
  # cast to DT for faster processing
  input_table <- data.table(input_data)
  # replace targetcol variables with targetm where sourcecol = sourcem
  input_table[input_table[[sourcecol]] == sourcem, targetcol] <- targetm
  # replace targetcol variables with targetf where sourcecol = sourcef
  input_table[input_table[[sourcecol]] == sourcef, targetcol] <- targetf

  # return table
  return(input_table)
}
