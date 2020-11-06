#' long_wide
#'
#' \code{long_wide} transforms data from long to wide format. Ideal for transforming output from growthcleanr::cleangrowth() into a format suitable for growthcleanr::ext_bmiz().
#'
#' @param long_df A data frame to be transformed. Expects columns: id, subjid, sex, agedays, param, measurement, and clean_value.
#' @param id name of observation ID column
#' @param subjid name of subject ID column
#' @param sex name of sex descriptor column
#' @param agedays name of age (in days) descriptor column
#' @param param name of parameter column to identify each type of measurement
#' @param measurement name of measurement column containing the actual measurement data
#' @param include_all Determines whether the function keeps all exclusion codes. If TRUE, all exclusion types are kept and the inclusion_types argument is ignored. Defaults to FALSE.
#' @param inclusion_types Vector indicating which exclusion codes from the cleaning algorithm should be included in the data, given that include_all is FALSE. For all options, see growthcleanr::cleangrowth(). Defaults to c("Include").
#' @param na.rm exclude rows where bli was not computed due to missing values.
#'
#' @return Returns a data frame transformed from long to wide. Includes only values flagged with indicated inclusion types. Note that, for each subject, heights without corresponding weights for a given age (and vice versa) will be dropped.
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#' # Run on a small subset of given data
#' df <- as.data.frame(syngrowth)
#' df <- df[df$subjid %in% unique(df[, "subjid"])[1:5], ]
#' df <- cbind(df,
#'   "clean_value" = cleangrowth(
#'     df$subjid,
#'     df$param,
#'     df$agedays,
#'     df$sex,
#'     df$measurement
#'   )
#' )
#' # Convert to wide format
#' long_df <- long_wide(df)
#'
#' # Include all inclusion types
#' long_df <- long_wide(df, include_all = TRUE)
#'
#' # Specify all inclusion codes
#' long_df <- long_wide(df, inclusion_types = c("Include", "Exclude-Carried-Forward"))
long_wide <- function(long_df,
                  id = "id",
                  subjid = "subjid",
                  sex = "sex",
                  agedays = "agedays",
                  param = "param",
                  measurement = "measurement",
                  clean_value = "clean_value",
                  include_all = FALSE,
                  inclusion_types = c("Include"),
                  na.rm = TRUE) {
  # ==== Dealing with "undefined global functions or variables" ==== #
  ## Only for variable which couldn't be quoted everywhere
  agey <- ht <- wt <- bmi <- sex_recoded <- agem <- NULL
  # ==== Dealing with "undefined global functions or variables" ==== #

  if (!is.logical(include_all)) {
   stop(paste0(
      "include_all is not a logical of length 1. It is a ",
      typeof(include_all), " of length ", length(include_all), "."
    ))
  }

  cols <- c(
    id, subjid, sex, agedays,
    param, measurement, clean_value
  )
  if (!all(cols %in% names(long_df))) {
    stop(paste0(
      "Not all needed columns were present\n",
      paste(paste0("  + '", cols[!cols %in% names(long_df)], "'"), collapse = "\n")
    ))
  }

  # selects each column with specified / default variable name
  setDT(long_df)
  obs_df <- long_df[, .SD, .SDcols = cols]
  setnames(
    x = obs_df,
    old = cols,
    new = c(
      "id", "subjid", "sex", "agedays",
      "param", "measurement", "clean_value"
    )
  )

  # extract values flagged with indicated inclusion types:
  if (include_all) {
    obs_df <- obs_df
  } else {
    obs_df <- obs_df[clean_value %in% inclusion_types]
  }


  # only include observations at least 24 months old
  obs_df <- obs_df[agedays >= 730]

  # calculate age in years
  obs_df[, "agey" := round(agedays / 365.25, 4)]

  # calculate age in months
  obs_df[, "agem" := round(agey * 12, 4)]

  # recode sex to expected ext_bmiz() format
  obs_df <- recode_sex(
    input_data = obs_df,
    sourcecol = "sex",
    sourcem = "0",
    sourcef = "1",
    targetcol = "sex_recoded",
    targetm = 1L,
    targetf = 2L
  )

  clean_df <- obs_df[,
    `:=`("sex" = sex_recoded, param = as.character(param))
  ][, c("subjid", "id", "agey", "agem", "agedays", "sex", "param", "measurement")]

  # check for unique weight and height ids
  if (anyDuplicated(clean_df$id)) {
    stop("Duplicate IDs in data.")
  }

  # Best approach but because of duplicated measurements ..
  # param_separated <- dcast(
  #   data = clean_df,
  #   formula = ... ~ param,
  #   value.var = c("measurement", "id"),
  #   fun.aggregate = mean
  # )
  # setnames(
  #   x = param_separated,
  #   old = c("measurement_HEIGHTCM", "measurement_WEIGHTKG", "id_HEIGHTCM", "id_WEIGHTKG"),
  #   new = c("ht", "wt", "ht_id", "wt_id")
  # )
  # it will be as below
  param_separated <- Reduce(
    f = function(x, y) merge(x, y, by = c("subjid", "agey", "agem", "agedays", "sex")),
    x = lapply(
      X = split(clean_df, clean_df[["param"]]),
        FUN = function(dt) {
        short_name <- c("HEIGHTCM" = "ht", "WEIGHTKG" = "wt")[unique(dt[["param"]])]
        setnames(dt, c("id", "measurement"), c(paste0(short_name, "_id"), short_name))
        dt[, -c("param")]
      }
    )
  )

  param_separated[!is.na(ht) & !is.na(wt), "bmi" := wt / ((ht / 100)^2)]

  out <- param_separated[,
   c("subjid", "agey", "agem", "bmi", "sex", "wt", "wt_id", "ht", "ht_id", "agedays")
  ]

  if (na.rm) out <- out[!is.na(bmi)]

  return(out)
}
