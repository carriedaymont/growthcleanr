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
#'
#' @return Returns a data frame transformed from long to wide. Includes only values flagged with indicated inclusion types. Note that, for each subject, heights without corresponding weights for a given age (and vice versa) will be dropped.
#'
#' @export
#' @rawNamespace import(tidyr, except = extract)
#' @rawNamespace import(dplyr, except = c(last, first, summarize, src, between))
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
                  inclusion_types = c("Include")) {
  # selects each column with specified / default variable name
  obs_df <- dplyr::select(
    long_df,
    id, subjid, sex, agedays,
    param, measurement, clean_value
  )

  # if all columns could be found,
  # 7 columns will be present in the correct order. Thus, rename
  if (ncol(obs_df) == 7) {
    names(obs_df) <- c(
      "id",
      "subjid",
      "sex",
      "agedays",
      "param",
      "measurement",
      "clean_value"
    )
  } else {
    # catch error if any variables were not found
    stop("not all needed columns were present")
  }

  # extract values flagged with indicated inclusion types:
  if (include_all == TRUE) {
    obs_df <- obs_df
  } else if (include_all == FALSE) {
    obs_df <- obs_df[obs_df$clean_value %in% inclusion_types, ]
  } else {
    stop(paste0(
      "include_all is not a logical of length 1. It is a ",
      typeof(include_all), " of length ", length(include_all)
    ))
  }


  # only include observations at least 24 months old
  obs_df <- obs_df[obs_df$agedays >= 730, ]

  # calculate age in years
  obs_df$agey <- round(obs_df$agedays / 365.25, 4)

  # calculate age in months
  obs_df$agem <- round((obs_df$agey * 12), 4)

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

  clean_df <- obs_df %>%
    dplyr::mutate(sex = sex_recoded) %>%
    dplyr::mutate(param = as.character(param)) %>%
    dplyr::select(subjid, id, agey, agem, agedays, sex, param, measurement)


  # check for unique weight and height ids
  if (any(duplicated(clean_df$id))) {
    stop("duplicate IDs in long_df")
  }

  # separate heights and weights using unique ids
  param_separated <- tidyr::pivot_wider(clean_df, names_from = param, values_from = measurement)

  # extract heights and weights attached to ids
  height <- param_separated %>%
    dplyr::filter(!is.na(HEIGHTCM)) %>%
    dplyr::filter(is.na(WEIGHTKG)) %>%
    dplyr::mutate(ht_id = id) %>%
    dplyr::select(-id) %>%
    dplyr::select(-WEIGHTKG)

  weight <- param_separated %>%
    dplyr::filter(is.na(HEIGHTCM)) %>%
    dplyr::filter(!is.na(WEIGHTKG)) %>%
    dplyr::mutate(wt_id = id) %>%
    dplyr::select(-id) %>%
    dplyr::select(-HEIGHTCM)


  # join based on subjid, age, and sex
  wide_df <- merge(height, weight, by = c("subjid", "agey", "agem", "agedays", "sex")) %>%
    dplyr::mutate(bmi = WEIGHTKG / ((HEIGHTCM * .01)^2)) %>% # calculate bmi
    dplyr::mutate(wt = WEIGHTKG, ht = HEIGHTCM) %>% # rename height and weight
    dplyr::select(subjid, agey, agem, bmi, sex, wt, wt_id, ht, ht_id, agedays)

  return(wide_df)
}


