#' splitinput
#'
#' \code{splitinput} Splits input based on keepcol specified, yielding csv files each with at least the mininum
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
#' splitinput(df)
#'
#' # Specifying the name, directory and minimum row size
#' splitinput(df, fname = "syngrowth", fdir = tempdir(), min_nrow = 5000)
#'
#' # Specifying a different subject ID column
#' colnames(df)[colnames(df) == "subjid"] <- "sub_id"
#' splitinput(df, keepcol = "sub_id")
#' }
splitinput <- function(df,
                       fname = deparse(substitute(df)),
                       fdir = "",
                       min_nrow = 10000,
                       keepcol = "subjid") {
    # first, check if the given directory exists
    if (fdir != "" & is.character(fdir) & !dir.exists(fdir)) {
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
        utils::write.csv(
          split_df,
          file = file.path(fdir, paste(fname, fname_counter_str, "csv", sep = ".")),
          row.names = FALSE
        )

        split_df <- data.frame() # reset split_df
        fname_counter <- fname_counter + 1
      } else if (name == utils::tail(split_sample_names, 1)) {
        # for last part, just write
        fname_counter_str <- sprintf("%05d", fname_counter) # pad 0s
        utils::write.csv(
          split_df,
          file = file.path(fdir, paste(fname, fname_counter_str, "csv", sep = ".")),
          row.names = FALSE,
          na = ""
        )
      }
    }

    return(fname_counter)
  }


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
  input_table <- data.table::data.table(input_data)
  # replace targetcol variables with targetm where sourcecol = sourcem
  input_table[input_table[[sourcecol]] == sourcem, targetcol] <- targetm
  # replace targetcol variables with targetf where sourcecol = sourcef
  input_table[input_table[[sourcecol]] == sourcef, targetcol] <- targetf

  # return table
  return(input_table)
}


#' longwide
#'
#' \code{longwide} transforms data from long to wide format. Ideal for transforming output from growthcleanr::cleangrowth() into a format suitable for growthcleanr::ext_bmiz().
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
#' long_df <- longwide(df)
#'
#' # Include all inclusion types
#' long_df <- longwide(df, include_all = TRUE)
#'
#' # Specify all inclusion codes
#' long_df <- longwide(df, inclusion_types = c("Include", "Exclude-Carried-Forward"))
longwide <- function(long_df,
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
