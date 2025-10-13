# Supporting pediatric growthcleanr functions
# Supporting functions for pediatric piece of algorithm

#' Helper function for cleanbatch to identify subset of observations that are either "included" or a "temporary extraneous"
#'
#' @keywords internal
#' @noRd
valid <- function(df,
                  include.temporary.extraneous = FALSE,
                  include.extraneous = FALSE,
                  include.carryforward = FALSE) {
  exclude <- if (is.data.frame(df))
    df$exclude
  else
    df
  return(
    exclude < 'Exclude'
    |
      include.temporary.extraneous &
      exclude == 'Exclude-Temporary-Extraneous-Same-Day'
    | include.extraneous & exclude == 'Exclude-Extraneous-Same-Day'
    |
      include.carryforward &
      exclude == 'Exclude-Carried-Forward'
  )
}

#' Helper function to treat NA values as FALSE
#'
#' @keywords internal
#' @noRd
na_as_false <- function(v) {
  v[is.na(v)] <- FALSE
  v
}

#' Function for temporary extraneous (step 5):
#' 5.  Temporary extraneous: I use extraneous to refer to more than more than one recorded value for a parameter on the same day,
#'     and we need to select which one to include in our analysis. The overall strategy will be to select a measurement using a simple
#'     strategy that will be used temporarily, and select permanently in a later step after we have a somewhat cleaner dataset that
#'     can help us identify the best extraneous.
#' a.  For subjects/parameters with extraneous: Determine median_tbc*sd for both parameters: the median tbc*sd for each subject
#'     and parameter including only non-extraneous values with exc_*==0. The median of the same parameter as the extraneous will
#'     be referred to as median_tbc*sd, the median of the other parameter will be referred to as median_tbcOsd.
#' b.	For each subject/parameter with extraneous and at least one value for the subject/parameter on a day with no extraneous,
#'     select the value closest to the median_tbc*sd for temporary inclusion, and assign all other extraneous exc_*=2.
#'  i.	For each subject/parameter with extraneous and no values for the subject/parameter on a day with no extraneous, select the
#'     value closest to the median_tbcOsd for temporary inclusion, and assign all other extraneous exc_*=2.
#'     If median_tbcOsd is missing because there are no values for the other parameter, randomly choose one extraneous value for
#'     each subject/parameter/age to keep as exc_*=0 and replace exc_*=2 for all other extraneous for that subject/parameter/age.
#'
#' @keywords internal
#' @noRd
temporary_extraneous <- function(df) {
  # avoid "no visible binding" warnings
  agedays <- delta.median.sd <- extraneous <- extraneous.this.day <- NULL
  index <- median.sd <- param <- subjid <- tbc.sd <- NULL

  # add subjid and param if needed (may be missing depending on where this is called from)
  if (is.null(df$subjid))
    df[, subjid := NA]
  if (is.null(df$param))
    df[, param := NA]
  # only operate on valid rows (but include rows that may have previously been flagged as a "temporary extraneous")
  valid.rows <- valid(df, include.temporary.extraneous = TRUE)
  # make a small copy of df with only the fields we need for efficiency
  df <- df[j = .(tbc.sd), keyby = .(subjid, param, agedays, index)]
  # initialize some useful fields
  df[, `:=`(
    median.sd = as.double(NA),
    delta.median.sd = as.double(NA),
    extraneous.this.day = FALSE,
    extraneous = FALSE
  )]
  # determine on which days there are extraneous measurements (more than 1 valid measurement on that age day)
  df[valid.rows, extraneous.this.day := (.N > 1), by = .(subjid, param, agedays)]
  # calculate median of measurements on days where there is no extraneous
  df[valid.rows &
       !extraneous.this.day, median.sd := median(tbc.sd), by = .(subjid, param)]
  # distribute median to other valid or potential extraneous rows for same parameter for that subject
  df[valid.rows, median.sd := sort(median.sd)[1], by = .(subjid, param)]
  # take care of subject/parameters with more than one day with a valid observation
  # determine the absolute difference between the measurements sd score and the median for that parameter for each child
  df[valid.rows &
       extraneous.this.day, delta.median.sd := abs(tbc.sd - median.sd), by = .(subjid, param)]
  # identify subjects that have duplicates/extraneous on all days of observation for that parameter (i.e. delta.median.sd undefined)
  subj.all.dups <- df[valid.rows &
                        extraneous.this.day &
                        is.na(delta.median.sd), unique(subjid)]
  df[valid.rows &
       subjid %in% subj.all.dups, delta.median.sd := (function(subj.df) {
         # iterate over parameters where delta.median.sd is not yet defined
         # pass 1: take median from other parameter(s)
         for (p in subj.df[is.na(delta.median.sd) &
                           extraneous.this.day, unique(param)]) {
           median.other.sd <- subj.df[param != p &
                                        !extraneous.this.day, median(tbc.sd)]
           subj.df[param == p, delta.median.sd := abs(tbc.sd - median.other.sd)]
         }
         return(subj.df$delta.median.sd)
       })(copy(.SD)), by = .(subjid)]
  # Final pass: take median as zero (i.e. if no measurements from a different parameter)
  # NOTE -- this is not exactly the same as taking a random parameter
  df[valid.rows &
       extraneous.this.day &
       is.na(delta.median.sd), delta.median.sd := abs(tbc.sd)]
  # flag any extraneous value on the same day that is not the minimum distance from the median sd score
  # NOTE: in the case of exact ducplicates, "which.min" will pick the first
  df[valid.rows &
       extraneous.this.day, extraneous := seq_along(delta.median.sd) != which.min(delta.median.sd), by =
       .(subjid, param, agedays)]
  # return the extraneous valid rows (i.e. the ones that should be temporarily excluded)
  return(df$extraneous & valid.rows)
}

#' Function to identify switches (step 6):
#'  7.  Identify switches (weight and height each recorded as the other parameter)
#' Utility function to find the value associated with the opposite parameter
#' For example, if param.1='WEIGHTKG' and param.2='HEIGHTCM', then a vector of height values observed on the same day as the weight values is returned and vice versa.
#'
#' @keywords internal
#' @noRd
swap_parameters <- function(param.1 = 'WEIGHTKG',
                            param.2 = 'HEIGHTCM',
                            field.name = 'tbc.sd',
                            df) {
  # avoid "no visible binding" warnings
  agedays <- agedays.other <- param <- param.other <- subjid <- subjid.other <- swap <- NULL

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
  setkey(valid.other, subjid.other, param.other, agedays.other)
  # clear swap field (should retain original type)
  df[, swap := NA]
  # use indexed table for speed -- assign all values in one statement
  df[valid.rows, swap := valid.other[list(subjid, param, agedays), swap]]
  return(df$swap)
}

#' Function to convert vector of ages, in days, into matrix
#'
#' @keywords internal
#' @noRd
as_matrix_delta <- function(agedays) {
  n <- length(agedays)
  delta <- abs(matrix(rep(agedays, n), n, byrow = TRUE) - agedays)

  return(delta)
}
