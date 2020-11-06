#' Function for temporary duplicates (step 5):
#' 5.  Temporary duplicates: I use duplicates to refer to more than more than one recorded value for a parameter on the same day,
#'     and we need to select which one to include in our analysis. The overall strategy will be to select a measurement using a simple
#'     strategy that will be used temporarily, and select permanently in a later step after we have a somewhat cleaner dataset that
#'     can help us identify the best duplicate.
#' a.  For subjects/parameters with duplicates: Determine median_tbc*sd for both parameters: the median tbc*sd for each subject
#'     and parameter including only non-duplicate values with exc_*==0. The median of the same parameter as the duplicate will
#'     be referred to as median_tbc*sd, the median of the other parameter will be referred to as median_tbcOsd.
#' b.	For each subject/parameter with duplicates and at least one value for the subject/parameter on a day with no duplicates,
#'     select the value closest to the median_tbc*sd for temporary inclusion, and assign all other duplicates exc_*=2.
#'  i.	For each subject/parameter with duplicates and no values for the subject/parameter on a day with no duplicates, select the
#'     value closest to the median_tbcOsd for temporary inclusion, and assign all other duplicates exc_*=2.
#'     If median_tbcOsd is missing because there are no values for the other parameter, randomly choose one duplicate value for
#'     each subject/parameter/age to keep as exc_*=0 and replace exc_*=2 for all other duplicates for that subject/parameter/age.
#'
#' @import data.table
#' @importFrom stats median
#'
#' @keywords internal
#' @noRd
temporary_duplicates <- function(df) {
  # ==== Dealing with "undefined global functions or variables" ==== #
  ## Only for variable which couldn't be quoted everywhere
  duplicates.this.day <- tbc.sd <- median.sd <- delta.median.sd <- subjid <- param <- NULL
  # ==== Dealing with "undefined global functions or variables" ==== #

  # add subjid and param if needed (may be missing depending on where this is called from)
  if (is.null(df[["subjid"]])) {
    df[, "subjid" := NA]
  }
  if (is.null(df[["param"]])) {
    df[, "param" := NA]
  }
  # only operate on valid rows (but include rows that may have previously been flagged as a "temporary duplicate")
  valid.rows <- valid(df, include.temporary.duplicates = TRUE)
  # make a small copy of df with only the fields we need for efficiency
  df <- df[j = .(tbc.sd), keyby = c("subjid", "param", "agedays", "index")]
  # initialize some useful fields
  df[, `:=`(
    "median.sd" = NA_real_,
    "delta.median.sd" = NA_real_,
    "duplicates.this.day" = FALSE,
    "duplicate" = FALSE
  )]
  # determine on which days there are duplicate measurements (more than 1 valid measurement on that age day)
  df[valid.rows, "duplicates.this.day" := (.N > 1), by = c("subjid", "param", "agedays")]
  # calculate median of measurements on days where there is no duplicate
  df[valid.rows & !duplicates.this.day, "median.sd" := stats::median(tbc.sd), by = c("subjid", "param")]
  # distribute median to other valid or potential duplicate rows for same parameter for that subject
  df[valid.rows, "median.sd" := sort(median.sd)[1], by = c("subjid", "param")]
  # take care of subject/parameters with more than one day with a valid observation
  # determine the absolute difference between the measurements sd score and the median for that parameter for each child
  df[valid.rows & duplicates.this.day, "delta.median.sd" := abs(tbc.sd - median.sd), by = c("subjid", "param")]
  # identify subjects that have duplicates on all days of observation for that parameter (i.e. delta.median.sd undefined)
  subj.all.dups <- df[valid.rows & duplicates.this.day & is.na(delta.median.sd), unique(subjid)]
  df[
    valid.rows & subjid %in% subj.all.dups,
    "delta.median.sd" := (function(subj.df) {
      # iterate over parameters where delta.median.sd is not yet defined
      # pass 1: take median from other parameter(s)
      for (p in subj.df[is.na(delta.median.sd) & duplicates.this.day, unique(param)]) {
        median.other.sd <- subj.df[param != p & !duplicates.this.day, stats::median(tbc.sd)]
        subj.df[param == p, "delta.median.sd" := abs(tbc.sd - median.other.sd)]
      }
      return(subj.df$delta.median.sd)
    })(copy(.SD)),
    by = "subjid"
  ]
  # Final pass: take median as zero (i.e. if no measurements from a different parameter)
  # NOTE -- this is not exactly the same as taking a random parameter
  df[valid.rows & duplicates.this.day & is.na(delta.median.sd), "delta.median.sd" := abs(tbc.sd)]
  # flag any duplicate value on the same day that is not the minimum distance from the median sd score
  # NOTE: in the case of exact ducplicates, "which.min" will pick the first
  df[
    valid.rows & duplicates.this.day,
    "duplicate" := seq_along(delta.median.sd) != which.min(delta.median.sd),
    by = c("subjid", "param", "agedays")
  ]
  # return the duplicated valid rows (i.e. the ones that should be temporarily excluded)
  return(df$duplicate & valid.rows)
}
