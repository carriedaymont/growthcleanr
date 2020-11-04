#' Exponentially Weighted Moving Average (EWMA)
#'
#' \code{ewma} calculates the exponentially weighted moving average (EWMA) for a set of numeric observations over time.
#'
#' @param agedays Vector of age in days for each z score (potentially transformed to adjust weighting).
#' @param z Input vector of numeric z-score data.
#' @param ewma.exp Exponent to use for weighting.
#' @param ewma.adjacent Specify whether EWMA values excluding adjacent measurements should be calculated.
#'     Defaults to TRUE.
#'
#' @return Data frame with 3 variables:
#' * The first variable (ewma.all) contains the EWMA at observation time
#'   excluding only the actual observation for that time point.
#' * The second variable (ewma.before) contains the EWMA for each observation excluding both the actual observation
#'   and the immediate prior observation.
#' * The third variable (ewma.after) contains the EWMA for each observation excluding both the actual observation
#'   and the subsequent observation.
#' @md
#'
#' @export
#' @examples
#' # Run on 1 subject, 1 type of parameter
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid == df_stats$subjid[1] &
#'   df_stats$param == "HEIGHTCM", ]
#'
#' # Get the uncentered z-scores
#' measurement_to_z <- read_anthro(cdc.only = TRUE)
#' sd <- measurement_to_z(
#'   df_stats$param,
#'   df_stats$agedays,
#'   df_stats$sex,
#'   df_stats$measurement,
#'   TRUE
#' )
#'
#' # Calculate exponentially weighted moving average
#' e_df <- ewma(df_stats$agedays, sd, ewma.exp = -1.5)
ewma <- function(agedays, z, ewma.exp, ewma.adjacent = TRUE) {
  # 6.  EWMA calculation description: Most of the next steps will involve calculating
  #     the exponentially weighted moving average for each subject and parameter.
  #     I will describe how to calculate EWMASDs, and will describe
  #     how it needs to be varied in subsequent steps.
  # a.	The overall goal of the EWMASD calculation is to identify the difference
  #     between the SD-score and what we might predict that DS-score should be, in order to
  #     determine whether it should be excluded.
  # b.	Only nonmissing SD-scores for a parameter that are not designated for exclusion
  #     are included in the following calculations.
  # c.	For each SD-score SDi and associated agedaysi calculate the following for every
  #     other z-score (SDj…SDn) and associated agedays (agedaysj…agedaysn)  for the
  #     same subject and parameter
  #   i.	ΔAgej=agedaysj-agedaysi
  #   ii.	EWMAZ=SDi=[Σj→n(SDj*((5+ΔAgej)^-1.5))]/[ Σj→n((5+ΔAgej)^-1.5)]
  #   iii.	For most EWMASD calculations, there are 3 EWMASDs that need to be calculated.
  #     I will note if not all of these need to be done for a given step.
  #     1.	EWMASDall calculated as above
  #     2.	EWMAZbef calculated excluding the SD-score just before
  #     the SD-score of interest (sorted by agedays). For the first observation for a parameter for a
  #         subject, this should be identical to EWMASDall rather than missing.
  #     3.	EWMAZaft calculated excluding the z-score just after
  #     the SD-score of interest (sorted by agedays). For the lastobservation for a parameter for a subject,
  #         this should be identical to EWMASDall rather than missing.
  #   iv.	For each of the three EWMASDs, calculate the dewma_*=SD-EWMASD
  # d.	EWMASDs and ΔEWMASDs will change if a value is excluded or manipulated using one of
  #     the methods below, therefore EWMASDs and ΔEWMASDs be recalculated for each
  #     step where they are needed.
  # e.	For these calculations, use variables that allow for precise storage of numbers
  #     (in Stata this is called 'double') because otherwise rounding errors can cause
  #     problems in a few circumstances

  n <- length(agedays)
  # initialize response variables
  ewma.all <- ewma.before <- ewma.after <- vector("numeric", 0)
  if (n > 0) {
    # organize into data frame and sort into order of increasing age,
    # but retain original sort order information in index
    if (!all(agedays == cummax(agedays))) {
      warning("EWMA ordering is not sorted; double check")
    } # add in a check to make sure the inputs are already sorted (they should be)
    index <- order(agedays)


    # calculate matrix of differences in age, and add 5 to each delta per Daymont algorithm
    delta <- as_matrix_delta(agedays)
    delta <- ifelse(delta == 0, 0, (delta + 5)^ewma.exp)

    # calculate EWMAs, and return in order of original data
    ewma.all[index] <- delta %*% z / apply(delta, 1, sum)

    if (ewma.adjacent) {
      if (n > 2) {
        delta2 <- delta
        delta2[col(delta2) == row(delta2) - 1] <- 0
        ewma.before[index] <- delta2 %*% z / apply(delta2, 1, sum)
        delta3 <- delta
        delta3[col(delta3) == row(delta3) + 1] <- 0
        ewma.after[index] <- delta3 %*% z / apply(delta3, 1, sum)
      } else {
        ewma.before <- ewma.after <- ewma.all
      }
    }
  }
  # return all 3 EWMAs as a data frame
  if (ewma.adjacent) {
    data.frame(ewma.all, ewma.before, ewma.after)
  } else {
    data.frame(ewma.all)
  }
}
