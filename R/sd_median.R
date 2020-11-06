#' Calculate median SD score by age for each parameter.
#'
#' @param param Vector identifying each measurement, may be 'WEIGHTKG', or 'HEIGHTCM'.
#' @param agedays Numeric vector containing the age in days at each measurement.
#' @param sex Vector identifying the gender of the subject, may be 'M', 'm', or 0 for males, vs. 'F', 'f' or 1 for females.
#' @param sd.orig Vector of previously calculated standard deviation (SD) scores for each measurement before re-centering.
#'
#' @return Table of data with median SD-scores per day of life by gender and parameter.
#'
#' @export
#' @import data.table
#' @importFrom stats approx median
#'
#' @examples
#' # Run on 1 subject
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid == df_stats$subjid[1], ]
#'
#' # Get the original standard deviations
#' measurement_to_z <- read_anthro(cdc.only = TRUE)
#' sd.orig <- measurement_to_z(
#'   df_stats$param,
#'   df_stats$agedays,
#'   df_stats$sex,
#'   df_stats$measurement,
#'   TRUE
#' )
#'
#' # Calculate median standard deviations
#' sd.m <- sd_median(
#'   df_stats$param,
#'   df_stats$sex,
#'   df_stats$agedays,
#'   sd.orig
#' )
sd_median <- function(param, sex, agedays, sd.orig) {
  # ==== Dealing with "undefined global functions or variables" ==== #
  ## Only for variable which couldn't be quoted everywhere
  ageyears <- NULL
  # ==== Dealing with "undefined global functions or variables" ==== #

  # 3.  SD-score recentering: Because the basis of the method is comparing SD-scores over time,
  #     we need to account for the fact that
  #     the mean SD-score for the population changes with age.
  # a.  Determine the median cdc*sd for each parameter by year of age (with sexes combined): median*sd.
  # b.	The median*sd should be considered to apply to midyear-age,
  #     defined as the age in days with the same value as the integer
  #     portion of (365.25*year + 365.25/2).
  # c.	Linearly interpolate median*sd for each parameter between each midyear-age,
  #     naming the interpolated values rc*sd.
  # d.	For ages below the first midyear-age, let rc*sd equal the median*sd for the earliest year.
  #     For ages above the last midyear_age, let rc*sd equal the median*sd for the last year.
  # e.	Subtract rcsd_* from SDorig to create the recentered SD-score.
  #     This recentered SD-score, labeled tbc*sd
  #     (stands for "to be cleaned") will be used for most of the rest of the analyses.
  # f.	In future steps I will sometimes refer to measprev and measnext which
  #     refer to the previous or next wt or ht measurement
  #     for which exc_*==0 for the subject and parameter, when the data are sorted by subject,
  #     parameter, and agedays. SDprev and SDnext refer to the tbc*sd of the previous or next measurement.
  dt <- data.table(param, sex, agedays, ageyears = floor(agedays / 365.25), sd.orig)
  setkeyv(dt, c("param", "sex", "agedays"))
  # determine ages (in days) we need to cover from min to max age in years
  agedays.to.cover <- dt[, (floor(min(ageyears) * 365.25):floor(max(ageyears + 1) * 365.25))]
  # group all measurements above age 19 together (copied from CD stata code e-mailed 4/3/15)
  dt[ageyears > 19, "ageyears" := 19]
  dt.median <- dt[
    !is.na(sd.orig),
    list(sex = c(0, 1), sd.median = rep(stats::median(sd.orig), 2)),
    by = c("param", "ageyears")
  ]
  dt.median <- dt.median[,
    list(
      agedays = agedays.to.cover,
      sd.median = stats::approx(
        floor((ageyears + 0.5) * 365.25),
        sd.median,
        xout = agedays.to.cover,
        rule = 2
      )$y
    ),
    by = c("param", "sex")
  ]
  setkeyv(dt.median, c("param", "sex", "agedays"))
  return(dt.median)
}