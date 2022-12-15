# Main pediatric growthcleanr function -- cleanbatch
# internal supporting functions for pediatrics can be found in: pediatric_support.R

#' Function to clean data (optionally in batches):
#' 4.  Dataset split (optional)
#' a.	In Stata, many of the rest of the steps require you to repeat steps as many times as there are observations for the subject with the highest number of
#'     observations. Therefore, in order to speed things up I separated the dataset at this point into three groups: short (subjects with <=10 total
#'     observations), medium (subjects with >10 and <=30 observations) and long (subjects with >30 observations). I performed all the rest of the steps
#'     separately for each of the 3 datasets and then re-combined them at the end. Note that all observations for a subject stay together. This step is
#'     definitely not necessary if it will not speed things up in whatever software program you're using.
#'
#' NOTE: to use mutliple processes in R we will process patients in batches using the plyr ddply function
#' set up function to process one patient.  Function to parallelize batches is below.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats median
#' @noRd
cleanbatch <- function(data.df,
                       log.path,
                       quietly,
                       parallel,
                       measurement.to.z,
                       ewma.fields,
                       ewma.exp,
                       recover.unit.error,
                       include.carryforward,
                       sd.extreme,
                       z.extreme,
                       exclude.levels,
                       tanner.ht.vel,
                       who.ht.vel,
                       lt3.exclude.mode,
                       error.load.threshold,
                       error.load.mincount) {
  # avoid "no visible binding" warnings
  abs.2ndlast.sd <- abs.tbc.sd <- abs.tbc.sd.next <- abs.tbc.sd.prev <- abssum2 <- NULL
  aft.g.befp1 <- agedays <- agedays.other <- bef.g.aftm1 <- delta <- NULL
  delta.agedays.next <- delta.agedays.prev <- delta.next.ht <- delta.next.sd <- NULL
  delta.prev.ht <- delta.prev.sd <- dewma.after <- dewma.after.prev <- NULL
  dewma.all <- dewma.before <- dewma.before.next <- dnext.sd <- dnext.sd.minus <- NULL
  dnext.sd.plus <- dprev.sd <- dprev.sd.minus <- dprev.sd.plus <- dup <- NULL
  dup.ratio <- ewma.after <- ewma.all <- ewma.before <- exclude <- exclude.count <- NULL
  extraneous <- extraneous.this.day <- first.of.three.or.more <- ht.exp <- NULL
  include.count <- index <- last.of.three.or.more <- line <- max.ht.vel <- NULL
  max.whoinc.1.ht <- max.whoinc.2.ht <- max.whoinc.3.ht <- max.whoinc.4.ht <- NULL
  max.whoinc.6.ht <- maxdiff.next.ht <- maxdiff.prev.ht <- median.other.sd <- NULL
  min.ht.vel <- mindiff.next.ht <- mindiff.prev.ht <- pair <- pair.next <- NULL
  pair.prev <- param <- param.other <- prev.v <- sd.median <- sex <- subjid <- NULL
  swap.flag.1 <- swap.flag.2 <- tanner.months <- tbc.other.sd <- tbc.sd <- NULL
  tbc.sd.d <- tbc.sd.max <- tbc.sd.min <- tbc.sd.minus <- tbc.sd.next <- tbc.sd.plus <- NULL
  tbc.sd.prev <- tbc.sd.sw <- tbc.sd.t <- temp.diff <- temp.exclude <- v <- NULL
  v.d <- v.minus <- v.next <- v.orig <- v.plus <- v.prev <- v.sw <- v.t <- NULL
  valid.interior.measurement <- who.maxdiff.next.ht <- who.mindiff.next.ht <- NULL
  whoagegrp.ht <- whoinc.1.ht <- whoinc.2.ht <- whoinc.3.ht <- whoinc.4.ht <- NULL
  whoinc.6.ht <- whoinc.age.ht <- z.orig <- NULL

  data.df <- data.table(data.df, key = c('subjid', 'param', 'agedays', 'index'))

  if (parallel & !is.na(log.path)) {
    sink(
      sprintf(
        "%s/cleangrowth-%s-batch-%03d.log",
        log.path,
        Sys.Date(),
        data.df$batch[1]
      )
    )
  }

  if (!quietly)
    cat(sprintf(
      "[%s] Processing Batch #%s...\n",
      Sys.time(),
      data.df$batch[1]
    ))

  # NOTE: in each step, we redo temp SDEs

  # initial temp SDEs ----

  # save a copy of all original measurement values before any transformation
  data.df[, v.orig := v]

  if (!quietly)
    cat(sprintf(
      "[%s] Preliminarily identify potential extraneous...\n",
      Sys.time()
    ))
  data.df$exclude[temporary_extraneous_infants(data.df)] <- 'Exclude-Temporary-Extraneous-Same-Day'

  # capture a list of subjects with possible extraneous for efficiency later
  subj.dup <- data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', unique(subjid)]

  # carried forwards ----

  # 9.  Exclude values that are carried forward. For the purposes of this analysis, any value that is identical to the preceding value for the same parameter
  #     and subject is considered carried forward. The chances of having identical measurements, even at an age/interval when little or no growth would be expected,
  #     is fairly small, and when this is the case the carried forward measurements provide little new information.
  # a.	Calculate d_prev_wt=wt-wtprev and d_prev_ht=ht-htprev. Use original measurements rather than transformed
  #     measurements (unit errors and switches).
  # b.	Unlike most steps, do this step for all extraneous values (exc_*==2) in addition to included values (exc_*==0), comparing all values for one day to all
  #     values from the prior day – if there are any values with a d_prev*==0, the value on the latter day should be excluded.
  # c.	Replace exc_*=3 for all values with d_prev*==0 & (exc_*==0 OR exc_*==2)
  if (!include.carryforward) {
    if (!quietly)
      cat(sprintf(
        "[%s] Exclude measurements carried forward...\n",
        Sys.time()
      ))
    # for efficiency, bring get.prev and get.next inline here (working on valid rows within a single parameter for a single subject)
    # structure c(NA, field.name[-.N]) == get.prev
    data.df[, prev.v := as.double(NaN)]
    data.df[valid(data.df), prev.v := c(NA, v.orig[-.N]), by = .(subjid, param)]

    # optimize "carry forward" for children without extraneous.
    data.df[!(subjid %in% subj.dup) &
              v.orig == prev.v, exclude := 'Exclude-Carried-Forward']

    # need to handle children with extraneous measurements on same day separately
    data.df[subjid %in% subj.dup &
              valid(data.df, include.temporary.extraneous = TRUE), exclude := (function(df) {
                setkey(df, agedays)
                ages = unique(agedays)
                # no point in looking for measurements carried forward if all measurements are from a single day of life
                if (length(ages) > 1) {
                  # iterate over each age
                  for (i in 2:length(ages)) {
                    # find the set of measurements from the previous age in days
                    all.prev.v = df[agedays == ages[i - 1], v.orig]
                    # if a measurement for the current age is found in the set of measurements from the previous age, then mark it as carried forward
                    df[agedays == ages[i] &
                         v.orig %in% all.prev.v, exclude := 'Exclude-Carried-Forward']
                  }
                }
                return(df$exclude)
              })(copy(.SD)), .SDcols = c('agedays', 'exclude', 'v.orig'), by = .(subjid, param)]
  }

  # 9d.  Replace exc_*=0 if exc_*==2 & redo step 5 (temporary extraneous)
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # BIV ----

  # Wt limits from do-file Oct 10 2022, in dataset, highest plausible weight=29.5kg, 65 lbs
  # HC Limits based on analysis in do-file from Oct 11 2022: in dataset, lowest plausible 18.5, z=-13, highest plausible 63.5, z=11.2  (at birth, highest=42, z=6)


  exc_nam <- "Exclude-Absolute-BIV"

  # identify absolute cutoffs
  # Min/max weight at birth based on published births
  data.df[valid(data.df) & param == "WEIGHTKG" & measurement < 0.2,
          exclude := exc_nam]
  # Min weight after birth based on rules about who can go home with alowance for significant weight loss -- this is for outpatient data only so very low weight babies would still be in NICU
  data.df[valid(data.df) & param == "WEIGHTKG" & measurement > 10.5 &
            agedays == 0,
          exclude := exc_nam]
  # Max weight for <2y based on eval (see do-file from Oct 10 2022), highest plaus weight 29.5kg
  data.df[valid(data.df) & param == "WEIGHTKG" & measurement < 1 &
            agedays == 0,
          exclude := exc_nam]
  # Max weight for <2y based on eval (see do-file from Oct 10 2022), highest plaus weight 29.5kg
  data.df[valid(data.df) & param == "WEIGHTKG" & measurement > 35 &
            agedays < 2,
          exclude := exc_nam]
  # Max weight for all based on published data
  data.df[valid(data.df) & param == "WEIGHTKG" & measurement > 600,
          exclude := exc_nam]

  # Min/max HC based on analysis in do file from
  # Also, 18 is z=-6 for 22 0/7 in Fenton and 65 is z=6 for 40 0/7
  data.df[valid(data.df) & param == "HEIGHTCM" & measurement < 18,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEIGHTCM" & measurement > 65 &
            agedays == 0,
          exclude := exc_nam]

  # Min/max HC based on analysis in do file from Oct 11 2022
  # Also, 13 is z=-6 for 22 0/7 in Fenton and
  data.df[valid(data.df) & param == "HEADCM" & measurement < 13,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEADCM" & measurement > 75,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEADCM" & measurement > 50 &
            agedays == 0,
          exclude := exc_nam]

  exc_nam <- "Exclude-Z-BIV"

  # identify z cutoff
  # ***Note, using unrecentered values***
  #  *For weight only do after birth
  data.df[valid(data.df) & param == "WEIGHTKG" & z.orig < - 25 &
            agedays < 1,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "WEIGHTKG" & z.orig < -15 &
            agedays >= 1,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "WEIGHTKG" & z.orig > 15,
          exclude := exc_nam]

  # *Max z-score for height based on analysis of CHOP data because 15/25 too loose for upper limits
  data.df[valid(data.df) & param == "HEIGHTCM" & z.orig < - 25 &
            agedays < 1,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEIGHTCM" & z.orig < -15 &
            agedays >= 1,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEIGHTCM" & z.orig > 8,
          exclude := exc_nam]

  # head circumference
  data.df[valid(data.df) & param == "HEADCM" & z.orig < -15,
          exclude := exc_nam]
  data.df[valid(data.df) & param == "HEADCM" & z.orig > 15,
          exclude := exc_nam]


  # 9d.  Replace exc_*=0 if exc_*==2 & redo step 5 (temporary extraneous)
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # end ----

  if (!quietly)
    cat(sprintf("[%s] Completed Batch #%s...\n", Sys.time(), data.df$batch[1]))
  if (!quietly & parallel)
    sink()
  return(data.df[j = .(line, exclude, tbc.sd, tbc.other.sd, param)]) #debugging
}
