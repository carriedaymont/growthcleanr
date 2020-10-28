### Define the helper function and the carryforward adjustment function ###
# Segments labeled with ADJUSTCF EDIT in comments are areas where substantive
# changes were made to the original cleangrowth logic. There are other changes
# as well that are not specifically called out.

# helper function to treat NA values as FALSE
na.as.false = function(v) {
  v[is.na(v)] = F
  v
}

#' adjustcarryforward
#' \code{adjustcarryforward} Uses absolute height velocity to identify values
#' excluded as carried forward values for reinclusion.
#' @param subjid Vector of unique identifiers for each subject in the database.
#' @param param Vector identifying each measurement, may be 'WEIGHTKG', 'HEIGHTCM', or 'LENGTHCM'
#'   'HEIGHTCM' vs. 'LENGTHCM' only affects z-score calculations between ages 24 to 35 months (730 to 1095 days).
#'   All linear measurements below 731 days of life (age 0-23 months) are interpreted as supine length, and
#'   all linear measurements above 1095 days of life (age 36+ months) are interpreted as standing height.
#' @param agedays Numeric vector containing the age in days at each measurement.
#' @param sex Vector identifying the gender of the subject, may be 'M', 'm', or 0 for males, vs. 'F',
#'  'f' or 1 for females.
#' @param measurement Numeric vector containing the actual measurement data.  Weight must be in
#'   kilograms (kg), and linear measurements (height vs. length) in centimeters (cm).
#' @param orig.exclude Vector of exclusion assessment results from cleangrowth()
#' @param sd.recenter Data frame or table with median SD-scores per day of life
#' @param ewma.exp Exponent to use for weighting measurements in the exponentially weighted moving
#'  average calculations. Defaults to -1.5. This exponent should be negative in order to weight growth
#'  measurements closer to the measurement being evaluated more strongly. Exponents that are further from
#'  zero (e.g. -3) will increase the relative influence of measurements close in time to the measurement
#'  being evaluated compared to using the default exponent.
#' @param ref.data.path Path to reference data. If not supplied, the year 2000
#' Centers for Disease Control (CDC) reference data will be used.
#' @param quietly Determines if function messages are to be displayed and if log files (parallel only)
#' are to be generated. Defaults to TRUE.
#' @param minfactor Sweep variable for computing mindiff.next.ht in 15f, default 0.5
#' @param maxfactor Sweep variable for computing maxdiff.next.ht in 15f, default 2
#' @param banddiff Sweep variable for computing mindiff.next.ht in 15f, default 3
#' @param banddiff_plus Sweep variable for computing maxdiff.next.ht in 15, default 5.5
#' @param min_ht.exp_under Sweep variable for computing ht.exp in 15f, default 2
#' @param min_ht.exp_over Sweep variable for computing ht.exp in 15f, default 0
#' @param max_ht.exp_under Sweep variable for computing ht.exp in 15f, default 0.33
#' @param max_ht.exp_over Sweep variable for computing ht.exp in 15f, default 1.5
#'
#' @return Re-evaluated exclusion assessments based on height velocity.
#'
#' @export
#' @rawNamespace import(plyr, except = c(failwith, id, summarize, count, desc, mutate, arrange, rename, is.discrete, summarise, summarize))
#' @rawNamespace import(dplyr, except = c(last, first, summarize, src, between))
#' @import data.table
adjustcarryforward <- function(subjid,
                               param,
                               agedays,
                               sex,
                               measurement,
                               orig.exclude,
                               sd.recenter = NA,
                               ewma.exp = -1.5,
                               ref.data.path = "",
                               quietly = T,
                               minfactor = 0.5,
                               maxfactor = 2,
                               banddiff = 3,
                               banddiff_plus = 5.5,
                               min_ht.exp_under = 2,
                               min_ht.exp_over = 0,
                               max_ht.exp_under = 0.33,
                               max_ht.exp_over = 1.5) {
  # organize data into a dataframe along with a line "index" so the original data order can be recovered
  data.all = data.table(
    line = seq_along(measurement),
    subjid = as.factor(subjid),
    param,
    agedays = as.integer(agedays),
    v = ifelse(measurement == 0, NaN, measurement),
    sex = as.integer(ifelse(
      sex %in% c(0, 'm', 'M'), 0, ifelse(sex %in% c(1, 'f', 'F'), 1, NA)
    )),
    orig.exclude = as.factor(orig.exclude)
  )

  ### ADJUSTCF EDIT
  data.all = data.all[, n := 1:.N]
  ### ENDEDIT

  data.orig = data.all

  setkey(data.all, subjid)

  subjid.unique = unique(data.all$subjid)

  #### ADJUSTCF EDIT ####
  # for this purpose, want to subset dataset down to just "Exclude-Carried-Forward" and "Include" - assume all other measurements are invalid
  # want to remove any carried forward values whose non-carried forward is also excluded
  # data.all <- data.all %>%
  #   mutate(orig.exclude.lag = lag(orig.exclude, n = 1)) %>%
  #   filter(!(
  #     orig.exclude == "Exclude-Carried-Forward" &
  #       grepl("exclude", orig.exclude.lag, ignore.case = T)
  #   )) %>%
  #   filter(orig.exclude %in% c("Exclude-Carried-Forward", "Include")) %>%
  #   select(-orig.exclude.lag)

  # NEW EDIT --
  # remove all the weight measurements
  data.all <- data.all %>%
    filter(param %in% c("HEIGHTCM", "LENGTHCM"))

  # filter to only subjects with possible carried forwards - n is here to merge back
  # if they have all includes, filter them out
  data.all <- data.all %>%
    filter(subjid %in% data.all$subjid[data.all$orig.exclude == "Exclude-Carried-Forward"]) %>%
    as.data.table()

  # here's what we want to filter out -- anything that's not carried forward/include
  # we're also going to include strings of carried forward
  # but we also need to make sure they're not coming from an excluded value

  # start of string to remove: everything that isn't include/excl-cf
  st <- which(!data.all$orig.exclude %in% c("Exclude-Carried-Forward", "Include"))
  # end of string: include or the end of a subject
  subj_end <- length(data.all$subjid)-match(unique(data.all$subjid),rev(data.all$subjid))+1
  end <- c(which(data.all$orig.exclude == "Include"),subj_end)
  end <- unique(sort(end))

  # remove anything between start and ends (including start, not including end)
  to_rem <- unlist(
    lapply(st, function(x){
      to_rem <- c(x:(end[end > x][1]))
      if (to_rem[length(to_rem)] %in% subj_end){
        # if it's the last value, we want to get rid of that end
        return(to_rem)
      } else {
        # if it's an include, we want to keep it (don't remove)
        return(to_rem[-length(to_rem)])
      }
    })
  )
  to_rem <- unique(to_rem)

  data.all <- data.all[-to_rem,]

  # filter to only subjects with possible carried forwards again
  data.all <- data.all %>%
    filter(subjid %in% data.all$subjid[data.all$orig.exclude == "Exclude-Carried-Forward"]) %>%
    as.data.table()

  ### END EDIT ####

  # load tanner height velocity data. sex variable is defined such that 0=male and 1=female
  # recode column names to match syntactic style ("." rather than "_" in variable names)
  tanner_ht_vel_path <- ifelse(
    ref.data.path == "",
    system.file("extdata/tanner_ht_vel.csv", package = "growthcleanr"),
    paste(ref.data.path, "tanner_ht_vel.csv", sep =
            "")
  )

  tanner.ht.vel = fread(tanner_ht_vel_path)

  setnames(tanner.ht.vel,
           colnames(tanner.ht.vel),
           gsub('_', '.', colnames(tanner.ht.vel)))
  setkey(tanner.ht.vel, sex, tanner.months)
  # keep track of column names in the tanner data
  tanner.fields = colnames(tanner.ht.vel)
  tanner.fields = tanner.fields[!tanner.fields %in% c('sex', 'tanner.months')]

  who_max_ht_vel_path <- ifelse(
    ref.data.path == "",
    system.file("extdata/who_ht_maxvel_3sd.csv", package = "growthcleanr"),
    paste(ref.data.path, "who_ht_maxvel_3sd.csv", sep =
            "")
  )

  who_ht_vel_3sd_path <- ifelse(
    ref.data.path == "",
    system.file("extdata/who_ht_vel_3sd.csv", package = "growthcleanr"),
    paste(ref.data.path, "who_ht_vel_3sd.csv", sep =
            "")
  )
  who.max.ht.vel = fread(who_max_ht_vel_path)
  who.ht.vel = fread(who_ht_vel_3sd_path)
  setkey(who.max.ht.vel, sex, whoagegrp_ht)
  setkey(who.ht.vel, sex, whoagegrp_ht)
  who.ht.vel = as.data.table(dplyr::full_join(who.ht.vel, who.max.ht.vel, by =
                                                c('sex', 'whoagegrp_ht')))

  setnames(who.ht.vel, colnames(who.ht.vel), gsub('_', '.', colnames(who.ht.vel)))
  setkey(who.ht.vel, sex, whoagegrp.ht)
  # keep track of column names in the who growth velocity data
  who.fields = colnames(who.ht.vel)
  who.fields = who.fields[!who.fields %in% c('sex', 'whoagegrp.ht')]

  # 1.  General principles
  # a.	All steps are done separately for each parameter unless otherwise noted
  # b.	All steps are done sorted by subject, parameter, and age (in days) for nonexcluded and nonmissing values only unless otherwise noted. This is very important.
  #     Sorting needs to be redone with each step to account for excluded and transformed  values.
  # c.	The next value refers to the value with the next highest age for the same parameter and the same subject, and the previous value refers to the value with the
  #     next lowest age for the same parameter and the same subject.
  # d.	You will need to set up a method for keeping track of whether a value is missing or excluded (and in what step). I use variables called exc_* that are =0
  #     if a value is to be included, =1 if missing, and =2 or higher if it is to be excluded, with each number indicating a different step. I also set up
  #     parameter-specific subjid_* variables that are = to subjid for included values and are blank if the value is missing or should be excluded. These subjid_*
  #     variables need to be updated with each step.
  # e.  All steps assume that data are sorted by subjid_*, parameter, and age (in days) for nonexcluded and nonmissing values only
  #     unless otherwise noted. Sorting needs to be redone after any transformations or exclusions to account for excluded and
  #     transformed values.
  # f.  The next value refers to the nonexcluded nonmissing value with the next highest age for the same parameter and the same
  #     subject, and the previous value refers to the nonexcluded nonmissing value with the next lowest age for the same parameter
  #     and the same subject.
  # g.  exc_* should only be replaced with a  higher value if exc_*==0 at the time of replacement, unless otherwise specified.


  # NOTE: in the R code below exclusion is documented as a series of factor levels, where all levels occuring before 'Exclude' in the sequence are considered
  # to be valid measurements.  We use the built in sorting of the data.table object and subsets rather than re-sorting at each step
  # to ensure that only valid measurements are used at the beginning of each step.
  # Also, unlike the Stata code, the measurement parameter (weight vs. height) is recorded as a factor in the data frame, rather than as a variable name

  # 2.  Data set-up
  # a.	I always code sex as 0=Male, 1=Female, so I recoded the variable sex that way and left a variable sexorigcode the way the data was sent to me (1=Female 2=Male)
  # b.	Remove rows that are duplicates for subjid, param, and measurement from further analysis
  #     NOTE: this step is not needed -- handled automatically by "temporary duplicate" step.
  # c.  I generated separate variables for weight (wt) and height (ht), as well as exc_* and subjid_* variables. Set exc_*=0 if value is not missing
  #     and exc_*=1 if value is missing. In all future steps, exc_* should only be changed if it is 0. This helps to keep track of which step excluded a value.
  #     I also kept the measurement variable there and untouched because sometimes wt and ht got transformed to something else.
  # d.	I made tables based on CDC growth curve parameters that include data for each day that I will send separately. The LMS parameters for each day are
  #     cubically interpolated from the values by month available on the CDC website. Create wt and ht z-scores for each value of each parameter (Z: WtZ and HtZ).
  # e.	There are variables in the table labelled cdc_*_csd_pos and cdc_*_csd_neg. For each age and sex, these correspond to ½ of the absolute value of the
  #     median and the value with a z-score of +2 (csd_pos) and -2 (csd_neg). These can be created to generate a score similar to the z-score but with an
  #     important difference. The z-scores created using the LMS method account for skewness in the distribution of the parameters (particularly weight), which
  #     can lead to small changes in z-score with large changes in weight in subjects with very high weight, and relatively large changes in z-score for smaller
  #     changes in weight in subjects with low weights.  The score we will create can be called an SD-score (SDorig: WtSDorig and HtSDorig that is calculated by
  #     dividing the difference between the value and the median by the SD score (use csd_pos if the value is above the median, csd_neg if the value is below the
  #     median). These SD-scores, rather than z-scores, now form the basis for the algorithm.

  # calculate z scores
  if (!quietly)
    cat(sprintf("[%s] Calculating z-scores...\n", Sys.time()))
  measurement.to.z = read.anthro(ref.data.path, cdc.only = T)
  data.all[, z.orig := measurement.to.z(param, agedays, sex, v)]

  # calculate "standard deviation" scores
  if (!quietly)
    cat(sprintf("[%s] Calculating SD-scores...\n", Sys.time()))
  data.all[, sd.orig := measurement.to.z(param, agedays, sex, v, T)]

  # sort by subjid, param, agedays
  setkey(data.all, subjid, param, agedays)

  # add a new convenience index for bookkeeping
  data.all[, index := 1:.N]

  # enumerate the different exclusion levels
  exclude.levels = c(
    'Missing',
    'No Change',
    'Include',
    'Exclude-Min-Height-Change',
    'Exclude-Max-Height-Change'
  )

  # Mark missing values for exclusion
  data.all[, exclude := factor(with(data.all, ifelse(
    is.na(v) |
      agedays < 0, 'Missing', 'No Change'
  )),
  levels = exclude.levels,
  ordered = T)] # why is this ordered??

  # after calculating z scores, for convenience, recategorize linear parameters as 'HEIGHTCM'
  data.all[param == 'LENGTHCM', param := 'HEIGHTCM']

  # define field names needed by helper functions
  ewma.fields = c('ewma.all', 'ewma.before', 'ewma.after')

  # 3.  SD-score recentering: Because the basis of the method is comparing SD-scores over time, we need to account for the fact that
  #     the mean SD-score for the population changes with age.
  # a.	Determine the median cdc*sd for each parameter by year of age (with sexes combined): median*sd.
  # b.	The median*sd should be considered to apply to midyear-age, defined as the age in days with the same value as the integer
  #     portion of (365.25*year + 365.25/2).
  # c.	Linearly interpolate median*sd for each parameter between each midyear-age, naming the interpolated values rc*sd.
  # d.	For ages below the first midyear-age, let rc*sd equal the median*sd for the earliest year.
  #     For ages above the last midyear_age, let rc*sd equal the median*sd for the last year.
  # e.	Subtract rcsd_* from SDorig to create the recentered SD-score.  This recentered SD-score, labeled tbc*sd
  #     (stands for "to be cleaned") will be used for most of the rest of the analyses.
  # f.	In future steps I will sometimes refer to measprev and measnext which refer to the previous or next wt or ht measurement
  #     for which exc_*==0 for the subject and parameter, when the data are sorted by subject, parameter, and agedays. SDprev and SDnext refer to the tbc*sd of the previous or next measurement.

  if (!quietly)
    cat(sprintf("[%s] Re-centering data...\n", Sys.time()))

  # see function definition below for explanation of the re-centering process
  # returns a data table indexed by param, sex, agedays
  if (!is.data.table(sd.recenter)) {
    # ADJUSTCF EDIT - be explicit about levels to keep
    keep.levels <- c(
      "Include",
      "Unit-Error-High",
      "Unit-Error-Low",
      "Unit-Error-Possible",
      "Swapped-Measurements"
    )
    sd.recenter = data.all[orig.exclude %in% keep.levels, sd.median(param, sex, agedays, sd.orig)]
    # END EDIT
  }
  # add sd.recenter to data, and recenter
  setkey(data.all, param, sex, agedays)
  data.all = sd.recenter[data.all]
  setkey(data.all, subjid, param, agedays)
  data.all[, tbc.sd := sd.orig - sd.median]

  # safety check: treat observations where tbc.sd cannot be calculated as missing
  data.all[is.na(tbc.sd), exclude := 'Missing']



  ######### END DATA PROCESING #########

  ######### START FLAGGING #########

  # 15.  Exclude heights based on absolute differences in measurement. The key to this step is that once we identify pairs of measurements with implausible
  #      amounts of absolute difference between them, we have to determine which of the two values in the pair is less likely to be representative and should
  #      be excluded. For subjects/parameters with 3 or more measurements, this is done by looking at the dewma_* of each of the 2 values in a pair using a ewma that
  #      excludes the other value in the pair. For subjects/parameters with 2 or more measurements, this is done by looking at the absolute value of the tbc*sd.
  #      The Tanner height velocity reference is used for measurements taken at >2yo, WHO will be used for <2yo. For a few pairs of measurements either could be used;
  #      WHO will be used if difference between ages is < 9 months.
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude heights based on growth velocity...\n",
      Sys.time()
    ))
  data.all[param == 'HEIGHTCM', exclude := (function(subj.df) {
    # assign some book keeping variables
    #subj.df[, `:=`(subjid = subjid, param='HEIGHTCM',index=1:.N)]
    subj.df[, index := 1:.N]

    num.height.excluded = 0
    while (T) {
      # use a closure to discard all the extra fields added to df with each iteration
      subj.df[, exclude := (function (df) {
        # initialize fields
        df[, (ewma.fields) := as.double(NaN)]
        df[, `:=`(
          v.prev = as.double(NaN),
          v.next = as.double(NaN),
          dewma.after.prev = as.double(NaN),
          dewma.before.next = as.double(NaN),
          abs.tbc.sd.prev = as.double(NaN),
          abs.tbc.sd.next = as.double(NaN),
          agedays.next = as.integer(NaN),
          abs.2ndlast.sd = as.double(NaN),
          mindiff.prev.ht = as.double(NaN),
          mindiff.next.ht = as.double(NaN),
          maxdiff.prev.ht = as.double(NaN),
          maxdiff.next.ht = as.double(NaN),
          pair.prev = F,
          pair.next = F
        )]

        # ewma fields are needed later -- calculate now for efficiency
        df[, (ewma.fields) := ewma(agedays, tbc.sd, ewma.exp, T)]

        # calculate some usefule values (e.g. dewma values and tbc.sd) for use in later steps
        df[, `:=`(
          dewma.all = tbc.sd - ewma.all,
          dewma.before = tbc.sd - ewma.before,
          dewma.after = tbc.sd - ewma.after,
          abs.tbc.sd = abs(tbc.sd)
        )]

        # 15a.  As with steps 11 and 14, only one value will be excluded per round, and the step will be repeated until there are no more values to exclude
        # b.  For each height, calculate the d_age=agedays of next value-agedays of current value
        # NOTE: obtain next measurement, ewma.before and abs.tbc.sd as well since they are needed later

        # for efficiency, bring get.next inline here (working on valid rows within a single parameter for a single subject)
        # structure c(field.name[-1], NA) == get.next
        df[, `:=`(
          agedays.next = c(agedays[-1], NA),
          v.next = c(v[-1], NA),
          dewma.before.next = c(dewma.before[-1], NA),
          abs.tbc.sd.next = c(abs.tbc.sd[-1], NA)
        )]
        df[, delta.agedays.next := agedays.next - agedays]

        # 15c.	For each height, calculate mid_agedays=0.5*(agedays of next value + agedays of current value)
        df[, mid.agedays := 0.5 * (agedays.next + agedays)]

        # 15d.	Generate variable tanner_months= 6+12*(round(mid_agedays/365.25))
        # only calculate for rows that relate to height (may speed up subsequent processing)
        df[, tanner.months := 6 + 12 * (round(mid.agedays / 365.25))]

        # 15e.	Merge with dataset tanner_ht_vel using sex and tanner_months – this will give you min_ht_vel and max_ht_vel
        setkey(df, sex, tanner.months)
        df = tanner.ht.vel[df]

        # 15f.	Calculate the following:
        #   i.	mindiff_ht=0.5*min_ht_vel*(d_agedays/365.25)^2-3 if d_agedays<365.25
        #   ii.	replace mindiff_ht=0.5*min_ht_vel-3 if d_agedays>365.25
        df[, ht.exp := ifelse(delta.agedays.next < 365.25,
                              min_ht.exp_under,
                              min_ht.exp_over)]
        df[, `:=`(maxdiff.next.ht = as.double(NA),
                  mindiff.next.ht = as.double(NaN))]
        df[, mindiff.next.ht := minfactor * min.ht.vel * (delta.agedays.next /
                                                            365.25) ^ ht.exp - banddiff]

        # 15f.iii.	maxdiff_ht=2*max_ht_vel*(d_agedays/365.25)^1.5+5.5 if d_agedays>365.25
        #   iv.	replace maxdiff_ht=2*max_ht_vel*(d_agedays/365.25)^0.33+5.5 if d_agedays<365.25
        df[, ht.exp := ifelse(delta.agedays.next < 365.25,
                              max_ht.exp_under,
                              max_ht.exp_over)]
        df[, maxdiff.next.ht := maxfactor * max.ht.vel * (delta.agedays.next /
                                                            365.25) ^ ht.exp + banddiff_plus]

        # 15g.	Generate variable whoagegrp_ht=agedays/30.4375 rounded to the nearest integer
        df[, whoagegrp.ht := round(agedays / 30.4375)]

        # 15h.	Generate variable whoinc_age_ht based on values of d_agedays_ht according the the following table
        #   d_agedays_ht	whoinc_age_ht
        #   20-45	        1
        #   46-75	        2
        #   76-106	      3
        #   107-152	      4
        #   153-198	      6
        #   All others	  missing
        df[, whoinc.age.ht :=
             ifelse(delta.agedays.next < 20 ,
                    NA,
                    ifelse(
                      delta.agedays.next <= 45,
                      1,
                      ifelse(
                        delta.agedays.next <= 75,
                        2,
                        ifelse(
                          delta.agedays.next <= 106,
                          3,
                          ifelse(
                            delta.agedays.next <= 152,
                            4,
                            ifelse(delta.agedays.next <= 198, 6, NA)
                          )
                        )
                      )
                    ))]

        # i.	Merge using sex and whoagegrp_ht using who_ht_vel_3sd and who_ht_maxvel_3sd; this will give you varaibles whoinc_i_ht and maxwhoinc_i_ht
        #     for various intervals where i is 1,2, 3,4, 6 and corresponds to whoinc_age_ht.
        setkey(df, sex, whoagegrp.ht)
        df = who.ht.vel[df]

        # restore original sort order (ensures valid.rows variable applies to correct rows)
        setkey(df, index)

        # 15j.	Generate variable who_mindiff_ht=whoinc_i_ht according to the value if whoinc_age_ht; make who_mindiff_ht missing if whoinc_i_ht or whoinc_age_ht is missing.
        df[, who.mindiff.next.ht := ifelse(
          delta.agedays.next < 20 ,
          NA,
          ifelse(
            delta.agedays.next <= 45,
            whoinc.1.ht,
            ifelse(
              delta.agedays.next <= 75,
              whoinc.2.ht,
              ifelse(
                delta.agedays.next <= 106,
                whoinc.3.ht,
                ifelse(
                  delta.agedays.next <= 152,
                  whoinc.4.ht,
                  ifelse(delta.agedays.next <= 198, whoinc.6.ht, NA)
                )
              )
            )
          )
        )]

        # 15k.	Generate variable who_maxdiff_ht=max_whoinc_i_ht according to the value if whoinc_age_ht; make who_maxdiff_ht missing if max_whoinc_i_ht or
        #     whoinc_age_ht is missing.
        df[, who.maxdiff.next.ht := ifelse(
          delta.agedays.next < 20 ,
          NA,
          ifelse(
            delta.agedays.next <= 45,
            max.whoinc.1.ht,
            ifelse(
              delta.agedays.next <= 75,
              max.whoinc.2.ht,
              ifelse(
                delta.agedays.next <= 106,
                max.whoinc.3.ht,
                ifelse(
                  delta.agedays.next <= 152,
                  max.whoinc.4.ht,
                  ifelse(delta.agedays.next <= 198, max.whoinc.6.ht, NA)
                )
              )
            )
          )
        )]

        # 15l.	Scale allowed value based on d_agedays_ht:
        #   1.	replace who_mindiff_`p'=who_mindiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'<(whoinc_age_`p'*30.4375)
        #   2.	replace who_maxdiff_`p'=who_maxdiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'>(whoinc_age_`p'*30.4375)
        df[delta.agedays.next < whoinc.age.ht * 30.4375,
           `:=`(
             who.mindiff.next.ht = who.mindiff.next.ht * delta.agedays.next / (whoinc.age.ht *
                                                                                 30.4375),
             who.maxdiff.next.ht = who.maxdiff.next.ht * delta.agedays.next /
               (whoinc.age.ht * 30.4375)
           )]

        # 15m.	Replace mindiff_ht/maxdiff_ht with adjusted WHO value if Tanner value is missing or if both are present and age difference is < 9 months:
        #   1.	replace mindiff_`p'=0.5*who_mindiff_`p'-3 if who_mindiff_`p' is not missing & d_agedays_`p'<(9*30.4375)
        #   2.	replace maxdiff_`p'=2*who_maxdiff_`p'+3 if who_maxdiff_`p' is not missing & d_agedays_`p'<(9*30.4375)
        #   3.	replace mindiff_`p'=0.5*who_mindiff_`p'-3 if mindiff_`p' is missing & who_mindiff_`p' is not missing
        #   4.	replace maxdiff_`p'=2*who_maxdiff_`p'+3 if maxdiff_`p is missing & who_maxdiff_`p' is not missing

        # refactored logic slightly for efficiency
        df[!is.na(who.mindiff.next.ht) &
             (delta.agedays.next < 9 * 30.4375 |
                is.na(mindiff.next.ht)),
           `:=`(
             mindiff.next.ht = minfactor * who.mindiff.next.ht - banddiff,
             maxdiff.next.ht = maxfactor * who.maxdiff.next.ht + banddiff
           )]

        # 15m.5.  replace mindiff_`p'=-3 if mindiff_`p' is missing
        df[is.na(mindiff.next.ht), mindiff.next.ht := -3]

        # 15n.	Determine the min/maxdiffs for the previous age: mindiff_prev_ht, maxdiff_prev_ht
        # NOTE: obtain previous height, ewma.after value and abs.tbc.sd as well since they are needed in next steps

        # for efficiency, bring get.prev inline here (working on valid rows within a single parameter for a single subject)
        # structure c(NA, tbc.sd[-.N]) == get.prev
        df[, `:=`(
          v.prev = c(NA, v[-.N]),
          dewma.after.prev = c(NA, dewma.after[-.N]),
          abs.tbc.sd.prev = c(NA, abs.tbc.sd[-.N]),
          mindiff.prev.ht = c(NA, mindiff.next.ht[-.N]),
          maxdiff.prev.ht = c(NA, maxdiff.next.ht[-.N])
        )]

        # 15o.	Determine d_prev_ht=ht-htprev (set to missing for the first value for a subject) and d_next_ht=htnext-ht (set to missing for the last value for a subject)
        df[, `:=`(delta.prev.ht = v - v.prev,
                  delta.next.ht = v.next - v)]

        # 15p.  Perform a EWMA calculation with the following modifications:
        #  i.	  Generate a variable pair=1 if (d_prev_ht<mindiff_prev_ht OR d_ht<mindiff_ht OR d_prev_ht>maxdiff_prev_ht  OR d_ht>maxdiff_ht) AND exc_ht==0
        df[, pair := na.as.false(
          delta.prev.ht < mindiff.prev.ht |
            delta.next.ht < mindiff.next.ht |
            delta.prev.ht > maxdiff.prev.ht |
            delta.next.ht > maxdiff.next.ht
        )]

        # for efficiency, bring get.prev and get.next inline here (working on valid rows within a single parameter for a single subject)
        # structure c(NA, field.name[-.N]) == get.prev
        # structure c(field.name[-1], NA) == get.next
        df[, `:=`(pair.prev = c(F, pair[-.N]),
                  pair.next = c(pair[-1], F))]

        #  ii.	Generate bef_g_aftm1=1 if |Δewma_htbef| for the value of interest is greater than |Δewma_htaft| for the previous value
        #       AND the value of interest is not the first height value for that subject AND pair==1 AND pair for the previous value==1

        #  iii.	Generate aft_g_befp1=1 if |Δewma_htaft| for the value of interest is greater than |Δewma_htbef| for the next value
        #       AND the value of interest is not the last height value for that subject AND pair==1 AND pair for the next value==1
        # NOTE: pair.next will be NA last height, which will result in a FALSE value below
        df[, `:=`(
          bef.g.aftm1 = na.as.false(
            abs(dewma.before) > abs(dewma.after.prev)  & pair & pair.prev
          ),
          aft.g.befp1 = na.as.false(
            abs(dewma.after)  > abs(dewma.before.next) & pair & pair.next
          )
        )]

        #  iv.	Determine tbchtsd for each value as well as the one before prev_tbchtsd and after next_tbchtsd it
        # NOTE: done previously for efficiency

        # 15p.v.  Determine the total number of ht values for each subject (tot_ht)
        # NOTE: all rows are valid due to constraint in subj.df[...] statement
        num.valid = .N

        # 15q.	Identify a value for possible exclusion if one of the following sets of criteria are met. For values identified by each set of criteria determine
        #       the value of temp_diff using the formula given
        #   i.	d_prev_ht<mindiff_prev_ht & bef_g_aftm1_ht==1 & exc_ht==0 & mindiff_prev_ht is not missing
        #     a.  (temp_diff=|dewma_ht_bef|)
        df[, temp.diff := as.double(NaN)]
        df[, temp.exclude := factor(NA, levels = exclude.levels, ordered =
                                      T)]
        df[delta.prev.ht < mindiff.prev.ht & bef.g.aftm1,
           `:=`(temp.diff = abs(dewma.before),
                temp.exclude = 'Exclude-Min-Height-Change')]

        #   ii.	d_ht<mindiff_ht & aft_g_befp1_ht==1 & exc_ht==0 & mindiff_ht is not missing
        #     a.	(temp_diff=|dewma_ht_aft|)
        df[delta.next.ht < mindiff.next.ht & aft.g.befp1,
           `:=`(temp.diff = abs(dewma.after),
                temp.exclude = 'Exclude-Min-Height-Change')]

        #   iii.	d_prev_ht>maxdiff_prev_ht & bef_g_aftm1_ht==1 & exc_ht==0 & mindiff_prev_ht is not missing
        #     a.  (temp_diff=|dewma_ht_bef|)
        df[delta.prev.ht > maxdiff.prev.ht & bef.g.aftm1,
           `:=`(temp.diff = abs(dewma.before),
                temp.exclude = 'Exclude-Max-Height-Change')]

        #   iv.	d_ht>maxdiff_ht & aft_g_befp1_ht==1 & exc_ht==0 & mindiff_ht is not missing
        #     a.  (temp_diff=|dewma_ht_aft|)
        df[delta.next.ht > maxdiff.next.ht & aft.g.befp1,
           `:=`(temp.diff = abs(dewma.after),
                temp.exclude = 'Exclude-Max-Height-Change')]

        #   v.	d_prev_ht<mindiff_prev_ht & tot_ht==2 & |tbchtsd|>|prev_tbchtsd|
        #     a. for v-viii temp_diff is kept as missing
        #   vi. d_ht<mindiff_ht & tot_ht==2 & |tbchtsd|>|next_tbchtsd|
        df[delta.prev.ht < mindiff.prev.ht &
             num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev
           |
             delta.next.ht < mindiff.next.ht &
             num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next,
           temp.exclude := 'Exclude-Min-Height-Change']

        #   vii.	d_prev_ht>maxdiff_prev_ht & tot_ht==2 & |tbchtsd|>|prev_tbchtsd|
        #   viii. d_ht>maxdiff_ht & tot_ht==2 & |tbchtsd|>|next_tbchtsd|
        df[delta.prev.ht > maxdiff.prev.ht &
             num.valid == 2 & abs.tbc.sd > abs.tbc.sd.prev
           |
             delta.next.ht > maxdiff.next.ht &
             num.valid == 2 & abs.tbc.sd > abs.tbc.sd.next,
           temp.exclude := 'Exclude-Max-Height-Change']

        # r.  If there is only one potential exclusion identified in step 15j for a subject and parameter,
        #     replace exc_ht=15 for that value if it met criteria i, ii, v, or vi  and exc_ht=16 if it met criteria iii, iv, vii, or viii
        # NOTE: these exclusions are assigned in the code above as 'Exclude-Min-Height-Change' and 'Exclude-Max-Height-Change'

        ##### ADJUSTCF EDIT #####
        # if you're outside of the bands OR if you are not a carry forward then you have no change
        df$temp.exclude <- sapply(1:nrow(df), function(x) {
          ifelse(
            (
              df$temp.exclude[x] %in% c(
                'Exclude-Min-Height-Change',
                'Exclude-Max-Height-Change'
              )
            ) |
              (df$orig.exclude[x] != "Exclude-Carried-Forward"),
            'No Change',
            'Include'
          )
        })
        ##### END EDIT #####

        rep = df$temp.exclude == "Include"
        num.exclude = sum(rep)
        if (num.exclude == 1)
          df[rep, exclude := temp.exclude]

        # s.  If there is more than one potential exclusion identified in step 14h for a subject and parameter, determine which value has the largest temp_diff and
        #     replace exc_ht=15 for that value if it met criteria i, ii, v, or vi and exc_ht=16 for that value if it met criteria iii,  iv, vii, or viii.

        if (num.exclude > 1) {
          # first order by decreasing temp.diff (where rep=T)
          worst.row = order(rep, df$temp.diff, decreasing = T)[1]
          df[worst.row, exclude := temp.exclude]
        }

        return(df$exclude)


      })(copy(.SD))]



      # t.  If there was at least one subject who had a potential exclusion identified in step 15q, repeat steps 15b-15q. If there were no subjects with potential
      #     exclusions identified in step 15q, move on to step 16.
      newly.excluded = sum(subj.df$exclude %in% c('Include'))
      if (newly.excluded > num.height.excluded) {
        num.height.excluded = newly.excluded
      } else {
        break
      }
    }

    setkey(subj.df, index)
    return(subj.df$exclude)
  })(copy(.SD)), by = .(subjid), .SDcols = c('sex', 'agedays', 'v', 'tbc.sd', 'exclude', 'orig.exclude')]

  return(rbind(
    data.frame(adjustcarryforward = data.all$exclude, n = data.all$n),
    data.frame(
      filter(data.orig,!n %in% data.all$n) %>% mutate(adjustcarryforward = "Missing")  %>% select(adjustcarryforward, n)
    )
  ))
}
