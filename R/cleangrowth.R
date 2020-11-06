#' Clean growth measurements
#'
#' @param subjid Vector of unique identifiers for each subject in the database.
#' @param param Vector identifying each measurement, may be 'WEIGHTKG', 'HEIGHTCM', or 'LENGTHCM'
#'   'HEIGHTCM' vs. 'LENGTHCM' only affects z-score calculations between ages 24 to 35 months (730 to 1095 days).
#'   All linear measurements below 731 days of life (age 0-23 months) are interpreted as supine length, and
#'   all linear measurements above 1095 days of life (age 36+ months) are interpreted as standing height.
#' @param agedays Numeric vector containing the age in days at each measurement.
#' @param sex Vector identifying the gender of the subject, may be 'M', 'm', or 0 for males, vs. 'F', 'f' or 1 for females.
#' @param measurement Numeric vector containing the actual measurement data.  Weight must be in
#'   kilograms (kg), and linear measurements (height vs. length) in centimeters (cm).
#' @param recover.unit.error Indicates whether the cleaning algorithm should
#' attempt to identify unit errors (I.e. inches vs. cm, lbs vs. kg). If unit
#' errors are identified, the value will be corrected and retained within the
#' cleaning algorithm as a valid measurement.  Defaults to FALSE.
#' @param sd.extreme Measurements more than sd.extreme standard deviations from
#' the mean (either above or below) will be flagged as invalid. Defaults to 25.
#' @param z.extreme Measurements with an absolute z-score greater than
#' z.extreme will be flagged as invalid. Defaults to 25.
#' @param lt3.exclude.mode Determines type of exclusion procedure to use for 1 or 2 measurements of one type without
#' matching same ageday measurements for the other parameter. Options include "default" (standard growthcleanr approach),
#' and "flag.both" (in case of two measurements of one type without matching values for the other parameter, flag both
#' for exclusion if beyond threshold)
#' @param error.load.mincount minimum count of exclusions on parameter before
#' considering excluding all measurements. Defaults to 2.
#' @param error.load.threshold threshold of percentage of excluded measurement count to included measurement
#' count that must be exceeded before excluding all measurements of either parameter. Defaults to 0.5.
#' @param sd.recenter Data frame or table with median SD-scores per day of life
#' by gender and parameter. Columns in the table must include param, sex,
#' agedays, and sd.median.  If not supplied, the median values will be
#' calculated using the growth data that is being cleaned. Defaults to NA.
#' @param sdmedian.filename Name of file to save sd.median data calculated on the input dataset to as CSV.
#' Defaults to "", for which this data will not be saved. Use for extracting medians for parallel processing
#' scenarios other than the built-in parallel option.
#' @param sdrecentered.filename Name of file to save re-centered data to as CSV. Defaults to "", for which this
#' data will not be saved. Useful for post-processing and debugging.
#' @param include.carryforward Determines whether Carry-Forward values are kept in the output. Defaults to False.
#' @param ewma.exp Exponent to use for weighting measurements in the
#' exponentially weighted moving average calculations. Defaults to -1.5.
#' This exponent should be negative in order to weight growth measurements
#' closer to the measurement being evaluated more strongly. Exponents that are
#' further from zero (e.g. -3) will increase the relative influence of
#' measurements close in time to the measurement being evaluated compared to
#' using the default exponent.
#' @param ref.data.path Path to reference data. If not supplied, the year 2000
#' Centers for Disease Control (CDC) reference data will be used.
#' @param log.path Path to log file output when running in parallel (non-quiet mode). Default is ".". A new
#' directory will be created if necessary.
#' @param parallel Determines if function runs in parallel.  Defaults to FALSE.
#' @param num.batches Specify the number of batches to run in parallel. Only
#' applies if parallel is set to TRUE. Defaults to the number of workers
#' returned by the getDoParWorkers function in the foreach package.
#' @param quietly Determines if function messages are to be displayed and if log files (parallel only) are to be generated.
#' Defaults to TRUE
#'
#' @return Vector of exclusion codes for each of the input measurements.
#'
#'    Possible values for each code are:
#'
#' * 'Include', 'Unit-Error-High', 'Unit-Error-Low', 'Swapped-Measurements', 'Missing',
#' *  'Exclude-Carried-Forward', 'Exclude-SD-Cutoff', 'Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair',
#' *  'Exclude-Duplicate',
#' *  'Exclude-EWMA-8', 'Exclude-EWMA-9', 'Exclude-EWMA-10', 'Exclude-EWMA-11', 'Exclude-EWMA-12',
#'    'Exclude-EWMA-13', 'Exclude-EWMA-14',
#' *  'Exclude-Min-Height-Change', 'Exclude-Max-Height-Change',
#' *  'Exclude-Pair-Delta-17', 'Exclude-Pair-Delta-18', 'Exclude-Pair-Delta-19',
#' *  'Exclude-Single-Outlier', 'Exclude-Too-Many-Errors', 'Exclude-Too-Many-Errors-Other-Parameter'
#' @md
#'
#' @export
#'
#' @import data.table
#' @import foreach
#' @import doParallel
#' @importFrom plyr ddply
#'
#' @examples
#' # Run calculation using a small subset of given data
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid %in% unique(df_stats[, "subjid"])[1:5], ]
#'
#' clean_stats <- cleangrowth(
#'   subjid = df_stats$subjid,
#'   param = df_stats$param,
#'   agedays = df_stats$agedays,
#'   sex = df_stats$sex,
#'   measurement = df_stats$measurement
#' )
#'
#' # Once processed you can filter data based on result value
#' df_stats <- cbind(df_stats, "clean_result" == clean_stats)
#' clean_df_stats <- df_stats[, df_stats$clean_result == "Include"]
#'
#' # Parallel processing: run using 3 cores and batches
#' df_stats <- cleangrowth(
#'   subjid = df_stats$subjid,
#'   param = df_stats$param,
#'   agedays = df_stats$agedays,
#'   sex = df_stats$sex,
#'   measurement = df_stats$measurement,
#'   parallel = TRUE,
#'   num.batches = 2
#' )
cleangrowth <- function(subjid,
                        param,
                        agedays,
                        sex,
                        measurement,
                        recover.unit.error = FALSE,
                        sd.extreme = 25,
                        z.extreme = 25,
                        lt3.exclude.mode = "default",
                        height.tolerance.cm = 2.5,
                        error.load.mincount = 2,
                        error.load.threshold = 0.5,
                        sd.recenter = NA,
                        sdmedian.filename = "",
                        sdrecentered.filename = "",
                        include.carryforward = FALSE,
                        ewma.exp = -1.5,
                        ref.data.path = NULL,
                        log.path = ".",
                        parallel = FALSE,
                        num.batches = NA,
                        quietly = TRUE) {
  # ==== Dealing with "undefined global functions or variables" ==== #
  ## Only for variable which couldn't be quoted everywhere
  z.orig <- v <- sd.orig <- index <- exclude <- tbc.sd <- sd.median <- NULL
  # ==== Dealing with "undefined global functions or variables" ==== #

  # organize data into a dataframe along with a line "index" so the original data order can be recovered
  data.all <- data.table(
    line = seq_along(measurement),
    subjid = as.factor(subjid),
    param,
    agedays = as.integer(agedays),
    v = ifelse(measurement == 0, NA_real_, measurement),
    sex = as.integer(ifelse(
      test = sex %in% c(0, "m", "M"),
      yes = 0L,
      no = ifelse(sex %in% c(1, "f", "F"), 1L, NA_integer_)
    ))
  )

  # if parallel processing is desired, load additional modules
  if (parallel) {
    registerDoParallel(cores = num.batches)
    if (is.na(num.batches)) num.batches <- getDoParWorkers()
  } else {
    if (is.na(num.batches)) num.batches <- 1
  }

  setkeyv(data.all, "subjid")
  subjid.unique <- data.all[j = unique(subjid)]
  batches.all <- data.table(
    subjid = subjid.unique,
    batch = sample(num.batches, length(subjid.unique), replace = TRUE),
    key = "subjid"
  )
  data.all <- batches.all[data.all]

  # load tanner height velocity data. sex variable is defined such that 0=male and 1=female
  # recode column names to match syntactic style ("." rather than "_" in variable names)
  tanner_ht_vel_path <- ifelse(
    is.null(ref.data.path),
    system.file("extdata/tanner_ht_vel.csv", package = "growthcleanr"),
    paste0(ref.data.path, "tanner_ht_vel.csv")
  )

  tanner.ht.vel <- fread(tanner_ht_vel_path)

  setnames(
    tanner.ht.vel,
    colnames(tanner.ht.vel),
    gsub("_", ".", colnames(tanner.ht.vel))
  )
  setkeyv(tanner.ht.vel, c("sex", "tanner.months"))
  # keep track of column names in the tanner data
  tanner.fields <- colnames(tanner.ht.vel)
  tanner.fields <- tanner.fields[!tanner.fields %in% c("sex", "tanner.months")]

  who.ht.vel <- merge(
    x = fread(ifelse(
      is.null(ref.data.path),
      system.file("extdata/who_ht_vel_3sd.csv", package = "growthcleanr"),
      paste0(ref.data.path, "who_ht_vel_3sd.csv")
    )),
    y = fread(ifelse(
      is.null(ref.data.path),
      system.file("extdata/who_ht_maxvel_3sd.csv", package = "growthcleanr"),
      paste0(ref.data.path, "who_ht_maxvel_3sd.csv")
    )),
    by = c("sex", "whoagegrp_ht")
  )

  setnames(who.ht.vel, colnames(who.ht.vel), gsub("_", ".", colnames(who.ht.vel)))
  setkeyv(who.ht.vel, c("sex", "whoagegrp.ht"))
  # keep track of column names in the who growth velocity data
  who.fields <- colnames(who.ht.vel)
  who.fields <- who.fields[!who.fields %in% c("sex", "whoagegrp.ht")]

  # 1.  General principles
  # a.	All steps are done separately for each parameter unless otherwise noted
  # b.	All steps are done sorted by subject, parameter, and age (in days) for nonexcluded
  #     and nonmissing values only unless otherwise noted. This is very important.
  #     Sorting needs to be redone with each step to account for excluded and transformed values.
  # c.	The next value refers to the value with the next highest age for the same parameter
  #     and the same subject, and the previous value refers to the value with the
  #     next lowest age for the same parameter and the same subject.
  # d.	You will need to set up a method for keeping track of whether a value is missing or
  #     excluded (and in what step). I use variables called exc_* that are =0
  #     if a value is to be included, =1 if missing, and =2 or higher if it is to be excluded,
  #     with each number indicating a different step. I also set up
  #     parameter-specific subjid_* variables that are = to subjid for included values and are
  #     blank if the value is missing or should be excluded. These subjid_*
  #     variables need to be updated with each step.
  # e.  All steps assume that data are sorted by subjid_*, parameter, and age (in days)
  #     for nonexcluded and nonmissing values only
  #     unless otherwise noted. Sorting needs to be redone after any transformations or
  #     exclusions to account for excluded and
  #     transformed values.
  # f.  The next value refers to the nonexcluded nonmissing value with the next highest
  #     age for the same parameter and the same
  #     subject, and the previous value refers to the nonexcluded nonmissing value with
  #     the next lowest age for the same parameter
  #     and the same subject.
  # g.  exc_* should only be replaced with a  higher value if exc_*==0 at the time of
  #     replacement, unless otherwise specified.


  # NOTE: in the R code below exclusion is documented as a series of factor levels,
  # where all levels occuring before 'Exclude' in the sequence are considered
  # to be valid measurements.  We use the built in sorting of the data.table object
  # and subsets rather than re-sorting at each step
  # to ensure that only valid measurements are used at the beginning of each step.
  # Also, unlike the Stata code, the measurement parameter (weight vs. height)
  # is recorded as a factor in the data frame, rather than as a variable name

  # 2.  Data set-up
  # a.	I always code sex as 0=Male, 1=Female, so I recoded the variable sex that way
  #     and left a variable sexorigcode the way the data was sent to me (1=Female 2=Male)
  # b.	Remove rows that are duplicates for subjid, param, and measurement from further analysis
  #     NOTE: this step is not needed -- handled automatically by "temporary duplicate" step.
  # c.  I generated separate variables for weight (wt) and height (ht), as well as exc_*
  #     and subjid_* variables. Set exc_*=0 if value is not missing
  #     and exc_*=1 if value is missing. In all future steps, exc_* should only be changed
  #     if it is 0. This helps to keep track of which step excluded a value.
  #     I also kept the measurement variable there and untouched because sometimes wt
  #     and ht got transformed to something else.
  # d.	I made tables based on CDC growth curve parameters that include data for each day
  #     that I will send separately. The LMS parameters for each day are
  #     cubically interpolated from the values by month available on the CDC website.
  #     Create wt and ht z-scores for each value of each parameter (Z: WtZ and HtZ).
  # e.	There are variables in the table labelled cdc_*_csd_pos and cdc_*_csd_neg.
  #     For each age and sex, these correspond to ½ of the absolute value of the
  #     median and the value with a z-score of +2 (csd_pos) and -2 (csd_neg).
  #     These can be created to generate a score similar to the z-score but with an
  #     important difference. The z-scores created using the LMS method account for skewness
  #     in the distribution of the parameters (particularly weight), which
  #     can lead to small changes in z-score with large changes in weight in subjects
  #     with very high weight, and relatively large changes in z-score for smaller
  #     changes in weight in subjects with low weights.
  #     The score we will create can be called an SD-score (SDorig: WtSDorig and HtSDorig
  #     that is calculated by
  #     dividing the difference between the value and the median by the SD score
  #     (use csd_pos if the value is above the median, csd_neg if the value is below the
  #     median). These SD-scores, rather than z-scores, now form the basis for the algorithm.

  # calculate z scores
  if (!quietly) {
    cat(sprintf("[%s] Calculating z-scores...\n", Sys.time()))
  }
  measurement.to.z <- read_anthro(ref.data.path, cdc.only = TRUE)
  data.all[, z.orig := measurement.to.z(param, agedays, sex, v)]

  # calculate "standard deviation" scores
  if (!quietly) {
    cat(sprintf("[%s] Calculating SD-scores...\n", Sys.time()))
  }
  data.all[, sd.orig := measurement.to.z(param, agedays, sex, v, TRUE)]

  # sort by subjid, param, agedays
  setkeyv(data.all, c("subjid", "param", "agedays"))

  # add a new convenience index for bookkeeping
  data.all[, index := 1:.N]

  # enumerate the different exclusion levels
  exclude.levels <- c(
    "Include",
    "Unit-Error-High",
    "Unit-Error-Low",
    "Unit-Error-Possible",
    "Swapped-Measurements",
    "Exclude",
    "Missing",
    "Exclude-Temporary-Duplicate",
    "Exclude-Carried-Forward",
    "Exclude-SD-Cutoff",
    "Exclude-EWMA-Extreme",
    "Exclude-EWMA-Extreme-Pair",
    "Exclude-Duplicate",
    "Exclude-EWMA-8",
    "Exclude-EWMA-9",
    "Exclude-EWMA-10",
    "Exclude-EWMA-11",
    "Exclude-EWMA-12",
    "Exclude-EWMA-13",
    "Exclude-EWMA-14",
    "Exclude-Min-Height-Change",
    "Exclude-Max-Height-Change",
    "Exclude-Pair-Delta-17",
    "Exclude-Pair-Delta-18",
    "Exclude-Pair-Delta-19",
    "Exclude-Single-Outlier",
    "Exclude-Too-Many-Errors",
    "Exclude-Too-Many-Errors-Other-Parameter"
  )

  # Mark missing values for exclusion
  data.all[,
   exclude := factor(with(data.all, ifelse(
        is.na(v) |
          agedays < 0, "Missing", "Include"
      )),
      levels = exclude.levels,
      ordered = TRUE
    )
  ]

  # after calculating z scores, for convenience, recategorize linear parameters as 'HEIGHTCM'
  data.all[param == "LENGTHCM", param := "HEIGHTCM"]

  # define field names needed by helper functions
  ewma.fields <- c("ewma.all", "ewma.before", "ewma.after")

  # 3.  SD-score recentering: Because the basis of the method is comparing SD-scores over time,
  #     we need to account for the fact that
  #     the mean SD-score for the population changes with age.
  # a.	Determine the median cdc*sd for each parameter by year of age (with sexes combined): median*sd.
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
  # f.	In future steps I will sometimes refer to measprev
  #     and measnext which refer to the previous or next wt or ht measurement
  #     for which exc_*==0 for the subject and parameter, when the data are sorted by subject,
  #     parameter, and agedays.
  #     SDprev and SDnext refer to the tbc*sd of the previous or next measurement.

  if (!quietly) cat(sprintf("[%s] Re-centering data...\n", Sys.time()))

  # see function definition below for explanation of the re-centering process
  # returns a data table indexed by param, sex, agedays
  if (!is.data.table(sd.recenter)) {
    sd.recenter <- data.all[exclude < "Exclude", sd_median(param, sex, agedays, sd.orig)]
    if (sdmedian.filename != "") {
      fwrite(x = sd.recenter, file = sdmedian.filename, row.names = FALSE)
      if (!quietly) {
        cat(
          sprintf(
            "[%s] Wrote re-centering medians to %s...\n",
            Sys.time(),
            sdmedian.filename
          )
        )
      }
    }
  } else {
    # ensure passed-in medians are sorted correctly
    setkeyv(sd.recenter, c("param", "sex", "agedays"))
  }
  # add sd.recenter to data, and recenter
  setkeyv(data.all, c("param", "sex", "agedays"))
  data.all <- sd.recenter[data.all]
  setkeyv(data.all, c("subjid", "param", "agedays"))
  data.all[, tbc.sd := sd.orig - sd.median]

  if (sdrecentered.filename != "") {
    fwrite(x = data.all, file = sdrecentered.filename, row.names = FALSE)
    if (!quietly) {
      cat(
        sprintf(
          "[%s] Wrote re-centered data to %s...\n",
          Sys.time(),
          sdrecentered.filename
        )
      )
    }
  }

  # safety check: treat observations where tbc.sd cannot be calculated as missing
  data.all[is.na(tbc.sd), exclude := "Missing"]

  # NOTE: the rest of cleangrowth's steps are done through cleanbatch().

  # optionally process batches in parallel
  if (!quietly) {
    cat(sprintf(
      "[%s] Cleaning growth data in %d batch(es)...\n",
      Sys.time(),
      num.batches
    ))
  }
  if (num.batches == 1) {
    ret.df <- cleanbatch(data.all,
      log.path = log.path,
      quietly = quietly,
      parallel = parallel,
      measurement.to.z = measurement.to.z,
      ewma.fields = ewma.fields,
      ewma.exp = ewma.exp,
      recover.unit.error = recover.unit.error,
      include.carryforward = include.carryforward,
      sd.extreme = sd.extreme,
      z.extreme = z.extreme,
      exclude.levels = exclude.levels,
      tanner.ht.vel = tanner.ht.vel,
      who.ht.vel = who.ht.vel,
      lt3.exclude.mode = lt3.exclude.mode,
      error.load.threshold = error.load.threshold,
      error.load.mincount = error.load.mincount
    )
  } else {
    # create log directory if necessary
    if (!quietly) {
      cat(sprintf("[%s] Writing batch logs to '%s'...\n", Sys.time(), log.path))
    }
    dir.create(log.path, showWarnings = FALSE)

    ret.df <- plyr::ddply(
      .data = data.all,
      .variables = "batch",
      .fun = cleanbatch,
      .parallel = parallel,
      .paropts = list(.packages = "data.table"),
      log.path = log.path,
      quietly = quietly,
      parallel = parallel,
      measurement.to.z = measurement.to.z,
      ewma.fields = ewma.fields,
      ewma.exp = ewma.exp,
      recover.unit.error = recover.unit.error,
      include.carryforward = include.carryforward,
      sd.extreme = sd.extreme,
      z.extreme = z.extreme,
      exclude.levels = exclude.levels,
      tanner.ht.vel = tanner.ht.vel,
      who.ht.vel = who.ht.vel,
      lt3.exclude.mode = lt3.exclude.mode,
      error.load.threshold = error.load.threshold,
      error.load.mincount = error.load.mincount
    )
    stopImplicitCluster()
  }
  if (!quietly) cat(sprintf("[%s] Done!\n", Sys.time()))

  return(ret.df$exclude[order(ret.df$line)])
}
