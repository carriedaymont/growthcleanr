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

# Main growthcleanr algorithm and other exported functions
# Main growthcleanr algorithm (cleangrowth): adult and pediatric are in *_clean.R
#(main bulk of algorithm) and *_support.R (supporting functions)

#' Clean growth measurements
#'
#' @param subjid Vector of unique identifiers for each subject in the database.
#' @param param Vector identifying each measurement, may be 'WEIGHTKG', 'WEIGHTLBS', 'HEIGHTCM', 'HEIGHTIN', 'LENGTHCM', or 'HEADCM'.
#'   'HEIGHTCM'/'HEIGHTIN' vs. 'LENGTHCM' only affects z-score calculations between ages 24 to 35 months (730 to 1095 days).
#'   All linear measurements below 731 days of life (age 0-23 months) are interpreted as supine length, and
#'   all linear measurements above 1095 days of life (age 36+ months) are interpreted as standing height.
#'   Note: at the moment, all LENGTHCM will be converted to HEIGHTCM. In the future, the algorithm will be updated to consider this difference.
#'   Additionally, imperial 'HEIGHTIN' and 'WEIGHTLBS' measurements are converted to
#'   metric during algorithm calculations.
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
#' @param height.tolerance.cm maximum decrease in height tolerated for sequential measurements
#' @param error.load.mincount minimum count of exclusions on parameter before
#' considering excluding all measurements. Defaults to 2.
#' @param error.load.threshold threshold of percentage of excluded measurement count to included measurement
#' count that must be exceeded before excluding all measurements of either parameter. Defaults to 0.5.
#' @param sd.recenter specifies how to recenter medians. May be a data frame or
#' table w/median SD-scores per day of life by gender and parameter, or "NHANES"
#' or "derive" as a character vector.
#' \itemize{
#'   \item If `sd.recenter` is specified as a data set, use the data set
#'   \item If `sd.recenter` is specified as "`nhanes`", use NHANES reference medians
#'   \item If `sd.recenter` is specified as "`derive`", derive from input
#'   \item If `sd.recenter` is not specified or `NA`:
#'     \itemize{
#'       \item If the input set has at least 5,000 observations, derive medians from input
#'       \item If the input set has fewer than 5,000 observations, use NHANES
#'     }
#' }
#'
#' If specifying a data set, columns must include param, sex, agedays, and sd.median
#' (referred to elsewhere as "modified Z-score"), and those medians will be used
#' for recentering. A summary of how the NHANES reference medians were derived is
#' available in README.md. Defaults to NA.
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
#' @param log.path Path to log file output when running in parallel (non-quiet mode). Default is NA. A new
#' directory will be created if necessary. Set to NA to disable log files.
#' @param parallel Determines if function runs in parallel.  Defaults to FALSE.
#' @param num.batches Specify the number of batches to run in parallel. Only
#' applies if parallel is set to TRUE. Defaults to the number of workers
#' returned by the getDoParWorkers function in the foreach package.
#' @param quietly Determines if function messages are to be displayed and if log files (parallel only) are to be generated.
#' Defaults to TRUE
#' @param adult_cutpoint Number between 18 and 20, describing ages when the
#'  pediatric algorithm should not be applied (< adult_cutpoint), and the adult
#'   algorithm should apply (>= adult_cutpoint). Numbers outside this range will be
#'   changed to the closest number within the range. Defaults to 20.
#' @param weight_cap Positive number, describing a weight cap in kg (rounded to the
#' nearest .1, +/- .1) within the adult dataset. If there is no weight cap, set
#'  to Inf. Defaults to Inf.
#' @param adult_columns_filename Name of file to save original adult data, with additional output columns to
#' as CSV. Defaults to "", for which this data will not be saved. Useful
#' for post-analysis. For more information on this output, please see README.
#' @param prelim_infants TRUE/FALSE. Run the in-development release of the infants algorithm (expands pediatric algorithm to improve performance for children 0 – 2 years). Not recommended for use in research. For more information regarding the logic of the algorithm, see the vignette 'Preliminary Infants Algorithm.' Defaults to FALSE.
#'
#' @return Vector of exclusion codes for each of the input measurements.
#'
#'    Possible values for each code are:
#'
#' * 'Include', 'Unit-Error-High', 'Unit-Error-Low', 'Swapped-Measurements', 'Missing',
#' *  'Exclude-Carried-Forward', 'Exclude-SD-Cutoff', 'Exclude-EWMA-Extreme', 'Exclude-EWMA-Extreme-Pair',
#' *  'Exclude-Extraneous-Same-Day',
#' *  'Exclude-EWMA-8', 'Exclude-EWMA-9', 'Exclude-EWMA-10', 'Exclude-EWMA-11', 'Exclude-EWMA-12', 'Exclude-EWMA-13', 'Exclude-EWMA-14',
#' *  'Exclude-Min-Height-Change', 'Exclude-Max-Height-Change',
#' *  'Exclude-Pair-Delta-17', 'Exclude-Pair-Delta-18', 'Exclude-Pair-Delta-19',
#' *  'Exclude-Single-Outlier', 'Exclude-Too-Many-Errors', 'Exclude-Too-Many-Errors-Other-Parameter'
#' @md
#'
#' @export
#' @import data.table
#' @rawNamespace import(plyr, except = c(failwith, id, summarize, count, desc, mutate, arrange, rename, is.discrete, summarise, summarize))
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom utils write.csv
#' @rawNamespace import(R.utils, except = c(extract))
#' @examples
#' \donttest{
#' # Run calculation using a small subset of given data
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid %in% unique(df_stats[, "subjid"])[1:5], ]
#'
#' clean_stats <-cleangrowth(subjid = df_stats$subjid,
#'                          param = df_stats$param,
#'                          agedays = df_stats$agedays,
#'                          sex = df_stats$sex,
#'                          measurement = df_stats$measurement)
#'
#' # Once processed you can filter data based on result value
#' df_stats <- cbind(df_stats, "clean_result" = clean_stats)
#' clean_df_stats <- df_stats[df_stats$clean_result == "Include",]
#'
#' # Parallel processing: run using 2 cores and batches
#' clean_stats <- cleangrowth(subjid = df_stats$subjid,
#'                            param = df_stats$param,
#'                            agedays = df_stats$agedays,
#'                            sex = df_stats$sex,
#'                            measurement = df_stats$measurement,
#'                            parallel = TRUE,
#'                            num.batches = 2)
#' }
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
                        ref.data.path = "",
                        log.path = NA,
                        parallel = FALSE,
                        num.batches = NA,
                        quietly = TRUE,
                        adult_cutpoint = 20,
                        weight_cap = Inf,
                        adult_columns_filename = "",
                        prelim_infants = FALSE,
                        id = NULL) {
  # avoid "no visible binding" warnings
  N <- age_years <- batch <- exclude <- index <- line <- NULL
  newbatch <- sd.median <- sd.orig <- tanner.months <- tbc.sd <- NULL
  v <- v_adult <- whoagegrp.ht <- whoagegrp_ht <- z.orig <- NULL
  z.orig_cdc <- z.orig_who <- sd.orig_cdc <- sd.orig_who <- NULL
  result <- NULL
  
  sd.orig_uncorr <- agemonths <- intwt <- fengadays <- pmagedays <- cagedays <-
    unmod_zscore <- fen_wt_m <- fen_wt_l <- fen_wt_s <- cwho_cv <- ccdc_cv <-
    sd.c_cdc <- sd.c_who <- sd.c <- sd.corr <- seq_win <- sd.corr_abssumdiff <-
    sd.orig_abssumdiff <- orig_colnames <- ctbc.sd <- sum_sde <- no_sde <-
    sum_val <- no_dup_val <- no_outliers <- no_bigdiff <- nottoofar <- nnte <-
    nnte_full <- NULL
  
  # preprocessing ----
  
  # organize data into a dataframe along with a line "index" so the original data order can be recovered
  data.all.ages <- data.table(
    line = seq_along(measurement),
    subjid = as.factor(subjid),
    param,
    agedays = as.integer(agedays),
    v = ifelse(measurement == 0, NaN, measurement),
    v_adult = measurement,
    sex = as.integer(ifelse(
      sex %in% c(0, 'm', 'M'), 0, ifelse(sex %in% c(1, 'f', 'F'), 1, NA)
    ))
  )
  
  # CP ADD -
  if (!is.null(id)) data.all.ages[, id := id]
  
  # CP ADD UP
  
  # quality checks
  if (!is.numeric(adult_cutpoint)){
    stop("adult_cutpoint not numeric. Please enter a number between 18 and 20.")
  }
  if (!is.numeric(weight_cap) | weight_cap < 0){
    stop("weight_cap not numeric. Please enter a positive number.")
  }
  if (any(!param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHIN",
                        "WEIGHTLBS", "HEADCM"))){
    cat(sprintf("[%s] Parameters included that do not match 'param' specifications. Marking as missing...\n", Sys.time()))
    data.all.ages <-
      data.all.ages[
        !param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHIN",
                      "WEIGHTLBS", "HEADCM"), v := NA]
    data.all.ages <-
      data.all.ages[
        !param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHIN",
                      "WEIGHTLBS", "HEADCM"), v_adult := NA]
  }
  
  # rate limit cutpoint -- min 18, max 20
  cutpoint_update <-
    if (adult_cutpoint < 18){
      18
    } else if (adult_cutpoint > 20){
      20
    } else {
      adult_cutpoint
    }
  
  # split by cutpoint
  # for ease, data.all will refer to pediatric data; data.adult will refer to
  # adult data -- copy to make sure they're separate
  data.all <- copy(data.all.ages[agedays < adult_cutpoint*365.25,])
  data.adult <- copy(data.all.ages[agedays >= adult_cutpoint*365.25,])
  
  # TODO: ADD PARALLEL FOR ADULTS
  
  # constants for pediatric
  # enumerate the different exclusion levels
  if (prelim_infants){
    # different for infants
    exclude.levels.peds <- c(
      'Include',
      'Unit-Error-High',
      'Unit-Error-Low',
      'Unit-Error-Possible',
      'Swapped-Measurements',
      'Exclude',
      'Missing',
      'Not cleaned',
      'Exclude-Temporary-Extraneous-Same-Day',
      'Exclude-Carried-Forward',
      # added CF exclusions
      "Exclude-1-CF-deltaZ-<0.05",
      "Exclude-1-CF-deltaZ-<0.1-wholehalfimp",
      "Exclude-Teen-2-plus-CF-deltaZ-<0.05",
      "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp",
      'Exclude-EWMA-Extreme',
      'Exclude-EWMA-Extreme-Pair',
      'Exclude-SDE-Identical',
      'Exclude-SDE-All-Exclude',
      'Exclude-SDE-All-Extreme',
      'Exclude-SDE-EWMA',
      'Exclude-SDE-One-Day',
      "Exclude-EWMA2-middle",
      "Exclude-EWMA2-birth-WT",
      "Exclude-EWMA2-birth-WT-ext",
      "Exclude-EWMA2-first",
      "Exclude-EWMA2-first-ext",
      "Exclude-EWMA2-last",
      "Exclude-EWMA2-last-high",
      "Exclude-EWMA2-last-ext",
      "Exclude-EWMA2-last-ext-high",
      "Exclude-EWMA2-birth-HT-HC",
      "Exclude-EWMA2-birth-HT-HC-ext",
      "Exclude-Min-diff",
      "Exclude-Max-diff",
      "Exclude-2-meas->1-year",
      "Exclude-2-meas-<1-year",
      "Exclude-1-meas",
      "Exclude-Error-load",
      
      # old
      
      'Exclude-Extraneous-Same-Day',
      'Exclude-SD-Cutoff',
      'Exclude-EWMA-8',
      'Exclude-EWMA-9',
      'Exclude-EWMA-10',
      'Exclude-EWMA-11',
      'Exclude-EWMA-12',
      'Exclude-EWMA-13',
      'Exclude-EWMA-14',
      'Exclude-Min-Height-Change',
      'Exclude-Max-Height-Change',
      'Exclude-Pair-Delta-17',
      'Exclude-Pair-Delta-18',
      'Exclude-Pair-Delta-19',
      'Exclude-Single-Outlier',
      'Exclude-Too-Many-Errors',
      'Exclude-Too-Many-Errors-Other-Parameter',
      
      #new
      "Exclude-Absolute-BIV",
      "Exclude-Standardized-BIV",
      "Exclude-Evil-Twins",
      "Exclude-EWMA1-Extreme"
    )
  } else {
    exclude.levels.peds <- c(
      'Include',
      'Unit-Error-High',
      'Unit-Error-Low',
      'Unit-Error-Possible',
      'Swapped-Measurements',
      'Exclude',
      'Missing',
      'Exclude-Temporary-Extraneous-Same-Day',
      'Exclude-Carried-Forward',
      'Exclude-SD-Cutoff',
      'Exclude-EWMA-Extreme',
      'Exclude-EWMA-Extreme-Pair',
      'Exclude-Extraneous-Same-Day',
      'Exclude-EWMA-8',
      'Exclude-EWMA-9',
      'Exclude-EWMA-10',
      'Exclude-EWMA-11',
      'Exclude-EWMA-12',
      'Exclude-EWMA-13',
      'Exclude-EWMA-14',
      'Exclude-Min-Height-Change',
      'Exclude-Max-Height-Change',
      'Exclude-Pair-Delta-17',
      'Exclude-Pair-Delta-18',
      'Exclude-Pair-Delta-19',
      'Exclude-Single-Outlier',
      'Exclude-Too-Many-Errors',
      'Exclude-Too-Many-Errors-Other-Parameter'
    )
  }
  
  exclude.levels.adult <- c(
    "Include",
    "Exclude-Adult-BIV",
    "Exclude-Adult-Hundreds",
    "Exclude-Adult-Hundreds-RV",
    "Exclude-Adult-Unit-Errors",
    "Exclude-Adult-Unit-Errors-RV",
    "Exclude-Adult-Transpositions",
    "Exclude-Adult-Transpositions-RV",
    "Exclude-Adult-Weight-Cap-Identical",
    "Exclude-Adult-Weight-Cap",
    "Exclude-Adult-Swapped-Measurements",
    "Exclude-Adult-Identical-Same-Day",
    "Exclude-Adult-Extraneous-Same-Day",
    "Exclude-Adult-Distinct-Pairs",
    "Exclude-Adult-Distinct-3-Or-More",
    "Exclude-Adult-EWMA-Extreme",
    "Exclude-Adult-EWMA-Extreme-RV",
    "Exclude-Adult-Distinct-Ordered-Pairs",
    "Exclude-Adult-EWMA-Moderate",
    "Exclude-Adult-Possibly-Impacted-By-Weight-Cap",
    "Exclude-Adult-Distinct-Single",
    "Exclude-Adult-Too-Many-Errors"
  )
  
  exclude.levels <- base::union(exclude.levels.peds, exclude.levels.adult)
  
  # if there's no pediatric data, no need to go through this rigamarole
  if (nrow(data.all) > 0){
    
    # pediatric: height velocity calculations and preprocessing ----
    
    # for pediatric data, convert in and lbs to cm and kg (adult is done within algo)
    data.all[param == "HEIGHTIN", v := v*2.54]
    data.all[param == "HEIGHTIN", param := "HEIGHTCM"]
    data.all[param == "WEIGHTLBS", v := v/2.2046226]
    data.all[param == "WEIGHTLBS", param := "WEIGHTKG"]
    
    # if parallel processing is desired, load additional modules
    if (parallel) {
      if (is.na(num.batches)) {
        num.batches <- getDoParWorkers()
      }
      if (prelim_infants){
        # variables needed for parallel workers
        var_for_par <- c("temporary_extraneous", "valid", "swap_parameters",
                         "na_as_false", "ewma", "read_anthro", "as_matrix_delta",
                         "sd_median",
                         
                         "temporary_extraneous_infants",
                         "get_dop", "calc_oob_evil_twins",
                         "calc_and_recenter_z_scores")
      } else {
        var_for_par <- c("temporary_extraneous", "valid", "swap_parameters",
                         "na_as_false", "ewma", "read_anthro", "as_matrix_delta",
                         "sd_median")
      }
      
      cl <- makeCluster(num.batches)
      clusterExport(cl = cl, varlist = var_for_par, envir = environment())
      registerDoParallel(cl)
    } else {
      if (is.na(num.batches))
        num.batches <- 1
    }
    
    setkey(data.all, subjid)
    subjid.unique <- data.all[j = unique(subjid)]
    batches.all <- data.table(
      subjid = subjid.unique,
      batch = sample(num.batches, length(subjid.unique), replace = TRUE),
      key = 'subjid'
    )
    data.all <- batches.all[data.all]
    
    if (!quietly){
      cat(sprintf("[%s] Begin processing pediatric data...\n", Sys.time()))
    }
    
    # load tanner height velocity data. sex variable is defined such that 0=male and 1=female
    # recode column names to match syntactic style ("." rather than "_" in variable names)
    tanner_ht_vel_path <- ifelse(
      ref.data.path == "",
      system.file(file.path("extdata", "tanner_ht_vel.csv.gz"), package = "growthcleanr"),
      file.path(ref.data.path, "tanner_ht_vel.csv.gz")
    )
    
    tanner.ht.vel <- fread(tanner_ht_vel_path)
    
    setnames(tanner.ht.vel,
             colnames(tanner.ht.vel),
             gsub('_', '.', colnames(tanner.ht.vel)))
    setkey(tanner.ht.vel, sex, tanner.months)
    # keep track of column names in the tanner data
    tanner.fields <- colnames(tanner.ht.vel)
    tanner.fields <- tanner.fields[!tanner.fields %in% c('sex', 'tanner.months')]
    
    if (!prelim_infants){
      who_max_ht_vel_path <- ifelse(
        ref.data.path == "",
        system.file(file.path("extdata", "who_ht_maxvel_3sd.csv.gz"), package = "growthcleanr"),
        file.path(ref.data.path, "who_ht_maxvel_3sd.csv.gz")
      )
      
      who_ht_vel_3sd_path <- ifelse(
        ref.data.path == "",
        system.file(file.path("extdata", "who_ht_vel_3sd.csv.gz"), package = "growthcleanr"),
        file.path(ref.data.path, "who_ht_vel_3sd.csv.gz")
      )
    } else {
      who_max_ht_vel_path <- ifelse(
        ref.data.path == "",
        system.file(file.path("extdata", "who_hc_maxvel_3sd_infants.csv.gz"), package = "growthcleanr"),
        file.path(ref.data.path, "who_hc_maxvel_3sd_infants.csv.gz")
      )
      
      who_ht_vel_3sd_path <- ifelse(
        ref.data.path == "",
        system.file(file.path("extdata", "who_hc_vel_3sd_infants.csv.gz"), package = "growthcleanr"),
        file.path(ref.data.path, "who_hc_vel_3sd_infants.csv.gz")
      )
    }
    who.max.ht.vel <- fread(who_max_ht_vel_path)
    who.ht.vel <- fread(who_ht_vel_3sd_path)
    setkey(who.max.ht.vel, sex, whoagegrp_ht)
    setkey(who.ht.vel, sex, whoagegrp_ht)
    who.ht.vel <- as.data.table(dplyr::full_join(who.ht.vel, who.max.ht.vel, by =
                                                   c('sex', 'whoagegrp_ht')))
    
    setnames(who.ht.vel, colnames(who.ht.vel), gsub('_', '.', colnames(who.ht.vel)))
    setkey(who.ht.vel, sex, whoagegrp.ht)
    # keep track of column names in the who growth velocity data
    who.fields <- colnames(who.ht.vel)
    who.fields <- who.fields[!who.fields %in% c('sex', 'whoagegrp.ht')]
    
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
    # b.	Remove rows that are extraneous for subjid, param, and measurement from further analysis
    #     NOTE: this step is not needed -- handled automatically by "temporary extraneous" step.
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
    
    # recategorize linear parameters as 'HEIGHTCM'
    # NOTE: this will be changed in future to consider this difference
    data.all[param == 'LENGTHCM', param := 'HEIGHTCM']
    
    # calculate z/sd scores
    if(prelim_infants){
      if (!quietly)
        cat(sprintf("[%s] Calculating z-scores...\n", Sys.time()))
      # removing z calculations, as they are not used
      # for infants, use z and who
      measurement.to.z <- read_anthro(ref.data.path, cdc.only = TRUE,
                                      prelim_infants = TRUE)
      measurement.to.z_who <- read_anthro(ref.data.path, cdc.only = FALSE,
                                          prelim_infants = TRUE)
      
      # calculate "standard deviation" scores
      if (!quietly)
        cat(sprintf("[%s] Calculating SD-scores...\n", Sys.time()))
      data.all[, sd.orig_cdc := measurement.to.z(param, agedays, sex, v, TRUE)]
      data.all[, sd.orig_who := measurement.to.z_who(param, agedays, sex, v, TRUE)]
      
      # smooth z-scores/SD scores between ages 1 - 3yo using weighted scores
      # older uses cdc, younger uses who
      data.all$ageyears <- data.all$agedays/365.25
      
      who_weight <- 4 - (data.all$ageyears)
      cdc_weight <- (data.all$ageyears) - 2
      
      smooth_val <- data.all$ageyears >= 2 &
        data.all$ageyears <= 4 &
        data.all$param != "HEADCM"
      data.all[smooth_val,
               sd.orig := (data.all$sd.orig_cdc[smooth_val]*cdc_weight[smooth_val] +
                             data.all$sd.orig_who[smooth_val]*who_weight[smooth_val])/2]
      # <- 
      # otherwise use WHO and CDC for older and younger, respectively
      who_val <- data.all$param == "HEADCM" |
        data.all$ageyears < 2
      
      data.all[(who_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_cdc)),
               sd.orig := data.all$sd.orig_who[(who_val & !smooth_val)  | (smooth_val & is.na(data.all$sd.orig_cdc))]]
      
      cdc_val <- data.all$param != "HEADCM" &
        data.all$ageyears > 4
      
      data.all[(cdc_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_who)),
               sd.orig := data.all$sd.orig_cdc[(cdc_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_who))]]
      
      # NOTE: SD SCORES IN CODE ARE Z IN INFANT DOCS -- USE sd.orig ONLY
      
      # keep the original, uncorrected, unrecentered zscores
      data.all[,sd.orig_uncorr := sd.orig]
      
      # NOTE: MAY WANT TO SUBSET HERE
      
      # 2b: corrected z scores ----
      
      # keep the original column names -- we're adding a ton of columns that we
      # want to filter out after correction
      orig_colnames <- copy(colnames(data.all))
      
      # start by reading in fenton data
      fentlms_foraga <- fread(
        system.file(file.path("extdata", "fentlms_foraga.csv.gz"),
                    package = "growthcleanr"))
      fentlms_forz <- fread(
        system.file(file.path("extdata", "fentlms_forz.csv.gz"),
                    package = "growthcleanr"))
      ##### CP MOD DOWN #####
      
      # add age in months
      data.all[, agemonths := agedays / 30.4375]
      
      # initialize the column so it's always present
      data.all[, potcorr := FALSE]
      
      # mark infants needing correction: first weight per subject, Z<-2, age<10mo
      data.all[param == "WEIGHTKG",
               potcorr := (seq_len(.N) == 1L & sd.orig < -2 & agemonths < 10),
               by = subjid]
      
      # propagate that flag to all rows for that subject
      data.all[, potcorr := any(potcorr), by = subjid]
      
      # ensure it sticks around through merges
      setkey(data.all, subjid, param, agedays)
      

      ### CP MODIFY ###
      # integer weight is in grams, rounded to the nearest 10
      data.all[potcorr == TRUE, intwt := trunc(v*100)*10]
      # data.all[potcorr == TRUE, intwt := round(v*10)]
      write.csv(fentlms_foraga, "../data/fentlms_foraga.csv")
      # replace to facilitate merging with fenton curves
      data.all[intwt >= 250 & intwt <=560, intwt := 570]
      
      # merge with fenton curves
      data.all <- merge(
        data.all, fentlms_foraga, by = c("sex", "intwt"),
        all.x = TRUE)
      
      data.all[fengadays < 259, pmagedays := agedays + fengadays]
      ### CP BELOW
      data.all[fengadays >= 259 | is.na(fengadays), pmagedays := NA_real_]
      ### CP ABOVE
      
      data.all[fengadays < 259, cagedays := pmagedays - 280]
      # replace fengadays with pmagedays to facilitate merging
      data.all[, fengadays := pmagedays]
      

      
      # merge with fenton curves
      data.all <- merge(
        data.all, fentlms_forz, by = c("sex", "fengadays"),
        all.x = TRUE)
      
      print(with(data.all[potcorr == TRUE], table(is.na(fen_wt_m))))
      str(data.all[, .(sex, intwt, fengadays, potcorr)][potcorr == TRUE])
      print(min(data.all[,intwt], na.rm = TRUE))
      print(min(data.all[,intwt], na.rm = TRUE))
      str(fentlms_foraga)
      str(fentlms_forz)
      
      # add unmodified zscore using weight in unrounded grams
      data.all[, unmod_zscore :=
                 ((v*1000/fen_wt_m)^fen_wt_l - 1)/(fen_wt_l * fen_wt_s)]
      
      
      ### CP MPD ABOVE ### 
      ### CP MOD BELOW ###
      # 
      # data.all[potcorr == TRUE & !is.na(unmod_zscore),
      #          sd.corr := unmod_zscore]
      # data.all[is.na(sd.corr), sd.corr := sd.orig]
      # 
      
      
      # --- Assign Fenton- and WHO-corrected z-scores based on postmenstrual age --- #
      # Use Fenton z-score (unmod_zscore) if < 350 postmenstrual days
      data.all[potcorr == TRUE & !is.na(unmod_zscore),
               sd.corr := unmod_zscore]
      
  
      
      # For everyone else (no Fenton or WHO correction available), fall back to original z
      data.all[is.na(sd.corr), sd.corr := sd.orig]
      
      ### CP MOD BELOW ###
      ### ALL GOOD UP TO HERE ON CORRECTED Z SCORE
      # read in who/cdc data
      growthfile_who <- fread(
        system.file(file.path("extdata", "growthfile_who.csv.gz"),
                    package = "growthcleanr"))
      growthfile_cdc <- fread(
        system.file(file.path("extdata", "growthfile_cdc_ext_infants.csv.gz"),
                    package = "growthcleanr"))
      
      # merge in the who/cdc data
      data.all <- merge(data.all, growthfile_who,
                        by.x = c("sex", "cagedays"),
                        by.y = c("sex", "agedays"),
                        all.x = TRUE)
      data.all <- merge(data.all, growthfile_cdc,
                        by.x = c("sex", "cagedays"),
                        by.y = c("sex", "agedays"),
                        all.x = TRUE)
      
      # adjust WHO and CDC heights based on age
      data.all[, cwho_cv := v]
      data.all[, ccdc_cv := v]
      data.all[param == "HEIGHTCM" & agedays > 730 & cagedays <= 730,
               cwho_cv := cwho_cv + 0.8]
      data.all[param == "HEIGHTCM" & agedays > 730 & cagedays <= 730,
               ccdc_cv := ccdc_cv + 0.7]
      
      # create the corrected z scores
      # use read anthro, but pass in different arguments
      # pass in the corrected height
      data.all[, sd.c_cdc :=
                 measurement.to.z(param, cagedays, sex, ccdc_cv, TRUE)]
      data.all[, sd.c_who :=
                 measurement.to.z_who(param, cagedays, sex, cwho_cv, TRUE)]
      
      # Use WHO corrected z-score if pmagedays >= 350 and chronological age <= 730 days (~2 years)
      data.all[potcorr == TRUE & pmagedays >= 350 & agedays <= 730 & !is.na(sd.c_who),
               sd.corr := sd.c_who]
      
      table(data.all$potcorr, is.na(data.all$sd.corr))
      print("AAA")
      print(summary(data.all[potcorr == TRUE, sd.corr - sd.orig]))
      
      # smooth using weights as in original z score creation
      who_weight <- 4 - (data.all$agedays/365.25)
      cdc_weight <- (data.all$agedays/365.25) - 2
      
      smooth_val <- data.all$agedays/365.25 >= 2 &
        data.all$agedays/365.25 <= 4 &
        data.all$param != "HEADCM"
      data.all[smooth_val,
               sd.c := (sd.c_cdc[smooth_val]*cdc_weight[smooth_val] +
                          sd.c_who[smooth_val]*who_weight[smooth_val])/2]
      
      # otherwise use WHO and CDC for older and younger, respectively
      who_val <- data.all$param == "HEADCM" |
        data.all$agedays/365.25 < 2
      data.all[who_val | (smooth_val & is.na(data.all$sd.c_cdc)),
               sd.c := data.all$sd.c_who[who_val  | (smooth_val & is.na(data.all$sd.c_cdc))]]
      
      cdc_val <- data.all$param != "HEADCM" |
        data.all$agedays/365.25 > 4
      data.all[cdc_val | (smooth_val & is.na(data.all$sd.c_who)),
               sd.c := data.all$sd.c_cdc[cdc_val | (smooth_val & is.na(data.all$sd.c_who))]]
      ### CP REPLACCE BLOCK BELOW 
      # # smooth corrected and uncorrected z scores
      # uncorrweight <-  4 - (data.all$agedays/365.25)
      # corrweight <- (data.all$agedays/365.25) - 2
      # smooth_val <- data.all$agedays/365.25 >= 2 &
      #   data.all$agedays/365.25 <= 4
      # data.all[smooth_val,
      #          sd.corr := (sd.c[smooth_val]*corrweight[smooth_val] +
      #                        sd.orig[smooth_val]*uncorrweight[smooth_val])/2]
      # # for < 2 & potential correction, use fenton corrected score
      # data.all[agedays/365.25 <= 2 & potcorr == TRUE, sd.corr := sd.c]
      # # for > 4, use original z score
      # data.all[agedays/365.25 >= 4, sd.corr := sd.orig]
      # # for not potential corrections, use the original z score
      # data.all[!potcorr == TRUE, sd.corr := sd.orig]
      # # if the who/fenton score is not available for any reason, use the
      # # original
      # data.all[is.na(sd.corr), sd.corr := sd.orig]
      
      ### CP WORK BELOW
      
      # smooth corrected and uncorrected z scores (2–4y)
      uncorrweight <- 4 - (data.all$agedays/365.25)
      corrweight <- (data.all$agedays/365.25) - 2
      smooth_val <- data.all$agedays/365.25 >= 2 &
        data.all$agedays/365.25 <= 4

      data.all[smooth_val & !is.na(sd.c) & !is.na(sd.orig),
               sd.corr := (sd.orig[smooth_val]*uncorrweight[smooth_val] +
                             sd.c[smooth_val]*corrweight[smooth_val]) / 2]

      # <2y corrected (Fenton/WHO)
      data.all[agedays/365.25 <= 2 & potcorr == TRUE & !is.na(sd.c),
               sd.corr := sd.c]

      # >4y or not correctable
      data.all[agedays/365.25 >= 4 | potcorr == FALSE | is.na(sd.corr),
               sd.corr := sd.orig]
      
      ### CP WORK ABOVE
      
      # --- Sequential Stata-equivalent corrections ---
    
      
    
      
      ### CP BELOW 
      # examine_only <- data.all$param == "WEIGHTKG" &
      #   data.all$subjid %in% data.all$subjid[potcorr]
      examine_only <- data.all$param == "WEIGHTKG" & data.all$potcorr
      
      
      ### CP ABOVE
      tmp <- copy(data.all[examine_only,])
      tmp <- tmp[order(subjid, agedays),]
      tmp[, seq_win := sequence(.N), by = subjid]
      # we're only looking at the first 4 values, and they need to be < 2 years
      tmp <- tmp[seq_win <= 4 & (agedays/365.25) < 2,]
      # don't look at subjects where there is only 1 score and there is no
      # value for either
      tmp <- tmp[!(is.na(sd.corr) & is.na(sd.orig)),]
      tmp <- tmp[subjid %in% names(table(subjid) > 1),]
      
      # create differences, absolute sum them
      tmp[, sd.corr_abssumdiff := abs(sum(sd.corr[1] - sd.corr)), by = subjid]
      tmp[, sd.orig_abssumdiff := abs(sum(sd.orig[1] - sd.orig)), by = subjid]
      # find subjects where corrected value needs to be replaced
      
      #### CP MODIFY ####
      sub_replace <- unique(tmp[sd.corr_abssumdiff > sd.orig_abssumdiff,
                                subjid])
      #### CP MODIFY ####
      # replace accordingly in the main dataframe
      data.all[subjid %in% sub_replace, sd.corr := sd.orig]
      
      orig_colnames <- c(orig_colnames, "sd.corr")
      
      # remove many added columns
      write.csv(data.all, "../data/midpoint_check_corrected_zscores.csv")
      data.all <- data.all[, colnames(data.all) %in% c(orig_colnames, id),
                           with = FALSE]
      
    } else {
      # calculate z scores
      if (!quietly)
        cat(sprintf("[%s] Calculating z-scores...\n", Sys.time()))
      measurement.to.z <- read_anthro(ref.data.path, cdc.only = TRUE)
      data.all[, z.orig := measurement.to.z(param, agedays, sex, v)]
      
      # calculate "standard deviation" scores
      if (!quietly)
        cat(sprintf("[%s] Calculating SD-scores...\n", Sys.time()))
      data.all[, sd.orig := measurement.to.z(param, agedays, sex, v, TRUE)]
    }
    
    # sort by subjid, param, agedays
    setkey(data.all, subjid, param, agedays)
    
    # add a new convenience index for bookkeeping
    data.all[, index := 1:.N]
    
    # Mark missing values for exclusion
    data.all[, exclude := factor(with(data.all, ifelse(
      is.na(v) |
        agedays < 0, 'Missing', 'Include'
    )),
    levels = exclude.levels,
    ordered = TRUE)]
    # also mark certain measurements to not consider
    data.all[param == "HEADCM" & agedays > (3*365.25), exclude := "Not cleaned"]
    
    # define field names needed by helper functions
    ewma.fields <- c('ewma.all', 'ewma.before', 'ewma.after')
    
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
    # returns a data table indexed by param, sex, agedays. can use NHANES reference
    # data, derive from input, or use user-supplied data.
    if (!is.data.table(sd.recenter)) {
      # INFANTS CHANGES:
      # use recentering file derived from work, independent of sex
      if (prelim_infants){
        infants_reference_medians_path <- ifelse(
          ref.data.path == "",
          system.file(file.path("extdata",
                                "rcfile-2023-08-15_format.csv.gz"),
                      package = "growthcleanr"),
          file.path(ref.data.path, "rcfile-2023-08-15_format.csv.gz")
        )
        sd.recenter <- fread(infants_reference_medians_path)
        if (!quietly)
          cat(
            sprintf(
              "[%s] Using infants reference medians...\n",
              Sys.time()
            )
          )
      } else if ((is.character(sd.recenter) &
                  tolower(sd.recenter) == "nhanes") |
                 (!(is.character(sd.recenter) &
                    tolower(sd.recenter) == "derive") & (data.all[, .N] < 5000))) {
        
        # Use NHANES medians if the string "nhanes" is specified instead of a data.table
        # or if sd.recenter is not specified as "derive" and N < 5000.
        
        nhanes_reference_medians_path <- ifelse(
          ref.data.path == "",
          system.file(file.path("extdata", "nhanes-reference-medians.csv.gz"), package = "growthcleanr"),
          file.path(ref.data.path, "nhanes-reference-medians.csv.gz")
        )
        sd.recenter <- fread(nhanes_reference_medians_path)
        if (!quietly)
          cat(
            sprintf(
              "[%s] Using NHANES reference medians...\n",
              Sys.time()
            )
          )
      } else {
        # Derive medians from input data
        sd.recenter <- data.all[exclude < 'Exclude', sd_median(param, sex, agedays, sd.orig)]
        if (!quietly)
          cat(
            sprintf(
              "[%s] Using re-centering medians derived from input...\n",
              Sys.time()
            )
          )
        if (sdmedian.filename != "") {
          write.csv(sd.recenter, sdmedian.filename, row.names = FALSE)
          if (!quietly)
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
      # Use specified data
      if (!quietly)
        cat(
          sprintf(
            "[%s] Using specified re-centering medians...\n",
            Sys.time()
          )
        )
    }
    
    # ensure recentering medians are sorted correctly
    setkey(sd.recenter, param, sex, agedays)
    
    # add sd.recenter to data, and recenter
    setkey(data.all, param, sex, agedays)
    data.all <- sd.recenter[data.all]
    
    setkey(data.all, subjid, param, agedays)
    data.all[, tbc.sd := sd.orig - sd.median]
    if (prelim_infants){
      # separate out corrected and noncorrected values
      data.all[, ctbc.sd := sd.corr - sd.median]
    }
    
    if (sdrecentered.filename != "") {
      write.csv(data.all, sdrecentered.filename, row.names = FALSE)
      if (!quietly)
        cat(
          sprintf(
            "[%s] Wrote re-centered data to %s...\n",
            Sys.time(),
            sdrecentered.filename
          )
        )
    }
    
    # notification: ensure awareness of small subsets in data
    if (!quietly) {
      year.counts <- data.all[, .N, floor(agedays / 365.25)]
      if (year.counts[N < 100, .N] > 0) {
        cat(
          sprintf(
            "[%s] Note: input data has at least one age-year with < 100 subjects...\n",
            Sys.time()
          )
        )
      }
    }
    
    # safety check: treat observations where tbc.sd cannot be calculated as missing
    data.all[is.na(tbc.sd), exclude := 'Missing']
    
    if (prelim_infants){
      # 4: identify subset that don't need to be cleaned (nnte) ----
      
      # identify those meeting all subjects meeting these criteria as "no need
      # to ewma"
      # does the subject have sdes
      
      # keep the original column names -- we're adding a ton of columns that we
      # want to filter out after correction
      orig_colnames <- copy(colnames(data.all))
      # no SDEs
      data.all[, sum_sde := .N, by = c("subjid", "param", "agedays")]
      data.all[, no_sde := sum_sde == 1]
      # does the subject have identical values
      data.all[, sum_val := .N, by = c("subjid", "param", "v")]
      data.all[, no_dup_val := sum_val == 1]
      # all tbc.sd are within [-3,3] -- 0 is false
      data.all[, no_outliers := sum((tbc.sd > -3 & tbc.sd < 3) |
                                      is.na(tbc.sd)) == (.N),
               by = c("subjid", "param")]
      data.all[, no_outliers := no_outliers == 1]
      # all max - min tbd.sc < 2.5
      data.all[, no_bigdiff :=
                 rep((abs(max(tbc.sd, na.rm = TRUE) - min(tbc.sd, na.rm = TRUE)) < 2.5),
                     .N),
               by = c("subjid", "param")]
      # the previous value can't be too far from the current value
      data.all[, seq_win := sequence(.N), by = c("subjid", "param")]
      data.all[, nottoofar :=
                 (abs(tbc.sd - dplyr::lag(tbc.sd)) < 1 | seq_win == 1) &
                 (abs(tbc.sd - dplyr::lead(tbc.sd)) < 1 | seq_win == .N),
               by = c("subjid", "param")]
      data.all[is.na(nottoofar),  nottoofar :=  TRUE]
      
      # cumulative: no need to ewma -- needs to work for all within a subject &
      # parameter
      data.all[, nnte := no_sde & no_dup_val & no_outliers & no_bigdiff & nottoofar]
      # NNTE can be calculated by parameter -- but it's occasionally easier for
      # calculations to require all parameters to be nnte
      data.all[, nnte_full := sum(nnte) == .N, by = c("subjid", "param")]
      data.all[, nnte := sum(nnte) == .N, by = c("subjid")]
      
      
      # remove many added columns -- except for nnte
      orig_colnames <- c(orig_colnames, "nnte", "nnte_full")
      data.all <- data.all[, colnames(data.all) %in% orig_colnames,
                           with = FALSE]
    }
    # pediatric: cleanbatch (most of steps) ----
    
    # NOTE: the rest of cleangrowth's steps are done through cleanbatch().
    
    # Store off vars for later - CP
    # --- capture key checkpoint columns for later diagnostics ---
    checkpoint_cols <- c("id", "subjid", "param", "agedays",
                         "sd.orig", "sd.corr", "sd.orig_uncorr",
                         "z.orig", "tbc.sd", "ctbc.sd")
    checkpoint_cols <- intersect(checkpoint_cols, names(data.all))
    checkpoint_data <- data.all[, ..checkpoint_cols]
    # Store off vars for later - CP
    
    
    # optionally process batches in parallel
    if (!quietly)
      cat(sprintf(
        "[%s] Cleaning growth data in %d batch(es)...\n",
        Sys.time(),
        num.batches
      ))
    if (num.batches == 1) {
      if (!prelim_infants){
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
                             error.load.mincount = error.load.mincount)
      } else {
        ret.df <- cleanbatch_infants(
          data.all,
          log.path = log.path,
          quietly = quietly,
          parallel = parallel,
          measurement.to.z = measurement.to.z,
          ewma.fields = ewma.fields,
          recover.unit.error = recover.unit.error,
          include.carryforward = include.carryforward,
          sd.extreme = sd.extreme,
          z.extreme = z.extreme,
          exclude.levels = exclude.levels,
          tanner.ht.vel = tanner.ht.vel,
          who.ht.vel = who.ht.vel,
          lt3.exclude.mode = lt3.exclude.mode,
          error.load.threshold = error.load.threshold,
          error.load.mincount = error.load.mincount,
          ref.data.path = ref.data.path)
      }
    } else {
      # create log directory if necessary
      if (!is.na(log.path)) {
        cat(sprintf("[%s] Writing batch logs to '%s'...\n", Sys.time(), log.path))
        ifelse(!dir.exists(log.path), dir.create(log.path, recursive = TRUE), FALSE)
      }
      
      if (!prelim_infants){
        ret.df <- ddply(
          data.all,
          .(batch),
          cleanbatch,
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
      } else {
        ret.df <- ddply(
          data.all,
          .(batch),
          cleanbatch_infants,
          .parallel = parallel,
          .paropts = list(.packages = "data.table"),
          log.path = log.path,
          quietly = quietly,
          parallel = parallel,
          measurement.to.z = measurement.to.z,
          ewma.fields = ewma.fields,
          recover.unit.error = recover.unit.error,
          include.carryforward = include.carryforward,
          sd.extreme = sd.extreme,
          z.extreme = z.extreme,
          exclude.levels = exclude.levels,
          tanner.ht.vel = tanner.ht.vel,
          who.ht.vel = who.ht.vel,
          lt3.exclude.mode = lt3.exclude.mode,
          error.load.threshold = error.load.threshold,
          error.load.mincount = error.load.mincount,
          ref.data.path = ref.data.path
        )
      }
      stopCluster(cl)
    }
    
    
    if (!quietly)
      cat(sprintf("[%s] Done with pediatric data!\n", Sys.time()))
  } else {
    ret.df <- data.table()
    
    if (!quietly)
      cat(sprintf("[%s] No pediatric data. Moving to adult data...\n", Sys.time()))
  }
  
  # adult: send to cleanadult to do most of the work ----
  
  # no need to do this if there's no data
  if (nrow(data.adult) > 0){
    if (!quietly){
      cat(sprintf("[%s] Begin processing adult data...\n", Sys.time()))
    }
    
    # if parallel processing is desired, load additional modules
    if (parallel) {
      if (is.na(num.batches)) {
        num.batches <- getDoParWorkers()
      }
      # variables needed for parallel workers
      var_for_par <- c(
        "cleanadult", "check_between", "round_pt", "get_float_rem",
        "as.matrix.delta_dn", "ewma_dn", "remove_biv", "remove_biv_high",
        "remove_biv_low", "identify_rv", "temp_sde", "redo_identify_rv",
        "rem_hundreds", "rem_unit_errors", "get_num_places", "switch_tens_ones",
        "rem_transpositions", "ht_allow", "ht_change_groups",
        "ht_3d_growth_compare", "remove_ewma_wt", "remove_mod_ewma_wt"
      )
      
      cl <- makeCluster(num.batches)
      clusterExport(cl = cl, varlist = var_for_par, envir = environment())
      registerDoParallel(cl)
    } else {
      if (is.na(num.batches))
        num.batches <- 1
    }
    
    # this is where we do most of the adult work
    # is the randomness necessary here?
    subjid.unique <- data.adult[j = unique(subjid)]
    batches.adult <- data.table(
      subjid = subjid.unique,
      newbatch = sample(num.batches, length(subjid.unique), replace = TRUE)
    )
    data.adult <- merge(data.adult, batches.adult, by = "subjid")
    
    # add age in years
    data.adult[, age_years := agedays/365.25]
    # rename for ease of use
    data.adult[, measurement := v_adult]
    data.adult[, id := line]
    
    if (num.batches == 1) {
      # do the cleaning
      res <- cleanadult(data.adult, weight_cap = weight_cap)
    } else {
      res <- ddply(
        data.adult,
        .(newbatch),
        cleanadult,
        .parallel = parallel,
        .paropts = list(.packages = "data.table"),
        weight_cap = weight_cap
      )
      
      res <- as.data.table(res)
    }
    
    # replace result with missing if measurement or agedays are missing
    res[is.na(measurement) | agedays < 0, result := "Missing"]
    
    if (parallel){
      stopCluster(cl)
    }
    
    if (adult_columns_filename != "") {
      write.csv(res, adult_columns_filename, row.names = FALSE, na = "")
      if (!quietly){
        cat(
          sprintf(
            "[%s] Wrote adult data with additional columns to to %s...\n",
            Sys.time(),
            adult_columns_filename
          )
        )
      }
    }
    
    if (!quietly)
      cat(sprintf("[%s] Done with adult data!\n", Sys.time()))
  } else {
    res <- data.table()
    
    if (!quietly){
      cat(sprintf("[%s] No adult data. Moving to postprocessing...\n", Sys.time()))
    }
  }
  
  # if (any(nrow(data.all) > 0, nrow(data.adult) > 0)) {
  #   # join with pediatric data
  #   full_out <- data.table(
  #     line = c(ret.df$line, res$line),
  #     exclude = c(as.character(ret.df$exclude), res$result),
  #     mean_sde = c(rep(NA, nrow(ret.df)), res$mean_sde)
  #   )
  #   full_out[, exclude := factor(exclude, levels = exclude.levels)]
  #   full_out <- full_out[order(line),]
  #   # remove column added for keeping track
  #   full_out[, line := NULL]
  #   
  #   return(full_out$exclude)
  # } else {
  #   return(c())
  # }
  
  # --- FINAL RETURN (modified to return all columns) ---
  
  # combine pediatric and adult results if present
  if (any(nrow(data.all) > 0, nrow(data.adult) > 0)) {
    
    # join pediatric and adult outputs together
    full_out <- data.table(
      line = c(ret.df$line, res$line),
      exclude = c(as.character(ret.df$exclude), res$result)
    )
    
    # preserve original exclude levels
    full_out[, exclude := factor(exclude, levels = exclude.levels)]
    
    # join back to original input data to preserve all columns
    # this assumes data.all.ages was the original full input
    all_results <- merge(
      data.all.ages,           # original dataset (with id, param, etc.)
      full_out, 
      by = "line",
      all.x = TRUE,
      sort = FALSE
    )
    # Merge checkpoint diagnostics by ID
    all_results <- merge(all_results, checkpoint_data, by = "id", all.x = TRUE)
    all_results <- all_results[match(data.all.ages$id, all_results$id)]
    
    # restore original row order
    setorder(all_results, line)
    
    # if you have an 'id' field in data.work, it will now be preserved automatically
    # (since it existed in the original input)
    
    # Return *everything* for downstream debugging
    return(all_results)
    
  } else {
    # if there was no data at all, return an empty table with consistent structure
    return(data.table())
  }
  
  
}

#' Function to calculate z-scores and csd-scores based on anthro tables.
#'
#' @param path Path to supplied reference anthro data. Defaults to package anthro tables.
#' @param cdc.only Whether or not only CDC data should be used. Defaults to false.
#' @param prelim_infants TRUE/FALSE. Run the in-development release of the infants algorithm (expands pediatric algorithm to improve performance for children 0 – 2 years). Not recommended for use in research. For more information regarding the logic of the algorithm, see the vignette 'Preliminary Infants Algorithm.' Defaults to FALSE.
#'
#' @return Function for calculating BMI based on measurement, age in days, sex, and measurement value.
#' @export
#' @import data.table
#' @importFrom utils read.csv read.table
#' @examples
#' \donttest{
#' # Return calculating function with all defaults
#' afunc <- read_anthro()
#'
#' # Return calculating function while specifying a path and using only CDC data
#' afunc <- read_anthro(path = system.file("extdata", package = "growthcleanr"),
#'                      cdc.only = TRUE)
#' }
read_anthro <- function(path = "", cdc.only = FALSE, prelim_infants = FALSE) {
  # avoid "no visible bindings" warning
  src <- param <- sex <- age <- ret <- m <- NULL
  csdneg <- csdpos <- s <- NULL
  
  # set correct path based on input reference table path (if any)
  if (!prelim_infants){
    weianthro_path <- ifelse(
      path == "",
      system.file(file.path("extdata", "weianthro.txt.gz"), package = "growthcleanr"),
      file.path(path, "weianthro.txt.gz")
    )
    lenanthro_path <- ifelse(
      path == "",
      system.file(file.path("extdata", "lenanthro.txt.gz"), package = "growthcleanr"),
      file.path(path, "lenanthro.txt.gz")
    )
    bmianthro_path <- ifelse(
      path == "",
      system.file(file.path("extdata", "bmianthro.txt.gz"), package = "growthcleanr"),
      file.path(path, "bmianthro.txt.gz")
    )
    growth_cdc_ext_path <- ifelse(
      path == "",
      system.file(file.path("extdata", "growthfile_cdc_ext.csv.gz"), package = "growthcleanr"),
      file.path(path, "growthfile_cdc_ext.csv.gz")
    )
  } else {
    weianthro_path <- lenanthro_path <- bmianthro_path <-
      ifelse(
        path == "",
        system.file(file.path("extdata", "growthfile_who.csv.gz"), package = "growthcleanr"),
        file.path(path, "growthfile_who.csv.gz")
      )
    growth_cdc_ext_path <- ifelse(
      path == "",
      system.file(file.path("extdata", "growthfile_cdc_ext_infants.csv.gz"), package = "growthcleanr"),
      file.path(path, "growthfile_cdc_ext_infants.csv.gz")
    )
  }
  growth_cdc_ext <- read.csv(gzfile(growth_cdc_ext_path))
  
  l <- if (!prelim_infants){
    list(
      with(
        read.table(gzfile(weianthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'WEIGHTKG',
          sex = sex - 1,
          age,
          l,
          m,
          s,
          csdpos = as.double(NA),
          csdneg = as.double(NA)
        )
      ),
      with(
        read.table(gzfile(lenanthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'HEIGHTCM',
          sex = sex - 1,
          age,
          l,
          m,
          s,
          csdpos = as.double(NA),
          csdneg = as.double(NA)
        )
      ),
      with(
        read.table(gzfile(bmianthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'BMI',
          sex = sex - 1,
          age,
          l,
          m,
          s,
          csdpos = as.double(NA),
          csdneg = as.double(NA)
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'WEIGHTKG',
          sex,
          age = agedays,
          l = cdc_wt_l,
          m = cdc_wt_m,
          s = cdc_wt_s,
          csdpos = cdc_wt_csd_pos,
          csdneg = cdc_wt_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'HEIGHTCM',
          sex,
          age = agedays,
          l = cdc_ht_l,
          m = cdc_ht_m,
          s = cdc_ht_s,
          csdpos = cdc_ht_csd_pos,
          csdneg = cdc_ht_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'BMI',
          sex,
          age = agedays,
          l = cdc_bmi_l,
          m = cdc_bmi_m,
          s = cdc_bmi_s,
          csdpos = cdc_bmi_csd_pos,
          csdneg = cdc_bmi_csd_neg
        )
      )
    )
  } else {
    list(
      with(
        read.csv(gzfile(weianthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'WEIGHTKG',
          sex,
          age = agedays,
          l = who_wt_l,
          m = who_wt_m,
          s = who_wt_s,
          csdpos = who_wt_csd_pos,
          csdneg =  who_wt_csd_neg
        )
      ),
      with(
        read.csv(gzfile(lenanthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'HEIGHTCM',
          sex,
          age = agedays,
          l = who_ht_l,
          m = who_ht_m,
          s = who_ht_s,
          csdpos = who_ht_csd_pos,
          csdneg =  who_ht_csd_neg
        )
      ),
      with(
        read.csv(gzfile(lenanthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'HEADCM',
          sex,
          age = agedays,
          l = who_hc_l,
          m = who_hc_m,
          s = who_hc_s,
          csdpos = who_hc_csd_pos,
          csdneg =  who_hc_csd_neg
        )
      ),
      with(
        read.csv(gzfile(bmianthro_path), header = TRUE),
        data.frame(
          src = 'WHO',
          param = 'BMI',
          sex,
          age = agedays,
          l = who_bmi_l,
          m = who_bmi_m,
          s = who_bmi_s,
          csdpos = who_bmi_csd_pos,
          csdneg =  who_bmi_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'WEIGHTKG',
          sex,
          age = agedays,
          l = cdc_wt_l,
          m = cdc_wt_m,
          s = cdc_wt_s,
          csdpos = cdc_wt_csd_pos,
          csdneg = cdc_wt_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'HEIGHTCM',
          sex,
          age = agedays,
          l = cdc_ht_l,
          m = cdc_ht_m,
          s = cdc_ht_s,
          csdpos = cdc_ht_csd_pos,
          csdneg = cdc_ht_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'HEADCM',
          sex,
          age = agedays,
          l = cdc_hc_l,
          m = cdc_hc_m,
          s = cdc_hc_s,
          csdpos = cdc_hc_csd_pos,
          csdneg = cdc_hc_csd_neg
        )
      ),
      with(
        growth_cdc_ext,
        data.frame(
          src = 'CDC',
          param = 'BMI',
          sex,
          age = agedays,
          l = cdc_bmi_l,
          m = cdc_bmi_m,
          s = cdc_bmi_s,
          csdpos = cdc_bmi_csd_pos,
          csdneg = cdc_bmi_csd_neg
        )
      )
    )
  }
  
  anthro <- rbindlist(l)
  
  
  setkey(anthro, src, param, sex, age)
  
  return(function(param, agedays, sex, measurement, csd = FALSE) {
    # For now, we will only use CDC growth reference data, note that the cubically interpolated file
    # we are using has linear measurments derived from length data for children < 731 days, and height thereafter
    src <- ifelse(agedays < 3*365.25 & !cdc.only, 'WHO', 'CDC')
    
    # keep column sequence the same fo efficient join
    dt <- data.table(src, param, sex, agedays, measurement)
    dt <- anthro[dt]
    
    dt[, ret := as.double(NA)]
    if (csd) {
      dt[measurement < m, ret := (measurement - m) / csdneg]
      dt[measurement >= m, ret := (measurement - m) / csdpos]
    } else {
      dt[l == 0, ret := log(measurement / m) / s]
      dt[l != 0, ret := (((measurement / m) ^ l) - 1) / (l * s)]
    }
    
    return(dt$ret)
  })
}

#' Exponentially Weighted Moving Average (EWMA)
#'
#' \code{ewma} calculates the exponentially weighted moving average (EWMA) for a set of numeric observations over time.
#'
#' @param agedays Vector of age in days for each z score (potentially transformed to adjust weighting).
#'
#' @param z Input vector of numeric z-score data.
#'
#' @param ewma.exp Exponent to use for weighting.
#'
#' @param ewma.adjacent Specify whether EWMA values excluding adjacent measurements should be calculated.  Defaults to TRUE.
#'
#' @return Data frame with 3 variables:
#' * The first variable (ewma.all) contains the EWMA at observation time
#'   excluding only the actual observation for that time point.
#' * The second variable (ewma.before) contains the EWMA for each observation excluding both the actual observation
#'   and the immediate prior observation.
#' * The third variable (ewma.after) contains the EWMA for each observation excluding both the actual observation
#'   and the subsequent observation.
#'@md
#'
#' @export
#' @examples
#' # Run on 1 subject, 1 type of parameter
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid == df_stats$subjid[1] &
#'                        df_stats$param == "HEIGHTCM", ]
#'
#' # Get the uncentered z-scores
#' measurement_to_z <- read_anthro(cdc.only = TRUE)
#' sd <- measurement_to_z(df_stats$param,
#'                        df_stats$agedays,
#'                        df_stats$sex,
#'                        df_stats$measurement,
#'                        TRUE)
#'
#' # Calculate exponentially weighted moving average
#' e_df <- ewma(df_stats$agedays, sd, ewma.exp = -1.5)
ewma <- function(agedays, z, ewma.exp, ewma.adjacent = TRUE) {
  # 6.  EWMA calculation description: Most of the next steps will involve calculating the exponentially weighted moving average for each subject and parameter. I will
  #     describe how to calculate EWMASDs, and will describe how it needs to be varied in subsequent steps.
  # a.	The overall goal of the EWMASD calculation is to identify the difference between the SD-score and what we might predict that DS-score should be, in order to
  #     determine whether it should be excluded.
  # b.	Only nonmissing SD-scores for a parameter that are not designated for exclusion are included in the following calculations.
  # c.	For each SD-score SDi and associated agedaysi calculate the following for every other z-score (SDj…SDn) and associated agedays (agedaysj…agedaysn)  for the
  #     same subject and parameter
  #   i.	ΔAgej=agedaysj-agedaysi
  #   ii.	EWMAZ=SDi=[Σj→n(SDj*((5+ΔAgej)^-1.5))]/[ Σj→n((5+ΔAgej)^-1.5)]
  #   iii.	For most EWMASD calculations, there are 3 EWMASDs that need to be calculated. I will note if not all of these need to be done for a given step.
  #     1.	EWMASDall calculated as above
  #     2.	EWMAZbef calculated excluding the SD-score just before the SD-score of interest (sorted by agedays). For the first observation for a parameter for a
  #         subject, this should be identical to EWMASDall rather than missing.
  #     3.	EWMAZaft calculated excluding the z-score just after the SD-score of interest (sorted by agedays). For the lastobservation for a parameter for a subject,
  #         this should be identical to EWMASDall rather than missing.
  #   iv.	For each of the three EWMASDs, calculate the dewma_*=SD-EWMASD
  # d.	EWMASDs and ΔEWMASDs will change if a value is excluded or manipulated using one of the methods below, therefore EWMASDs and ΔEWMASDs be recalculated for each
  #     step where they are needed.
  # e.	For these calculations, use variables that allow for precise storage of numbers (in Stata this is called 'double') because otherwise rounding errors can cause
  #     problems in a few circumstances
  
  n <- length(agedays)
  # initialize response variables
  ewma.all <- ewma.before <- ewma.after <- vector('numeric', 0)
  if (n > 0) {
    # organize into data frame and sort into order of increasing age,
    # but retain original sort order information in index
    if (!all(agedays == cummax(agedays)))
      warning("EWMA ordering is not sorted; double check") #add in a check to make sure the inputs are already sorted (they should be)
    index <- order(agedays)
    
    
    # calculate matrix of differences in age, and add 5 to each delta per Daymont algorithm
    delta <- as_matrix_delta(agedays)
    delta <- ifelse(delta == 0, 0, (delta + 5) ^ ewma.exp)
    
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
  return(if (ewma.adjacent)
    data.frame(ewma.all, ewma.before, ewma.after)
    else
      data.frame(ewma.all))
}

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
#' @examples
#' # Run on 1 subject
#' df_stats <- as.data.frame(syngrowth)
#' df_stats <- df_stats[df_stats$subjid == df_stats$subjid[1], ]
#'
#' # Get the original standard deviations
#' measurement_to_z <- read_anthro(cdc.only = TRUE)
#' sd.orig <- measurement_to_z(df_stats$param,
#'                        df_stats$agedays,
#'                        df_stats$sex,
#'                        df_stats$measurement,
#'                        TRUE)
#'
#' # Calculate median standard deviations
#' sd.m <- sd_median(df_stats$param,
#'                   df_stats$sex,
#'                   df_stats$agedays,
#'                   sd.orig)
sd_median <- function(param, sex, agedays, sd.orig) {
  # 3.  SD-score recentering: Because the basis of the method is comparing SD-scores over time, we need to account for the fact that
  #     the mean SD-score for the population changes with age.
  # a.  Determine the median cdc*sd for each parameter by year of age (with sexes combined): median*sd.
  # b.	The median*sd should be considered to apply to midyear-age, defined as the age in days with the same value as the integer
  #     portion of (365.25*year + 365.25/2).
  # c.	Linearly interpolate median*sd for each parameter between each midyear-age, naming the interpolated values rc*sd.
  # d.	For ages below the first midyear-age, let rc*sd equal the median*sd for the earliest year.
  #     For ages above the last midyear_age, let rc*sd equal the median*sd for the last year.
  # e.	Subtract rcsd_* from SDorig to create the recentered SD-score.  This recentered SD-score, labeled tbc*sd
  #     (stands for "to be cleaned") will be used for most of the rest of the analyses.
  # f.	In future steps I will sometimes refer to measprev and measnext which refer to the previous or next wt or ht measurement
  #     for which exc_*==0 for the subject and parameter, when the data are sorted by subject, parameter, and agedays. SDprev and SDnext refer to the tbc*sd of the previous or next measurement.
  
  # avoid "no visible binding" warnings
  ageyears <- sd.median <- NULL
  
  dt <- data.table(param, sex, agedays, ageyears = floor(agedays / 365.25), sd.orig)
  setkey(dt, param, sex, agedays)
  # determine ages (in days) we need to cover from min to max age in years
  agedays.to.cover <- dt[, (floor(min(ageyears) * 365.25):floor(max(ageyears +
                                                                      1) * 365.25))]
  # group all measurements above age 19 together (copied from CD stata code e-mailed 4/3/15)
  dt[ageyears > 19, ageyears := 19]
  dt.median <- dt[!is.na(sd.orig), list(sex = c(0, 1), sd.median = rep(median(sd.orig), 2)), by =
                    .(param, ageyears)]
  dt.median <- dt.median[, list(
    agedays = agedays.to.cover,
    sd.median = approx(
      floor((ageyears + 0.5) * 365.25),
      sd.median,
      xout = agedays.to.cover,
      rule = 2
    )$y
  ), by = .(param, sex)]
  setkey(dt.median, param, sex, agedays)
  return(dt.median)
}



#### CP ####

# Infants support

# supporting (all) ---

# function that takes in a parameter name and returns the designated other
# parameter (DOP)
#' @keywords internal
#' @noRd
get_dop <- function(param_name){
  dop <- if (param_name == "WEIGHTKG"){
    "HEIGHTCM"
  } else if (param_name == "HEIGHTCM"){
    "WEIGHTKG"
  } else { # HEADCM
    "HEIGHTCM"
  }
  
  return(dop)
}

# temporary extraneous ----
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
temporary_extraneous_infants <- function(df) {
  # avoid "no visible binding" warnings
  agedays <- delta.median.sd <- extraneous <- extraneous.this.day <- NULL
  index <- median.sd <- param <- subjid <- tbc.sd <- NULL
  
  # add subjid and param if needed (may be missing depending on where this is called from)
  if (is.null(df$subjid))
    df[, subjid := NA]
  if (is.null(df$param))
    df[, param := NA]
  # only operate on valid rows (but include rows that may have previously been flagged as a "temporary extraneous")
  # also do not calculate temporary extraneous for NNTEs
  valid.rows <- valid(df, include.temporary.extraneous = TRUE) &
    !df$nnte
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
  
  # add the "other map" -- now including head circumference
  other_p_map <- c(
    "HEIGHTCM" = "WEIGHTKG",
    "WEIGHTKG" = "HEIGHTCM",
    "HEADCM" = "HEIGHTCM"
  )
  
  df[valid.rows &
       subjid %in% subj.all.dups, delta.median.sd := (function(subj.df) {
         # iterate over parameters where delta.median.sd is not yet defined
         # pass 1: take median from other parameter(s)
         for (p in subj.df[is.na(delta.median.sd) &
                           extraneous.this.day, unique(param)]) {
           median.other.sd <- subj.df[param == other_p_map[p] &
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

# evil twins ----

#' Function to calculate out of bounds (OOB) measurements for Evil Twins step
#'
#' @param df data table with all parameters
#'
#' @keywords internal
#' @noRd
calc_oob_evil_twins <- function(df){
  # # start by determining if a measurement is out of bounds (oob)
  
  # for differences for adjacent, we pad with a high number to default to true
  # we don't want to consider the comparison with the next subj/param
  oob <- ((
    (abs(c(df$tbc.sd[2:nrow(df)], Inf) - df$tbc.sd) > 5 &
       rev(duplicated(rev(paste(df$subjid, "_", df$param))))) |
      (abs(c(Inf, df$tbc.sd[1:(nrow(df)-1)]) - df$tbc.sd) > 5 &
         duplicated(paste(df$subjid, "_", df$param)))
  )) & ((
    (abs(c(df$ctbc.sd[2:nrow(df)], Inf) - df$ctbc.sd) > 5 &
       rev(duplicated(rev(paste(df$subjid, "_", df$param))))) |
      (abs(c(Inf, df$ctbc.sd[1:(nrow(df)-1)]) - df$ctbc.sd) > 5 &
         duplicated(paste(df$subjid, "_", df$param)))
  ))
  
  df[, "oob" := oob]
  
  return(df)
}

# moderate ewma ----

# function to calculate and recenter z scores for a given column
# df: data frame with parameter, agedays, sex, cn
# cn: column name to calculate, smooth, and recenter
# ref.data.path: reference data path
#
# returns df with additional column, tbc.(cn), which is the recentered z-score
# for the input
#' @keywords internal
#' @noRd
calc_and_recenter_z_scores <- function(df, cn, ref.data.path){
  # avoid no visible warning errors
  cn.orig_cdc <- param <- agedays <- sex <- cn.orig_who <- cn.orig <- subjid <-
    tbc.cn <- sd.median <- NULL
  
  # for infants, use z and who
  measurement.to.z <- read_anthro(ref.data.path, cdc.only = TRUE,
                                  prelim_infants = TRUE)
  measurement.to.z_who <- read_anthro(ref.data.path, cdc.only = FALSE,
                                      prelim_infants = TRUE)
  
  # calculate "standard deviation" scores
  df[, cn.orig_cdc := measurement.to.z(param, agedays, sex, get(cn), TRUE)]
  df[, cn.orig_who := measurement.to.z_who(param, agedays, sex, get(cn), TRUE)]
  
  # smooth z-scores/SD scores between ages 1 - 3yo using weighted scores
  # older uses cdc, younger uses who
  who_weight <- 4 - (df$agedays/365.25)
  cdc_weight <- (df$agedays/365.25) - 2
  
  smooth_val <- df$agedays/365.25 >= 2 &
    df$agedays/365.25 <= 4 &
    df$param != "HEADCM"
  df[smooth_val,
     cn.orig := (cn.orig_cdc[smooth_val]*cdc_weight[smooth_val] +
                   cn.orig_who[smooth_val]*who_weight[smooth_val])/2]
  
  # otherwise use WHO and CDC for older and younger, respectively
  who_val <- df$param == "HEADCM" |
    df$agedays/365.25 < 2
  df[who_val | (smooth_val & is.na(df$cn.orig_cdc)),
     cn.orig := df$cn.orig_who[who_val  | (smooth_val & is.na(df$cn.orig_cdc))]]
  
  cdc_val <- df$param != "HEADCM" |
    df$agedays/365.25 > 4
  df[cdc_val | (smooth_val & is.na(df$cn.orig_who)),
     cn.orig := df$cn.orig_cdc[cdc_val | (smooth_val & is.na(df$sd.orig_who))]]
  
  # now recenter -- already has the sd.median from the original recentering
  setkey(df, subjid, param, agedays)
  df[, tbc.cn := cn.orig - sd.median]
  
  # rename ending column
  setnames(df, "tbc.cn", paste0("tbc.", cn))
  
  return(df)
}

# Main pediatric growthcleanr function -- infants cleanbatch
# internal supporting functions for pediatrics can be found in: infants_support.R
# note

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
#' @importFrom stats median embed
#' @noRd
cleanbatch_infants <- function(data.df,
                               log.path,
                               quietly,
                               parallel,
                               measurement.to.z,
                               ewma.fields,
                               recover.unit.error,
                               include.carryforward,
                               sd.extreme,
                               z.extreme,
                               exclude.levels,
                               tanner.ht.vel,
                               who.ht.vel,
                               lt3.exclude.mode,
                               error.load.threshold,
                               error.load.mincount,
                               ref.data.path) {
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
  
  oob <- sd_med <- med_diff <- max_diff <- sum_oob <- i.exclude <- NULL
  
  
  # avoid no visible warning errors
  sum_sde <- no_sde <- cf <- wholehalfimp <- seq_win <- cs <- absdiff <-
    sd.orig_uncorr <- seq_win <- absdiff <- wholehalfimp <- ageyears <- ctbc.sd <-
    ..col_replace <- c.ewma.all <- pot_excl <- c.dewma.all <- p_plus <-
    p_minus <- ctbc.sd <- c.ewma.all <- tbc_diff_next <- tbc_diff_prior <-
    tbc_diff_plus_next <- tbc.p_plus <- tbc_diff_plus_prior <-
    tbc_diff_minus_next <- tbc.p_minus <- tbc_diff_minus_prior <- addcrithigh <-
    addcritlow <- tbc_dop <- i.tbc.sd <- rowind <- abssum <- c.dewma.all <-
    whoagegrp_ht <- d_agedays <- mindiff <- maxdiff <- who_mindiff_ht <-
    who_maxdiff_ht <- mindiff_prev <- maxdiff_prev <- whoinc.age.hc <-
    who_maxdiff_hc <- who_mindiff_hc <- diff_prev <-
    diff_next <- aft.g.aftm1 <- val_excl <-
    absval <- comp_diff <- err_ratio <-
    NULL
  
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
  
  # 5: temporary SDEs ----
  
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
  
  # 6:  carried forwards ----
  
  if (!include.carryforward) {
    if (!quietly)
      cat(sprintf(
        "[%s] Exclude measurements carried forward...\n",
        Sys.time()
      ))
    
    # save original column names
    orig_colnames <- colnames(data.df)
    
    # we only running carried forwards on valid values, non NNTE values,
    # and non single values
    tmp <- table(paste0(data.df$subjid, "_", data.df$param))
    not_single <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 1]
    valid_set <- valid(data.df, include.temporary.extraneous = TRUE) &
      !data.df$nnte_full &
      not_single
    data.sub <- data.df[valid_set,]
    
    # for ease, order by subject, parameter, and agedays
    data.sub <- data.sub[order(subjid, param, agedays)]
    
    # if no sdes, exclude carried forward values
    data.sub[, sum_sde := .N, by = c("subjid", "param", "agedays")]
    data.sub[, no_sde := !paste0(data.sub$subjid, "_", data.sub$param) %in%
               paste0(data.sub$subjid, "_", data.sub$param)[data.sub$sum_sde > 1]]
    data.sub[no_sde == TRUE, cf :=  (v.orig - shift(v.orig)) == 0, by = c("subjid", "param")]
    data.sub[is.na(cf), cf := FALSE]
    
    # for values with sdes, compare to all values from the prior day
    data.sub[no_sde == FALSE, cf := (function(df){
      ages <- unique(agedays)
      if (length(ages) > 1){
        for (i in 2:length(ages)) {
          # find the set of measurements from the previous age in days
          all.prev.v <- df[agedays == ages[i - 1], v.orig]
          # if a measurement for the current age is found in the set of measurements from the previous age, then mark it as carried forward
          df[agedays == ages[i] &
               v.orig %in% all.prev.v, cf := TRUE]
        }
      }
    })(copy(.SD)), .SDcols = c('agedays', 'cf', 'v.orig'), by = .(subjid, param)]
    
    # merge in the carried forwards
    cf_idx <- data.sub$index[data.sub$cf]
    data.df[index %in% cf_idx, exclude := "Exclude-Carried-Forward"]
    
    ### CP NOTE - Check closely - potential overwrite of carried forwar
    
    # redo temp sde
    # data.df$exclude[temporary_extraneous_infants(data.df)] <- 'Exclude-Temporary-Extraneous-Same-Day'
    
    ### CP NOTE - Check closely - potential overwrite of carried forwar
    
    
    # find out if values are whole or half imperial
    data.df[param == "WEIGHTKG",
            wholehalfimp := (v.orig * 2.20462262)%%1 < 0.01]
    data.df[param == "HEIGHTCM",
            wholehalfimp := (v.orig * 1.27)%%1 < 0.01]
    data.df[param == "HEIGHTCM",
            wholehalfimp := (v.orig * 1.27 / 2)%%1 < 0.01]
    data.df[param == "HEADCM",
            wholehalfimp := (v.orig * 1.27)%%1 < 0.01]
    data.df[param == "HEADCM",
            wholehalfimp := (v.orig * 1.27 / 2)%%1 < 0.01]
    # find string carried forward length
    data.df[, seq_win := sequence(rle(as.character(exclude))$lengths),
            by = c("subjid", "param")]
    data.df[exclude != "Exclude-Carried-Forward", seq_win := NA]
    # relabel all the initial values as "0" -- reorder
    data.df <- data.df[order(subjid, param, v.orig, agedays)]
    data.df[which(seq_win == 1)-1, seq_win :=  0]
    # reorder correctly
    data.df <- data.df[order(subjid, param, agedays)]
    # also create labels
    data.df[, cs := rep(1:length(rle(as.character(exclude))$lengths),
                        times = rle(as.character(exclude))$lengths),
            by = c("subjid", "param")]
    # fix the first one -- have to do this outside of data table
    data.df[which(seq_win == 0), "cs"] <-
      data.df[which(seq_win == 0)+1, "cs"]
    data.df[is.na(seq_win), cs := NA]
    
    # find th diff between initial and cf
    data.df[!is.na(seq_win), absdiff :=
              abs(sd.orig_uncorr - sd.orig_uncorr[1]),
            by = c("subjid", "param", "cs")]
    
    # handle CFs by case
    data.df[!is.na(seq_win), exclude := (function(df){
      # only 1 cf in string
      if (max(seq_win) == 1){
        df[seq_win != 0 & absdiff < 0.05,
           exclude := "Exclude-1-CF-deltaZ-<0.05"]
        
        # use short circuiting to have the lower end covered (>= 0.05)
        df[seq_win != 0 & absdiff < .1 & wholehalfimp,
           exclude := "Exclude-1-CF-deltaZ-<0.1-wholehalfimp"]
      } else if (max(seq_win) > 1){
        df[seq_win != 0 & agedays/365.25 > 16 & sex == 1 & absdiff < 0.05,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.05"]
        df[seq_win != 0 & agedays/365.25 > 17 & sex == 0 & absdiff < 0.05,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.05"]
        
        df[seq_win != 0 & agedays/365.25 > 16 & sex == 1 & absdiff < .1 &
             wholehalfimp,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp"]
        df[seq_win != 0 & agedays/365.25 > 17 & sex == 0 & absdiff < .1 &
             wholehalfimp,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp"]
      }
      
      return(df$exclude)
    })(copy(.SD)),
    .SDcols = c('agedays', "seq_win", 'absdiff', "sex", "cs", "wholehalfimp",
                "exclude"),
    by = c("subjid", "param", "cs")]
  }
  
  # BIV ----
  
  # Wt limits from do-file Oct 10 2022, in dataset, highest plausible weight=29.5kg, 65 lbs
  # HC Limits based on analysis in do-file from Oct 11 2022: in dataset, lowest plausible 18.5, z=-13, highest plausible 63.5, z=11.2  (at birth, highest=42, z=6)
  
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude BIVs...\n",
      Sys.time()
    ))
  
  # add age in years
  data.df[, ageyears := agedays/365.25]
  
  # calculate the valid step
  valid_set <- valid(data.df, include.temporary.extraneous = TRUE) &
    !data.df$nnte_full
  
  exc_nam <- "Exclude-Absolute-BIV"
  
  # identify absolute cutoffs
  # Min/max weight at birth based on published births
  data.df[valid_set & param == "WEIGHTKG" & v < 0.2 & agedays == 0,
          exclude := exc_nam]
  data.df[valid_set & param == "WEIGHTKG" & v < 1 & agedays != 0,
          exclude := exc_nam]
  # Min weight after birth based on rules about who can go home with alowance for significant weight loss -- this is for outpatient data only so very low weight babies would still be in NICU
  data.df[valid_set & param == "WEIGHTKG" & v > 10.5 &
            agedays == 0,
          exclude := exc_nam]
  # Max weight for <2y based on eval (see do-file from Oct 10 2022), highest plaus weight 29.5kg
  data.df[valid_set & param == "WEIGHTKG" & v < 1 &
            agedays == 0,
          exclude := exc_nam]
  # Max weight for <2y based on eval (see do-file from Oct 10 2022), highest plaus weight 29.5kg
  data.df[valid_set & param == "WEIGHTKG" & v > 35 &
            ageyears < 2,
          exclude := exc_nam]
  # Max weight for all based on published data
  data.df[valid_set & param == "WEIGHTKG" & v > 600,
          exclude := exc_nam]
  
  # Min/max HC based on analysis in do file from
  # Also, 18 is z=-6 for 22 0/7 in Fenton and 65 is z=6 for 40 0/7
  data.df[valid_set & param == "HEIGHTCM" & v < 18,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & v > 244,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & v > 65 &
            agedays == 0,
          exclude := exc_nam]
  
  # Min/max HC based on analysis in do file from Oct 11 2022
  # Also, 13 is z=-6 for 22 0/7 in Fenton and
  data.df[valid_set & param == "HEADCM" & v < 13,
          exclude := exc_nam]
  data.df[valid_set & param == "HEADCM" & v > 75,
          exclude := exc_nam]
  data.df[valid_set & param == "HEADCM" & v > 50 &
            agedays == 0,
          exclude := exc_nam]
  
  exc_nam <- "Exclude-Standardized-BIV"
  
  # identify z cutoff
  # ***Note, using unrecentered values***
  #  *For weight only do after birth
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < -25 &
            ageyears < 1,
          exclude := exc_nam]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < -15 &
            ageyears >= 1,
          exclude := exc_nam]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr > 22,
          exclude := exc_nam]
  
  # *Max z-score for height based on analysis of CHOP data because 15/25 too loose for upper limits
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < -25 &
            ageyears < 1,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < -15 &
            ageyears >= 1,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr > 8,
          exclude := exc_nam]
  
  # head circumference
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr < -15,
          exclude := exc_nam]
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr > 15,
          exclude := exc_nam]
  
  # then, do some explicit overwriting for the 0 case (otherwise will be
  # set as missing)
  data.df[v == 0, exclude := exc_nam]
  
  # 9d.  Replace exc_*=0 if exc_*==2 & redo step 5 (temporary extraneous)
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']
  
  # evil twins ----
  # Evil Twins: An important weakness in the original pediatric growthcleanr algorithm was that it often failed to identify two or more implausible measurements that occurred next to each other, even if they were extremely deviant from a child’s other measurements. This step is now added to identify these multiple extreme values, although it also identifies some single values.
  
  exc_nam <- "Exclude-Evil-Twins"
  
  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non pair
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single_pairs <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 2]
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    !data.df$nnte_full & # does not use the "other"
    not_single_pairs
  
  # make sure it's ordered by subjid and parameter
  data.df <- data.df[order(subjid, param),]
  
  # 9A/B/C
  # first, find out of any re possible evil twins in the first place
  start_df <- calc_oob_evil_twins(data.df[valid_set,])
  
  if (any(start_df$oob)) {
    if (!quietly)
      cat(sprintf(
        "[%s] Exclude evil twins...\n",
        Sys.time()
      ))
    
    # start to evaluate and remove evil twins
    data.df[valid_set, exclude := (function(df) {
      # where we are updating results
      upd.df <- copy(df)
      upd.df <- calc_oob_evil_twins(upd.df)
      # count the amount of oobs for each subject/param and distribute it out
      upd.df[, `:=` (sum_oob = sum(oob, na.rm = TRUE)), by =.(subjid, param)]
      
      any_oob <- any(upd.df$sum_oob >= 2)
      # while there are multiple oob, we want to remove
      while (any_oob){
        
        # 9D
        # now calculate the maximum difference from the median tbc.sd
        upd.df[, `:=` (sd_med = median(tbc.sd, na.rm = TRUE)), by =.(subjid, param)]
        upd.df[, `:=` (med_diff = abs(tbc.sd - sd_med)), by =.(subjid, param)]
        upd.df[, `:=` (max_diff = med_diff  == max(med_diff)), by =.(subjid, param)]
        # for ones with no tbc.sd, mark as false
        upd.df[is.na(max_diff), max_diff := FALSE]
        
        upd.df[sum_oob > 0 & max_diff, exclude := exc_nam]
        
        df[upd.df[exclude == exc_nam,], exclude := i.exclude, on = .(line)]
        
        #9E
        # reupdate valid (to recalculate OOB -- others are not included)
        upd.df <- calc_oob_evil_twins(df[valid(df),])
        upd.df[, `:=` (sum_oob = sum(oob, na.rm = TRUE)), by =.(subjid, param)]
        
        any_oob <- any(upd.df$sum_oob >= 2)
        
      }
      
      return(df$exclude)
    })(copy(.SD))]
    
  }
  
  # 9F.  redo temp extraneous
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']
  
  # Extreme EWMA ----
  
  # 11.  Exclude extreme errors with EWMA
  # a.	Erroneous measurements can distort the EWMA for measurements around them. Therefore, if the EWMA method identifies more than one value for a subject and
  #     parameter that meets criteria for exclusion, we will only exclude the value that deviates the most from expected in any given step. Then we will repeat the
  #     entire process until no measurements are identified that meet criteria for exclusion.
  # b.	Perform a EWMA calculation
  #   i.	Only use values where exc_*==0 to determine the EWMAs. However, calculate dewma_* variables for values where exc_*==0 or exc_*==2
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude extreme measurements based on EWMA...\n",
      Sys.time()
    ))
  data.df <- data.df[, exclude := (function(df) {
    num.ewma.excluded <- 0
    # optimization: determine whether this subject has any extraneous
    has.extraneous <- subjid %in% subj.dup
    
    # only enter if there's more than one valid value
    tmp <- table(paste0(df$subjid, "_", df$param))
    not_single_pairs <- paste0(df$subjid, "_", df$param) %in% names(tmp)[tmp > 2]
    valid_set <- valid(df, include.temporary.extraneous = TRUE) &
      !df$nnte_full & # does not use the "other"
      not_single_pairs
    
    if (sum(valid_set) > 1){
      while (TRUE) {
        tmp <- table(paste0(df$subjid, "_", df$param))
        not_single_pairs <- paste0(df$subjid, "_", df$param) %in% names(tmp)[tmp > 2]
        valid_set <- valid(df, include.temporary.extraneous = TRUE) &
          !df$nnte_full & # does not use the "other"
          not_single_pairs
        
        
        df[, (ewma.fields) := as.double(NaN)]
        
        # first, calculate which exponent we want to put through (pass a different
        # on for each exp)
        # subset df to only valid rows
        df_sub <- df[valid_set,]
        tmp <- data.frame(
          "before" = abs(df_sub$agedays - c(NA, df_sub$agedays[1:(nrow(df_sub)-1)])),
          "after" = abs(df_sub$agedays - c(df_sub$agedays[2:(nrow(df_sub))], NA))
        )
        maxdiff <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
        exp_vals <- rep(-1.5, nrow(tmp))
        exp_vals[maxdiff > 365.25] <- -2.5
        exp_vals[maxdiff > 730.5] <- -3.5
        df[valid_set, exp_vals := exp_vals]
        
        # calculate ewma -- need to do it in a special way to not include
        # temp extaneous
        
        # first calculate ewma for rows without temp extraneous
        no_tmp_extr <- valid(df, include.temporary.extraneous = FALSE)
        
        df[valid_set & no_tmp_extr, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
        df[valid_set & no_tmp_extr, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]
        
        # now go through the temp extraneous and do the ewmas just for those
        # rows
        if (sum(df$exclude == "Exclude-Temporary-Extraneous-Same-Day") > 0){
          # go through each temp extraneous
          for (j in
               which(df$exclude == "Exclude-Temporary-Extraneous-Same-Day")){
            # subset to only the temp extraneous and valid values
            df_sub <- df[valid_set & (no_tmp_extr | (1:nrow(df) == j)),]
            
            # run ewma on those values
            df_sub[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
            df_sub[, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]
            
            # merge the value back in
            col_replace <- c(ewma.fields, paste0("c.",ewma.fields))
            df_val <-
              df_sub[df_sub$exclude == "Exclude-Temporary-Extraneous-Same-Day",
                     ..col_replace]
            df[j, (col_replace) := df_val]
          }
          
        }
        
        # include the corrected
        df[, `:=`(
          dewma.all = tbc.sd - ewma.all,
          dewma.before = tbc.sd - ewma.before,
          dewma.after = tbc.sd - ewma.after,
          
          c.dewma.all = ctbc.sd - c.ewma.all
        )]
        
        
        # calculate potential exclusions
        df[valid_set, pot_excl :=
             (dewma.all > 3.5 & dewma.before > 3 & dewma.after > 3 & tbc.sd > 3.5 &
                ((!is.na(ctbc.sd) & c.dewma.all > 3.5) | is.na(ctbc.sd))
             ) |
             # for both
             (
               dewma.all < -3.5 &
                 dewma.before < -3 & dewma.after < -3 & tbc.sd < -3.5 &
                 ((!is.na(ctbc.sd) & c.dewma.all > 3.5) | is.na(ctbc.sd))
             )
        ]
        df[valid_set,][c(1, nrow(df)), pot_excl := FALSE]
        df[!valid_set, pot_excl := FALSE] # to compensate for later
        
        
        # 11c.  Identify all values that meet all of the following criteria as potential exclusions:
        #   i.	There are 3 or more measurements for that subject and parameter
        #   ii.	(dewma_*>3.5 & dewma_*_bef>3 & dewma_*_aft>3 & tbc*sd>3.5) OR (dewma_*<-3.5 & dewma_*_bef<-3 & d & dewma_*_aft<-3 & tbc*sd<-3.5)
        #   iii.	exc_*==0
        
        num.exclude <- sum(df$pot_excl)
        # 11d.  If there is only one potential exclusion identified in step 11c for a subject and parameter, replace exc_*=5 for that value
        if (num.exclude == 1)
          df[pot_excl, exclude := 'Exclude-EWMA1-Extreme']
        # 11e.  If there is more than one potential exclusion identified in step 11c for a subject and parameter, calculate abssum_*=|tbc*sd+dewma_*| for each exclusion and
        #     replace exc_*=5 for the value with the highest abssum_*
        if (num.exclude > 1) {
          # first order by decreasing abssum
          worst.row <- with(df, order(pot_excl, abs(tbc.sd + dewma.all),
                                      decreasing = TRUE))[1]
          df[worst.row, exclude := 'Exclude-EWMA1-Extreme']
        }
        
        # 11h.  Recalculate temporary extraneous as in step 5
        # optimize: only perform these steps if this subject is known to have extraneous measurements
        if (has.extraneous) {
          df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
          df[temporary_extraneous_infants(df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']
        }
        
        # 11i.  If there was at least one subject who had a potential exclusion identified in step 11c, repeat steps 11b-11g. If there were no subjects with potential
        #     exclusions identified in step 11c, move on to step 12.
        newly.excluded <- sum(df$exclude %in% c('Exclude-EWMA1-Extreme'))
        if (newly.excluded > num.ewma.excluded) {
          num.ewma.excluded <- newly.excluded
        } else {
          break
        }
      }
    }
    return(df$exclude)
  })(copy(.SD)), by = .(subjid, param), .SDcols = c('index', 'sex', 'agedays', 'tbc.sd', 'ctbc.sd', 'nnte_full', 'exclude')]
  
  # 9d.  Replace exc_*=0 if exc_*==2 & redo step 5 (temporary extraneous)
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']
  
  # 13: SDEs ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude same day extraneous...\n",
      Sys.time()
    ))
  
  # prepare a list of valid rows and initialize variables for convenience
  valid.rows <- valid(data.df, include.temporary.extraneous = TRUE) &
    !data.df$nnte_full  # does not use the "other"
  
  data.df <- data.df[valid.rows, exclude := (function(subj_df) {
    # keep the original indices and exclude to overwrite
    orig_idx <- copy(subj_df$index)
    exclude_all <- copy(subj_df$exclude)
    if (nrow(subj_df) > 0 & any(subj_df$exclude == "Exclude-Temporary-Extraneous-Same-Day")){
      # for convenience
      subj_df[, extraneous := exclude == "Exclude-Temporary-Extraneous-Same-Day"]
      
      step <- "Exclude-SDE-Identical"
      
      # identify duplicate days
      dup_days <- unique(subj_df$agedays[subj_df$extraneous])
      
      # first, exclude identical measurements on the same day
      ide_ids <- c() # where we keep the identical ones
      for (dd in dup_days){
        # count amount of unique values
        s_df <- copy(subj_df[subj_df$agedays == dd,])
        ide_tab <- table(s_df$v)
        if (any(ide_tab > 1)){
          # for each identical, keep only the first one, by id
          # (ordered initially)
          ide_ids <- c(ide_ids, s_df$index[
            as.character(s_df$v) %in% names(ide_tab[ide_tab > 1])
          ][duplicated(
            s_df$v[
              as.character(s_df$v) %in% names(ide_tab[ide_tab > 1])
            ]
          )])
        }
      }
      criteria <- subj_df$index %in% ide_ids
      
      # we're going to update excludes before moving on to the rest of this
      # subject
      exclude_all[criteria] <- step
      
      subj_df <- subj_df[!criteria,]
      
      if (any(criteria)){
        # reevaluate temp same day
        subj_df[, extraneous := temporary_extraneous_infants(subj_df)]
        
        # identify duplicate days
        dup_days <- unique(subj_df$agedays[subj_df$extraneous])
      }
      
      # 13B: identify similar groups
      similar_ids <- c()
      for (dd in dup_days){
        # first, calculate if all SDEs on a single day are similar
        s_df <- copy(subj_df[subj_df$agedays == dd,])
        
        similar <- if (subj_df$param[1] == "WEIGHTKG"){
          max(s_df$v)/min(s_df$v) <= 1.03 & max(s_df$v) - min(s_df$v) <= 2.5
        } else if (subj_df$param[1] == "HEIGHTCM"){
          (max(s_df$v) - min(s_df$v) <  2.541 & min(s_df$v) < 127) |
            (max(s_df$v) - min(s_df$v) <  5.081 & min(s_df$v) >= 127)
        } else { # head circumference
          max(s_df$v) - min(s_df$v) <  1.271
        }
        
        if (similar){
          similar_ids <- c(similar_ids, s_df$index)
        }
      }
      
      # 13C: exclude non similar groups
      
      step <- "Exclude-SDE-All-Exclude"
      # now the rest!
      
      # next, calculate the duplicate ratio -- what proportion of days are
      # duplicated
      dup_ratio <-
        length(unique(subj_df$agedays[duplicated(subj_df$agedays)]))/
        length(unique(subj_df$agedays))
      # also check whether or not any same-days are adjacent -- need 4 at least
      # rolling windows of day differences -- we are looking for 0,x,0
      if (nrow(subj_df) > 3){
        roll <- embed(diff(subj_df$agedays), 3)
        adjacent <- sapply(1:nrow(roll), function(x){
          all(c(roll[x,1] == 0, roll[x,2] != 0, roll[x,3] == 0))
        })
        
        # get age days & rows that had adjacent
        idx_roll <- c(embed(1:nrow(subj_df),4)[adjacent,, drop = FALSE])
        idx_adj <- which(subj_df$agedays %in% subj_df$agedays[idx_roll])
      } else {
        adjacent <- FALSE
        idx_adj <- idx_roll <- c()
      }
      
      # if dup ratio is too high, we exclude all
      # if adjacent, we exclude only those relevant sdes
      # same day extraneous -- for strict (default)
      criteria <-
        # looser criteria:
        # 1.	Exclude for ratio if ratio >1/3
        # 2.	Exclude for adjacent if 2 are adjacent AND ratio >1/4
        # 3.	Exclude for adjacent if 2 are adjacent AND those 2 are the first 2 or last 2 measurements for the subject/parameter
        #   note: if the first or last is sde and it's adjacent, by default
        #   that means it's the first or last two that are adjacent
        # 4.	Exclude for adjacent if 3 or more are adjacent
        if (dup_ratio > (1/3)){
          subj_df$extraneous
        } else if (
          (sum(adjacent) == 1 & dup_ratio > (1/4)) |
          (sum(adjacent) == 1 & c(1) %in% idx_adj) |
          (sum(adjacent) == 1 &
           c(nrow(subj_df)) %in% idx_adj) |
          (sum(adjacent) >= 2)
        ){
          # adjacent sum == 1 means 2 are adjacent
          # only exclude sdes that are adjacent
          subj_df$agedays %in% subj_df$agedays[idx_roll]
        } else {
          rep(FALSE, nrow(subj_df))
        }
      # re-include similar groups
      criteria[subj_df$index %in% similar_ids] <- FALSE
      
      # we're going to update excludes before moving on to the rest of this
      # subject
      exclude_all[orig_idx %in% subj_df$index[criteria]] <- step
      
      subj_df <- subj_df[!criteria,]
      
      if (any(criteria)){
        # reevaluate temp same day
        subj_df[, extraneous := temporary_extraneous_infants(subj_df)]
        
        # identify duplicate days
        dup_days <- unique(subj_df$agedays[subj_df$extraneous])
      }
      
      step_extreme <- "Exclude-SDE-All-Extreme"
      step <- "Exclude-SDE-EWMA"
      # next, check for trivial differences for SDEs on the same day
      # this only works if there are non-sde days
      
      if (any(subj_df$extraneous)){
        # 13D: calculating ewmas
        
        # check for SDEs by EWMA -- alternate calculate excluding other SDEs
        all_sdes <- duplicated(subj_df$agedays) |
          duplicated(subj_df$agedays, fromLast = TRUE)
        
        # first, calculate which exponent we want to put through (pass a different
        # on for each exp)
        tmp <- data.frame(
          "before" = abs(subj_df$agedays - c(NA, subj_df$agedays[1:(nrow(subj_df)-1)])),
          "after" = abs(subj_df$agedays - c(subj_df$agedays[2:(nrow(subj_df))], NA))
        )
        maxdiff <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
        exp_vals <- rep(-1.5, nrow(tmp))
        exp_vals[maxdiff > 365.25] <- -2.5
        exp_vals[maxdiff > 730.5] <- -3.5
        subj_df[, exp_vals := exp_vals]
        
        # calculate the time difference for all values, as well as exponential
        delta <- as_matrix_delta(subj_df$agedays)
        delta <- ifelse(delta == 0, 0, (delta) ^ subj_df$exp_vals)
        
        rem_ids_extreme <- c()
        rem_ids <- c()
        for (dd in dup_days){
          # first, calculate if all SDEs on a single day are similar
          s_df <- copy(subj_df[subj_df$agedays == dd,])
          
          similar <- if (subj_df$param[1] == "WEIGHTKG"){
            max(s_df$v)/min(s_df$v) <= 1.03 & max(s_df$v) - min(s_df$v) <= 2.5
          } else if (subj_df$param[1] == "HEIGHTCM"){
            (max(s_df$v) - min(s_df$v) <  2.541 & min(s_df$v) < 127) |
              (max(s_df$v) - min(s_df$v) <  5.081 & min(s_df$v) >= 127)
          } else { # head circumference
            max(s_df$v) - min(s_df$v) <  1.271
          }
          
          if (!similar){
            # 13E
            
            # calculate dewma for each SDE on this day
            dewma <- sapply(1:nrow(s_df), function(x){
              ind <- subj_df$index == s_df[[x, "index"]]
              
              ewma_res <- sum(subj_df$tbc.sd[!all_sdes]*delta[ind,!all_sdes])/
                sum(delta[ind, !all_sdes])
              
              # delta ewma
              return(s_df$tbc.sd[x] - ewma_res)
            })
            absdewma <- abs(dewma)
            
            if (min(absdewma) > 1){
              rem_ids_extreme <- c(rem_ids_extreme, s_df$index)
            } else {
              de_day <- which.min(abs(dewma))
              
              # keep the value with the lowest EWMA -- do not keep rest
              rem_ids <- c(rem_ids, s_df$index[-de_day])
            }
          }
          # re-include similar groups
          criteria_extreme <- subj_df$index %in% rem_ids_extreme
          criteria <- subj_df$index %in% rem_ids
        }
        # we're going to update excludes before moving on to the rest of this
        # subject
        exclude_all[orig_idx %in% subj_df$index[criteria_extreme]] <- step_extreme
        exclude_all[orig_idx %in% subj_df$index[criteria]] <- step
        
        subj_df <- subj_df[!(criteria | criteria_extreme),]
        
        if (any(c(criteria, criteria_extreme))){
          # reevaluate temp same day
          subj_df[, extraneous := temporary_extraneous_infants(subj_df)]
          
          # identify duplicate days
          dup_days <- unique(subj_df$agedays[subj_df$extraneous])
        }
      }
      
      # all remaining are similar
      step_extreme <- "Exclude-SDE-All-Extreme"
      step <- "Exclude-SDE-One-Day"
      if (any(subj_df$extraneous)){
        # check for SDEs
        all_sdes <- duplicated(subj_df$agedays) |
          duplicated(subj_df$agedays, fromLast = TRUE)
        
        rem_ids <- c()
        rem_ids_extreme <- c()
        for (dd in dup_days){
          # first, calculate if all SDEs on a single day are similar
          s_df <- copy(subj_df[subj_df$agedays == dd,])
          
          # find the DOP (designated other parameter)
          dop <-
            data.df[subjid == s_df$subjid[1] & param == get_dop(s_df$param[1]) &
                      agedays == s_df$agedays[1],]
          
          if (nrow(dop) > 0){
            # 13F
            med_diff <- abs(median(dop$tbc.sd)-s_df$tbc.sd)
            
            if (min(med_diff) > 2){
              rem_ids_extreme <- c(rem_ids_extreme, s_df$index)
            } else {
              de_day <- which.min(med_diff)
              
              # keep the value with the lowest median diff -- do not keep rest
              rem_ids <- c(rem_ids, s_df$index[-de_day])
            }
          } else if (nrow(s_df) > 3){
            #13G
            med_diff <- abs(median(s_df$tbc.sd)-s_df$tbc.sd)
            
            # if even and the middle two values are equidistant
            if (nrow(s_df)%%2 == 0 &
                diff(sort(med_diff)[1:2]) < 1e-8){
              # if the ageday is even, keep the lower z score; otherwise,
              # keep higher
              keep_val <-
                if (s_df$agedays[1]%%2 == 0){
                  which(s_df$v == min(s_df$v[order(med_diff)[1:2]]))
                } else {
                  which(s_df$v == max(s_df$v[order(med_diff)[1:2]]))
                }
              rem_ids <- c(rem_ids, s_df$index[-keep_val])
            }
            
            de_day <- which.min(med_diff)
            # keep the value with the lowest median diff -- do not keep rest
            rem_ids <- c(rem_ids, s_df$index[-de_day])
          } else {
            # 13H
            
            # 2 values, keep based on the age days
            keep_val <-
              if (s_df$agedays[1]%%2 == 0){
                which.min(s_df$v)
              } else {
                which.max(s_df$v)
              }
            rem_ids <- c(rem_ids, s_df$index[-keep_val])
          }
        }
        
        # re-include similar groups
        criteria_extreme <- subj_df$index %in% rem_ids_extreme
        criteria <- subj_df$index %in% rem_ids
        # we're going to update excludes before moving on to the rest of this
        # subject
        exclude_all[orig_idx %in% subj_df$index[criteria_extreme]] <- step_extreme
        exclude_all[orig_idx %in% subj_df$index[criteria]] <- step
      }
    }
    
    # replace all the "exclude temp extraneous" with includes -- other SDEs
    # were removed
    exclude_all[exclude_all == "Exclude-Temporary-Extraneous-Same-Day"] <-
      "Include"
    
    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param), .SDcols = colnames(data.df)]
  
  # 15: moderate EWMA ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude moderate EWMA...\n",
      Sys.time()
    ))
  
  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non pair
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single_pairs <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 2]
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    !data.df$nnte_full & # does not use the "other"
    not_single_pairs
  
  # work in an a function order to encapsulate and not keep all the additional
  # columns
  
  # calculate plus/minus values
  
  # order just for ease later
  data.df <- data.df[order(subjid, param, agedays),]
  data.df <- data.df[valid_set, exclude := (function(df) {
    # 15A: calc plus/minus values
    df[param == "WEIGHTKG", p_plus := 1.05*v]
    df[param == "WEIGHTKG", p_minus := .95*v]
    
    df[param == "HEIGHTCM", p_plus := v+1]
    df[param == "HEIGHTCM", p_minus := v-1]
    
    df[param == "HEADCM", p_plus := v+1]
    df[param == "HEADCM", p_minus := v-1]
    
    #15B: smooth and recenter
    df <- calc_and_recenter_z_scores(df, "p_plus", ref.data.path)
    df <- calc_and_recenter_z_scores(df, "p_minus", ref.data.path)
    
    #15c: exclude agedays = 0 for ht, hc
    #15d: calculate ewma -- run within a subject/parameter
    df[(param != "HEIGHTCM" | param != "HEADCM") & agedays != 0,
       exclude := (function(df_sub) {
         # save initial exclusions to keep track
         ind_all <- copy(df_sub$index)
         exclude_all <- copy(df_sub$exclude)
         
         testing <- TRUE
         
         while (testing & nrow(df_sub) >= 3){
           df_sub[, (ewma.fields) := as.double(NaN)]
           
           # first, calculate which exponent we want to put through (pass a different
           # on for each exp)
           # subset df to only valid rows
           
           tmp <- data.frame(
             "before" = abs(df_sub$agedays - c(NA, df_sub$agedays[1:(nrow(df_sub)-1)])),
             "after" = abs(df_sub$agedays - c(df_sub$agedays[2:(nrow(df_sub))], NA))
           )
           maxdiff <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
           exp_vals <- rep(-1.5, nrow(tmp))
           exp_vals[maxdiff > 365.25] <- -2.5
           exp_vals[maxdiff > 730.5] <- -3.5
           df_sub[, exp_vals := exp_vals]
           
           # calculate ewma
           df_sub[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
           df_sub[, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]
           
           # calculate dewma
           df_sub[, `:=`(
             dewma.all = tbc.sd - ewma.all,
             dewma.before = tbc.sd - ewma.before,
             dewma.after = tbc.sd - ewma.after,
             
             c.dewma.all = ctbc.sd - c.ewma.all
           )]
           
           # 15E: calculate prior and next differences
           df_sub[, tbc_diff_next := tbc.sd - c(tbc.sd[2:nrow(df_sub)], NA)]
           df_sub[, tbc_diff_prior := tbc.sd - c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
           
           df_sub[, tbc_diff_plus_next := tbc.p_plus - c(tbc.sd[2:nrow(df_sub)], NA)]
           df_sub[, tbc_diff_plus_prior :=
                    tbc.p_plus- c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
           
           df_sub[, tbc_diff_minus_next := tbc.p_minus - c(tbc.sd[2:nrow(df_sub)], NA)]
           df_sub[, tbc_diff_minus_prior :=
                    tbc.p_minus- c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
           
           # 15F: all exclusion criteria
           df_sub[,
                  addcrithigh :=
                    dewma.before > 1 & dewma.after > 1 &
                    ((tbc_diff_next > 1 & tbc_diff_plus_next > 1 & tbc_diff_minus_next > 1) | is.na(tbc_diff_next)) &
                    ((tbc_diff_prior > 1 & tbc_diff_plus_prior > 1 & tbc_diff_minus_prior > 1) | is.na(tbc_diff_prior))
           ]
           df_sub[,
                  addcritlow :=
                    dewma.before < -1 & dewma.after < -1 &
                    ((tbc_diff_next < -1 & tbc_diff_plus_next < -1 & tbc_diff_minus_next < -1) | is.na(tbc_diff_next)) &
                    ((tbc_diff_prior < -1 & tbc_diff_plus_prior < -1 & tbc_diff_minus_prior < -1) | is.na(tbc_diff_prior))
           ]
           
           #15G: find all the values for the DOP
           # find the comparison -- making sure to keep only valid rows
           compare_df <- data.df[subjid == df_sub$subjid[1] &
                                   param == get_dop(df_sub$param[1]) &
                                   exclude == "Include",]
           if (nrow(compare_df) > 0){
             df_sub[compare_df, tbc_dop := i.tbc.sd,
                    on = c("agedays")]
             df_sub[is.na(tbc_dop), tbc_dop := median(compare_df$tbc.sd)]
           } else {
             df_sub[, tbc_dop := NA]
           }
           
           #15H: all exclusions
           df_sub$rowind <- 1:nrow(df_sub)
           
           #a
           df_sub[-c(1, nrow(df_sub)),][
             dewma.all > 1 & (c.dewma.all > 1 | is.na(c.dewma.all)) & addcrithigh,
             exclude := "Exclude-EWMA2-middle"]
           df_sub[-c(1, nrow(df_sub)),][
             dewma.all < -1 & (c.dewma.all < -1 | is.na(c.dewma.all)) & addcritlow,
             exclude := "Exclude-EWMA2-middle"]
           
           #b
           df_sub[agedays == 0 & agedays[2] < 365.25 &
                    dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-birth-WT"]
           df_sub[agedays == 0 & agedays[2] < 365.25 &
                    dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-birth-WT"]
           
           #c
           df_sub[agedays == 0 & agedays[2] >= 365.25 &
                    dewma.all > 4 & (c.dewma.all > 4 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-birth-WT-ext"]
           df_sub[agedays == 0 & agedays[2] >= 365.25 &
                    dewma.all < -4 & (c.dewma.all < -4 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-birth-WT-ext"]
           
           #d
           df_sub[rowind == 1 & agedays[1] != 0 & agedays[2]-agedays[1] < 365.25 &
                    dewma.all > 2 & (c.dewma.all > 2 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-first"]
           df_sub[rowind == 1 & agedays[1] != 0 & agedays[2]-agedays[1] < 365.25 &
                    dewma.all < -2 & (c.dewma.all < -2 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-first"]
           
           # e
           df_sub[rowind == 1 & agedays[1] != 0 & agedays[2]-agedays[1] >= 365.25 &
                    dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-first-ext"]
           df_sub[rowind == 1 & agedays[1] != 0 & agedays[2]-agedays[1] >= 365.25 &
                    dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-first-ext"]
           
           # f
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] < 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) < 2 &
                    dewma.all > 2 & (c.dewma.all > 2 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-last"]
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] < 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) < 2 &
                    dewma.all < -2 & (c.dewma.all < -2 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-last"]
           
           # g
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] < 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) >= 2 &
                    dewma.all > abs(tbc.sd[nrow(df_sub)-1]) &
                    (c.dewma.all > 3 | is.na(c.dewma.all)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-last-high"]
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] < 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) >= 2 &
                    dewma.all < abs(tbc.sd[nrow(df_sub)-1]) &
                    (c.dewma.all < -3 | is.na(c.dewma.all)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-last-high"]
           
           # h
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] >= 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) < 2 &
                    dewma.all > 3 &
                    (c.dewma.all > 3 | is.na(c.dewma.all)) &
                    (tbc.sd - tbc_dop > 4 | is.na(tbc_dop)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-last-ext"]
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] >= 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) < 2 &
                    dewma.all < -3 &
                    (c.dewma.all <- 3 | is.na(c.dewma.all)) &
                    (tbc.sd - tbc_dop < -4 | is.na(tbc_dop)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-last-ext"]
           
           # i
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] >= 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) >= 2 &
                    dewma.all > (1+abs(tbc.sd[nrow(df_sub)-1])) &
                    (c.dewma.all > 3 | is.na(c.dewma.all)) &
                    (tbc.sd - tbc_dop > 4 | is.na(tbc_dop)) &
                    addcrithigh
                  , exclude := "Exclude-EWMA2-last-ext-high"]
           df_sub[rowind == nrow(df_sub) &
                    agedays[nrow(df_sub)]-agedays[nrow(df_sub)-1] >= 365.25*2 &
                    abs(tbc.sd[nrow(df_sub)-1]) >= 2 &
                    dewma.all < (-1-abs(tbc.sd[nrow(df_sub)-1])) &
                    (c.dewma.all <- 3 | is.na(c.dewma.all)) &
                    (tbc.sd - tbc_dop < -4 | is.na(tbc_dop)) &
                    addcritlow
                  , exclude := "Exclude-EWMA2-last-ext-high"]
           
           # figure out if any of the exclusions hit
           count_exclude <- sum(df_sub$exclude != "Include")
           if (count_exclude > 0){
             df_sub[, abssum := abs(tbc.sd + dewma.all)]
             
             # choose the highest abssum for exclusion
             idx <- df_sub$index[which.max(df_sub[df_sub$exclude != "Include",
                                                  abssum])]
             
             exclude_all[ind_all == idx] <- df_sub[index == idx, exclude]
             
             #set up to continue on
             testing <- TRUE
             
             df_sub <- df_sub[index != idx, ]
           } else {
             testing <- FALSE
           }
         }
         
         return(exclude_all)
       })(copy(.SD)), , by = .(subjid, param), .SDcols = colnames(df)]
    
    return(df$exclude)
  })(copy(.SD)), .SDcols = colnames(data.df)]
  
  # 16: moderate EWMA for birth HT and HC ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude moderate EWMA for birth height and head circumference...\n",
      Sys.time()
    ))
  
  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non weight
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 1]
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    !data.df$nnte_full & # does not use the "other"
    not_single &
    data.df$param != "WEIGHTKG" # do not run on weight
  
  # work in an a function order to encapsulate and not keep all the additional
  # columns
  
  # NOTE: check NNTE na's
  
  # calculate plus/minus values
  
  # order just for ease later
  data.df <- data.df[order(subjid, param, agedays),]
  data.df <- data.df[valid_set, exclude := (function(df) {
    # =calc plus/minus values
    df[param == "HEIGHTCM", p_plus := v+1]
    df[param == "HEIGHTCM", p_minus := v-1]
    
    df[param == "HEADCM", p_plus := v+1]
    df[param == "HEADCM", p_minus := v-1]
    
    #smooth and recenter
    df <- calc_and_recenter_z_scores(df, "p_plus", ref.data.path)
    df <- calc_and_recenter_z_scores(df, "p_minus", ref.data.path)
    
    #16a: calculate ewma -- run within a subject/parameter
    # inly evaluate on subject/parameters with a birth ageday
    df <- df[subjid %in% subjid[agedays == 0],
             exclude := (function(df_sub) {
               # save initial exclusions to keep track
               ind_all <- copy(df_sub$index)
               exclude_all <- copy(df_sub$exclude)
               
               testing <- TRUE
               
               while (testing & nrow(df_sub) >= 3){
                 df_sub[, (ewma.fields) := as.double(NaN)]
                 
                 # first, calculate which exponent we want to put through (pass a different
                 # on for each exp)
                 # subset df to only valid rows
                 
                 tmp <- data.frame(
                   "before" = abs(df_sub$agedays - c(NA, df_sub$agedays[1:(nrow(df_sub)-1)])),
                   "after" = abs(df_sub$agedays - c(df_sub$agedays[2:(nrow(df_sub))], NA))
                 )
                 maxdiff <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
                 exp_vals <- rep(-1.5, nrow(tmp))
                 exp_vals[maxdiff > 365.25] <- -2.5
                 exp_vals[maxdiff > 730.5] <- -3.5
                 df_sub[, exp_vals := exp_vals]
                 
                 # calculate ewma
                 df_sub[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
                 df_sub[, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]
                 
                 # calculate dewma
                 df_sub[, `:=`(
                   dewma.all = tbc.sd - ewma.all,
                   dewma.before = tbc.sd - ewma.before,
                   dewma.after = tbc.sd - ewma.after,
                   
                   c.dewma.all = ctbc.sd - c.ewma.all
                 )]
                 
                 # 16B: calculate prior and next differences
                 df_sub[, tbc_diff_next := tbc.sd - c(tbc.sd[2:nrow(df_sub)], NA)]
                 df_sub[, tbc_diff_prior := tbc.sd - c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
                 
                 df_sub[, tbc_diff_plus_next := tbc.p_plus - c(tbc.sd[2:nrow(df_sub)], NA)]
                 df_sub[, tbc_diff_plus_prior :=
                          tbc.p_plus- c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
                 
                 df_sub[, tbc_diff_minus_next := tbc.p_minus - c(tbc.sd[2:nrow(df_sub)], NA)]
                 df_sub[, tbc_diff_minus_prior :=
                          tbc.p_minus- c(NA, tbc.sd[1:(nrow(df_sub)-1)])]
                 
                 # 16C: all exclusion criteria
                 df_sub[,
                        addcrithigh :=
                          dewma.before > 1 & dewma.after > 1 &
                          ((tbc_diff_next > 1 & tbc_diff_plus_next > 1 & tbc_diff_minus_next > 1) | is.na(tbc_diff_next)) &
                          ((tbc_diff_prior > 1 & tbc_diff_plus_prior > 1 & tbc_diff_minus_prior > 1) | is.na(tbc_diff_prior))
                 ]
                 df_sub[,
                        addcritlow :=
                          dewma.before < -1 & dewma.after < -1 &
                          ((tbc_diff_next < -1 & tbc_diff_plus_next < -1 & tbc_diff_minus_next < -1) | is.na(tbc_diff_next)) &
                          ((tbc_diff_prior < -1 & tbc_diff_plus_prior < -1 & tbc_diff_minus_prior < -1) | is.na(tbc_diff_prior))
                 ]
                 
                 #16D: all exclusions
                 df_sub$rowind <- 1:nrow(df_sub)
                 
                 #a
                 df_sub[agedays == 0 & agedays[2] < 365.25 &
                          dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) &
                          addcrithigh
                        , exclude := "Exclude-EWMA2-birth-HT-HC"]
                 df_sub[agedays == 0 & agedays[2] < 365.25 &
                          dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) &
                          addcritlow
                        , exclude := "Exclude-EWMA2-birth-HT-HC"]
                 
                 #b
                 df_sub[agedays == 0 & agedays[2] >= 365.25 &
                          dewma.all > 4 & (c.dewma.all > 4 | is.na(c.dewma.all)) &
                          addcrithigh
                        , exclude := "Exclude-EWMA2-birth-HT-HC-ext"]
                 df_sub[agedays == 0 & agedays[2] >= 365.25 &
                          dewma.all < -4 & (c.dewma.all < -4 | is.na(c.dewma.all)) &
                          addcritlow
                        , exclude := "Exclude-EWMA2-birth-HT-HC-ext"]
                 
                 # figure out if any of the exclusions hit
                 count_exclude <- sum(df_sub$exclude != "Include")
                 if (count_exclude > 0){
                   df_sub[, abssum := abs(tbc.sd + dewma.all)]
                   
                   # choose the highest abssum for exclusion
                   idx <- df_sub$index[which.max(df_sub[df_sub$exclude != "Include",
                                                        abssum])]
                   
                   exclude_all[ind_all == idx] <- df_sub[index == idx, exclude]
                   
                   #set up to continue on
                   testing <- TRUE
                   
                   df_sub <- df_sub[index != idx, ]
                 } else {
                   testing <- FALSE
                 }
               }
               return(exclude_all)
             })(copy(.SD)), , by = .(subjid, param), .SDcols = colnames(df)]
    
    return(df$exclude)
  })(copy(.SD)), .SDcols = colnames(data.df)]
  
  # 17: raw differences ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude raw differences...\n",
      Sys.time()
    ))
  
  # read in tanner data
  tanner_ht_vel_rev_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "tanner_ht_vel.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "tanner_ht_vel.csv.gz")
  )
  
  tanner.ht.vel.rev <- fread(tanner_ht_vel_rev_path)
  
  setnames(tanner.ht.vel.rev,
           colnames(tanner.ht.vel.rev),
           gsub('_', '.', colnames(tanner.ht.vel.rev)))
  setkey(tanner.ht.vel.rev, sex, tanner.months)
  
  # read in the who height data
  who_max_ht_vel_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "who_ht_maxvel_3sd.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "who_ht_maxvel_3sd.csv.gz")
  )
  
  who_ht_vel_3sd_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "who_ht_vel_3sd.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "who_ht_vel_3sd.csv.gz")
  )
  
  who.max.ht.vel <- fread(who_max_ht_vel_path)
  who.ht.vel <- fread(who_ht_vel_3sd_path)
  setkey(who.max.ht.vel, sex, whoagegrp_ht)
  setkey(who.ht.vel, sex, whoagegrp_ht)
  who.ht.vel <- as.data.table(dplyr::full_join(who.ht.vel, who.max.ht.vel, by =
                                                 c('sex', 'whoagegrp_ht')))
  
  setnames(who.ht.vel, colnames(who.ht.vel), gsub('_', '.', colnames(who.ht.vel)))
  setkey(who.ht.vel, sex, whoagegrp.ht)
  
  # read in who hc data
  who_max_hc_vel_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "who_hc_maxvel_3sd_infants.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "who_hc_maxvel_3sd_infants.csv.gz")
  )
  
  who_hc_vel_3sd_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "who_hc_vel_3sd_infants.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "who_hc_vel_3sd_infants.csv.gz")
  )
  
  who.max.hc.vel <- fread(who_max_hc_vel_path)
  who.hc.vel <- fread(who_hc_vel_3sd_path)
  setkey(who.max.hc.vel, sex, whoagegrp_ht)
  setkey(who.hc.vel, sex, whoagegrp_ht)
  who.hc.vel <- as.data.table(dplyr::full_join(who.hc.vel, who.max.hc.vel, by =
                                                 c('sex', 'whoagegrp_ht')))
  
  setnames(who.hc.vel, colnames(who.hc.vel), gsub('_', '.', colnames(who.hc.vel)))
  setkey(who.hc.vel, sex, whoagegrp.ht)
  
  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non weight
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 1]
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    !data.df$nnte_full & # does not use the "other"
    not_single &
    data.df$param != "WEIGHTKG" # do not run on weight
  
  # order just for ease later
  data.df <- data.df[order(subjid, param, agedays),]
  data.df <- data.df[valid_set, exclude := (function(df) {
    # work inside a closure to drop added column values
    
    # save initial exclusions to keep track
    ind_all <- copy(df$index)
    exclude_all <- copy(df$exclude)
    
    testing <- TRUE
    
    while (testing & nrow(df) > 1){
      # sort df since it got reordered with keys
      df <- df[order(agedays),]
      
      # 17A
      df[, d_agedays := dplyr::lead(agedays) - agedays]
      
      # 17E -- only applies to height
      if (df$param[1] == "HEIGHTCM"){
        # 17B
        df[, tanner.months := 6+12*(round(.5*(agedays + dplyr::lead(agedays))/365.25))]
        # set tanner to missing for smaller values
        df[(agedays/30.4375) < 30, tanner.months := NA]
        
        # merge with the tanner info
        setkey(df, sex, tanner.months)
        df <- tanner.ht.vel.rev[df]
        
        # 17C
        
        # a
        df[max.ht.vel < 2.54, max.ht.vel := 2.54]
        df[d_agedays > 2*30.4375 & max.ht.vel < 2*2.54, max.ht.vel := 2*2.54]
        df[d_agedays > .5*365.25 & max.ht.vel < 4*2.54, max.ht.vel := 4*2.54]
        df[d_agedays > 365.25 & max.ht.vel < 8*2.54, max.ht.vel := 8*2.54]
        
        # b
        df[d_agedays < 365.25, mindiff := .5*min.ht.vel*(d_agedays/365.25)-3 ]
        df[d_agedays > 365.25, mindiff := .5*min.ht.vel-3 ]
        df[d_agedays < 365.25,
           maxdiff := 2*min.ht.vel*(d_agedays/365.25)^1.5 + 5.5 ]
        df[d_agedays > 365.25,
           maxdiff := 2*min.ht.vel*(d_agedays/365.25)^0.33 + 5.5 ]
        
        # 17D
        # generate the who age group variable
        df[, whoagegrp.ht := round(agedays/30.4375)]
        df[whoagegrp.ht > 24 | dplyr::lead(whoagegrp.ht) > 24,
           whoagegrp.ht := NA]
        
        
        # 17E
        df[d_agedays >= 20 & d_agedays < 46, whoinc.age.ht := 1]
        df[d_agedays >= 46 & d_agedays < 76, whoinc.age.ht := 2]
        df[d_agedays >= 76 & d_agedays < 107, whoinc.age.ht := 3]
        df[d_agedays >= 107 & d_agedays < 153, whoinc.age.ht := 4]
        df[d_agedays >= 153 & d_agedays < 199, whoinc.age.ht := 6]
        
        # update the edge intervals
        df[d_agedays < 20, whoinc.age.ht := 1]
        df[d_agedays > 199, d_agedays := 200]
        df[d_agedays == 200, whoinc.age.ht := 6]
        
        # 17F
        # merge with WHO
        # add the column name we want to grab
        m_who_ht_vel <- merge(df, who.ht.vel, by = c("sex", "whoagegrp.ht"),
                              all.x = TRUE, sort = FALSE)
        for (i in unique(df$whoinc.age.ht[!is.na(df$whoinc.age.ht) &
                                          !is.na(df$whoagegrp.ht)])){
          sub_m_who_ht_vel <- m_who_ht_vel[whoinc.age.ht == i,]
          
          cn <- paste0("whoinc.", i, ".ht")
          df[whoinc.age.ht == i,
             who_mindiff_ht := as.numeric(sub_m_who_ht_vel[, get(cn)])]
          cn <- paste0("max.whoinc.", i, ".ht")
          df[whoinc.age.ht == i,
             who_maxdiff_ht := as.numeric(sub_m_who_ht_vel[, get(cn)])]
        }
        # if there are none, preallocate for ease
        if (length(unique(df$whoinc.age.ht[!is.na(df$whoinc.age.ht) &
                                           !is.na(df$whoagegrp.ht)])) < 1){
          df[, who_mindiff_ht := NA_real_]
          df[, who_maxdiff_ht := NA_real_]
        }
        df[, who_mindiff_ht := as.numeric(who_mindiff_ht)]
        df[, who_maxdiff_ht := as.numeric(who_maxdiff_ht)]
        
        # 17G
        df[d_agedays < whoinc.age.ht*30.4375, who_mindiff_ht :=
             who_mindiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        df[d_agedays > whoinc.age.ht*30.4375, who_maxdiff_ht :=
             who_maxdiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        
        df[d_agedays < 9*30.4375, who_mindiff_ht := who_mindiff_ht*.5-3]
        df[d_agedays < 9*30.4375, who_maxdiff_ht := who_maxdiff_ht*2+3]
        
        # 17H
        # tanner is implicit
        # greater than 9 months, use tanner if available, otherwise who
        df[(d_agedays < 9*30.4375 | is.na(min.ht.vel)) &
             !is.na(who_mindiff_ht), mindiff := who_mindiff_ht]
        df[(d_agedays < 9*30.4375 | is.na(min.ht.vel)) &
             !is.na(who_mindiff_ht), maxdiff := who_maxdiff_ht]
        # otherwise, fill in
        df[is.na(mindiff), mindiff := 3]
        # for birth measurements, add allowance of 1.5cm
        df[agedays == 0, mindiff := mindiff - 1.5]
        df[agedays == 0, maxdiff := maxdiff + 1.5]
        
        # 17I
        # sort df since it got reordered with keys
        df <- df[order(agedays),]
        df[, mindiff_prev := dplyr::lag(mindiff)]
        df[, maxdiff_prev := dplyr::lag(maxdiff)]
      } else { # head circumference
        # 17D
        # generate the who age group variable
        df[, whoagegrp.ht := round(agedays/30.4375)]
        df[whoagegrp.ht > 24 | dplyr::lead(whoagegrp.ht) > 24,
           whoagegrp.ht := NA]
        
        # 17J
        df[d_agedays >= 46 & d_agedays < 76, whoinc.age.hc := 2]
        df[d_agedays >= 76 & d_agedays < 107, whoinc.age.hc := 3]
        df[d_agedays >= 107 & d_agedays < 153, whoinc.age.hc := 4]
        df[d_agedays >= 153 & d_agedays < 199, whoinc.age.hc := 6]
        
        # update the edge intervals to missing
        df[d_agedays < 46 | d_agedays > 199, whoinc.age.hc := NA]
        
        # 17K
        # merge with WHO
        # add the column name we want to grab
        
        m_who_hc_vel <- merge(df, who.hc.vel, by = c("sex", "whoagegrp.ht"),
                              all.x = TRUE, sort = FALSE)
        for (i in unique(df$whoinc.age.hc[!is.na(df$whoinc.age.hc) &
                                          !is.na(df$whoagegrp.hc)])){
          sub_m_who_hc_vel <- m_who_hc_vel[whoinc.age.ht == i,]
          
          cn <- paste0("whoinc.", i, ".ht")
          df[whoinc.age.hc == i,
             who_mindiff_hc := as.numeric(sub_m_who_ht_vel[, get(cn)])]
          cn <- paste0("max.whoinc.", i, ".ht")
          df[whoinc.age.hc == i,
             who_maxdiff_hc := as.numeric(sub_m_who_ht_vel[, get(cn)])]
        }
        # if there are none, preallocate for ease
        if (length(unique(df$whoinc.age.hc[!is.na(df$whoinc.age.hc) &
                                           !is.na(df$whoagegrp.hc)])) < 1){
          df[, who_mindiff_hc := NA_real_]
          df[, who_maxdiff_hc := NA_real_]
        }
        df[, who_mindiff_hc := as.numeric(who_mindiff_hc)]
        df[, who_maxdiff_hc := as.numeric(who_maxdiff_hc)]
        
        # 17L
        df[d_agedays < whoinc.age.hc*30.4375, who_mindiff_hc :=
             who_mindiff_hc * d_agedays/(whoinc.age.hc*30.4375)]
        df[d_agedays > whoinc.age.hc*30.4375, who_maxdiff_hc :=
             who_maxdiff_hc * d_agedays/(whoinc.age.hc*30.4375)]
        
        df[d_agedays < 9*30.4375, who_mindiff_hc := who_mindiff_hc*.5-1.5]
        df[d_agedays < 9*30.4375, who_maxdiff_hc := who_maxdiff_hc*2+1.5]
        
        # 17M
        # use who
        df[, mindiff := who_mindiff_hc]
        df[, maxdiff := who_maxdiff_hc]
        # otherwise, fill in
        df[is.na(mindiff), mindiff := -1.5]
        # for birth measurements, add allowance of .5cm
        df[agedays == 0, mindiff := mindiff - .5]
        df[agedays == 0, maxdiff := maxdiff + .5]
        
        # 17N
        # sort df since it got reordered with keys
        df <- df[order(agedays),]
        df[, mindiff_prev := dplyr::lag(mindiff)]
        df[, maxdiff_prev := dplyr::lag(maxdiff)]
      }
      
      # 17O: generate ewma
      df[, (ewma.fields) := as.double(NaN)]
      
      # first, calculate which exponent we want to put through (pass a different
      # on for each exp)
      # subset df to only valid rows
      tmp <- data.frame(
        "before" = abs(df$agedays - c(NA, df$agedays[1:(nrow(df)-1)])),
        "after" = abs(df$agedays - c(df$agedays[2:(nrow(df))], NA))
      )
      maxdiff_e <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
      exp_vals <- rep(-1.5, nrow(tmp))
      exp_vals[maxdiff_e > 365.25] <- -2.5
      exp_vals[maxdiff_e > 730.5] <- -3.5
      df[, exp_vals := exp_vals]
      
      # calculate ewma
      df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
      
      # calculate dewma
      df[, `:=`(
        dewma.all = tbc.sd - ewma.all,
        dewma.before = tbc.sd - ewma.before,
        dewma.after = tbc.sd - ewma.after
      )]
      
      # add differences for convenience
      df[, diff_prev := (v-dplyr::lag(v))]
      df[, diff_next := (dplyr::lead(v)-v)]
      
      if (nrow(df) > 2){
        # 17P/R: identify pairs and calculate exclusions
        df[, pair := (v-dplyr::lag(v)) < mindiff_prev |
             (dplyr::lead(v)-v) < mindiff |
             (v-dplyr::lag(v)) > maxdiff_prev | (dplyr::lead(v)-v) > maxdiff
        ]
        df[is.na(pair), pair := FALSE]
        df[pair & abs(dewma.before) > dplyr::lag(abs(dewma.after)),
           bef.g.aftm1 := TRUE]
        df[pair & abs(dewma.after) > dplyr::lead(abs(dewma.before)),
           aft.g.aftm1 := TRUE]
        
        # Q
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prev & bef.g.aftm1, val_excl := "Exclude-Min-diff"]
        df[diff_next < mindiff & aft.g.aftm1, val_excl := "Exclude-Min-diff"]
        df[diff_prev > maxdiff_prev & bef.g.aftm1, val_excl := "Exclude-Max-diff"]
        df[diff_next > maxdiff & aft.g.aftm1, val_excl := "Exclude-Max-diff"]
        
        
        
      } else { # only 2 values
        # 17Q/R -- exclusions for pairs
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prev & abs(tbc.sd) > dplyr::lag(abs(tbc.sd)),
           val_excl := "Exclude-Min-diff"]
        df[diff_next < mindiff & abs(tbc.sd) > dplyr::lead(abs(tbc.sd)),
           val_excl := "Exclude-Min-diff"]
        df[diff_prev > maxdiff_prev & abs(tbc.sd) > dplyr::lag(abs(tbc.sd)),
           val_excl := "Exclude-Max-diff"]
        df[diff_next > maxdiff &  abs(tbc.sd) > dplyr::lead(abs(tbc.sd)),
           val_excl := "Exclude-Max-diff"]
      }
      
      # figure out if any of the exclusions hit
      count_exclude <- sum(df$val_excl != "Include")
      if (count_exclude > 0){
        if (nrow(df) > 2){
          df[, absval := abs(dewma.all)]
        } else {
          df[, absval := abs(tbc.sd)]
        }
        
        # choose the highest abssum for exclusion
        idx <- df$index[which.max(df[df$val_excl != "Include",
                                     absval])]
        
        exclude_all[ind_all == idx] <- df[index == idx, val_excl]
        
        #set up to continue on
        if (count_exclude > 1){
          testing <- TRUE
          
          df <- df[index != idx, ]
        } else {
          testing <- FALSE
        }
      } else {
        testing <- FALSE
      }
    }
    
    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param), .SDcols = colnames(data.df)]
  
  # 19: 1 or 2 measurements ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude 1 or 2 measurements...\n",
      Sys.time()
    ))
  
  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non weight
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  only_single_pairs <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp <= 2]
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    !data.df$nnte_full & # does not use the "other"
    only_single_pairs
  
  # order just for ease later
  data.df <- data.df[order(subjid, param, agedays),]
  data.df <- data.df[valid_set, exclude := (function(df) {
    # save initial exclusions to keep track
    ind_all <- copy(df$index)
    exclude_all <- copy(df$exclude)
    
    # 19A: is it a single or a pair?
    is_single <- nrow(df) == 1
    
    # find the DOP (designated other parameter)
    dop <-
      data.df[subjid == df$subjid[1] & param == get_dop(df$param[1]) &
                exclude == "Include",]
    
    # 19D: calculate the voi comparison
    if (nrow(dop) > 0){
      for (i in 1:nrow(df)){
        comp_val <-
          if (df$agedays[i] %in% dop$agedays){
            abs(dop[dop$agedays == df$agedays[i], tbc.sd] - df[i, tbc.sd])
          } else {
            abs(median(dop$tbc.sd) - df[i, tbc.sd])
          }
        
        df[i, comp_diff := comp_val]
      }
    } else {
      df[, comp_diff := rep(NA, nrow(df))]
    }
    
    # 19B/C: for pairs, calculate info
    if (!is_single){
      diff_tbc.sd <- df[2, tbc.sd] - df[1, tbc.sd]
      diff_ctbc.sd <- df[2, ctbc.sd] - df[1, ctbc.sd]
      diff_agedays <- df[2, agedays] - df[1, agedays]
      
      abs_tbd.sd <- abs(df[1, tbc.sd])
      abs_ctbd.sd <- abs(df[1, ctbc.sd])
      
      med_dop <-
        if (nrow(dop) > 0){
          median(dop$tbc.sd)
        } else {
          NA
        }
      med_cdop <-
        if (nrow(dop) > 0){
          median(dop$ctbc.sd)
        } else {
          NA
        }
      
      # 19E: which is larger
      max_ind <- if (!all(is.na(df$comp_diff))){
        which.max(df$comp_diff)
      } else {
        which.max(abs(df$tbc.sd))
      }
      
      # 19F/G: which exclusion
      if (diff_tbc.sd > 4 & (diff_ctbc.sd > 4 | is.na(diff_ctbc.sd)) &
          diff_agedays >=365.25){
        df[max_ind, exclude := "Exclude-2-meas->1-year"]
      } else if (diff_tbc.sd > 2.5 & (diff_ctbc.sd > 2.5 | is.na(diff_ctbc.sd)) &
                 diff_agedays < 365.25){
        df[max_ind, exclude := "Exclude-2-meas-<1-year"]
      }
      
      # save the results
      exclude_all <- df$exclude
      
      # 19H
      # if one needs to get removed, we want to reevaluate as a single
      df <- df[exclude == "Include",]
    }
    
    if (nrow(df) == 1){
      df[(abs(tbc.sd) > 3 & !is.na(comp_diff) & comp_diff > 5) |
           (abs(tbc.sd) > 5 & is.na(comp_diff)),
         exclude := "Exclude-1-meas"]
    }
    
    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param), .SDcols = colnames(data.df)]
  
  # 21: error load ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude error load...\n",
      Sys.time()
    ))
  
  valid_set <- !data.df$nnte_full
  
  data.df[valid_set,
          err_ratio := sum(!exclude %in%
                             c("Include",
                               "Exclude-SDE-Identical",
                               "Exclude-SDE-All-Exclude",
                               "Exclude-SDE-All-Extreme",
                               "Exclude-SDE-EWMA",
                               "Exclude-SDE-All-Extreme",
                               "Exclude-SDE-One-Day",
                               "Exclude-Carried-Forward",
                               "Exclude-1-CF-deltaZ-<0.05",
                               "Exclude-1-CF-deltaZ-<0.1-wholehalfimp",
                               "Exclude-Teen-2-plus-CF-deltaZ-<0.05",
                               "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp"
                             ))/.N,
          by = c("subjid", "param")]
  
  data.df[valid_set & err_ratio > .4 & exclude == "Include",
          exclude := "Exclude-Error-load"]
  
  
  # end ----
  
  if (!quietly)
    cat(sprintf("[%s] Completed Batch #%s...\n", Sys.time(), data.df$batch[1]))
  if (!quietly & parallel)
    sink()
  
  return(data.df[j = .(line, exclude, tbc.sd, param)]) #debugging
}
