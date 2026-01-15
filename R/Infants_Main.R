################################################################################
# GROWTHCLEANR INFANT ALGORITHM - R IMPLEMENTATION
################################################################################
#
# PURPOSE:
#   Identifies and flags implausible pediatric growth measurements (weight,
#   height, head circumference) from electronic health records for children
#   ages 0-20 years, with enhanced methods for infants 0-2 years.
#
# VERSION: 2026-01-09
#
# AUTHOR: Carrie Daymont, Penn State College of Medicine
#
# PARALLEL IMPLEMENTATION:
#   This R implementation is algorithmically equivalent to the Stata version
#   but optimized for parallel processing of large datasets.
#

################################################################################
# CITATION
################################################################################
#
# If using this algorithm, please cite:
#   Daymont C, et al. (2017). Automated identification of implausible values
#   in growth data from pediatric electronic health records. JAMIA, 24(6):1080-1087.
#
################################################################################

################################################################################
# ALGORITHM STRUCTURE
################################################################################
#
# The algorithm proceeds through 22 steps:
#
#   STEP 0:  Dataset preparation and variable setup
#   STEP 1:  Missing data and initialization
#   STEP 2:  Z-score calculation (WHO/CDC) and recentering
#   STEP 2b: Gestational age correction for premature/small infants
#   STEP 3:  Unit error recovery (optional)
#   STEP 4:  Swapped parameters detection
#   STEP 5:  Temporary same-day extraneous (SDE) resolution
#   STEP 6:  Carried forward identification
#   STEP 7:  Absolute biologically implausible values (BIV)
#   STEP 9:  Evil Twins (adjacent extreme values)
#   STEP 11: EWMA-based extreme value detection (EWMA1)
#   STEP 13: Final SDE resolution (one-day and EWMA-based)
#   STEP 15: Moderate EWMA exclusions (EWMA2)
#   STEP 17: Height velocity checks
#   STEP 19: Pairs and singles evaluation
#   STEP 21: Error load assessment
#   STEP 22: Output preparation
#
# Note: Step numbers are not consecutive due to alignment with parallel
# Stata algorithm development.
#
# KEY CONCEPTS:
#   - SDE (Same-Day Extraneous): Multiple measurements of same parameter on
#     same day. Algorithm selects most plausible value.
#   - Carried Forward: Identical consecutive values, likely copied from prior
#     visit rather than independently measured.
#   - EWMA (Exponentially Weighted Moving Average): Uses weighted average of
#     surrounding measurements to identify outliers, giving more weight to
#     temporally closer values.
#   - tbc.sd (To-Be-Cleaned SD score): Recentered z-score used throughout
#   - DOP (Designated Other Parameter): For weight->height, for height/HC->weight
#   - Evil Twins: Multiple consecutive extreme values that fooled earlier steps
#
################################################################################
# DATA STRUCTURE NOTES
################################################################################
#
# R IMPLEMENTATION SPECIFICS:
#   - Uses data.table for efficient processing of large datasets
#   - Parallel processing via foreach/doParallel when parallel=TRUE
#   - Subjects are split into batches for parallel computation
#   - Each batch processes independently, then results are combined
#
# VARIABLE NAMING DIFFERENCES FROM STATA:
#   - Stata: swtz, shtz, shcz -> R: sd.orig (smoothed original z-score)
#   - Stata: tbcwtz, tbchtz, tbchcz -> R: tbc.sd (to-be-cleaned SD score)
#   - Stata: exc_wt, exc_ht, exc_hc -> R: exclude (single character variable)
#   - Stata: subjid_`p' -> R: param-based subsetting
#
# SORT ORDER:
#   Critical for deterministic results. Primary sort:
#     data.table::setkey(data.df, subjid, param, agedays, id)
#
#   The 'id' field is required and should be unique per row. For datasets
#   with multiple measurements on the same day, 'id' serves as a tiebreaker
#   to ensure consistent ordering across runs. If you know the timing of values
#   on individual days, you may wish to assign higher ids to later values.
#
################################################################################
# SUPPORTING FUNCTIONS
################################################################################
#
# This file contains the main cleanbatch() function and supporting utilities:
#   - calc_oob_evil_twins(): Identifies Evil Twins (Step 9)
#   - temporary_extraneous_infants(): Temporary SDE resolution (Step 5)
#   - clean_infants_with_sde(): Final SDE resolution (Step 13)
#   - Various EWMA calculation helpers
#
# Additional support functions are in pediatric_support.R:
#   - valid(): Identifies included/partially-included rows
#   - swap_parameters(): Finds DOP values for comparison
#   - temporary_extraneous(): Generic SDE resolution (used in Step 5)
#
################################################################################
# USAGE NOTES
################################################################################
#
# 1. PARALLEL PROCESSING: Set parallel=TRUE and num.batches for large datasets.
#    Each batch processes independently, so results are deterministic.
#
# 2. MEMORY: Large datasets may require increasing R memory limits.
#    Use data.table operations for efficiency.
#
# 3. DEBUGGING: Results from intermediate steps can be captured by inserting
#    saveRDS() or write.csv() calls at key points.
#
# 4. DATA.TABLE MODIFICATION: Many operations modify data.df by reference
#    (using := operator) for efficiency. Be aware of this when debugging.
#

#
#
################################################################################
# INFORMATION ABOUT IDENTIFIERS AND PRECISION HANDLING
################################################################################
#
# A substantial amount of effort is made to handle floating point precision
# and to order by row id. This serves two purposes
# 1) Ensures consistent handling of a dataset across multiple runs
# 2) Facilitates debugging of R and Stata by eliminating sources of variation
#
# Row IDENTIFIERS are discussed further above.
#
# FLOATING POINT PRECISION / ROUNDING SCHEME:
#   To ensure R and Stata implementations produce identical results, this
#   algorithm uses a consistent rounding strategy:
#
#   1. Z-SCORE ROUNDING: All z-scores are rounded to 3 decimal places
#      immediately after calculation using a custom round_stata() function
#      or janitor::round_half_up(). This cleans up initial calculation noise.
#
#   2. THRESHOLD DOUBLE-ROUNDING: When comparing values to decision thresholds,
#      we use "double-rounding" - round to 3 decimals, then to 2:
#
#      Example: janitor::round_half_up(janitor::round_half_up(variable, 3), 2) > threshold
#
#      Why? When subtracting two already-rounded values, floating-point
#      arithmetic can introduce noise at finer precision levels. The first
#      rounding (3 decimals = 0.001) cleans up subtraction noise, and the
#      second rounding (2 decimals = 0.01) standardizes the comparison point.
#
#   3. ROUNDING DIRECTION: R's default round() uses "round half to even"
#      (banker's rounding), but Stata rounds "half away from zero".
#      This implementation uses janitor::round_half_up() to match Stata's
#      behavior exactly.
#
#   CRITICAL: Do NOT use base R's round() function for any precision-sensitive
#   operations in this algorithm. Always use janitor::round_half_up() or the
#   custom round_stata() function.
#
#   IMPORTANT: Do NOT round raw measurements (weight, height, head circumference).
#   Rounding only applies to derived z-scores and threshold comparisons.
#


  #### R-Oxygen Markup (Hidden)


# Main growthcleanr algorithm and other exported functions
# Main growthcleanr algorithm (cleangrowth): adult and pediatric are in *_clean.R
#(main bulk of algorithm) and *_support.R (supporting functions)

#' Clean growth measurements
#'
#' @param subjid Vector of unique identifiers for each subject in the database.
#' @param id Vector of unique ids for each row in the database.
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
#' @param prelim_infants TRUE/FALSE. Run the in-development release of the infants algorithm (expands pediatric algorithm to improve performance for children 0 - 2 years). Not recommended for use in research. For more information regarding the logic of the algorithm, see the vignette 'Preliminary Infants Algorithm.' Defaults to FALSE.
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
                        id)
                        { # id parameter is required
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

  # Stata-style rounding (round half away from zero)
  # R's default round() uses banker's rounding (round half to even)
  # which can cause discrepancies at boundary cases
  round_stata <- function(x, digits = 3) {
    multiplier <- 10^digits
    sign(x) * floor(abs(x) * multiplier + 0.5) / multiplier
  }

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

  # Keep user's id as 'id' (matches Stata), create 'internal_id' for R processing
  # id is required in input data
  if (length(id) != length(measurement)) {
    stop("id must be provided and have the same length as measurement")
  }
  data.all.ages[, id := id]
  data.all.ages[, internal_id := seq_len(.N)]  # Sequential IDs for internal use

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
      "Exclude-EWMA1-Extreme",
      "Exclude-SDE-EWMA-All-Extreme"
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
      # Round z-scores to 0.001
      # Use Stata-style rounding (round half away from zero)
      data.all[, sd.orig_cdc := round_stata(sd.orig_cdc, 3)]
      data.all[, sd.orig_who := round_stata(sd.orig_who, 3)]

      # Changed smoothing from 2-4y to 2-5y to match Stata
      # smooth z-scores/SD scores between ages 2-5yo using weighted scores
      # older uses cdc, younger uses who
      data.all$ageyears <- data.all$agedays/365.25

      who_weight <- 5 - (data.all$ageyears)
      cdc_weight <- (data.all$ageyears) - 2

      smooth_val <- data.all$ageyears >= 2 &
        data.all$ageyears <= 5 &
        data.all$param != "HEADCM"
      data.all[smooth_val,
               sd.orig := (data.all$sd.orig_cdc[smooth_val]*cdc_weight[smooth_val] +
                             data.all$sd.orig_who[smooth_val]*who_weight[smooth_val])/3]
      # <- 
      # otherwise use WHO and CDC for older and younger, respectively
      who_val <- data.all$param == "HEADCM" |
        data.all$ageyears < 2
      
      data.all[(who_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_cdc)),
               sd.orig := data.all$sd.orig_who[(who_val & !smooth_val)  | (smooth_val & is.na(data.all$sd.orig_cdc))]]
      
      cdc_val <- data.all$param != "HEADCM" &
        data.all$ageyears > 5
      
      data.all[(cdc_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_who)),
               sd.orig := data.all$sd.orig_cdc[(cdc_val & !smooth_val) | (smooth_val & is.na(data.all$sd.orig_who))]]

      # Round z-scores to 0.001
      # Use Stata-style rounding
      data.all[, sd.orig := round_stata(sd.orig, 3)]

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
      
      # add age in months
      data.all[, agemonths := agedays / 30.4375]

      # Sort by id, age-dependent direction
      # Age-dependent ID selection for SDEs:
      # At age 0: prefer LOWEST id (earliest measurement, before fluid/interventions)
      # At age > 0: prefer HIGHEST id (consistent with other SDE handling)
      # Create sort key: for age 0 use id ascending, for age > 0 use id descending
      data.all[, id_sort := ifelse(agedays == 0, id, -id)]
      setorder(data.all, subjid, param, agedays, id_sort)

      # initialize the column so it's always present
      # Create row-level potcorr_wt flag
      # Equivalent to Stata line 293: gen potcorr_wt=sn_wt==1 & swtz<-2 & agemonths<10
      data.all[, potcorr := FALSE]
      data.all[, potcorr_wt := FALSE]

      # Mark infants needing Fenton correction
      # potcorr_wt: row-level flag for FIRST weight per subject (earliest agedays, not necessarily birth)
      #             that qualifies: Z<-2 AND age<10 months
      # After sort by agedays + id_sort, seq_len(.N)==1L selects first weight with age-dependent ID
      # Double-round threshold comparisons (3→2 decimals)
      data.all[param == "WEIGHTKG",
               potcorr_wt := (seq_len(.N) == 1L & janitor::round_half_up(janitor::round_half_up(sd.orig, 3), 2) < -2 & agemonths < 10),
               by = subjid]

      # Clean up sort key
      data.all[, id_sort := NULL]

      # propagate that flag to all rows for that subject (subject-level)
      # Equivalent to Stata line 294: bysort subjid: egen potcorr=max(potcorr_wt)
      data.all[, potcorr := any(potcorr_wt), by = subjid]
      
      # ensure it sticks around through merges
      # Include id for deterministic order
      setkey(data.all, subjid, param, agedays, id)
      
      # integer weight is in grams, rounded to the nearest 10
      data.all[potcorr == TRUE, intwt := trunc(v*100)*10]
      # data.all[potcorr == TRUE, intwt := round(v*10)]
      # write.csv(fentlms_foraga, "../data/fentlms_foraga.csv")
      # replace to facilitate merging with fenton curves
      data.all[intwt >= 250 & intwt <=560, intwt := 570]
      
      # merge with fenton curves
      data.all <- merge(
        data.all, fentlms_foraga, by = c("sex", "intwt"),
        all.x = TRUE)

      # Propagate fengadays ONLY from potcorr_wt rows
      # Equivalent to Stata lines 309-310: gen aga_i=fengadays if potcorr_wt==1; bysort subjid: egen aga=min(aga_i)
      # Only use fengadays from the qualifying first weight (potcorr_wt==TRUE), not from any row
      data.all[, fengadays_subj := ifelse(any(!is.na(fengadays) & potcorr_wt),
                                           min(fengadays[potcorr_wt], na.rm = TRUE),
                                           NA_real_),
               by = subjid]
      # Use subject-level fengadays for all calculations
      data.all[!is.na(fengadays_subj), fengadays := fengadays_subj]
      data.all[, fengadays_subj := NULL]  # Clean up temporary variable

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

      # Reset potcorr_wt when Fenton merge fails
      # Equivalent to Stata lines 306-308
      # If merge failed (fen_wt_m is NA), reset potcorr_wt for that row
      # Note: only potcorr_wt rows contribute to fengadays, so this prevents failed merges from being used
      # Debug statements removed
      data.all[potcorr_wt == TRUE & is.na(fen_wt_m), potcorr_wt := FALSE]
      # Recalculate subject-level potcorr: if ANY row for subject has potcorr_wt==TRUE, ALL rows get potcorr=TRUE
      # Equivalent to Stata lines 307-308: drop potcorr; bysort subjid: egen potcorr=max(potcorr_wt)
      data.all[, potcorr := any(potcorr_wt), by = subjid]
      
      # Fix Fenton z-score calculation for HT/HC
      # Convert measurement to Fenton units: WT needs kg→grams (*1000), HT/HC already in cm
      # Calculate Fenton z-score using appropriate LMS parameters for each measurement type
      data.all[, v_fenton := ifelse(param == "WEIGHTKG", v * 1000, v)]

      # Get the correct L, M, S parameters based on param type
      data.all[, fen_l := ifelse(param == "WEIGHTKG", fen_wt_l,
                                  ifelse(param == "HEIGHTCM", fen_ht_l, fen_hc_l))]
      data.all[, fen_m := ifelse(param == "WEIGHTKG", fen_wt_m,
                                  ifelse(param == "HEIGHTCM", fen_ht_m, fen_hc_m))]
      data.all[, fen_s := ifelse(param == "WEIGHTKG", fen_wt_s,
                                  ifelse(param == "HEIGHTCM", fen_ht_s, fen_hc_s))]

      # Calculate z-score using LMS formula
      data.all[, unmod_zscore := ((v_fenton/fen_m)^fen_l - 1)/(fen_l * fen_s)]

      # Clean up temporary variables
      data.all[, c("v_fenton", "fen_l", "fen_m", "fen_s") := NULL]
      
      
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
      
      # table(data.all$potcorr, is.na(data.all$sd.corr))

      # Fix WHO↔CDC smoothing to match Stata 2bH
      # Stata smooths ages 2-5 (not 2-4), divides by 3 (not 2), handles HC separately
      # Equivalent to Stata lines: gen cwhoweight=5-ageyears; gen ccdcweight=ageyears-2
      # foreach p in wt ht { replace c`p'z=(cwho`p'z * cwhoweight + ccdc`p'z * ccdcweight)/3 if ageyears>2 }

      # For WT and HT: smooth WHO↔CDC for ages 2-5
      data.all[
        (param %in% c("WEIGHTKG", "HEIGHTCM")) &
          (agedays/365.25 > 2) &
          (agedays/365.25 < 5) &
          !is.na(sd.c_who) &
          !is.na(sd.c_cdc),

        sd.c := (sd.c_who * (5 - agedays/365.25) +
                   sd.c_cdc * ((agedays/365.25) - 2)) / 3
      ]

      # For WT and HT: ages ≥5 use CDC only
      data.all[
        (param %in% c("WEIGHTKG", "HEIGHTCM")) &
          (agedays/365.25 >= 5) &
          !is.na(sd.c_cdc),
        sd.c := sd.c_cdc
      ]

      # For WT and HT: ages ≤2 use WHO only
      data.all[
        (param %in% c("WEIGHTKG", "HEIGHTCM")) &
          (agedays/365.25 <= 2) &
          !is.na(sd.c_who),
        sd.c := sd.c_who
      ]

      # For HC: use WHO for all ages >2 (no smoothing for HC)
      # Equivalent to Stata: replace chcz=cwhohcz if ageyears>2
      data.all[
        (param == "HEADCM") &
          (agedays/365.25 > 2) &
          !is.na(sd.c_who),
        sd.c := sd.c_who
      ]

      # For HC: ages ≤2 use WHO
      data.all[
        (param == "HEADCM") &
          (agedays/365.25 <= 2) &
          !is.na(sd.c_who),
        sd.c := sd.c_who
      ]

      # Prefer Fenton over corrected WHO for ages ≤2
      # Equivalent to Stata 2bG lines 375-376:
      #   gen c`p'z=fen`p'z if potcorr==1
      #   replace c`p'z=cwho`p'z if potcorr==1 & c`p'z==. & ageyears<=2
      # Stata logic: Use Fenton when available, fall back to corrected WHO only if Fenton missing

      # Save corrected WHO/CDC in case needed as fallback
      data.all[, sd.c_temp := sd.c]

      # For ages ≤2 and potcorr: prefer Fenton
      data.all[potcorr == TRUE & agedays/365.25 <= 2 & !is.na(unmod_zscore),
               sd.c := unmod_zscore]

      # For ages ≤2 and potcorr: use corrected WHO only if Fenton missing
      data.all[potcorr == TRUE & agedays/365.25 <= 2 & is.na(sd.c) & !is.na(sd.c_temp),
               sd.c := sd.c_temp]

      # Clean up temp column
      data.all[, sd.c_temp := NULL]

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
      # uncorrweight <- 4 - (data.all$agedays/365.25)
      # corrweight <- (data.all$agedays/365.25) - 2
      # smooth_val <- data.all$agedays/365.25 >= 2 &
      #   data.all$agedays/365.25 <= 4
      # 
      # data.all[smooth_val & !is.na(sd.c) & !is.na(sd.orig),
      #          sd.corr := (sd.orig[smooth_val]*uncorrweight[smooth_val] +
      #                        sd.c[smooth_val]*corrweight[smooth_val]) / 2]
      
      data.all[
        (agedays/365.25 >= 2) &
          (agedays/365.25 <= 4) &
          !is.na(sd.c) &
          !is.na(sd.orig),
        
        sd.corr :=
          (sd.orig * (4 - agedays/365.25) +
             sd.c     * ((agedays/365.25) - 2)) / 2
      ]
      
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

      # Filter SDE-Identicals BEFORE uncorr calculation
      # Matches Stata's Early Step 13 (lines 172-180) which removes identicals before section 2bK
      # At birth (age 0): Keep LOWEST id (earliest, before fluid/interventions)
      # At age > 0: Keep HIGHEST id (consistent with other SDE handling)

      # Identify same-day identical values
      tmp[, `:=`(
        n_on_day = .N,
        n_unique_vals = uniqueN(v)
      ), by = .(subjid, agedays)]

      tmp[, all_identical := n_on_day > 1 & n_unique_vals == 1]

      # Age-dependent ID selection for identicals
      tmp[, keep_id := ifelse(agedays == 0,
                              min(id, na.rm = TRUE),
                              max(id, na.rm = TRUE)),
          by = .(subjid, agedays)]

      # Filter out non-selected identicals
      tmp <- tmp[!(all_identical & id != keep_id)]

      # Clean up temporary columns
      tmp[, `:=`(n_on_day = NULL, n_unique_vals = NULL, all_identical = NULL, keep_id = NULL)]

      # Now create seq_win on the filtered data (without SDE-Identicals)
      # Age-dependent ID sorting for consistent sequence numbering:
      # At age 0: sort by id ascending (lowest first)
      # At age > 0: sort by id descending (highest first)
      tmp[, id_sort := ifelse(agedays == 0, id, -id)]
      tmp <- tmp[order(subjid, agedays, id_sort),]
      tmp[, id_sort := NULL]

      # Create seq_win (matches Stata's sn_wt creation at line 303, after identicals removed)
      tmp[, seq_win := sequence(.N), by = subjid]
      # we're only looking at the first 4 values, and they need to be < 2 years
      tmp <- tmp[seq_win <= 4 & (agedays/365.25) < 2,]
      # don't look at subjects where there is only 1 score and there is no
      # value for either
      tmp <- tmp[!(is.na(sd.corr) & is.na(sd.orig)),]
      tmp <- tmp[subjid %in% names(table(subjid) > 1),]
      
      # abssumdiff was actually correct
      # Stata line 429: abs(dswtz2 + dswtz3 + dswtz4) = abs(sum(diffs))
      # R original: abs(sum(sd.corr[1] - sd.corr)) = abs(sum(diffs)) - CORRECT!
      # Previous "fix" changed to sum(abs(diffs)) which was WRONG
      # create differences, absolute sum them
      tmp[, sd.corr_abssumdiff := abs(sum(sd.corr[1] - sd.corr)), by = subjid]
      tmp[, sd.orig_abssumdiff := abs(sum(sd.orig[1] - sd.orig)), by = subjid]
      # find subjects where corrected value needs to be replaced

      # Only evaluate uncorr for first weight (potcorr_wt==TRUE)
      # Equivalent to Stata: gen uncorr_i=abssumdiff_sc>abssumdiff_s & abssumdiff_s!=. & potcorr_wt==1
      # Filter to first measurement only (equivalent to potcorr_wt==1, earliest agedays per subject)
      tmp[, is_first := seq_win == 1]

      sub_replace <- unique(tmp[sd.corr_abssumdiff > sd.orig_abssumdiff &
                                !is.na(sd.orig_abssumdiff) &
                                is_first == TRUE,
                                subjid])

      # Debug output for uncorr
      if (!quietly) {
        cat(sprintf("DEBUG uncorr: %d subjects flagged for replacement\n", length(sub_replace)))
        if (length(sub_replace) > 0 && length(sub_replace) <= 10) {
          cat("Subjects:", paste(sub_replace, collapse=", "), "\n")
        }
      }

      # replace accordingly in the main dataframe
      data.all[subjid %in% sub_replace, sd.corr := sd.orig]

      # Round z-scores to 0.001
      # Use Stata-style rounding
      data.all[, sd.corr := round_stata(sd.corr, 3)]

      # Add uncorr flag for comparison with Stata
      # uncorr=1 means correction was reverted (subject in sub_replace list)
      data.all[, uncorr := as.integer(subjid %in% sub_replace)]

      # Add sd.corr, potcorr, uncorr to orig_colnames
      # so they don't get dropped during column filtering at line 964
      # Added uncorr to list
      orig_colnames <- c(orig_colnames, "sd.corr", "potcorr", "uncorr")

      # Mid-correction debug save (before columns dropped)
#      if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 2) {
#        saveRDS(data.all[, .(index, subjid, param, agedays, ageyears, sex, v,
#                             potcorr, pmagedays, cagedays,
#                             sd.orig, sd.corr, sd.c, sd.c_who, sd.c_cdc,
#                             unmod_zscore, sd.orig_uncorr)],
#                paste0("R-step2-mid-correction-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#        message("DEBUG: Saved R-step2-mid-correction.rds (BEFORE column filtering)")
#      }

      # remove many added columns
      # write.csv(data.all, "../data/midpoint_check_corrected_zscores.csv")
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
    # Include id in setkey for deterministic SDE order
    # Without id, SDEs (same subjid/param/agedays) have undefined order, causing different
    # index values in sequential vs parallel mode
    setkey(data.all, subjid, param, agedays, id)

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
    
    # Include id for deterministic order
    setkey(data.all, subjid, param, agedays, id)
    data.all[, tbc.sd := sd.orig - sd.median]
    if (prelim_infants){
      # separate out corrected and noncorrected values
      data.all[, ctbc.sd := sd.corr - sd.median]
    }

    # Round z-scores to 0.001
    # Use Stata-style rounding
    data.all[, tbc.sd := round_stata(tbc.sd, 3)]
    if (prelim_infants) {
      data.all[, ctbc.sd := round_stata(ctbc.sd, 3)]
    }

    # Debug save point for Step 2/3 (after z-score calc and recentering)
#    if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 2) {
      # Save key intermediate variables for debugging correction logic
      # Create final_tbc before saving
      # For potcorr subjects, use ctbc.sd (corrected); for others, use tbc.sd (uncorrected)
#      if ("potcorr" %in% colnames(data.all)) {
#       data.all[, final_tbc := ifelse(!is.na(potcorr) & potcorr == TRUE & !is.na(ctbc.sd), ctbc.sd, tbc.sd)]
#        if (!quietly) {
#          potcorr_count <- sum(data.all[, !is.na(potcorr) & potcorr == TRUE & !is.na(ctbc.sd)])
#          cat(sprintf("DEBUG Issue 16: %d rows with potcorr==TRUE using ctbc.sd\n", potcorr_count))
#        }
#      } else {
#        data.all[, final_tbc := tbc.sd]
#      }

      # Added uncorr to debug output
      # Added id to debug output for comparison with Stata
      # Commented out - was running unconditionally after debug_step removal
#      debug_cols <- intersect(colnames(data.all),
#                              c("id", "index", "subjid", "param", "agedays", "sex", "v", "ageyears",
#                                "sd.orig", "sd.corr", "tbc.sd", "ctbc.sd", "final_tbc",
#                                "sd.orig_cdc", "sd.orig_who", "sd.c_cdc", "sd.c_who", "sd.c",
#                                "pmagedays", "cagedays", "potcorr", "uncorr", "unmod_zscore"))
#      saveRDS(data.all[, ..debug_cols],
#              paste0("R-step2-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#      message("DEBUG: Saved R-step2-output.rds (after z-score calc and recentering) - STOPPING")
#
#      return(data.all[, .(line, exclude, tbc.sd = final_tbc, param)])
#    }

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

    # Set HC >= 5 years to Missing
    # WHO reference for HC only goes up to 5 years
    data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Missing']

    if (prelim_infants){
      # Removed NNTE calculation - was not helpful for efficiency
      # NNTE tried to predict "won't need EWMA" based on trajectory smoothness, but:
      # 1. Prediction-based filtering is less effective than deterministic filtering
      # 2. Better approach: pre-filter based on what's PHYSICALLY POSSIBLE (no SDEs, no duplicate values)
      # 3. Those deterministic filters are now in CF Step 6 and SDE Step 13
      # Set nnte=FALSE for all rows so existing filters become no-ops
      data.all[, nnte := FALSE]
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
    all_results <- merge(all_results, checkpoint_data %>% select(-c(subjid, param, agedays)), by = "id", all.x = TRUE)
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
    # Changed WHO cutoff from 3y to 5y for smoothing
    # WHO data available up to 5 years for height/weight, needed for 2-5y smoothing
    src <- ifelse(agedays < 5*365.25 & !cdc.only, 'WHO', 'CDC')
    
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
ewma <- function(agedays, z, ewma.exp, ewma.adjacent = TRUE, window = 25) {
  # Added window parameter to limit EWMA to max window values on each side
  # Changed default to 25 for better accuracy with minimal efficiency loss
  # Set window = Inf to disable windowing
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
    # Apply exponent ROW-wise, not column-wise
    # Each observation's EWMA should use that observation's exponent for all weights
    # sweep(..., 1, ...) applies the vector to each ROW (margin 1)
    # Old (buggy): delta <- ifelse(delta == 0, 0, (delta + 5) ^ ewma.exp)  # column-wise recycling
    delta <- ifelse(delta == 0, 0, sweep(delta + 5, 1, ewma.exp, FUN = "^"))

    # Apply windowing: zero out weights beyond window positions
    # Matches Stata: abs(sn_`ep' - `n') <= `window'
    if (!is.null(window) && is.finite(window)) {
      pos_diff <- abs(outer(1:n, 1:n, "-"))  # position difference matrix
      delta[pos_diff > window] <- 0
    }

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
# Revised temp SDE logic to match Stata ewmacode-2025-12-11.do
# Key changes:
#   1. SP median calculated from ALL included values (not just non-SDE days)
#   2. DOP median distance used as secondary tiebreaker (always, not just fallback)
#   3. Sort order: absdmedian_spz -> absdmedian_dopz -> index (spz/dopz ascending, index descending)
# Note: Identical values are removed earlier (Step 13b SDE-Identicals), so by the time
# tempsde runs, all same-day values are dissimilar.

temporary_extraneous_infants <- function(df, exclude_from_dop_ids = NULL) {
  # Added exclude_from_dop_ids parameter
  # In Step 13, temp SDEs should be excluded from DOP median calculation (but included in SP median)
  # Pass the ids of temp SDEs to exclude them from DOP median only

  # avoid "no visible binding" warnings
  agedays <- absdmedian.spz <- absdmedian.dopz <- extraneous <- NULL
  extraneous.this.day <- index <- median.spz <- median.dopz <- NULL
  param <- subjid <- tbc.sd <- NULL

  # Make copy before modifying to avoid data.table alloccol error
  df <- copy(df)

  # add subjid and param if needed (may be missing depending on where this is called from)
  if (is.null(df$subjid))
    df[, subjid := NA]
  if (is.null(df$param))
    df[, param := NA]

  # Fix keyby reordering bug
  # The old code computed valid.rows BEFORE keyby, but keyby reorders rows.
  # Also, the returned boolean vector was in keyby order, not original order.
  # Fix: Store original order, compute valid.rows AFTER keyby, return in original order

  # Use id not index for deterministic SDE order
  # index is assigned based on input order after sorting, which is not deterministic for SDEs
  # id is user-provided and deterministic

  # Store original row positions for mapping result back to original order
  df[, orig_row := .I]
  original_rows <- df$orig_row

  # make a small copy of df with fields we need
  # Include id for age-dependent tiebreaker
  # Use id in keyby for deterministic SDE order
  df <- df[j = .(tbc.sd, exclude, id, orig_row), keyby = .(subjid, param, agedays, id)]

  # Now compute valid.rows on the keyby-sorted data
  # Removed nnte filter (nnte calculation removed)
  valid.rows <- valid(df, include.temporary.extraneous = TRUE)

  # initialize fields
  df[, `:=`(
    median.spz = as.double(NA),
    median.dopz = as.double(NA),
    absdmedian.spz = as.double(NA),
    absdmedian.dopz = as.double(NA),
    extraneous.this.day = FALSE,
    extraneous = FALSE
  )]

  # determine on which days there are extraneous measurements (more than 1 valid measurement on that age day)
  df[valid.rows, extraneous.this.day := (.N > 1), by = .(subjid, param, agedays)]

  # DOP (designated other parameter) mapping
  dop_map <- c(
    "HEIGHTCM" = "WEIGHTKG",
    "WEIGHTKG" = "HEIGHTCM",
    "HEADCM" = "HEIGHTCM"
  )

  # ----- STEP 1: Calculate SP medians for each parameter -----
  # Stata: bysort subjid_`sdep': egen median_spz_`sdep'_i = median(tbc`sdep'z)
  # Calculate median for each subject-param (all valid values, not just non-SDE days)
  df[valid.rows, median.spz := median(tbc.sd, na.rm = TRUE), by = .(subjid, param)]

  # DEBUG 2025-12-18: Check Subject 725 median calculation - REMOVED (medians match)

  # ----- STEP 2: Distribute SP medians across all rows for each subject -----
  # This enables cross-parameter lookup (DOP medians)
  # Stata uses two-step: first by subjid_p, then max() by subjid to distribute
  # In R, we create a lookup table and merge

  # Create lookup table of SP medians by subject
  sp_medians <- df[valid.rows & !is.na(median.spz),
                   .(median.spz = median.spz[1]),
                   by = .(subjid, param)]

  # ----- STEP 3: Calculate DOP medians -----
  # DOP median logic differs between Step 5 and Step 13
  # Step 5: Use ALL measurements (no temp SDEs exist yet)
  # Step 13: Exclude temp SDEs from DOP median (passed via exclude_from_dop_ids)

  # For each parameter, look up the median of its designated other parameter
  df[, median.dopz := as.double(NA)]

  if (is.null(exclude_from_dop_ids)) {
    # Step 5: Use sp_medians calculated from all values
    for (p in c("WEIGHTKG", "HEIGHTCM", "HEADCM")) {
      dop <- dop_map[p]
      dop_medians <- sp_medians[param == dop, .(subjid, dop_median = median.spz)]
      dop_lookup <- setNames(dop_medians$dop_median, dop_medians$subjid)
      df[param == p, median.dopz := dop_lookup[as.character(subjid)]]
    }
  } else {
    # Step 13: Calculate DOP median excluding temp SDEs (matches Stata line 1516)
    # Stata: bysort subjid_`o': egen median_`o'z_i=median(tbc`o'z) if exc_`o'==0
    dop_valid_rows <- valid.rows & !(df$id %in% exclude_from_dop_ids)

    # Calculate DOP-specific medians excluding temp SDEs
    df[, median.dopz.calc := as.double(NA)]
    df[dop_valid_rows, median.dopz.calc := median(tbc.sd, na.rm = TRUE), by = .(subjid, param)]

    # Create lookup table for DOP medians (excludes temp SDEs)
    dop_sp_medians <- df[dop_valid_rows & !is.na(median.dopz.calc),
                         .(median.spz = median.dopz.calc[1]),
                         by = .(subjid, param)]

    for (p in c("WEIGHTKG", "HEIGHTCM", "HEADCM")) {
      dop <- dop_map[p]
      dop_medians <- dop_sp_medians[param == dop, .(subjid, dop_median = median.spz)]
      dop_lookup <- setNames(dop_medians$dop_median, dop_medians$subjid)
      df[param == p, median.dopz := dop_lookup[as.character(subjid)]]
    }

    # Clean up temporary column
    df[, median.dopz.calc := NULL]
  }

  # ----- STEP 4: Calculate absolute distances from medians -----
  # Only for SDE days
  df[valid.rows & extraneous.this.day, `:=`(
    absdmedian.spz = abs(tbc.sd - median.spz),
    absdmedian.dopz = abs(tbc.sd - median.dopz)
  )]

  # Handle cases where SP median is NA (all days are SDE days for this param)
  # Fall back to using DOP median as the primary comparison
  df[valid.rows & extraneous.this.day & is.na(absdmedian.spz) & !is.na(median.dopz),
     absdmedian.spz := abs(tbc.sd - median.dopz)]

  # Final fallback: if both medians unavailable, use abs(tbc.sd) (treat 0 as median)
  df[valid.rows & extraneous.this.day & is.na(absdmedian.spz),
     absdmedian.spz := abs(tbc.sd)]

  # Set NA DOP distances to Inf so they sort last
  df[valid.rows & extraneous.this.day & is.na(absdmedian.dopz),
     absdmedian.dopz := Inf]

  # ----- STEP 5: Select which value to keep -----
  # Sort by: absdmedian.spz (asc), absdmedian.dopz (asc), index (desc = keep highest)
  # Keep first after sort, mark others as extraneous

  df[valid.rows & extraneous.this.day, extraneous := {
    # Sort by absdmedian.spz (asc), then absdmedian.dopz (asc), then id tiebreaker
    # Keep the first after sort, mark others as extraneous
    # Round to 3 decimals to avoid floating point precision issues
    # Age-dependent id tiebreaker to match Stata
    # Stata uses obsid (observation ID) for tiebreaker, not row index
    # At agedays=0: pick lowest id (sort ascending)
    # At agedays>0: pick highest id (sort descending via -id)
    # Use id for selection, not index
    tiebreaker <- if(agedays[1] == 0) id else -id
    # Restored rounding with round_half_up to match Stata
    # Double-round to handle floating point noise from z-score subtraction
    ord <- order(janitor::round_half_up(janitor::round_half_up(absdmedian.spz, 3), 2),
                 janitor::round_half_up(janitor::round_half_up(absdmedian.dopz, 3), 2), tiebreaker)
    keep_id <- id[ord[1]]
    id != keep_id
  }, by = .(subjid, param, agedays)]

  # Debug save point for mid-Step 5
#  if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 5) {
#    saveRDS(df, paste0("R-step5-mid-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step5-mid.rds (with median vars)")
#    stop("DEBUG: Stopping at mid-Step 5 as requested")
#  }

  # Map result back to original row order
  # The df is in keyby order, but caller expects result in original order.
  # Use orig_row instead of index for deterministic ordering
  extraneous_result <- df$extraneous & valid.rows

  # Create a named vector: orig_row -> is_extraneous
  extraneous_lookup <- setNames(extraneous_result, df$orig_row)

  # Return in original order by looking up each original row position
  return(as.logical(extraneous_lookup[as.character(original_rows)]))
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

  # Handle minimal datasets with < 2 rows
  # When there are 0 or 1 rows, there are no adjacent pairs to compare for evil twins
  if (nrow(df) < 2) {
    df[, "oob" := FALSE]
    return(df)
  }

  # for differences for adjacent, we pad with a high number to default to true
  # we don't want to consider the comparison with the next subj/param
  # Fix logic - SAME pair must exceed both thresholds
  # Old logic allowed mixing: tbc from one neighbor, ctbc from another
  # Correct: (tbc_next > 5 AND ctbc_next > 5) OR (tbc_prev > 5 AND ctbc_prev > 5)
  same_sp_next <- rev(duplicated(rev(paste(df$subjid, "_", df$param))))
  same_sp_prev <- duplicated(paste(df$subjid, "_", df$param))

  tbc_next_diff <- abs(c(df$tbc.sd[2:nrow(df)], Inf) - df$tbc.sd)
  tbc_prev_diff <- abs(c(Inf, df$tbc.sd[1:(nrow(df)-1)]) - df$tbc.sd)
  ctbc_next_diff <- abs(c(df$ctbc.sd[2:nrow(df)], Inf) - df$ctbc.sd)
  ctbc_prev_diff <- abs(c(Inf, df$ctbc.sd[1:(nrow(df)-1)]) - df$ctbc.sd)

  # Double-round threshold comparisons (3→2 decimals)
  # Handles floating-point noise from subtraction of rounded values
  tbc_next_diff_rounded <- janitor::round_half_up(janitor::round_half_up(tbc_next_diff, 3), 2)
  tbc_prev_diff_rounded <- janitor::round_half_up(janitor::round_half_up(tbc_prev_diff, 3), 2)
  ctbc_next_diff_rounded <- janitor::round_half_up(janitor::round_half_up(ctbc_next_diff, 3), 2)
  ctbc_prev_diff_rounded <- janitor::round_half_up(janitor::round_half_up(ctbc_prev_diff, 3), 2)

  oob <- (tbc_next_diff_rounded > 5 & ctbc_next_diff_rounded > 5 & same_sp_next) |
         (tbc_prev_diff_rounded > 5 & ctbc_prev_diff_rounded > 5 & same_sp_prev)

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
  
  # Fix smoothing formula
  # Match base z-score smoothing: ages 2-5, who_weight = 5-age, divide by 3
  who_weight <- 5 - (df$agedays/365.25)
  cdc_weight <- (df$agedays/365.25) - 2

  smooth_val <- df$agedays/365.25 >= 2 & df$agedays/365.25 <= 5 & df$param != "HEADCM"

  # Fix data.table double-indexing bug
  # When using df[smooth_val, ...], columns cn.orig_cdc/cn.orig_who are already subset
  # Do NOT use [smooth_val] on column references inside the assignment
  # Only use [smooth_val] on external vectors (cdc_weight, who_weight)
  df[smooth_val,
     cn.orig := (cn.orig_cdc * cdc_weight[smooth_val] +
                   cn.orig_who * who_weight[smooth_val])/3]
  
  # otherwise use WHO and CDC for older and younger, respectively
  ### CP REPLACE D
  # who_val <- df$param == "HEADCM" |
  #   df$agedays/365.25 < 2
  # df[who_val | (smooth_val & is.na(df$cn.orig_cdc)),
  #    cn.orig := df$cn.orig_who[who_val  | (smooth_val & is.na(df$cn.orig_cdc))]]
  
  who_val <- df$param == "HEADCM" | df$agedays/365.25 < 2
  df[who_val, cn.orig := df$cn.orig_who[who_val]]
  
  # Fix CDC z-score assignment for ages >= 4
  # Match Stata behavior: use CDC z-scores for WT/HT at ages >= 4 years
  # Changed | to & and cn.orig_who to cn.orig_cdc
  cdc_val <- df$param != "HEADCM" & df$agedays/365.25 >= 4
  df[cdc_val, cn.orig := df$cn.orig_cdc[cdc_val]]
  
  
  # df[cdc_val | (smooth_val & is.na(df$cn.orig_who)), cn.orig := df$cn.orig_cdc[cdc_val | (smooth_val & is.na(df$sd.orig_who))]]
  
  
  ### CP REPLACE U
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
    originator_row <- cf_string_id <- ageday_has_include <-
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
  
  # Use id instead of index for deterministic SDE order
  # index depends on input order, which can vary between batches
  # id is user-provided and deterministic
  data.df <- data.table(data.df, key = c('subjid', 'param', 'agedays', 'id'))
  # Recreate index for batch processing
  # index was created on full dataset before batching (line 1022), so batches have
  # non-contiguous indices. This breaks merge/subset operations that use index.
  data.df[, index := 1:.N]

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

  # EARLY STEP 13: SDE-Identicals (before Steps 5/6) ----
  # Updated to use age-dependent id preference
  # Same-day identical values should be excluded before CF detection to match Stata
  # Stata removes identicals early (lines 179-215 in do file) before Steps 5/6
  # At birth (age 0): Keep LOWEST id (earliest, before fluid/interventions)
  # At age > 0: Keep HIGHEST id (consistent with other SDE handling)

  if (!quietly)
    cat(sprintf(
      "[%s] Exclude same-day identical measurements (early SDE-Identical)...\n",
      Sys.time()
    ))

  # SDE-Identical handles partial identicals - duplicates are marked even when
  # mixed with non-duplicates. Groups by subjid/param/agedays/v to find
  # duplicates of each specific value.
  data.df <- data.df %>%
    arrange(subjid, param, agedays, id) %>%
    # Group by subjid-param-agedays-measurement to find identical values
    group_by(subjid, param, agedays, v) %>%
    mutate(
      n_same_value = sum(exclude == "Include"),
      # Has duplicates if more than 1 included value with same measurement
      has_dup = (n_same_value > 1),
      # Age-dependent: keep lowest id at birth, highest id otherwise
      keep_id = ifelse(agedays == 0,
                       min(id[exclude == "Include"], na.rm = TRUE),
                       max(id[exclude == "Include"], na.rm = TRUE)),
      exclude = case_when(
        has_dup & exclude == "Include" & id != keep_id ~ "Exclude-SDE-Identical",
        TRUE ~ exclude
      )
    ) %>%
    ungroup() %>%
    select(-n_same_value, -has_dup, -keep_id)

  data.df <- data.table(data.df)
  # Include id for deterministic SDE order
  setkey(data.df, subjid, param, agedays, id)

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

  # Debug check after Step 5 - stop if duplicates found
  step5_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(step5_dupe_check) > 0) {
    message("DEBUG Step5: Found ", nrow(step5_dupe_check), " Include duplicates after Step 5 temp SDE")
    fwrite(step5_dupe_check, "debug-step5-duplicates.csv")
    # Save full rows for duplicates
    dupe_rows <- data.df[exclude == "Include"][step5_dupe_check, on = .(subjid, param, agedays)]
    fwrite(dupe_rows, "debug-step5-duplicate-rows.csv")
    # Save ALL data at this point for full diagnostics
    fwrite(data.df, "debug-step5-full-data.csv")
    stop("DEBUG: Stopping after Step 5 - duplicate Includes found. See debug-step5-*.csv files")
  }

  # 6:  carried forwards ----

  # Initialize sde_identical_rows before CF block so it exists when referenced later
  sde_identical_rows <- data.df[0]  # Empty data.table with same structure

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
    # CP TOGGLED TO FALSE but then back to true because nothing changed?
    # Removed nnte filter (nnte calculation removed)
    valid_set <- valid(data.df, include.temporary.extraneous = TRUE) &
      not_single
    data.sub <- data.df[valid_set,]

    # Pre-filter to subject-params with potential CFs
    # CF requires duplicate values - if all values unique, skip CF processing entirely
    # This excludes SDE-Identicals from consideration (already handled in early Step 13)
    data.sub[, sp_key := paste0(subjid, "_", param)]
    sp_has_dups <- data.sub[, .(has_dup_vals = uniqueN(v.orig) < .N), by = sp_key]
    sp_with_potential_cf <- sp_has_dups[has_dup_vals == TRUE, sp_key]
    n_total_sp <- uniqueN(data.sub$sp_key)
    n_with_cf <- length(sp_with_potential_cf)
    if (!quietly)
      cat(sprintf("  CF pre-filter: %d/%d subject-params have potential CFs (%.1f%%)\n",
                  n_with_cf, n_total_sp, 100*n_with_cf/n_total_sp))

    # Initialize cf=FALSE for all, only process those with potential CFs
    data.sub[, cf := FALSE]

    # Only process subject-params with duplicate values
    if (length(sp_with_potential_cf) > 0) {
      cf_subset <- data.sub[sp_key %in% sp_with_potential_cf]
      # Add id for consistent SDE order
      setorder(cf_subset, subjid, param, agedays, id)

      # CF logic updated to match Stata
      # Replaced dplyr/map_lgl with data.table for CF detection
      # Key logic fix (2025-12-11): Only compare if prior ageday has exactly ONE value
      #   - If prior ageday has multiple values (SDEs), skip CF comparison
      # Performance: data.table vectorized ops vs row-by-row map_lgl

      # Step 1: Count values per (subjid, param, ageday)
      day_counts <- cf_subset[, .(n_on_day = .N), by = .(subjid, param, agedays)]

      # Step 2: Get unique agedays per subject-param, then find prior ageday for each
      unique_days <- unique(cf_subset[, .(subjid, param, agedays)])
      setorder(unique_days, subjid, param, agedays)
      unique_days[, prior_ageday := shift(agedays, type = "lag"), by = .(subjid, param)]

      # Step 3: Get the value from each unique (subjid, param, ageday) where count == 1
      # For days with >1 value, we set prior_val to NA so CF check fails
      single_val_days <- cf_subset[, .(
        single_val = if (.N == 1) v.orig[1] else NA_real_
      ), by = .(subjid, param, agedays)]

      # Step 4: Merge prior_ageday into cf_subset
      cf_subset[unique_days, prior_ageday := i.prior_ageday, on = .(subjid, param, agedays)]

      # Step 5: Look up the single value from prior ageday
      cf_subset[single_val_days, prior_single_val := i.single_val,
               on = .(subjid, param, prior_ageday = agedays)]

      # Step 6: CF = TRUE if prior day had exactly 1 value AND current matches it
      # Use exact equality to match Stata (no tolerance)
      cf_subset[, cf := !is.na(prior_single_val) & v.orig == prior_single_val]

      # Merge cf results back to data.sub by index
      data.sub[cf_subset, cf := i.cf, on = "index"]
    }
    # End of pre-filter if block - subject-params without potential CFs already have cf=FALSE

    # Cleanup temp columns
    data.sub[, sp_key := NULL]

    #### CP Untested, replicate of above code in data.table() if that sequence of mods works. This will be way faster than dplyr.
    
    # data.temp <- as.data.frame(data.sub) %>%
    #   arrange(subjid, param, agedays) %>%
    #   group_by(subjid, param) %>%
    #   mutate(
    #     cf = map_lgl(
    #       seq_along(agedays),
    #       function(i) {
    #         # first observation can't be carried forward
    #         if (i == 1) return(FALSE)
    # 
    #         # find the immediately previous unique ageday for this subject/param
    #         prev_ages <- unique(agedays[agedays < agedays[i]])
    #         if (length(prev_ages) == 0) return(FALSE)
    #         prev_age <- max(prev_ages)
    # 
    #         # get all values from that previous ageday
    #         prev_vals <- v.orig[agedays == prev_age]
    # 
    #         # if current value matches any of those (within tolerance), flag as CF
    #         any(abs(v.orig[i] - prev_vals) < 1e-6)
    #       }
    #     )
    #   ) %>%
    #   ungroup()
    
    
    
    
    
    
    
    # # for values with sdes, compare to all values from the prior day
    # data.sub[no_sde == FALSE, cf := (function(df){
    #   ages <- unique(agedays)
    #   if (length(ages) > 1){
    #     for (i in 2:length(ages)) {
    #       # find the set of measurements from the previous age in days
    #       all.prev.v <- df[agedays == ages[i - 1], v.orig]
    #       # if a measurement for the current age is found in the set of measurements from the previous age, then mark it as carried forward
    #       df[agedays == ages[i] &
    #            v.orig %in% all.prev.v, cf := TRUE]
    #     }
    #   }
    # })(copy(.SD)), .SDcols = c('agedays', 'cf', 'v.orig'), by = .(subjid, param)]
    # 
    # merge in the carried forwards
    cf_idx <- data.sub$index[data.sub$cf]
    data.df[index %in% cf_idx, exclude := "Exclude-Carried-Forward"]

    # Check if CFs exist before rescue processing
    # This optimization skips rescue logic if no CFs are present (matches Stata lines 775-780)
    any_cf <- any(data.df$exclude == "Exclude-Carried-Forward")
    if (!quietly)
      cat(sprintf("  CF rescue pre-filter: CFs exist = %s\n", any_cf))

    # Only process CF rescue if CFs exist
    if (any_cf) {

    # Redo temporary SDEs after CF identification
    # Temp SDEs should be re-evaluated now that CFs are excluded
    # This matches Stata lines 773-780

    # Check if any SDEs exist (matches Stata's anysde_chunk check)
    any_sde <- any(data.df$exclude == "Exclude-Temporary-Extraneous-Same-Day")

    if (any_sde) {
      # Temporarily convert Temp SDEs to potential includes for re-evaluation
      # (Stata does this by resetting exc==2 to exc==0 and restoring subjid)
      data.df[exclude == "Exclude-Temporary-Extraneous-Same-Day", exclude := "Include"]

      # Re-run temp SDE logic (now CFs are excluded, so SDE evaluation will differ)
      data.df$exclude[temporary_extraneous_infants(data.df)] <- 'Exclude-Temporary-Extraneous-Same-Day'
    }

    ### CP DROP IN
    # find out if values are whole or half imperial
    # data.df[param == "WEIGHTKG",
    #         wholehalfimp := (v.orig * 2.20462262)%%1 < 0.01]
    # data.df[param == "HEIGHTCM",
    #         wholehalfimp := (v.orig * 1.27)%%1 < 0.01]
    # data.df[param == "HEIGHTCM",
    #         wholehalfimp := (v.orig * 1.27 / 2)%%1 < 0.01]
    # data.df[param == "HEADCM",
    #         wholehalfimp := (v.orig * 1.27)%%1 < 0.01]
    # data.df[param == "HEADCM",
    #         wholehalfimp := (v.orig * 1.27 / 2)%%1 < 0.01]
    #
    # determine if measurements are in whole or half imperial units
    
    # wholehalfimp is ROW-LEVEL (not subject-level)
    # Stata calculates wholehalfimp per-row based on that row's own parameter:
    # - WEIGHTKG rows: TRUE if measurement is whole pounds
    # - HEIGHTCM rows: TRUE if measurement is whole or half inches
    # - HEADCM rows: TRUE if measurement is whole or half inches
    # The flag does NOT propagate across parameters or across rows within a subject

    # Initialize wholehalfimp to FALSE for all rows
    data.df[, wholehalfimp := FALSE]

    # WEIGHTKG: check if whole pounds (measurement * 2.20462262 is near integer)
    data.df[param == "WEIGHTKG",
            wholehalfimp := abs((v.orig * 2.20462262) %% 1) < 0.01]

    # HEIGHTCM: check if whole or half inches (measurement / 2.54 is near integer or .5)
    data.df[param == "HEIGHTCM",
            wholehalfimp := abs((v.orig / 2.54) %% 1) < 0.01 |
                            abs((v.orig / 2.54) %% 0.5) < 0.01]

    # HEADCM: check if whole or half inches (measurement / 2.54 is near integer or .5)
    data.df[param == "HEADCM",
            wholehalfimp := abs((v.orig / 2.54) %% 1) < 0.01 |
                            abs((v.orig / 2.54) %% 0.5) < 0.01]
    
    ### CP UP

    # Fix: Exclude SDE-Identicals from entire CF string calculation
    # SDE-Identicals break multiple parts of the CF logic (rle, originator detection, cs assignment)
    # Temporarily remove them, do CF calculations, then add them back (like Stata's subjidresc approach)

    # Save SDE-Identical rows and remove from data.df temporarily
    sde_identical_rows <- data.df[exclude == "Exclude-SDE-Identical"]
    data.df <- data.df[exclude != "Exclude-SDE-Identical"]

    # Complete rewrite of CF string detection
    # Previous approach used v.orig sort + position indexing which failed when
    # multiple CF strings had same measurement value (see Subject 746 bug)
    # New approach: directly match each CF to its most recent prior Include with same value
    # Calculate ageday_has_include BEFORE string detection
    # CFs on days with Includes are not eligible for rescue and should be excluded from
    # string detection (like Stata's subjidresc approach)

    # First, identify CFs on days with Includes - these are NOT eligible for rescue
    # ageday_has_include only checks for Includes
    # Removed Temp SDE check - only Includes should block CF rescue eligibility
    # This matches Stata line 817: egen ageday_include = max(exc==0)
    data.df[, ageday_has_include := any(as.character(exclude) == "Include"),
            by = c("subjid", "param", "agedays")]

    # POSITIONAL STRING DETECTION
    # Changed from value-based to positional approach to match Stata lines 849-883
    # Originators are Includes where the NEXT value is a CF
    # Strings propagate forward until interrupted by non-CF
    # Originators are ANY Include where next is CF
    # The ageday_has_include check only applies to CFs (for rescue eligibility), not originators

    # Initialize variables
    data.df[, cf_binary := exclude == "Exclude-Carried-Forward"]

    # Process by subject-param to maintain ordering
    data.df[, ':=' (
      nextcf = shift(cf_binary, type = "lead", fill = FALSE),
      priorcf = shift(cf_binary, type = "lag", fill = FALSE),
      originator = FALSE,
      cf_string_num = NA_integer_,
      originator_z = NA_real_
    ), by = c("subjid", "param")]

    # Identify originators: Include where next value is a CF
    # Any Include can be an originator - ageday_has_include only restricts CF rescue eligibility
    data.df[exclude == "Include" & nextcf == TRUE, originator := TRUE]

    # Assign sequential string numbers to originators (per subject-param)
    data.df[, originator_seq := cumsum(originator), by = c("subjid", "param")]
    data.df[originator == TRUE, cf_string_num := originator_seq]

    # Store originator z-scores (use sd.orig_uncorr to match Stata's s<param>z)
    data.df[originator == TRUE, originator_z := sd.orig_uncorr]

    # Get max number of CFs to determine loop iterations
    max_cf_count <- data.df[cf_binary == TRUE, .N, by = c("subjid", "param")]
    max_iterations <- if (nrow(max_cf_count) > 0) max(max_cf_count$N) else 0

    # Propagate string numbers and originator z-scores forward to consecutive CFs
    # Only propagate to CFs that are eligible (not on days with includes)
    if (max_iterations > 0) {
      for (i in 1:max_iterations) {
        data.df[, ':=' (
          cf_string_num = ifelse(
            cf_binary == TRUE &
            (ageday_has_include == FALSE | is.na(ageday_has_include)) &
            is.na(cf_string_num),
            shift(cf_string_num, type = "lag"),
            cf_string_num
          ),
          originator_z = ifelse(
            cf_binary == TRUE &
            (ageday_has_include == FALSE | is.na(ageday_has_include)) &
            is.na(originator_z),
            shift(originator_z, type = "lag"),
            originator_z
          )
        ), by = c("subjid", "param")]
      }
    }

    # Calculate string length (number of CFs in each string, excluding originator)
    data.df[, cf_string_length := sum(cf_binary == TRUE), by = c("subjid", "param", "cf_string_num")]

    # Create seq_win variable (position in string: 0 for originator, 1/2/3... for CFs)
    data.df[, seq_win := NA_integer_]
    data.df[originator == TRUE, seq_win := 0]
    data.df[!is.na(cf_string_num) & cf_binary == TRUE,
            seq_win := seq_len(.N),
            by = c("subjid", "param", "cf_string_num")]

    # Create cs variable for compatibility with rescue code logic
    data.df[, cs := cf_string_num]

    # Calculate absdiff (absolute z-score difference from originator)
    # Double-round to handle floating point - first to 0.001 (clean up noise),
    # then to 0.01 for threshold comparison. This matches Stata where subtraction of 0.001-rounded
    # values produces a clean 0.001 multiple before rounding to 0.01.
    data.df[!is.na(seq_win), absdiff := janitor::round_half_up(
      janitor::round_half_up(abs(sd.orig_uncorr - originator_z), 3), 2)]

    # Clean up temporary variables
    data.df[, c("cf_binary", "is_eligible_include", "nextcf", "priorcf",
                "originator", "originator_seq", "cf_string_num", "originator_z", "cf_string_length") := NULL]
    
    # Fix CF rescue logic
    # Fix 1: Add same-day include check - CFs on days with includes are not eligible for rescue
    # Fix 2: Changed > to >= for adolescent thresholds
    # ageday_has_include now calculated earlier and used in string detection
    # CFs with ageday_has_include already excluded from cs/seq_win assignment

    # Only apply rescue codes to CFs with valid cs assignment (already filtered for ageday_has_include)
    data.df[!is.na(seq_win) & (ageday_has_include == FALSE | is.na(ageday_has_include)),
            exclude := (function(df){
      # only 1 cf in string - codes 4 and 5
      if (max(seq_win) == 1){

        # Code 4: only 1 CF AND |diff| < 0.05
        # Reverted to original threshold
        df[seq_win != 0 & absdiff < 0.05,
           exclude := "Exclude-1-CF-deltaZ-<0.05"]

        # Code 5: only 1 CF AND |diff| >= 0.05 AND < 0.1 AND wholehalfimp
        # Reverted to original threshold
        df[seq_win != 0 & absdiff >= 0.05 & absdiff < .1 & wholehalfimp,
           exclude := "Exclude-1-CF-deltaZ-<0.1-wholehalfimp"]

      } else if (max(seq_win) > 1){
        # 2+ CFs in string - only teens get rescue (codes 6 and 7)

        # Code 6: any # of consecutive CFs AND adol AND |diff| < 0.05
        # Changed > to >= for adolescent thresholds
        # Reverted to original threshold
        df[seq_win != 0 & agedays/365.25 >= 16 & sex == 1 & absdiff < 0.05,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.05"]
        df[seq_win != 0 & agedays/365.25 >= 17 & sex == 0 & absdiff < 0.05,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.05"]

        # Code 7: any # of consecutive CFs AND adol AND |diff| >= 0.05 AND < 0.1 AND wholehalfimp
        # Reverted to original threshold
        df[seq_win != 0 & agedays/365.25 >= 16 & sex == 1 &
             absdiff >= 0.05 & absdiff < .1 & wholehalfimp,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp"]
        df[seq_win != 0 & agedays/365.25 >= 17 & sex == 0 &
             absdiff >= 0.05 & absdiff < .1 & wholehalfimp,
           exclude := "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp"]
        # REMOVED: "CP ADD" section - codes 4/5 should NOT apply to multi-CF strings
      }

      return(df$exclude)
    })(copy(.SD)),
    .SDcols = c('agedays', "seq_win", 'absdiff', "sex", "cs", "wholehalfimp",
                "exclude"),
    by = c("subjid", "param", "cs")]

    } # End if (any_cf)
  }

  # Add SDE-Identical rows back after CF rescue
  if (nrow(sde_identical_rows) > 0) {
    data.df <- rbind(data.df, sde_identical_rows, fill = TRUE)
    data.df <- data.df[order(subjid, param, agedays, id)]
  }

  # Debug save point for Step 6 output
#  if (!is.null(debug_step) && debug_step == 6) {
#    saveRDS(data.df, paste0("R-step6-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step6-output.rds - STOPPING EXECUTION")
#    return(data.df[, .(line, exclude, tbc.sd, param)])
#  }

  # Step 7: BIV ----
  # NOTE: Step numbers aligned with Stata 2025-12-10: BIV is Step 7, Evil Twins is Step 9

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
  # Removed nnte filter (nnte calculation removed)
  valid_set <- valid(data.df, include.temporary.extraneous = TRUE)

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
  # data.df[valid_set & param == "WEIGHTKG" & v < 1 &
  #           agedays == 0,
  #         exclude := exc_nam]
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
  # Added check to skip rows already marked Absolute-BIV
  # This prevents Standardized-BIV from overwriting Absolute-BIV when both conditions are true.
  # NOTE: This is the only step where an exclusion code can overwrite a non-temporary exclusion
  # code (e.g., this could overwrite CF codes). No other steps should overwrite non-temporary codes.
  abs_biv <- "Exclude-Absolute-BIV"

  # identify z cutoff
  # ***Note, using unrecentered values***
  #  *For weight only do after birth
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < -25 &
            ageyears < 1 & exclude != abs_biv,
          exclude := exc_nam]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < -15 &
            ageyears >= 1 & exclude != abs_biv,
          exclude := exc_nam]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr > 22 &
            exclude != abs_biv,
          exclude := exc_nam]

  # *Max z-score for height based on analysis of CHOP data because 15/25 too loose for upper limits
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < -25 &
            ageyears < 1 & exclude != abs_biv,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < -15 &
            ageyears >= 1 & exclude != abs_biv,
          exclude := exc_nam]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr > 8 &
            exclude != abs_biv,
          exclude := exc_nam]

  # head circumference
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr < -15 &
            exclude != abs_biv,
          exclude := exc_nam]
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr > 15 &
            exclude != abs_biv,
          exclude := exc_nam]

  # then, do some explicit overwriting for the 0 case (otherwise will be
  # set as missing)
  data.df[v == 0, exclude := exc_nam]

  # 7d.  Replace exc_*=0 if exc_*==2 & redo step 5 (temporary extraneous)
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # Step 9: Evil Twins ----
  # Evil Twins: An important weakness in the original pediatric growthcleanr algorithm was that it often failed to identify two or more implausible measurements that occurred next to each other, even if they were extremely deviant from a child's other measurements. This step is now added to identify these multiple extreme values, although it also identifies some single values.

  exc_nam <- "Exclude-Evil-Twins"

  # MUST reorder BEFORE computing valid_set
  # The valid_set boolean vector is computed based on row positions.
  # If we reorder AFTER computing valid_set, the boolean indices no longer match
  # the correct rows, causing temp SDEs to be included and non-temp-SDEs excluded.
  # Must include agedays and id for consistent order
  # Without agedays, calc_oob_evil_twins compares non-adjacent-in-time values
  # Without id, SDE rows (same ageday) have undefined order, causing parallel inconsistency
  data.df <- data.df[order(subjid, param, agedays, id),]

  # create the valid set
  # we only running carried forwards on valid values, non NNTE values,
  # and non single values, and non pair
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single_pairs <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 2]
  # Removed nnte filter (nnte calculation removed)
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    not_single_pairs

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

      # Fix iteration condition
      # Changed >= 2 to >= 1 to allow cascading exclusions
      # After excluding 1 value, recalculating may reveal new OOB values
      any_oob <- any(upd.df$sum_oob >= 1)
      # while there are oob values, we want to remove
      while (any_oob){

        # 9D
        # now calculate the maximum difference from the median tbc.sd
        upd.df[, `:=` (sd_med = median(tbc.sd, na.rm = TRUE)), by =.(subjid, param)]
        upd.df[, `:=` (med_diff = abs(tbc.sd - sd_med)), by =.(subjid, param)]
        # Replace max_diff logic with tiebreaker hierarchy
        # Old logic marked ALL tied values, causing multiple exclusions per iteration
        # New logic: Select exactly ONE value per subject/param using tiebreaker:
        #   1. Furthest from median (highest med_diff)
        #   2. Most extreme overall (highest abs(tbc.sd))
        #   3. Lowest id (deterministic)

        # Get unique subject/param combinations with OOB values
        sp_combos <- unique(upd.df[oob == TRUE, .(subjid, param)])

        if (nrow(sp_combos) > 0) {
          # For each subject/param, find the one value to exclude
          # Use data.table by-group processing with .SD
          worst_lines <- upd.df[oob == TRUE, {
            # Sort by tiebreaker hierarchy and take first
            ord <- order(-med_diff, -abs(tbc.sd), id)
            .(worst_line = line[ord[1]])
          }, by = .(subjid, param)]

          # Mark the selected values for exclusion
          upd.df[line %in% worst_lines$worst_line, exclude := exc_nam]
        }

        df[upd.df[exclude == exc_nam,], exclude := i.exclude, on = .(line)]

        #10E
        # reupdate valid (to recalculate OOB -- others are not included)
        upd.df <- calc_oob_evil_twins(df[valid(df),])
        upd.df[, `:=` (sum_oob = sum(oob, na.rm = TRUE)), by =.(subjid, param)]

        # Changed >= 2 to >= 1
        any_oob <- any(upd.df$sum_oob >= 1)

      }

      return(df$exclude)
    })(copy(.SD))]

  }

  # 9F.  redo temp extraneous
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # Step 11: Extreme EWMA ----
  # Restructured to use global iterations for efficiency
  # Key changes:
  #   1. Global while loop instead of per-subject-param while loops
  #   2. Only Include values participate in EWMA (temp SDEs do NOT get EWMA calculated)
  #   3. After each iteration, only re-process subject-params that had new exclusions
  #   4. Temp SDE recalculation only for subjects with new exclusions
  # This provides major speedup when most subjects are clean after early iterations.

  # 11.  Exclude extreme errors with EWMA
  # a.  Erroneous measurements can distort the EWMA for measurements around them. Therefore, if the EWMA
  #     method identifies more than one value for a subject and parameter that meets criteria for exclusion,
  #     we will only exclude the value that deviates the most from expected in any given step. Then we will
  #     repeat the entire process until no measurements are identified that meet criteria for exclusion.
  # b.  Perform a EWMA calculation
  #   i.  Only use values where exc_*==0 to determine the EWMAs and for exclusion candidates.
  #       Temp SDEs (exc_*==2) do NOT participate in EWMA calculation or exclusion in this step.

  if (!quietly)
    cat(sprintf("[%s] Exclude extreme measurements based on EWMA...\n", Sys.time()))

  # Pre-identify subject-params that need processing:
  # Count INCLUDE values only (not temp SDEs) - need >2 for EWMA processing
  data.df[, sp_key := paste0(subjid, "_", param)]
  # Removed nnte filter (nnte calculation removed)
  include_counts <- data.df[valid(data.df, include.temporary.extraneous = FALSE),
                            .(n_include = .N), by = sp_key]
  sp_to_process <- include_counts[n_include > 2, sp_key]

  # Track subjects with SDEs for targeted temp SDE recalculation
  subj_with_sde <- unique(data.df[exclude == "Exclude-Temporary-Extraneous-Same-Day", subjid])

  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      cat(sprintf("  EWMA1 iteration %d: %d subject-params\n", iteration, length(sp_to_process)))

    # Track which rows had EWMA1 exclusions before this iteration
    data.df[sp_key %in% sp_to_process, had_ewma1_before := (exclude == 'Exclude-EWMA1-Extreme')]

    # Process each subject-param - ONE exclusion per subject-param per iteration
    data.df[sp_key %in% sp_to_process,
            exclude := (function(df) {
              # Only Include values participate (no temp SDEs)
              # Removed nnte filter
              include_set <- valid(df, include.temporary.extraneous = FALSE)

              if (sum(include_set) > 2) {
                # Initialize EWMA fields
                df[, (ewma.fields) := as.double(NaN)]

                # Calculate exp_vals based on age gaps (Include values only)
                # Use linear interpolation for exponent
                # Use Stata's linear formula
                df_sub <- df[include_set,]
                tmp_ages <- data.frame(
                  "before" = abs(df_sub$agedays - c(NA, df_sub$agedays[1:(nrow(df_sub)-1)])),
                  "after" = abs(df_sub$agedays - c(df_sub$agedays[2:(nrow(df_sub))], NA))
                )
                maxdiff <- sapply(1:nrow(tmp_ages), function(x){max(tmp_ages[x,], na.rm = TRUE)})
                ageyears <- maxdiff / 365.25
                exp_vals <- ifelse(ageyears <= 1, -1.5,
                            ifelse(ageyears >= 3, -3.5,
                                   -1.5 - (ageyears - 1)))
                df[include_set, exp_vals := exp_vals]

                # Calculate EWMA for Include values only (no temp SDEs)
                df[include_set, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
                df[include_set, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]

                # Calculate dewma for Include values only
                # No intermediate rounding - double-round at comparison points only
                df[include_set, `:=`(
                  dewma.all = tbc.sd - ewma.all,
                  dewma.before = tbc.sd - ewma.before,
                  dewma.after = tbc.sd - ewma.after,
                  c.dewma.all = ctbc.sd - c.ewma.all
                )]

                # Identify potential exclusions - Include values only
                df[, pot_excl := FALSE]
                # Fixed cdewma sign for negative outliers
                # Was: c.dewma.all > 3.5 for negative case (wrong - should be < -3.5)
                # Reverted to rounding with round_half_up
                # Changed threshold comparisons from 3 to 2 decimal places
                # Double-round (3 then 2 decimals) to handle floating-point noise
                df[include_set, pot_excl :=
                     (janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 3.5 & janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) > 3 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) > 3 & janitor::round_half_up(janitor::round_half_up(tbc.sd, 3), 2) > 3.5 &
                        ((!is.na(ctbc.sd) & janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3.5) | is.na(ctbc.sd))
                     ) |
                     (janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -3.5 & janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) < -3 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) < -3 & janitor::round_half_up(janitor::round_half_up(tbc.sd, 3), 2) < -3.5 &
                        ((!is.na(ctbc.sd) & janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3.5) | is.na(ctbc.sd))
                     )
                ]

                # First and last Include values can't be excluded in this step
                include_indices <- which(include_set)
                if (length(include_indices) > 0) {
                  df[include_indices[c(1, length(include_indices))], pot_excl := FALSE]
                }

                # Exclude at most ONE value (the worst)
                num.exclude <- sum(df$pot_excl, na.rm = TRUE)
                if (num.exclude == 1) {
                  df[pot_excl == TRUE, exclude := 'Exclude-EWMA1-Extreme']
                } else if (num.exclude > 1) {
                  # Round to 3 decimals to avoid floating point precision issues
                  # Use id tie-breaker for deterministic selection
                  # Double-round (3 then 2 decimals) for consistency
                  worst.row <- with(df, order(pot_excl, janitor::round_half_up(janitor::round_half_up(abs(tbc.sd + dewma.all), 3), 2), -id, decreasing = TRUE))[1]
                  df[worst.row, exclude := 'Exclude-EWMA1-Extreme']
                }
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('index', 'id', 'sex', 'agedays', 'tbc.sd', 'ctbc.sd', 'nnte', 'exclude')]

    # Capture EWMA1 iteration 1 values for debugging
    if (iteration == 1) {
      # Copy EWMA values from first iteration to permanent columns
      # Naming: ewma1_it1.ewma_all, ewma1_it1.ewma_before, ewma1_it1.dewma_all, etc.
      ewma1_cols <- c("ewma.all", "ewma.before", "ewma.after", "dewma.all", "dewma.before", "dewma.after")
      for (col in ewma1_cols) {
        new_col <- paste0("ewma1_it1.", gsub("\\.", "_", col))
        if (col %in% names(data.df)) {
          data.df[, (new_col) := get(col)]
        }
      }
    }

    # Find subject-params with NEW exclusions this iteration
    data.df[sp_key %in% sp_to_process, has_new_excl :=
              (exclude == 'Exclude-EWMA1-Extreme') & !had_ewma1_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])
    subjects_with_new_excl <- unique(data.df[has_new_excl == TRUE, subjid])

    # Recalculate temp SDEs only for subjects that had exclusions AND have SDEs
    affected_sde_subj <- intersect(subjects_with_new_excl, subj_with_sde)
    if (length(affected_sde_subj) > 0) {
      # Reset temp SDEs to Include for affected subjects
      data.df[subjid %in% affected_sde_subj & exclude == 'Exclude-Temporary-Extraneous-Same-Day',
              exclude := 'Include']
      # Recalculate temp SDEs for subset
      sde_subset <- data.df[subjid %in% affected_sde_subj]
      sde_result <- temporary_extraneous_infants(sde_subset)
      # Apply results back using index
      sde_indices_to_mark <- sde_subset$index[sde_result]
      data.df[index %in% sde_indices_to_mark, exclude := 'Exclude-Temporary-Extraneous-Same-Day']
    }

    # Next iteration: only subject-params with new exclusions
    sp_to_process <- sp_with_new_excl

    # Cleanup iteration tracking columns
    data.df[, c("had_ewma1_before", "has_new_excl") := NULL]
  }

  # Cleanup
  data.df[, sp_key := NULL]

  # Final temp SDE recalculation (end of Step 11, before Step 13)
  # Note: Does NOT use exclude_from_dop_ids - that's only for Step 13
  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df), exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # Debug exit for step 11
  # Fixed to check global environment like step 2
#  if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 11) {
#    saveRDS(data.df, paste0("R-step11-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    if (!quietly) cat(sprintf("[%s] DEBUG: Stopping after Step 11 (EWMA1 Extreme)\n", Sys.time()))
#    if (!quietly) cat(sprintf("[%s] DEBUG: Output saved to R-step11-output-*.rds\n", Sys.time()))
#    return(data.df[, .(id, exclude)])
#  }

  # 13: SDEs ----

  # Debug save point for Step 13 input
  # Fixed to check global environment like step 11
#  if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 13) {
#    saveRDS(data.df, paste0("R-step13-input-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step13-input.rds")
#  }

  if (!quietly)
    cat(sprintf(
      "[%s] Exclude same day extraneous...\n",
      Sys.time()
    ))
  
  # Pass temp SDE ids to exclude from DOP median
  # This matches Stata Step 13 line 1516: DOP median uses only fully included values (exc==0)
  temp_sde_ids_step13 <- data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', id]

  # Debug - check state BEFORE restoring temp SDEs
  pre_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(pre_dupe_check) > 0) {
    message("DEBUG Step13-PRE: Found ", nrow(pre_dupe_check), " Include duplicates BEFORE restoring temp SDEs")
    fwrite(pre_dupe_check, "debug-step13-PRE-duplicates.csv")
    # Save full rows for duplicates
    dupe_rows <- data.df[exclude == "Include"][pre_dupe_check, on = .(subjid, param, agedays)]
    fwrite(dupe_rows, "debug-step13-PRE-duplicate-rows.csv")
    # Save ALL data at this point
    fwrite(data.df, "debug-step13-PRE-full-data.csv")
    stop("DEBUG: Stopping at Step 13-PRE - duplicate Includes found. See debug-step13-PRE-*.csv files")
  }

  data.df[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.df[temporary_extraneous_infants(data.df, exclude_from_dop_ids = temp_sde_ids_step13),
          exclude := 'Exclude-Temporary-Extraneous-Same-Day']

  # Debug - check state AFTER temp SDE marking
  post_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(post_dupe_check) > 0) {
    message("DEBUG Step13-POST temp SDE: Found ", nrow(post_dupe_check), " Include duplicates AFTER temp SDE marking")
    fwrite(post_dupe_check, "debug-step13-POST-tempSDE-duplicates.csv")
    # Save full rows for duplicates
    dupe_rows <- data.df[exclude == "Include"][post_dupe_check, on = .(subjid, param, agedays)]
    fwrite(dupe_rows, "debug-step13-POST-tempSDE-duplicate-rows.csv")
    # Save ALL data at this point
    fwrite(data.df, "debug-step13-POST-tempSDE-full-data.csv")
    stop("DEBUG: Stopping at Step 13-POST temp SDE - duplicate Includes found. See debug-step13-POST-tempSDE-*.csv files")
  }
  keep_cols_sde <- names(data.df)

  # Pre-filter to subjects with same-day measurements (potential SDEs)
  # If a subject has no same-day measurements, SDE steps can do nothing - skip entirely
  # Exclude SDE-Identicals from count - they're already handled and shouldn't trigger SDE processing
  sde_day_counts <- data.df[exclude %in% c("Include", "Exclude-Temporary-Extraneous-Same-Day"),
                            .(n_on_day = .N), by = .(subjid, param, agedays)]
  subj_with_sde_days <- unique(sde_day_counts[n_on_day > 1, subjid])
  n_total_subj <- uniqueN(data.df$subjid)
  n_with_sde <- length(subj_with_sde_days)
  if (!quietly)
    cat(sprintf("  SDE pre-filter: %d/%d subjects have same-day measurements (%.1f%%)\n",
                n_with_sde, n_total_subj, 100*n_with_sde/n_total_subj))

  # Only process subjects with same-day measurements
  if (n_with_sde == 0) {
    # No SDEs to process - skip entire SDE section
    if (!quietly) cat("  No SDEs to process, skipping SDE step\n")
    data.sde <- data.df[0, ]  # Empty data.table with same structure
  } else {
    # Filter to subjects with potential SDEs before dplyr chain
    data.df_sde_subset <- data.df[subjid %in% subj_with_sde_days]

    data.df_sde_subset <- data.df_sde_subset %>%
      as.data.frame()  %>%
      # subset(valid(data.df, include.temporary.extraneous = TRUE) & nnte_full == FALSE) %>% mutate(id = id) %>%
      # select(subjid, id = id, agedays, param, v, tbc.sd, exclude, nnte_full) %>%
      group_by(subjid, param, agedays) %>%
      mutate(sde_this = case_when(
        any(exclude == "Exclude-Temporary-Extraneous-Same-Day") ~ TRUE,
        TRUE ~ FALSE
      ))
    data.sde <- data.df_sde_subset %>%
      group_by(subjid) %>%
      # Removed nnte filter (nnte calculation removed)
      filter(valid(cur_data(), include.temporary.extraneous = TRUE)) %>%
      mutate(id = id) %>%
      filter(any(sde_this == TRUE)) %>%
    arrange(subjid, param, agedays, id) %>%
    group_by(subjid, param, agedays) %>%
    mutate(
      # Age-dependent id selection
      # Age 0: keep lowest id, age > 0: keep highest id
      keep_id = ifelse(agedays == 0, min(id, na.rm = TRUE), max(id, na.rm = TRUE)),
      exclude = case_when(
        length(unique(v)) == 1 &
          id != keep_id ~ "Exclude-SDE-Identical",
        TRUE ~ exclude
      )
    ) %>%
    group_by(subjid, param, agedays, v) %>%
    mutate(
      # Age-dependent id selection for duplicate values
      keep_id_dup = ifelse(agedays == 0, min(id, na.rm = TRUE), max(id, na.rm = TRUE)),
      dup_count = n(),
      # Fix SDE-Identical bug: was exclude != "Include" (backwards)
      # Should mark currently-included duplicates as SDE-Identical (matches Stata line 206: exc==0)
      exclude = case_when(
        dup_count > 1 &
          id != keep_id_dup &
          exclude == "Include" ~ "Exclude-SDE-Identical",
        TRUE ~ exclude
      )
    ) %>%
    select(-keep_id, -keep_id_dup) 
  data.sde <- data.sde %>% data.table()
  # Include id for deterministic SDE order
  setkey(data.sde, subjid, param, agedays, id)


  # Track temp SDE status before restoration
  data.sde[, was_temp_sde := exclude == 'Exclude-Temporary-Extraneous-Same-Day']
  data.sde[exclude == 'Exclude-Temporary-Extraneous-Same-Day', exclude := 'Include']
  data.sde <- data.sde %>% as.data.frame() %>% arrange(subjid, param, agedays, id)

  data.sde <- data.sde %>%
    # -------------------------------------------
  # Identify subject–parameter pairs eligible for "One-Day SDE" logic
  # -------------------------------------------
  group_by(subjid, param) %>%
    mutate(
      # Use was_temp_sde instead of checking exclude
      # (temp SDEs were already converted to Include at line 3120)
      n_days_with_data = n_distinct(agedays[exclude == "Include"]),
      has_sde_day = any(was_temp_sde == TRUE),
      one_day_sde_flag = (n_days_with_data == 1 & has_sde_day)
    ) %>%
    
    # -------------------------------------------
  # Only run the following steps for one-day SDE subjects/params
  # -------------------------------------------
  group_by(subjid, param, agedays) %>%
    mutate(
      # Include temp SDEs in median for SDE-All-Extreme check
      # For one-day SDEs, Stata uses ALL same-day values to calculate median (including temp SDEs)
      # This ensures both values in a 2-value SDE have the same absdiff_rel_to_median
      median_tbc = median(tbc.sd[exclude == "Include"], na.rm = TRUE),

      absdiff_rel_to_median = case_when(
        exclude == "Include" ~ abs(tbc.sd - median_tbc),
        TRUE ~ NA
      ),

      min_absdiff_rel_to_median = min(absdiff_rel_to_median, na.rm = TRUE),

      # Double-round threshold comparisons (3→2 decimals)
      median_tbc_gt2 = ifelse(
        one_day_sde_flag &
          janitor::round_half_up(janitor::round_half_up(min(absdiff_rel_to_median, na.rm = TRUE), 3), 2) > 2,
        TRUE,
        FALSE
      )
    ) %>%
    ungroup() %>%
    mutate(
      # SDE-All-Extreme must check one_day_sde_flag
      # Double-round threshold comparisons (3→2 decimals)
      exclude = ifelse(
        one_day_sde_flag & !is.na(absdiff_rel_to_median) & janitor::round_half_up(janitor::round_half_up(min_absdiff_rel_to_median, 3), 2) > 2,
        "Exclude-SDE-All-Extreme",
        exclude
      ),
      # Fixed condition from ==1 to >=2
      # Need >=2 Includes to trigger selection logic (select one to keep, exclude others)
      # Must group by subjid/param/agedays for sum/min
      # Without grouping, sum/min operate on entire dataset, causing incorrect exclusions
    )

  # Refactored SDE-One-Day to use sort logic (simpler, consistent with temp SDE)
  # Calculate DOP medians first, then use single sort to select one value to keep
  data.sde <- data.sde %>%
    # Calculate DOP medians (needs cross-parameter grouping)
    # DOP median uses FULLY included values only
    group_by(subjid, agedays) %>%
    mutate(
      HT_dop_med = ifelse(one_day_sde_flag & param == "HEIGHTCM",
                          median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "WEIGHTKG"], na.rm = TRUE), NA),
      WT_dop_med = ifelse(one_day_sde_flag & param == "WEIGHTKG",
                          median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "HEIGHTCM"], na.rm = TRUE), NA),
      HC_dop_med = ifelse(one_day_sde_flag & param == "HEADCM",
                          median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "HEIGHTCM"], na.rm = TRUE), NA),
      absdiff_dop_med = case_when(
        one_day_sde_flag & exclude == "Include" & param == "HEIGHTCM" ~ abs(tbc.sd - HT_dop_med),
        one_day_sde_flag & exclude == "Include" & param == "WEIGHTKG" ~ abs(tbc.sd - WT_dop_med),
        one_day_sde_flag & exclude == "Include" & param == "HEADCM" ~ abs(tbc.sd - HC_dop_med),
        TRUE ~ NA
      )
    ) %>%
    ungroup() %>%
    # Now use sort logic to select one value per SDE group
    group_by(subjid, param, agedays) %>%
    mutate(
      # Round values for sorting (consistent with previous logic)
      # Double-round to handle floating point noise from z-score subtraction
      absdiff_median_rounded = janitor::round_half_up(janitor::round_half_up(absdiff_rel_to_median, 3), 2),
      absdiff_dop_rounded = janitor::round_half_up(janitor::round_half_up(absdiff_dop_med, 3), 2),
      # For sorting: NA DOP values should sort last (use Inf)
      absdiff_dop_for_sort = ifelse(is.na(absdiff_dop_rounded), Inf, absdiff_dop_rounded),
      # Age-dependent id tiebreaker: agedays==0 picks lowest internal_id, otherwise highest
      # Use internal_id (sequential 1:n) to match Stata's obsid
      tiebreaker_oneday = if(agedays[1] == 0) internal_id else -internal_id,
      # Find the id to keep: sort by absdiff_median (asc), absdiff_dop (asc), tiebreaker (asc)
      keep_id_oneday = {
        eligible_mask <- one_day_sde_flag & exclude == "Include"
        if(sum(eligible_mask, na.rm = TRUE) < 2) NA_integer_
        else {
          eligible_ids <- id[eligible_mask]
          eligible_median <- absdiff_median_rounded[eligible_mask]
          eligible_dop <- absdiff_dop_for_sort[eligible_mask]
          eligible_tiebreaker <- tiebreaker_oneday[eligible_mask]
          ord <- order(eligible_median, eligible_dop, eligible_tiebreaker)
          eligible_ids[ord[1]]
        }
      },
      # Simple logic: keep_id gets Include, all other eligible values get Exclude-SDE-One-Day
      exclude = case_when(
        one_day_sde_flag & exclude == "Include" & !is.na(keep_id_oneday) &
          id == keep_id_oneday ~ "Include",
        one_day_sde_flag & exclude == "Include" & !is.na(keep_id_oneday) &
          id != keep_id_oneday ~ "Exclude-SDE-One-Day",
        TRUE ~ exclude
      )
    ) %>%
    ungroup()

  # h. no need to calculate for R.

  # D. SDE-EWMAs

  # a. Calculate EWMA for full included vals
  # Exclude temp SDEs from EWMA calculation
  # Stata uses only Include values (not temp SDEs) for EWMA. This ensures:
  # 1. Only one value per same-day contributes to EWMA (matches Stata behavior)
  # 2. The diff calculation for exponent is correct (no same-day 0 diffs)
  ewma_df <- as.data.table(data.sde %>% filter(exclude == "Include" & !was_temp_sde))
  # Include id for deterministic order
  setkey(ewma_df, subjid, param, agedays, id)
  
  # tmp <- data.frame(
  #           "before" = abs(ewma_df$agedays - c(NA, ewma_df$agedays[1:(nrow(ewma_df)-1)])),
  #           "after" = abs(ewma_df$agedays - c(ewma_df$agedays[2:(nrow(ewma_df))], NA))
  #         )
  # maxdiff <- pmax(tmp$before, tmp$after, na.rm = TRUE)
  # # maxdiff <- sapply(1:nrow(tmp), function(x){max(tmp[x,], na.rm = TRUE)})
  # # exp_vals <- rep(-1.5, nrow(tmp))
  # # exp_vals[maxdiff > 365.25] <- -2.5
  # # exp_vals[maxdiff > 730.5] <- -3.5
  # # approximate continuous exponential decay equivalent to Stata's add(5)
  # # convert each observation's mean time gap to a smooth exponent
  # For each subject and parameter, calculate before/after gaps and assign exponent
  ewma_df[, c("diff_before", "diff_after") :=
            .(c(NA, diff(agedays)), c(diff(agedays), NA)),
          by = .(subjid, param)]
  
  # Take the larger gap (the "Ba" value) and convert to years
  ewma_df[, maxdiff := pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)]
  ewma_df[, ageyears := maxdiff / 365.25]
  
  # Apply the exact linear exponent rule from the specification
  ewma_df[, exp_val := fcase(
    ageyears <= 1, -1.5,
    ageyears >= 3, -3.5,
    default = -1.5 - (ageyears - 1)
  )]
  
  # Compute EWMA using these exponents
  ewma.fields <- c("ewma.all", "ewma.before", "ewma.after")
  ewma_df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_val, TRUE), by = .(subjid, param)]
  # avg_gap <- rowMeans(tmp, na.rm = TRUE)
  # exp_vals <- -1.5 * exp(-avg_gap / 365.25)
  # ewma.fields <- c('ewma.all', 'ewma.before', 'ewma.after')
  # ewma_df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)] %>% as.data.frame()
  data.sde <- merge(data.sde, ewma_df %>%
                      select(id, which(grepl("ewma", names(ewma_df)))), by = "id", all.x = TRUE)
  data.sde <- data.sde %>%
    group_by(subjid, param, agedays) %>%
    # b. Assign EWMAs to partially included vals (temp SDEs).
    # Use was_temp_sde marker instead of exclude
    # At this point, temp SDEs have been restored to Include, so we use the marker
    # Temp SDEs need EWMA assigned so they can participate in SDE-EWMA resolution
    mutate(ewma.all = case_when(is.na(ewma.all) & was_temp_sde ~ max(ewma.all[!was_temp_sde], na.rm = TRUE),
                                TRUE ~ ewma.all),
           # Use max(ewma) as reference like Stata
           # Stata: bysort subjidsde_`p' agedays: egen spa_ewma_`p'=max(ewma_`p')
           # Stata: gen absdewma_`p'=abs(tbc`p'z-spa_ewma_`p')
           spa_ewma = max(ewma.all, na.rm = TRUE),
           # Restored rounding with round_half_up
           # Double-round to handle floating point noise from z-score subtraction
           absdewma = janitor::round_half_up(janitor::round_half_up(abs(tbc.sd - spa_ewma), 3), 2))

  # Refactored SDE-EWMA to use sort logic (simpler, no edge case bugs)
  # Previously used complex Rule set 1/2 logic that missed non-min values when ties existed
  data.sde <- data.sde %>%
    group_by(subjid, param, agedays) %>%
    mutate(
      n_available = sum(exclude %in% c("Exclude-Temporary-Extraneous-Same-Day",
                                       "Include"), na.rm = TRUE),
      min_absdewma = suppressWarnings(min(absdewma[exclude %in% c("Exclude-Temporary-Extraneous-Same-Day",
                                                                  "Include")],
                                          na.rm = TRUE)),
      min_absdewma = ifelse(is.infinite(min_absdewma), NA_real_, min_absdewma),
      # Double-round min_absdewma before threshold comparison
      min_absdewma = janitor::round_half_up(janitor::round_half_up(min_absdewma, 3), 2),
      # First check: SDE-All-Extreme if min_absdewma > 1
      exclude = case_when(
        n_available >= 2 &
          exclude %in% c("Exclude-Temporary-Extraneous-Same-Day", "Include") &
          !is.na(min_absdewma) & min_absdewma > 1 ~ "Exclude-SDE-All-Extreme",
        TRUE ~ exclude
      )
    ) %>%
    # Now use sort logic: sort by absdewma, then id tiebreaker; keep first, exclude others
    group_by(subjid, param, agedays) %>%
    mutate(
      # Age-dependent id tiebreaker: agedays==0 picks lowest internal_id, otherwise highest
      # Use internal_id (sequential 1:n) to match Stata's obsid
      tiebreaker_ewma = if(agedays[1] == 0) internal_id else -internal_id,
      # Find the id to keep: sort eligible values by absdewma (asc), then tiebreaker (asc)
      keep_id_ewma = {
        eligible_mask <- exclude %in% c("Exclude-Temporary-Extraneous-Same-Day", "Include")
        if(sum(eligible_mask) == 0) NA_integer_
        else {
          eligible_ids <- id[eligible_mask]
          eligible_absdewma <- absdewma[eligible_mask]
          eligible_tiebreaker <- tiebreaker_ewma[eligible_mask]
          ord <- order(eligible_absdewma, eligible_tiebreaker)
          eligible_ids[ord[1]]
        }
      },
      # Simple logic: keep_id gets Include, all other eligible values get Exclude-SDE-EWMA
      exclude = case_when(
        exclude %in% c("Exclude-Temporary-Extraneous-Same-Day", "Include") &
          id == keep_id_ewma ~ "Include",
        exclude %in% c("Exclude-Temporary-Extraneous-Same-Day", "Include") &
          id != keep_id_ewma ~ "Exclude-SDE-EWMA",
        TRUE ~ exclude
      )
    ) %>%
    ungroup()
  }  # End of SDE pre-filter else block

  #
  # print(data.sde)
  data.df <- merge(data.df, data.sde %>% select(id, sde_exclude = exclude),by = "id", all.x = TRUE) %>% mutate(exclude = case_when(!is.na(sde_exclude) ~ sde_exclude,
                                                                                                                                   TRUE ~ exclude))
  # write.csv(data.sde, "../data/testoutsdeneww.csv")
  data.df <- data.df %>% select(which(names(data.df) %in% keep_cols_sde))
  data.df <- as.data.table(data.df)
  # Include id for deterministic SDE order
  setkey(data.df, subjid, param, agedays, id)

  # 15-16: moderate EWMA ----
  # Restructured to use global iterations for efficiency
  # Key changes:
  #   1. Pre-calculate p_plus/p_minus and their z-scores ONCE for all valid rows
  #   2. Global while loop instead of per-subject-param while loops
  #   3. After each iteration, only re-process subject-params that had new exclusions
  # This provides major speedup when most subjects are clean after early iterations.

  if (!quietly)
    cat(sprintf("[%s] Exclude moderate EWMA...\n", Sys.time()))

  # Order data for processing
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, id),]

  # Pre-identify subject-params that need processing: Include values, >2 values
  # Include temp SDEs in count so they get p_plus/p_minus calculated
  data.df[, sp_key := paste0(subjid, "_", param)]
  include_counts <- data.df[valid(data.df, include.temporary.extraneous = TRUE),
                            .(n_include = .N), by = sp_key]
  sp_to_process_15 <- include_counts[n_include > 2, sp_key]

  # 15A: Pre-calculate plus/minus values for ALL valid rows
  data.df[, p_plus := NA_real_]
  data.df[, p_minus := NA_real_]
  data.df[param == "WEIGHTKG" & sp_key %in% sp_to_process_15, `:=`(p_plus = 1.05*v, p_minus = 0.95*v)]
  data.df[param == "HEIGHTCM" & sp_key %in% sp_to_process_15, `:=`(p_plus = v+1, p_minus = v-1)]
  data.df[param == "HEADCM" & sp_key %in% sp_to_process_15, `:=`(p_plus = v+1, p_minus = v-1)]

  # Debug save point for Step 13 output
  # Fixed to check global environment like step 11
#  if (exists("debug_step", envir = .GlobalEnv) && !is.null(get("debug_step", envir = .GlobalEnv)) && get("debug_step", envir = .GlobalEnv) == 13) {
#    saveRDS(data.df, paste0("R-step13-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step13-output.rds - STOPPING EXECUTION")
#    return(data.df[, .(id, exclude)])
#  }

  # 15B: Pre-calculate z-scores for p_plus/p_minus (once for all data)
  # This is expensive, so we do it once upfront
  valid_for_zscore <- data.df$sp_key %in% sp_to_process_15 & !is.na(data.df$p_plus)
  if (sum(valid_for_zscore) > 0) {
    zscore_subset <- data.df[valid_for_zscore]
    zscore_subset <- calc_and_recenter_z_scores(zscore_subset, "p_plus", ref.data.path)
    zscore_subset <- calc_and_recenter_z_scores(zscore_subset, "p_minus", ref.data.path)
    # Merge z-scores back
    data.df[zscore_subset, `:=`(tbc.p_plus = i.tbc.p_plus, tbc.p_minus = i.tbc.p_minus), on = "index"]
  }

  # Pre-calculate first_meas indicator
  # Fix first_meas logic
  # Handle HT/HC differently from WT
  # "Non-birth first" means: the first value, when that first value is not at birth
  # For WEIGHTKG: birth IS included in Step 15, so first_meas only if position 1 has agedays > 0
  # For HEIGHTCM/HEADCM: birth is EXCLUDED from Step 15, so calculate position among non-birth only
  data.df[, first_meas := FALSE]
  # WEIGHTKG: first_meas = TRUE only if position 1 AND agedays > 0
  data.df[sp_key %in% sp_to_process_15 & exclude == "Include" & param == "WEIGHTKG",
          first_meas := (seq_len(.N) == 1 & agedays > 0),
          by = .(subjid, param)]
  # HEIGHTCM/HEADCM: birth excluded from Step 15, so first among non-birth IS the first
  data.df[sp_key %in% sp_to_process_15 & exclude == "Include" & param %in% c("HEIGHTCM", "HEADCM") & agedays > 0,
          first_meas := (seq_len(.N) == 1),
          by = .(subjid, param)]

  # Step 15: EWMA2 global iterations (excluding agedays=0 for HT/HC)

  # Debug save point for Step 15 input
#  if (!is.null(debug_step) && debug_step == 15) {
#    saveRDS(data.df, paste0("R-step15-input-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step15-input.rds")
#  }

  # Get initial subject-params to process
  step15_filter <- data.df$sp_key %in% sp_to_process_15 &
    valid(data.df, include.temporary.extraneous = FALSE) &
    !((data.df$param == "HEIGHTCM" | data.df$param == "HEADCM") & data.df$agedays == 0)
  sp_counts_15 <- data.df[step15_filter, .(n = .N), by = sp_key]
  sp_to_process <- sp_counts_15[n > 2, sp_key]

  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      cat(sprintf("  EWMA2 iteration %d: %d subject-params\n", iteration, length(sp_to_process)))

    # **DEBUGGING CODE** - Added 2025-12-25-1100: Save Step 15 WT intermediate values (first iteration only)
    # Commented out debug_subjid block
#    if (iteration == 1) {
#      debug_subjid <- "4-2ff6188c-182c-6347-1692-7b188bd5be34"
#      debug_filter <- data.df$subjid == debug_subjid & data.df$param == "WEIGHTKG" &
#        valid(data.df, include.temporary.extraneous = FALSE)
#      if (sum(debug_filter) > 0) {
#        debug_df <- copy(data.df[debug_filter])
#        setorder(debug_df, agedays)
#        # Calculate EWMA values
#        ewma_result <- ewma(debug_df$agedays, debug_df$tbc.sd,
#                            fcase(pmax(abs(c(NA, diff(debug_df$agedays))), abs(c(diff(debug_df$agedays), NA)), na.rm=TRUE)/365.25 <= 1, -1.5,
#                                  pmax(abs(c(NA, diff(debug_df$agedays))), abs(c(diff(debug_df$agedays), NA)), na.rm=TRUE)/365.25 >= 3, -3.5,
#                                  default = -1.5 - (pmax(abs(c(NA, diff(debug_df$agedays))), abs(c(diff(debug_df$agedays), NA)), na.rm=TRUE)/365.25 - 1)),
#                            TRUE)
#        debug_df[, `:=`(ewma.all = ewma_result$ewma.all, ewma.before = ewma_result$ewma.before, ewma.after = ewma_result$ewma.after)]
#        debug_df[, `:=`(dewma.all = tbc.sd - ewma.all, dewma.before = tbc.sd - ewma.before, dewma.after = tbc.sd - ewma.after)]
#        debug_df[, `:=`(tbc_diff_next = tbc.sd - c(tbc.sd[2:.N], NA), tbc_diff_prior = tbc.sd - c(NA, tbc.sd[1:(.N-1)]))]
#        debug_df[, `:=`(tbc_diff_plus_next = tbc.p_plus - c(tbc.sd[2:.N], NA), tbc_diff_plus_prior = tbc.p_plus - c(NA, tbc.sd[1:(.N-1)]))]
#        debug_df[, `:=`(tbc_diff_minus_next = tbc.p_minus - c(tbc.sd[2:.N], NA), tbc_diff_minus_prior = tbc.p_minus - c(NA, tbc.sd[1:(.N-1)]))]
#        # Epsilon approach - debug code updated to match main code
#        debug_df[, addcrithigh := dewma.before > 0.99 & dewma.after > 0.99 &
#             ((tbc_diff_next > 0.99 & tbc_diff_plus_next > 0.99 & tbc_diff_minus_next > 0.99) | is.na(tbc_diff_next)) &
#             ((tbc_diff_prior > 0.99 & tbc_diff_plus_prior > 0.99 & tbc_diff_minus_prior > 0.99) | is.na(tbc_diff_prior))]
#        debug_df[, addcritlow := dewma.before < -0.99 & dewma.after < -0.99 &
#             ((tbc_diff_next < -0.99 & tbc_diff_plus_next < -0.99 & tbc_diff_minus_next < -0.99) | is.na(tbc_diff_next)) &
#             ((tbc_diff_prior < -0.99 & tbc_diff_plus_prior < -0.99 & tbc_diff_minus_prior < -0.99) | is.na(tbc_diff_prior))]
#        fwrite(debug_df, "R-step15-intermediate-wt.csv")
#        message("DEBUG: Saved R-step15-intermediate-wt.csv for subject ", debug_subjid)
#      }
#    }

    # RECALCULATE filter each iteration - valid() depends on current exclusions
    # Must recalculate filter inside loop
    step15_filter <- data.df$sp_key %in% sp_to_process &
      valid(data.df, include.temporary.extraneous = FALSE) &
      !((data.df$param == "HEIGHTCM" | data.df$param == "HEADCM") & data.df$agedays == 0)

    # Recalculate first_meas each iteration
    # Handle HT/HC differently from WT (same as initial calc)
    data.df[sp_key %in% sp_to_process, first_meas := FALSE]
    # WEIGHTKG: first_meas = TRUE only if position 1 AND agedays > 0
    data.df[sp_key %in% sp_to_process & exclude == "Include" & param == "WEIGHTKG",
            first_meas := (seq_len(.N) == 1 & agedays > 0),
            by = .(subjid, param)]
    # HEIGHTCM/HEADCM: birth excluded from Step 15, so first among non-birth IS the first
    data.df[sp_key %in% sp_to_process & exclude == "Include" & param %in% c("HEIGHTCM", "HEADCM") & agedays > 0,
            first_meas := (seq_len(.N) == 1),
            by = .(subjid, param)]

    # Track exclusions before this iteration
    ewma2_codes <- c("Exclude-EWMA2-middle", "Exclude-EWMA2-birth-WT", "Exclude-EWMA2-birth-WT-ext",
                     "Exclude-EWMA2-first", "Exclude-EWMA2-first-ext", "Exclude-EWMA2-last",
                     "Exclude-EWMA2-last-high", "Exclude-EWMA2-last-ext", "Exclude-EWMA2-last-ext-high")
    data.df[sp_key %in% sp_to_process, had_ewma2_before := (exclude %in% ewma2_codes)]

    # Process each subject-param - ONE exclusion per subject-param per iteration
    data.df[step15_filter,
            exclude := (function(df) {
              if (nrow(df) < 3) return(df$exclude)

              # Initialize EWMA fields
              df[, (ewma.fields) := as.double(NaN)]

              # Calculate exp_vals based on age gaps
              df[, c("diff_before", "diff_after") := .(c(NA, diff(agedays)), c(diff(agedays), NA))]
              df[, maxdiff := pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)]
              df[, ageyears := maxdiff / 365.25]
              df[, exp_vals := fcase(ageyears <= 1, -1.5, ageyears >= 3, -3.5, default = -1.5 - (ageyears - 1))]

              # Calculate EWMA
              df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
              df[, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]

              # Calculate dewma
              # No intermediate rounding - double-round at comparison points only
              df[, `:=`(dewma.all = tbc.sd - ewma.all,
                        dewma.before = tbc.sd - ewma.before,
                        dewma.after = tbc.sd - ewma.after,
                        c.dewma.all = ctbc.sd - c.ewma.all)]

              # Calculate prior/next differences
              # No intermediate rounding - double-round at comparison points only
              df[, `:=`(tbc_diff_next = tbc.sd - c(tbc.sd[2:.N], NA),
                        tbc_diff_prior = tbc.sd - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_plus_next = tbc.p_plus - c(tbc.sd[2:.N], NA),
                        tbc_diff_plus_prior = tbc.p_plus - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_minus_next = tbc.p_minus - c(tbc.sd[2:.N], NA),
                        tbc_diff_minus_prior = tbc.p_minus - c(NA, tbc.sd[1:(.N-1)]))]

              # Additional criteria
              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              df[, addcrithigh := janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) > 1 &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_next, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_next, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_next, 3), 2) > 1) | is.na(tbc_diff_next)) &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_prior, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_prior, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_prior, 3), 2) > 1) | is.na(tbc_diff_prior))]
              df[, addcritlow := janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) < -1 &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_next, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_next, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_next, 3), 2) < -1) | is.na(tbc_diff_next)) &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_prior, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_prior, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_prior, 3), 2) < -1) | is.na(tbc_diff_prior))]

              # DOP lookup
              compare_df <- data.df[subjid == df$subjid[1] & param == get_dop(df$param[1]) & exclude == "Include",]
              if (nrow(compare_df) > 0) {
                df[compare_df, tbc_dop := i.tbc.sd, on = "agedays"]
                df[is.na(tbc_dop), tbc_dop := median(compare_df$tbc.sd)]
              } else {
                df[, tbc_dop := NA]
              }

              df$rowind <- 1:nrow(df)
              n <- nrow(df)

              # All exclusion rules (mark candidates)
              df[, pot_excl := ""]

              # Middle
              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              df[janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 1 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 1 | is.na(c.dewma.all)) & addcrithigh & !(rowind %in% c(1, n)),
                 pot_excl := "Exclude-EWMA2-middle"]
              df[janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -1 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -1 | is.na(c.dewma.all)) & addcritlow & !(rowind %in% c(1, n)),
                 pot_excl := "Exclude-EWMA2-middle"]

              # Birth WT
              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              df[agedays == 0 & c(agedays[2:.N], NA) < 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-birth-WT"]
              df[agedays == 0 & c(agedays[2:.N], NA) < 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-birth-WT"]
              df[agedays == 0 & c(agedays[2:.N], NA) >= 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 4 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 4 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-birth-WT-ext"]
              df[agedays == 0 & c(agedays[2:.N], NA) >= 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -4 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -4 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-birth-WT-ext"]

              # First
              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              df[first_meas & (c(agedays[2:.N], NA) - agedays < 365.25) & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 2 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 2 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-first"]
              df[first_meas & (c(agedays[2:.N], NA) - agedays < 365.25) & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -2 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -2 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-first"]

              df[first_meas & (c(agedays[2:.N], NA) - agedays >= 365.25) & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-first-ext"]
              df[first_meas & (c(agedays[2:.N], NA) - agedays >= 365.25) & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-first-ext"]

              # Last
              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              if (n >= 2) {
                gap_last <- df$agedays[n] - df$agedays[n-1]
                tbc_prev <- df$tbc.sd[n-1]
                df[rowind == n & gap_last < 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) < 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 2 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 2 | is.na(c.dewma.all)) & addcrithigh,
                   pot_excl := "Exclude-EWMA2-last"]
                df[rowind == n & gap_last < 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) < 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -2 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -2 | is.na(c.dewma.all)) & addcritlow,
                   pot_excl := "Exclude-EWMA2-last"]
                df[rowind == n & gap_last < 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) >= 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) & addcrithigh,
                   pot_excl := "Exclude-EWMA2-last-high"]
                df[rowind == n & gap_last < 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) >= 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) & addcritlow,
                   pot_excl := "Exclude-EWMA2-last-high"]
                df[rowind == n & gap_last >= 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) < 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) &
                     (janitor::round_half_up(janitor::round_half_up(tbc.sd - tbc_dop, 3), 2) > 4 | is.na(tbc_dop)) & addcrithigh, pot_excl := "Exclude-EWMA2-last-ext"]
                df[rowind == n & gap_last >= 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) < 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) &
                     (janitor::round_half_up(janitor::round_half_up(tbc.sd - tbc_dop, 3), 2) < -4 | is.na(tbc_dop)) & addcritlow, pot_excl := "Exclude-EWMA2-last-ext"]
                df[rowind == n & gap_last >= 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) >= 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > janitor::round_half_up(janitor::round_half_up(1+abs(tbc_prev), 3), 2) & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) &
                     (janitor::round_half_up(janitor::round_half_up(tbc.sd - tbc_dop, 3), 2) > 4 | is.na(tbc_dop)) & addcrithigh, pot_excl := "Exclude-EWMA2-last-ext-high"]
                df[rowind == n & gap_last >= 365.25*2 & janitor::round_half_up(janitor::round_half_up(abs(tbc_prev), 3), 2) >= 2 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < janitor::round_half_up(janitor::round_half_up(-1-abs(tbc_prev), 3), 2) & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) &
                     (janitor::round_half_up(janitor::round_half_up(tbc.sd - tbc_dop, 3), 2) < -4 | is.na(tbc_dop)) & addcritlow, pot_excl := "Exclude-EWMA2-last-ext-high"]
              }

              # Exclude the worst candidate (highest abssum)
              # All pot_excl candidates are eligible - birth HT/HC is already excluded
              # from Step 15 processing (temporarily excluded before iteration starts)
              candidates <- df[pot_excl != ""]
              if (nrow(candidates) > 0) {
                candidates[, abssum := abs(tbc.sd + dewma.all)]
                # Use id tie-breaker for deterministic selection
                # which.max returns first tie, which depends on data order
                # Use order() with id as tie-breaker for reproducibility
                ord <- order(-candidates$abssum, candidates$id)
                worst_idx <- candidates$index[ord[1]]
                df[index == worst_idx, exclude := pot_excl]
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('index', 'id', 'subjid', 'param', 'agedays', 'v', 'sex', 'tbc.sd', 'ctbc.sd',
                        'tbc.p_plus', 'tbc.p_minus', 'first_meas', 'exclude')]

    # Find subject-params with NEW exclusions this iteration
    data.df[sp_key %in% sp_to_process, has_new_excl := (exclude %in% ewma2_codes) & !had_ewma2_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])

    # Next iteration: only subject-params with new exclusions
    sp_to_process <- sp_with_new_excl

    # Cleanup iteration tracking
    data.df[, c("had_ewma2_before", "has_new_excl") := NULL]
  }

  # Debug save point for Step 15 output
#  if (!is.null(debug_step) && debug_step == 15) {
#    saveRDS(data.df, paste0("R-step15-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step15-output.rds - STOPPING EXECUTION")
#    return(data.df[, .(line, exclude, tbc.sd, param)])
#  }

  # Step 16: Birth HT/HC ----

  if (!quietly)
    cat(sprintf("[%s] Exclude moderate EWMA for birth height and head circumference...\n", Sys.time()))

  # Filter for Step 16: HT/HC only, subjects with birth measurement, >2 values
  # Get subjects with birth measurements (HT/HC at agedays=0)
  subj_with_birth <- unique(data.df[param %in% c("HEIGHTCM", "HEADCM") &
                                     valid(data.df, include.temporary.extraneous = FALSE) &
                                     agedays == 0, subjid])

  # Initial filter for counting
  step16_filter <- data.df$param %in% c("HEIGHTCM", "HEADCM") &
    valid(data.df, include.temporary.extraneous = FALSE) &
    data.df$subjid %in% subj_with_birth
  sp_counts_16 <- data.df[step16_filter, .(n = .N), by = sp_key]
  sp_to_process <- sp_counts_16[n > 2, sp_key]

  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      cat(sprintf("  EWMA2-birth-HT-HC iteration %d: %d subject-params\n", iteration, length(sp_to_process)))

    # RECALCULATE filter each iteration - valid() depends on current exclusions
    # Must recalculate filter inside loop
    step16_filter <- data.df$sp_key %in% sp_to_process &
      data.df$param %in% c("HEIGHTCM", "HEADCM") &
      valid(data.df, include.temporary.extraneous = FALSE) &
      data.df$subjid %in% subj_with_birth

    ewma2_hthc_codes <- c("Exclude-EWMA2-birth-HT-HC", "Exclude-EWMA2-birth-HT-HC-ext")
    data.df[sp_key %in% sp_to_process, had_ewma2_before := (exclude %in% ewma2_hthc_codes)]

    data.df[step16_filter,
            exclude := (function(df) {
              if (nrow(df) < 3) return(df$exclude)

              df[, (ewma.fields) := as.double(NaN)]
              df[, c("diff_before", "diff_after") := .(c(NA, diff(agedays)), c(diff(agedays), NA))]
              df[, maxdiff := pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)]
              df[, ageyears := maxdiff / 365.25]
              df[, exp_vals := fcase(ageyears <= 1, -1.5, ageyears >= 3, -3.5, default = -1.5 - (ageyears - 1))]

              df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE)]
              df[, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE)]
              df[, `:=`(dewma.all = tbc.sd - ewma.all, dewma.before = tbc.sd - ewma.before,
                        dewma.after = tbc.sd - ewma.after, c.dewma.all = ctbc.sd - c.ewma.all)]

              df[, `:=`(tbc_diff_next = tbc.sd - c(tbc.sd[2:.N], NA),
                        tbc_diff_prior = tbc.sd - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_plus_next = tbc.p_plus - c(tbc.sd[2:.N], NA),
                        tbc_diff_plus_prior = tbc.p_plus - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_minus_next = tbc.p_minus - c(tbc.sd[2:.N], NA),
                        tbc_diff_minus_prior = tbc.p_minus - c(NA, tbc.sd[1:(.N-1)]))]

              # Reverted epsilon approach to rounding with round_half_up
              # Changed threshold comparisons from 3 to 2 decimal places
              # Double-round (3 then 2 decimals) to handle floating-point noise
              df[, addcrithigh := janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) > 1 &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_next, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_next, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_next, 3), 2) > 1) | is.na(tbc_diff_next)) &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_prior, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_prior, 3), 2) > 1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_prior, 3), 2) > 1) | is.na(tbc_diff_prior))]
              df[, addcritlow := janitor::round_half_up(janitor::round_half_up(dewma.before, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(dewma.after, 3), 2) < -1 &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_next, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_next, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_next, 3), 2) < -1) | is.na(tbc_diff_next)) &
                   ((janitor::round_half_up(janitor::round_half_up(tbc_diff_prior, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_plus_prior, 3), 2) < -1 & janitor::round_half_up(janitor::round_half_up(tbc_diff_minus_prior, 3), 2) < -1) | is.na(tbc_diff_prior))]

              df[, pot_excl := ""]
              next_age <- if (nrow(df) > 1) df$agedays[2] else NA

              df[agedays == 0 & !is.na(next_age) & next_age < 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-birth-HT-HC"]
              df[agedays == 0 & !is.na(next_age) & next_age < 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -3 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-birth-HT-HC"]
              df[agedays == 0 & !is.na(next_age) & next_age >= 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) > 4 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) > 4 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := "Exclude-EWMA2-birth-HT-HC-ext"]
              df[agedays == 0 & !is.na(next_age) & next_age >= 365.25 & janitor::round_half_up(janitor::round_half_up(dewma.all, 3), 2) < -4 & (janitor::round_half_up(janitor::round_half_up(c.dewma.all, 3), 2) < -4 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := "Exclude-EWMA2-birth-HT-HC-ext"]

              candidates <- df[pot_excl != ""]
              if (nrow(candidates) > 0) {
                candidates[, abssum := abs(tbc.sd + dewma.all)]
                # Use id tie-breaker for deterministic selection
                ord <- order(-candidates$abssum, candidates$id)
                worst_idx <- candidates$index[ord[1]]
                df[index == worst_idx, exclude := pot_excl]
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('index', 'id', 'subjid', 'param', 'agedays', 'v', 'sex', 'tbc.sd', 'ctbc.sd',
                        'tbc.p_plus', 'tbc.p_minus', 'exclude')]

    data.df[sp_key %in% sp_to_process, has_new_excl := (exclude %in% ewma2_hthc_codes) & !had_ewma2_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])
    sp_to_process <- sp_with_new_excl
    data.df[, c("had_ewma2_before", "has_new_excl") := NULL]
  }

  # Cleanup
  data.df[, sp_key := NULL]

  # Debug save point for Step 16 output
#  if (!is.null(debug_step) && debug_step == 16) {
#    saveRDS(data.df, paste0("R-step16-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step16-output.rds - STOPPING EXECUTION")
#    return(data.df[, .(line, exclude, tbc.sd, param)])
#  }

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
  
  ### CP ADD SECTION ###
  # ---- FIX: enforce Stata-style WHO binning ----
  # Convert ages to WHO monthly bins using floor() to match Stata
  # who.ht.vel[, whoagegrp_ht := as.integer(whoagegrp_ht)]
  # setkey(who.ht.vel, sex, whoagegrp_ht)
  
  # When merging into data.df later, derive the same floor-binned variable
  # data.df[, whoagegrp_ht := floor(agedays / 30.4375)]
  # setkey(data.df, sex, whoagegrp_ht)
  
  # Merge using roll=TRUE so values carry forward correctly
  # data.df <- who.ht.vel[data.df, roll = TRUE]
  # ---- END FIX ----
  ### CP ADD SECTION ###
  
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
  
  # order just for ease later
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, id),]

  # Compute valid_set AFTER sort, not before
  # valid_set is a logical vector indexed by row position, so it must be computed
  # on the sorted data, not before sorting (which would cause index misalignment)
  # create the valid set
  # we only run on valid values, non single values, and non weight
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 1]
  # Removed nnte_full filter to match Stata (NNTE appended before Step 17)
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    not_single &
    data.df$param != "WEIGHTKG" # do not run on weight

  data.df <- data.df[valid_set, exclude := (function(df) {
    # work inside a closure to drop added column values
    
    # save initial exclusions to keep track
    ind_all <- copy(df$index)
    id_all <- copy(df$id)
    
    exclude_all <- copy(df$exclude)
    
    testing <- TRUE
    iter_count <- 0
    
    while (testing & nrow(df) > 1){
      iter_count <- iter_count + 1
      
      # sort df since it got reordered with keys
      # Include id for deterministic SDE order
      df <- df[order(agedays, id),]
      
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
        
        ### CP ADJUST THIS BLOCK D
        
        # b
        df[d_agedays < 365.25, mindiff := .5*min.ht.vel*(d_agedays/365.25)^2-3 ]
        df[d_agedays > 365.25, mindiff := .5*min.ht.vel-3 ]
        # Fix maxdiff exponents to match Stata
        # Was: ^1.5 for <365.25, ^0.33 for >365.25 (swapped)
        # Fixed: ^0.33 for <365.25 (cube root), ^1.5 for >365.25
        df[d_agedays < 365.25,
           maxdiff := 2*max.ht.vel*(d_agedays/365.25)^0.33 + 5.5 ]
        df[d_agedays > 365.25,
           maxdiff := 2*max.ht.vel*(d_agedays/365.25)^1.5 + 5.5 ]
        
        #### CP ADJUST UP
        
        # 17D
        # generate the who age group variable
        # CP change to floor from round
        ### CP MODIFY THIS to conditional change
        # df[, whoagegrp.ht := round(agedays/30.4375)]
        df[agedays/30.4375 <= 24, whoagegrp.ht := round(agedays/30.4375)]
        ### CP MOdify up to conditional change
        
        # CP change to floor from round
        df[whoagegrp.ht > 24 | dplyr::lead(whoagegrp.ht) > 24,
           whoagegrp.ht := NA]
        
        
        # 17E
        df[d_agedays >= 20 & d_agedays < 46, whoinc.age.ht := 1]
        df[d_agedays >= 46 & d_agedays < 76, whoinc.age.ht := 2]
        df[d_agedays >= 76 & d_agedays < 107, whoinc.age.ht := 3]
        df[d_agedays >= 107 & d_agedays < 153, whoinc.age.ht := 4]
        df[d_agedays >= 153 & d_agedays < 199, whoinc.age.ht := 6]
        
        
        # Chris updated this to >= from ==
        # update the edge intervals
        df[d_agedays < 20, whoinc.age.ht := 1]
        df[d_agedays == 200, d_agedays := 200]
        df[d_agedays >= 200, whoinc.age.ht := 6]
        # Chris updated this to >= from ==
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
        
        # CP ADD IN OR EQUAL TO FROM < or > 
        # 17G
        df[d_agedays < whoinc.age.ht*30.4375, who_mindiff_ht :=
             who_mindiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        df[d_agedays > whoinc.age.ht*30.4375, who_maxdiff_ht :=
             who_maxdiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        # CP ADD IN OR EQUAL TO FROM < or > 
        df[d_agedays < 9*30.4375, who_mindiff_ht := who_mindiff_ht*.5-3]
        df[d_agedays < 9*30.4375, who_maxdiff_ht := who_maxdiff_ht*2+3]
        
        # 17H
        ### WAS an error CP below where it said whomindiff for the second on einstead of maxdiff##
        # tanner is implicit
        # For gap < 9 months: use WHO (already transformed at lines 4030-4031)
        df[d_agedays < 9*30.4375 & !is.na(who_mindiff_ht), mindiff := who_mindiff_ht]
        df[d_agedays < 9*30.4375 & !is.na(who_maxdiff_ht), maxdiff := who_maxdiff_ht]
        # Apply transformation for >= 9 month gaps
        # When gap >= 9 months, Tanner is preferred. But if Tanner not available (min.ht.vel is NA),
        # use WHO as fallback WITH the transformation (matches Stata lines 2601-2602)
        df[d_agedays >= 9*30.4375 & is.na(min.ht.vel) & !is.na(who_mindiff_ht),
           mindiff := who_mindiff_ht*.5-3]
        df[d_agedays >= 9*30.4375 & is.na(min.ht.vel) & !is.na(who_maxdiff_ht),
           maxdiff := who_maxdiff_ht*2+3]
        
        ### CP change to -3 from 3
        # otherwise, fill in
        df[is.na(mindiff), mindiff := -3]
        ### CP change up to -3 from 3
        # for birth measurements, add allowance of 1.5cm
        df[agedays == 0, mindiff := mindiff - 1.5]
        df[agedays == 0, maxdiff := maxdiff + 1.5]
        
        # 17I
        # sort df since it got reordered with keys
        # Include id for deterministic SDE order
        df <- df[order(agedays, id),]
        # Reverted epsilon approach to rounding with round_half_up
        # Changed from 3 to 2 decimals for threshold comparisons
        # Double-round (3 then 2 decimals) to handle floating-point noise
        df[, mindiff := janitor::round_half_up(janitor::round_half_up(mindiff, 3), 2)]
        df[, maxdiff := janitor::round_half_up(janitor::round_half_up(maxdiff, 3), 2)]
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
        # Include id for deterministic SDE order
        df <- df[order(agedays, id),]
        # Reverted epsilon approach to rounding with round_half_up
        # Changed from 3 to 2 decimals for threshold comparisons
        # Double-round (3 then 2 decimals) to handle floating-point noise
        df[, mindiff := janitor::round_half_up(janitor::round_half_up(mindiff, 3), 2)]
        df[, maxdiff := janitor::round_half_up(janitor::round_half_up(maxdiff, 3), 2)]
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
      # Use linear interpolation for exponent (Stata's linear formula)
      # Same approach as Step 11 EWMA1
      exp_vals <- rep(-1.5, nrow(tmp))
      maxdiff_years <- maxdiff_e / 365.25
      exp_vals[maxdiff_years > 1 & maxdiff_years < 3] <-
        -1.5 - (maxdiff_years[maxdiff_years > 1 & maxdiff_years < 3] - 1)
      exp_vals[maxdiff_years >= 3] <- -3.5
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
      # Reverted epsilon approach to rounding with round_half_up
      # Double-round (3 then 2 decimals) to match mindiff/maxdiff
      df[, diff_prev := janitor::round_half_up(janitor::round_half_up(v-dplyr::lag(v), 3), 2)]
      df[, diff_next := janitor::round_half_up(janitor::round_half_up(dplyr::lead(v)-v, 3), 2)]



      if (nrow(df) > 2){
        # 17P/R: identify pairs and calculate exclusions
        # Reverted epsilon approach to rounding with round_half_up
        df[, pair := diff_prev < mindiff_prev |
             diff_next < mindiff |
             diff_prev > maxdiff_prev |
             diff_next > maxdiff
        ]
        df[is.na(pair), pair := FALSE]
        # Reset tie-breaker flags each iteration (they persist from prior)
        # Double-round dewma comparisons (3 then 2 decimals)
        df[, bef.g.aftm1 := NA]
        df[, aft.g.aftm1 := NA]
        df[(pair & dplyr::lag(pair)) & janitor::round_half_up(janitor::round_half_up(abs(dewma.before), 3), 2) > dplyr::lag(janitor::round_half_up(janitor::round_half_up(abs(dewma.after), 3), 2)),
           bef.g.aftm1 := TRUE]
        df[(pair & dplyr::lead(pair)) & janitor::round_half_up(janitor::round_half_up(abs(dewma.after), 3), 2) > dplyr::lead(janitor::round_half_up(janitor::round_half_up(abs(dewma.before), 3), 2)),
           aft.g.aftm1 := TRUE]

        # Q
        # Reverted epsilon approach to rounding
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prev & bef.g.aftm1, val_excl := "Exclude-Min-diff"]
        df[diff_next < mindiff & aft.g.aftm1, val_excl := "Exclude-Min-diff"]
        df[diff_prev > maxdiff_prev & bef.g.aftm1, val_excl := "Exclude-Max-diff"]
        df[diff_next > maxdiff & aft.g.aftm1, val_excl := "Exclude-Max-diff"]
        df[diff_prev < mindiff_prev & bef.g.aftm1, val_excl_code := "1"]
        df[diff_next < mindiff & aft.g.aftm1, val_excl_code := "2"]
        df[diff_prev > maxdiff_prev & bef.g.aftm1, val_excl_code := "3"]
        df[diff_next > maxdiff & aft.g.aftm1, val_excl_code := "4"]
      }else { # only 2 values
        # 17Q/R -- exclusions for pairs
        # Reverted epsilon approach to rounding
        # Double-round tbc.sd comparisons (3 then 2 decimals)
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prev & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lag(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl := "Exclude-Min-diff"]
        df[diff_next < mindiff & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lead(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl := "Exclude-Min-diff"]
        df[diff_prev > maxdiff_prev & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lag(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl := "Exclude-Max-diff"]
        df[diff_next > maxdiff & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lead(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl := "Exclude-Max-diff"]
        df[diff_prev < mindiff_prev & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lag(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl_code := "5"]
        df[diff_next < mindiff & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lead(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl_code := "6"]
        df[diff_prev > maxdiff_prev & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lag(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl_code := "7"]
        df[diff_next > maxdiff & janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2) > dplyr::lead(janitor::round_half_up(janitor::round_half_up(abs(tbc.sd), 3), 2)),
           val_excl_code := "8"]
      } 
      
      
      
      
      # figure out if any of the exclusions hit
      count_exclude <- sum(df$val_excl != "Include")

      if (count_exclude > 0){
        if (nrow(df) > 2){
          # Fix Min-diff tie-breaking
          # Use the appropriate DEWMA based on exclusion type:
          # - When comparing to previous value (codes 1,3): use dewma.before (excludes previous)
          # - When comparing to next value (codes 2,4): use dewma.after (excludes next)
          df[, absval := NA_real_]
          df[val_excl_code %in% c("1", "3"), absval := abs(dewma.before)]
          df[val_excl_code %in% c("2", "4"), absval := abs(dewma.after)]
        } else {
          df[, absval := abs(tbc.sd)]
        }

        # choose the highest absval for exclusion
        # Use id tie-breaker for deterministic selection
        # which.max returns first tie, which depends on data order
        candidates <- df[val_excl != "Include"]
        ord <- order(-candidates$absval, candidates$id)
        idx <- candidates$index[ord[1]]
        
        
        exclude_all[ind_all == idx] <- df[index == idx, val_excl]
        if (unique(df$subjid) == "76234") {
          print(paste("idx:", idx))
          print(paste("index in df:", paste(df$index, collapse=", ")))
          print(paste("ind_all:", paste(ind_all, collapse=", ")))
          print(paste("id_all", id_all))
        }
        
        # if (unique(df$subjid) == "76234") {
        # print("AAAA")
        
        # Add loop ID column (optional, helps track iterations)
        # df[, loop_iter := iter_count]
        
        # Define output path
        # out_path <- "../data/mindiff_subjid_76234_id_20275.csv"
        
        # Append if file exists, otherwise write with header
        #   if (!file.exists(out_path)) {
        #     write.table(df, out_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        #   } else {
        #     write.table(df, out_path, sep = ",", row.names = FALSE, col.names = TRUE)
        #   }
        # }
        #set up to continue on
        # Fix iteration condition
        # Continue iterating when count_exclude >= 1 (not > 1)
        # After excluding 1 candidate and re-evaluating, NEW candidates may emerge
        if (count_exclude >= 1){
          testing <- TRUE

          df <- df[index != idx, ]
        } else {
          testing <- FALSE
        }
      } else {
        testing <- FALSE
      }
    }
    
    #       if (count_exclude > 0){
    #   if (nrow(df) > 2){
    #     # Use the appropriate DEWMA based on exclusion type
    #     df[, absval := NA_real_]
    #     df[val_excl == "Exclude-Min-diff" & diff_prev < mindiff_prev, 
    #        absval := abs(dewma.before)]
    #     df[val_excl == "Exclude-Min-diff" & diff_next < mindiff, 
    #        absval := abs(dewma.after)]
    #     df[val_excl == "Exclude-Max-diff" & diff_prev > maxdiff_prev, 
    #        absval := abs(dewma.before)]
    #     df[val_excl == "Exclude-Max-diff" & diff_next > maxdiff, 
    #        absval := abs(dewma.after)]
    #   } else {
    #     df[, absval := abs(tbc.sd)]
    #   }
    #   
    #   # Choose the highest absval for exclusion
    #   idx <- df$index[which.max(df[df$val_excl != "Include", absval])]
    #   exclude_all[ind_all == idx] <- df[index == idx, val_excl]
    #   
    #   # Set up to continue on
    #   if (count_exclude > 1){
    #     testing <- TRUE
    #     df <- df[index != idx, ]
    #   } else {
    #     testing <- FALSE
    #   }
    # } else {
    #   testing <- FALSE
    # }}
    
    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param), .SDcols = colnames(data.df)]


  # Debug save point for Step 17 output
#  if (!is.null(debug_step) && debug_step == 17) {
#    saveRDS(data.df, paste0("R-step17-output-", format(Sys.time(), "%Y-%m-%d-%H%M"), ".rds"))
#    message("DEBUG: Saved R-step17-output.rds - STOPPING EXECUTION")
#    return(data.df[, .(line, exclude, tbc.sd, param)])
#  }

  # 19: 1 or 2 measurements ----
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude 1 or 2 measurements...\n",
      Sys.time()
    ))
  
  # order just for ease later
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, id),]

  # Compute valid_set AFTER sort, not before
  # create the valid set
  # we only run on valid values and single/pair values
  # Count only Include rows for singles/pairs determination
  # Stata line 2769: gen vistot_`p'=_N if exc_`p'==0 - counts only remaining included measurements
  # Old code counted ALL rows, so subjects with 4 heights (2 excluded + 2 included) were skipped
  include_df <- data.df[exclude == "Include"]
  tmp <- table(paste0(include_df$subjid, "_", include_df$param))
  only_single_pairs <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp <= 2]
  # Removed nnte_full filter to match Stata (NNTE included in Step 19)
  valid_set <- valid(data.df, include.temporary.extraneous = FALSE) &
    only_single_pairs

  # Snapshot DOP data BEFORE by-group processing
  # Problem: data.df is modified in-place during by=(subjid,param) processing.
  # When HEIGHT is processed first and excluded, WEIGHT's DOP lookup won't find it.
  # Solution: Snapshot all Include rows with z-scores before processing.
  # This ensures consistent DOP lookups regardless of processing order.
  dop_snapshot <- data.df[exclude == "Include", .(subjid, param, agedays, tbc.sd, ctbc.sd)]
  setkey(dop_snapshot, subjid, param, agedays)

  data.df <- data.df[valid_set, exclude := (function(df) {
    # save initial exclusions to keep track
    ind_all <- copy(df$index)
    exclude_all <- copy(df$exclude)

    # 19A: is it a single or a pair?
    is_single <- nrow(df) == 1

    # find the DOP (designated other parameter) from pre-Step-19 snapshot
    # Use dop_snapshot instead of data.df for consistent results
    dop <- dop_snapshot[subjid == df$subjid[1] & param == get_dop(df$param[1])]
    
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
      # Use id tie-breaker for deterministic selection
      max_ind <- if (!all(is.na(df$comp_diff))){
        ord <- order(-df$comp_diff, df$id)
        ord[1]
      } else {
        ord <- order(-abs(df$tbc.sd), df$id)
        ord[1]
      }
      
      # 19F/G: which exclusion
      # Use absolute difference to match Stata
      # Stata line 2778: absd_tbc`p'z=abs(tbc`p'z-tbc`p'z_other)
      # Stata line 2830: absd_tbc`p'z>2.5 (uses absolute value)
      # Round for thresholds
      # Double-round to handle floating point noise from z-score subtraction
      if (janitor::round_half_up(janitor::round_half_up(abs(diff_tbc.sd), 3), 2) > 4 &
          (janitor::round_half_up(janitor::round_half_up(abs(diff_ctbc.sd), 3), 2) > 4 | is.na(diff_ctbc.sd)) &
          diff_agedays >=365.25){
        df[max_ind, exclude := "Exclude-2-meas->1-year"]
      } else if (janitor::round_half_up(janitor::round_half_up(abs(diff_tbc.sd), 3), 2) > 2.5 &
                 (janitor::round_half_up(janitor::round_half_up(abs(diff_ctbc.sd), 3), 2) > 2.5 | is.na(diff_ctbc.sd)) &
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
      # Check if 1-meas exclusion applies
      # Round for thresholds=3 and 5
      # Changed from 3 to 2 decimals, use round_half_up to match Stata
      # Double-round (3 then 2 decimals) to handle floating-point noise
      one_meas_cond <- (janitor::round_half_up(janitor::round_half_up(abs(df$tbc.sd), 3), 2) > 3 & !is.na(df$comp_diff) & janitor::round_half_up(janitor::round_half_up(df$comp_diff, 3), 2) > 5) |
                       (janitor::round_half_up(janitor::round_half_up(abs(df$tbc.sd), 3), 2) > 5 & is.na(df$comp_diff))
      if (one_meas_cond) {
        df[, exclude := "Exclude-1-meas"]
        # Only update exclude_all if 1-meas exclusion happens
        # Previous code at line 4529 unconditionally overwrote exclude_all, losing 2-meas results
        exclude_all[ind_all == df$index] <- "Exclude-1-meas"
      }
    }

    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param), .SDcols = colnames(data.df)]
  
  # 21: error load ----
  
  if (!quietly)
    cat(sprintf(
      "[%s] Exclude error load...\n",
      Sys.time()
    ))
  
  # Removed nnte_full filter to match Stata (NNTE included in Step 21)
  # Fix error-load denominator to exclude SDE/CF
  # Stata formula: tot_exc / (tot_exc + tot_inc) - excludes SDE/CF from denominator
  # R was using: errors / total - incorrectly included SDE/CF in denominator
  valid_set <- rep(TRUE, nrow(data.df))

  # Non-error codes that should be excluded from both numerator AND denominator
  non_error_codes <- c("Exclude-SDE-Identical",
                       "Exclude-SDE-All-Exclude",
                       "Exclude-SDE-All-Extreme",
                       "Exclude-SDE-EWMA",
                       "Exclude-SDE-One-Day",
                       "Exclude-Carried-Forward",
                       "Exclude-1-CF-deltaZ-<0.05",
                       "Exclude-1-CF-deltaZ-<0.1-wholehalfimp",
                       "Exclude-Teen-2-plus-CF-deltaZ-<0.05",
                       "Exclude-Teen-2-plus-CF-deltaZ-<0.1-wholehalfimp",
                       "Missing")

  data.df[valid_set,
          c("err_ratio", "n_errors") := {
            # Count errors (not Include and not non-error codes)
            n_errors <- sum(!exclude %in% c("Include", non_error_codes))
            # Count includes
            n_includes <- sum(exclude == "Include")
            # Denominator is errors + includes (excludes SDE/CF/Missing)
            denom <- n_errors + n_includes
            err_ratio <- if (denom > 0) n_errors / denom else 0
            list(err_ratio, n_errors)
          },
          by = c("subjid", "param")]

  # Add error.load.mincount check
  # Stata requires at least 2 errors before applying Error-load
  # Round for threshold=0.4
  # Changed from 3 to 2 decimals, use round_half_up to match Stata
  # Double-round (3 then 2 decimals) to handle floating-point noise
  data.df[valid_set & janitor::round_half_up(janitor::round_half_up(err_ratio, 3), 2) > .4 & n_errors >= error.load.mincount & exclude == "Include",
          exclude := "Exclude-Error-load"]

  # Cleanup
  data.df[, c("err_ratio", "n_errors") := NULL]
  
  
  # end ----
  
  if (!quietly)
    cat(sprintf("[%s] Completed Batch #%s...\n", Sys.time(), data.df$batch[1]))
  if (!quietly & parallel)
    sink()

  # Return appropriate z-score based on potcorr status
  # For potcorr subjects, return ctbc.sd (corrected); for others, return tbc.sd (uncorrected)
  # First check if we have potcorr and ctbc.sd columns
  has_potcorr <- "potcorr" %in% colnames(data.df)
  has_ctbc <- "ctbc.sd" %in% colnames(data.df)

  if (has_potcorr && has_ctbc) {
    # Merge potcorr from checkpoint_data if needed
    if (!"potcorr" %in% colnames(data.df) && exists("checkpoint_data") && "potcorr" %in% colnames(checkpoint_data)) {
      data.df <- merge(data.df, checkpoint_data[, .(id, potcorr)], by = "id", all.x = TRUE)
    }

    # Create final_tbc column: use ctbc.sd for potcorr subjects, tbc.sd for others
    data.df[, final_tbc := ifelse(!is.na(potcorr) & potcorr == TRUE & !is.na(ctbc.sd), ctbc.sd, tbc.sd)]
  } else {
    # Fallback to tbc.sd if columns not available
    data.df[, final_tbc := tbc.sd]
  }

  # Return z-scores and EWMA1 iteration 1 values for comparison
  # Build list of columns to return - start with essential columns
  return_cols <- c("id", "line", "exclude", "param")

  # Add z-score columns if they exist
  zscore_cols <- c("sd.orig_who", "sd.orig_cdc", "sd.orig", "tbc.sd", "ctbc.sd")
  for (col in zscore_cols) {
    if (col %in% names(data.df)) {
      return_cols <- c(return_cols, col)
    }
  }

  # Add EWMA1 iteration 1 columns if they exist (pattern: ewma1_it1.ewma_*, ewma1_it1.dewma_*)
  ewma1_it1_cols <- grep("^ewma1_it1\\.", names(data.df), value = TRUE)
  for (col in ewma1_it1_cols) {
    if (!(col %in% return_cols)) {
      return_cols <- c(return_cols, col)
    }
  }

  # Also add final_tbc for reference
  return_cols <- c(return_cols, "final_tbc")

  return(data.df[, ..return_cols])
}

# Oriignal Valid Toggled off
# Supporting pediatric growthcleanr functions
# Supporting functions for pediatric piece of algorithm

#' Helper function for cleanbatch to identify subset of observations that are either "included" or a "temporary extraneous"
#'
#' @keywords internal
#' @noRd
# valid <- function(df,
#                   include.temporary.extraneous = FALSE,
#                   include.extraneous = FALSE,
#                   include.carryforward = FALSE) {
#   exclude <- if (is.data.frame(df))
#     df$exclude
#   else
#     df
#   return(
#     exclude < 'Exclude'
#     |
#       include.temporary.extraneous &
#       exclude == 'Exclude-Temporary-Extraneous-Same-Day'
#     | include.extraneous & exclude == 'Exclude-Extraneous-Same-Day'
#     |
#       include.carryforward &
#       exclude == 'Exclude-Carried-Forward'
#   )
# }

valid <- function(df,
                  include.temporary.extraneous = FALSE,
                  include.extraneous = FALSE,
                  include.carryforward = FALSE) {
  exclude <- if (is.data.frame(df)) df$exclude else df
  exclude <- as.character(exclude)

  # Use regex for string comparison (prefix matching)
  # The string 'Include' sorts lexicographically after 'Exclude', so use grepl instead
  keep <- !grepl("^Exclude", exclude)
  
  if (include.temporary.extraneous)
    keep <- keep | exclude == "Exclude-Temporary-Extraneous-Same-Day"
  if (include.extraneous)
    keep <- keep | exclude == "Exclude-Extraneous-Same-Day"
  if (include.carryforward)
    keep <- keep | exclude == "Exclude-Carried-Forward"
  
  return(keep)
}

