################################################################################
# GROWTHCLEANR CHILD ALGORITHM - R IMPLEMENTATION
################################################################################
#
# PURPOSE:
#   Identifies and flags implausible pediatric growth measurements (weight,
#   height, head circumference) from electronic health records for children
#   ages 0-20 years.
#
# AUTHOR: Carrie Daymont, Penn State College of Medicine
#
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
# Preprocessing (handled in cleangrowth() before dispatch to cleanchild()):
#   - Input validation, internal_id assignment, imperial->metric conversion
#   - Z-score calculation (WHO/CDC blend via CSD method)
#   - Step 2b: Gestational age correction for potcorr (preterm) subjects
#   - Recentering -> tbc.sd and ctbc.sd
#   - Exclude-Missing and Exclude-Not-Cleaned assignment
#
# Algorithm steps (within cleanchild()):
#   Early 13: SDE-Identicals (pre-CF identical-value removal)
#   STEP 5:   Temporary same-day extraneous (SDE) flagging
#   STEP 6:   Carried-forward identification + optional rescue
#   STEP 7:   Biologically implausible values (BIV)
#   STEP 9:   Evil Twins (adjacent extreme values)
#   STEP 11:  EWMA1 - extreme EWMA outliers
#   STEP 13:  Final SDE resolution
#   STEP 15:  EWMA2 - moderate EWMA outliers
#   STEP 16:  Birth HT/HC EWMA2 variant
#   STEP 17:  Height/HC velocity checks
#   STEP 19:  Pairs and singles evaluation
#   STEP 21:  Error load assessment
#   STEP 22:  Output assembly
#
# Note: Step numbers are not consecutive. Steps 1-4, 8, 10, 12, 14, 18, 20 do
# not exist as top-level steps in the algorithm.
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
#   - DOP (Designated Other Parameter): For weight->height, for height->weight, for HC->height
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
# SORT ORDER:
#   Critical for deterministic results. Primary sort:
#     data.table::setkey(data.df, subjid, param, agedays, internal_id)
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
# This file contains the main cleanchild() function and supporting utilities:
#   - identify_temp_sde(): Temporary SDE resolution (Step 5)
#   - calc_otl_evil_twins(): Identifies Evil Twins (Step 9)
#   - calc_and_recenter_z_scores(): Recomputes recentered z-scores for the
#       p_plus / p_minus perturbation columns (Child Step 15/16 pre-loop)
#   - ewma() / ewma_cache_init() / ewma_cache_update(): EWMA machinery
#   - .child_valid(): Identifies included/partially-included rows for each step
#   - .child_exc(): Generates Exclude-C-{Reason} codes
#   - .cf_rescue_lookup(): Age/interval/param/rounding CF rescue thresholds
#   - get_dop(): Maps a parameter to its Designated Other Parameter
#   - read_anthro(): Builds z-score closures from reference tables
#   - gc_preload_refs(): Pre-loads reference closures for repeated calls
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
# and to order by row id. This ensures consistent handling of a dataset
# across multiple runs.
#
# Row IDENTIFIERS are discussed further above.
#


#' Generate a child-algorithm exclusion code.
#'
#' Convenience wrapper that returns `paste0("Exclude-C-", suffix)`. The
#' `param_val` argument is accepted but ignored — child exclusion codes
#' are no longer param-specific (the row's `param` column carries that
#' information) — and is retained only to avoid changing ~50 call sites.
#' The function is vectorised over `suffix` (and `param_val`), so it
#' works inside data.table `j` expressions.
#'
#' @param param_val ignored. Kept for call-site compatibility.
#' @param suffix reason code (e.g. `"BIV"`, `"CF"`, `"Traj"`). Must
#'   correspond to a level in `exclude.levels` or downstream factor
#'   assignment will silently produce NA.
#'
#' @return Character scalar or vector: `"Exclude-C-"` prepended to
#'   `suffix`. E.g. `.child_exc("WEIGHTKG", "BIV")` returns
#'   `"Exclude-C-BIV"`; the first arg's value is irrelevant.
#'
#' @keywords internal
#' @noRd
.child_exc <- function(param_val, suffix) {
  paste0("Exclude-C-", suffix)
}

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
#' @param biv.z.wt.low.young Lower unrecentered CSD z-score cutoff for
#' WEIGHTKG under one year of age (`ageyears < 1`). Weights with
#' `sd.orig_uncorr < biv.z.wt.low.young` are flagged as standardized BIV.
#' Defaults to -25.
#' @param biv.z.wt.low.old Lower unrecentered CSD z-score cutoff for
#' WEIGHTKG at or above one year (`ageyears >= 1`). Weights with
#' `sd.orig_uncorr < biv.z.wt.low.old` are flagged as standardized BIV.
#' Defaults to -15.
#' @param biv.z.wt.high Upper unrecentered CSD z-score cutoff for
#' WEIGHTKG (all ages). Weights with `sd.orig_uncorr > biv.z.wt.high` are
#' flagged as standardized BIV. Defaults to 22.
#' @param biv.z.ht.low.young Lower unrecentered CSD z-score cutoff for
#' HEIGHTCM under one year of age. Heights with
#' `sd.orig_uncorr < biv.z.ht.low.young` are flagged as standardized BIV.
#' Defaults to -25.
#' @param biv.z.ht.low.old Lower unrecentered CSD z-score cutoff for
#' HEIGHTCM at or above one year. Heights with
#' `sd.orig_uncorr < biv.z.ht.low.old` are flagged as standardized BIV.
#' Defaults to -15.
#' @param biv.z.ht.high Upper unrecentered CSD z-score cutoff for HEIGHTCM
#' (all ages). Heights with `sd.orig_uncorr > biv.z.ht.high` are flagged
#' as standardized BIV. Defaults to 8 (tighter than the weight upper limit,
#' based on analysis of CHOP data showing the ±15/±25 range was too loose
#' for heights).
#' @param biv.z.hc.low Lower unrecentered CSD z-score cutoff for HEADCM
#' (all ages). Head-circumference values with
#' `sd.orig_uncorr < biv.z.hc.low` are flagged as standardized BIV.
#' Defaults to -15.
#' @param biv.z.hc.high Upper unrecentered CSD z-score cutoff for HEADCM
#' (all ages). Head-circumference values with
#' `sd.orig_uncorr > biv.z.hc.high` are flagged as standardized BIV.
#' Defaults to 15.
#' @param error.load.mincount minimum count of exclusions on parameter before
#' considering excluding all measurements. Defaults to 2.
#' @param error.load.threshold threshold of percentage of excluded measurement count to included measurement
#' count that must be exceeded before excluding all measurements of either parameter. Defaults to 0.5.
#' @param sd.recenter Optional. A data.table with columns \code{param}, \code{sex},
#'   \code{agedays}, and \code{sd.median} to use for recentering instead of the
#'   built-in reference file (\code{rcfile-2023-08-15_format.csv.gz}). Defaults to
#'   NA, which uses the built-in file.
#' @param cf_rescue CF rescue mode for the child algorithm. One of:
#'   \describe{
#'     \item{"standard"}{(Default) Use age/interval/param-specific lookup thresholds.
#'       Identical-to-prior values are rescued (included) if |deltaZ| < threshold,
#'       excluded as CF if |deltaZ| >= threshold. NR cells get no rescue (all excluded).}
#'     \item{"none"}{No rescue. All identical-to-prior values are excluded as CF.}
#'     \item{"all"}{All identical-to-prior values are rescued (no CF exclusions).}
#'   }
#' @param cf_detail Logical. If TRUE, include \code{cf_status} and \code{cf_deltaZ}
#'   columns in the output. \code{cf_status} is one of: NA (not a CF candidate),
#'   "CF-NR" (CF, does not meet rescue criteria), or "CF-Resc" (CF, meets rescue
#'   criteria). \code{cf_deltaZ} is the absolute z-score difference used for rescue
#'   evaluation. Defaults to FALSE.
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
#' @param adult_permissiveness Permissiveness level for the adult algorithm:
#'   "loosest", "looser" (default), "tighter", or "tightest". Sets defaults
#'   for all adult BIV limits, EWMA caps, height bands, etc. See
#'   `adult_clean.R` for sub-parameter details.
#' @param adult_scale_max_lbs Physical scale upper limit in pounds for adult
#'   weight data. Weights at or above this value are excluded. Defaults to Inf.
#' @param ref_tables Optional. Pre-loaded reference closures from \code{\link{gc_preload_refs}}.
#'   When provided, skips all \code{read_anthro()} file reads (~0.93 sec per call on 200
#'   subjects; ~13 hours saved over 50K calls). Recommended for repeated calls (e.g.,
#'   simulation loops). Defaults to NULL (reads files each call).
#' @param cached_results Optional. A data.table of GC results from a prior \code{cleangrowth()}
#'   call on the same or similar dataset (e.g., a baseline run before error injection).
#'   Two usage modes:
#'   \itemize{
#'     \item \strong{Auto-detect mode} (\code{changed_subjids = NULL}): \code{cleangrowth()}
#'       compares the current input data against the cached results to automatically identify
#'       which subjects have any differences (added, removed, or modified rows). Only those
#'       subjects are re-processed; unchanged subjects receive their cached results. This is
#'       the recommended mode for error-injection workflows where most subjects are unmodified.
#'     \item \strong{Explicit mode} (\code{changed_subjids} provided): Only subjects in
#'       \code{changed_subjids} are re-processed. No auto-detection is performed. Use this
#'       when you already know which subjects changed (e.g., from an external comparison).
#'   }
#'   Defaults to NULL (process all subjects).
#' @param changed_subjids Optional. Vector of subject IDs whose measurements have changed
#'   since the cached baseline run. When provided with \code{cached_results}, only these
#'   subjects are re-processed. When NULL and \code{cached_results} is provided, changed
#'   subjects are auto-detected by comparing input data to cached results.
#'   Defaults to NULL.
#' @param tri_exclude Logical. If TRUE, include a \code{tri_exclude} column in the output
#'   that classifies each row as \code{"Include"}, \code{"Same-Day"} (same-day extraneous
#'   values), or \code{"Exclude"} (all other exclusions). Defaults to FALSE.
#' @param batch_size Number of subjects per processing batch. Subjects are never split
#'   across batches. Defaults to 2000.
#'
#' @return A data.table with columns: \code{id} (user-provided row identifier),
#'   \code{exclude} (exclusion code or \code{"Include"}), \code{param},
#'   \code{cf_rescued}, \code{sd.orig_who}, \code{sd.orig_cdc}, \code{sd.orig},
#'   \code{tbc.sd}, \code{ctbc.sd}, \code{final_tbc}. Adult rows additionally
#'   include \code{mean_ht} and \code{bin_result} (NA for child rows).
#'
#'   Child exclusion codes: \code{Exclude-C-CF}, \code{Exclude-C-BIV},
#'   \code{Exclude-C-Evil-Twins}, \code{Exclude-C-Traj-Extreme},
#'   \code{Exclude-C-Identical}, \code{Exclude-C-Extraneous},
#'   \code{Exclude-C-Traj}, \code{Exclude-C-Abs-Diff},
#'   \code{Exclude-C-Pair}, \code{Exclude-C-Single},
#'   \code{Exclude-C-Too-Many-Errors}, \code{Exclude-Missing},
#'   \code{Exclude-Not-Cleaned}.
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
                        biv.z.wt.low.young = -25,
                        biv.z.wt.low.old = -15,
                        biv.z.wt.high = 22,
                        biv.z.ht.low.young = -25,
                        biv.z.ht.low.old = -15,
                        biv.z.ht.high = 8,
                        biv.z.hc.low = -15,
                        biv.z.hc.high = 15,
                        error.load.mincount = 2,
                        error.load.threshold = 0.5,
                        sd.recenter = NA,
                        cf_rescue = "standard",
                        cf_detail = FALSE,
                        ref.data.path = "",
                        log.path = NA,
                        parallel = FALSE,
                        num.batches = NA,
                        quietly = TRUE,
                        adult_cutpoint = 20,
                        adult_permissiveness = "looser",
                        adult_scale_max_lbs = Inf,
                        ewma_window = 15,
                        id = NULL,
                        ref_tables = NULL,
                        cached_results = NULL,
                        changed_subjids = NULL,
                        tri_exclude = FALSE,
                        batch_size = 2000)
                        {
  # ewma_window: number of neighbors on each side for EWMA weighting.

  # Validate cf_rescue
  cf_rescue <- match.arg(cf_rescue, c("standard", "none", "all"))

  # avoid "no visible binding" warnings
  N <- age_years <- batch <- exclude <- index <- line <- NULL
  newbatch <- sd.median <- sd.orig <- tanner.months <- tbc.sd <- NULL
  v <- v_adult <- whoagegrp.ht <- whoagegrp_ht <- z.orig <- NULL
  z.orig_cdc <- z.orig_who <- sd.orig_cdc <- sd.orig_who <- NULL
  result <- NULL

  sd.orig_uncorr <- agemonths <- ageyears_2b <- intwt <- fengadays <- pmagedays <- cagedays <-
    unmod_zscore <- fen_param <- M <- S_upper <- S_lower <- cwho_cv <- ccdc_cv <-
    sd.c_cdc <- sd.c_who <- sd.c <- sd.corr <- seq_win <- sd.corr_abssumdiff <-
    sd.orig_abssumdiff <- orig_colnames <- ctbc.sd <- sum_sde <- no_sde <-
    NULL

  # preprocessing ----

  # Partial-run mode: filter input to changed subjects only.
  # This must happen before data.all.ages is built.
  #
  # Two modes:
  #   1. Auto-detect: cached_results provided, changed_subjids NULL.
  #      Compare input data to cached_results to find subjects with any
  #      differences (added/removed/modified rows). Only those subjects
  #      are re-processed.
  #   2. Explicit: both cached_results and changed_subjids provided.
  #      Only subjects in changed_subjids are re-processed.
  do_partial <- !is.null(cached_results)
  if (do_partial && is.null(changed_subjids)) {
    # Auto-detect changed subjects by comparing input to cached results.
    # Build a minimal data.table from input vectors for comparison.
    # cached_results contains subjid, param, agedays, sex, v (measurement
    # with 0→NaN applied). Compare on the original measurement values
    # before that transformation, using the same columns cached_results has.
    input_dt <- data.table(
      subjid = as.factor(subjid),
      param  = param,
      agedays = as.integer(agedays),
      sex    = as.integer(ifelse(
        sex %in% c(0, "m", "M"), 0L, ifelse(sex %in% c(1, "f", "F"), 1L, NA_integer_)
      )),
      v      = ifelse(measurement == 0, NaN, measurement)
    )
    # Sort both datasets identically for per-subject comparison
    compare_cols <- c("subjid", "param", "agedays", "sex", "v")
    setkeyv(input_dt, compare_cols)
    cached_compare <- cached_results[, ..compare_cols]
    setkeyv(cached_compare, compare_cols)

    # Subjects only in input (added) or only in cache (removed)
    input_subj  <- unique(as.character(input_dt$subjid))
    cached_subj <- unique(as.character(cached_results$subjid))
    added   <- setdiff(input_subj, cached_subj)
    removed <- setdiff(cached_subj, input_subj)

    # For common subjects: compare row-by-row after sorting
    common <- intersect(input_subj, cached_subj)
    modified <- character(0)
    if (length(common) > 0) {
      # Per-subject comparison: concatenate all values into a single string
      # and compare. This is fast and handles row count differences.
      input_hash <- input_dt[as.character(subjid) %in% common,
                             .(hash = paste(param, agedays, sex, v,
                                            sep = "|", collapse = "\n")),
                             by = subjid]
      cached_hash <- cached_compare[as.character(subjid) %in% common,
                                    .(hash = paste(param, agedays, sex, v,
                                                   sep = "|", collapse = "\n")),
                                    by = subjid]
      merged <- merge(input_hash, cached_hash, by = "subjid",
                      suffixes = c("_new", "_old"))
      modified <- as.character(merged[hash_new != hash_old, subjid])
    }

    changed_subjids <- c(added, modified)
    # Note: removed subjects are not in input data, so they simply won't
    # appear in output. No special handling needed.

    if (!quietly) {
      message(sprintf(
        "[%s] Auto-detected %d changed subjects (%d added, %d modified, %d removed, %d unchanged)",
        Sys.time(), length(changed_subjids), length(added), length(modified),
        length(removed), length(common) - length(modified)
      ))
    }
  }
  if (do_partial) {
    if (length(changed_subjids) == 0) return(cached_results)
    keep <- as.character(subjid) %in% as.character(changed_subjids)
    if (sum(keep) == 0) return(cached_results)  # nothing changed — return cache as-is
    subjid      <- subjid[keep]
    param       <- param[keep]
    agedays     <- agedays[keep]
    sex         <- sex[keep]
    measurement <- measurement[keep]
    id          <- id[keep]
  }

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

  # Keep user's id as 'id', create 'internal_id' for R processing.
  # internal_id is a sequential integer assigned AFTER sorting by the canonical key
  # (subjid, param, agedays, id). This ensures internal_id reflects id-sorted order,
  # so results are deterministic regardless of input row order.
  # id is required in input data.
  if (length(id) != length(measurement)) {
    stop("id must be provided and have the same length as measurement")
  }
  data.all.ages[, id := id]
  setkey(data.all.ages, subjid, param, agedays, id)
  data.all.ages[, internal_id := seq_len(.N)]

  # quality checks
  if (!is.numeric(adult_cutpoint)){
    stop("adult_cutpoint not numeric. Please enter a number between 18 and 20.")
  }
  if (!is.numeric(adult_scale_max_lbs) || adult_scale_max_lbs < 0){
    stop("adult_scale_max_lbs not numeric. Please enter a positive number.")
  }
  if (any(!param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHTIN",
                        "WEIGHTLBS", "HEADCM"))){
    message(sprintf("[%s] Parameters included that do not match 'param' specifications. Marking as missing...", Sys.time()))
    data.all.ages <-
      data.all.ages[
        !param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHTIN",
                      "WEIGHTLBS", "HEADCM"), v := NA]
    data.all.ages <-
      data.all.ages[
        !param %in% c("LENGTHCM", "HEIGHTCM", "WEIGHTKG", "HEIGHTIN",
                      "WEIGHTLBS", "HEADCM"), v_adult := NA]
  }

  # rate limit cutpoint -- min 18, max 20
  if (adult_cutpoint < 18) {
    adult_cutpoint <- 18
  } else if (adult_cutpoint > 20) {
    adult_cutpoint <- 20
  }
  ### BATCHING ###
  if (!quietly) message("Batching subjects...")
  results_list <- list()
  i <- 1
  patients <- unique(data.all.ages[, .(subjid)])
  patients[, batch := (seq_len(.N) - 1L) %/% batch_size + 1L]
  batches <- patients

  # --- Batch-invariant constants (computed once, not per batch) ---

  # .child_exc() is defined at file level (before cleangrowth) so it's
  # visible in data.table j expressions within cleanchild()

  # enumerate the different exclusion levels
  exclude.levels.peds <- c(
    'Include',
    'Exclude-Missing',
    'Exclude-Not-Cleaned',
    'Exclude-C-Temp-Same-Day',
    # child exclusion codes (not param-specific; param is in the data)
    'Exclude-C-CF',
    'Exclude-C-Traj-Extreme',
    'Exclude-C-Identical',
    'Exclude-C-Extraneous',
    'Exclude-C-Traj',
    'Exclude-C-Abs-Diff',
    'Exclude-C-Pair',
    'Exclude-C-Single',
    'Exclude-C-Too-Many-Errors',
    'Exclude-C-BIV',
    'Exclude-C-Evil-Twins'
  )

  # Adult exclusion codes used by cleanadult()
  # Not param-specific; param is in the data
  exclude.levels.adult <- c(
    "Exclude-A-BIV",
    "Exclude-A-Scale-Max",
    "Exclude-A-Scale-Max-Identical",
    "Exclude-A-Scale-Max-RV-Propagated",
    "Exclude-A-Evil-Twins",
    "Exclude-A-Identical",
    "Exclude-A-Extraneous",
    "Exclude-A-Traj-Extreme",
    "Exclude-A-Traj-Extreme-firstRV",
    "Exclude-A-Traj-Extreme-allRV",
    "Exclude-A-Traj-Extreme-firstRV-RV-Propagated",
    "Exclude-A-Ord-Pair",
    "Exclude-A-Ord-Pair-All",
    "Exclude-A-Window",
    "Exclude-A-Window-All",
    "Exclude-A-2D-Ordered",
    "Exclude-A-2D-Non-Ordered",
    "Exclude-A-Traj-Moderate",
    "Exclude-A-Traj-Moderate-allRV",
    "Exclude-A-Traj-Moderate-RV-Propagated",
    "Exclude-A-Traj-Moderate-Error-Load",
    "Exclude-A-Traj-Moderate-Error-Load-RV",
    "Exclude-A-Single",
    "Exclude-A-Too-Many-Errors"
  )
  exclude.levels <- c(exclude.levels.peds, exclude.levels.adult)

  # Load velocity reference tables (same for all batches, read once here)
  # These are only used by the pediatric algorithm (Step 17). For all-adult
  # datasets they are loaded but unused — negligible cost (~0.01 sec).
  # Tanner height velocity (used by Step 17 for ages 2.5+ years)
  tanner_ht_vel_path <- ifelse(
    ref.data.path == "",
    system.file(file.path("extdata", "tanner_ht_vel.csv.gz"), package = "growthcleanr"),
    file.path(ref.data.path, "tanner_ht_vel.csv.gz")
  )
  tanner.ht.vel <- fread(tanner_ht_vel_path)
  setnames(tanner.ht.vel, colnames(tanner.ht.vel), gsub('_', '.', colnames(tanner.ht.vel)))
  setkey(tanner.ht.vel, sex, tanner.months)

  # WHO velocity tables are loaded inside cleanchild() (HT and HC separately).

  # --- Parallel setup (once, before outer batch loop) ---
  # Create a single cluster for both child and adult dispatch across all
  # outer batches. This avoids the overhead of creating/destroying a cluster
  # per batch (~2-5 sec each). All internal functions needed by both
  # cleanchild() and cleanadult() are exported to workers.
  if (parallel) {
    if (is.na(num.batches)) {
      num.batches <- getDoParWorkers()
    }
    var_for_par <- c(
      # Child algorithm functions
      ".child_valid",
      "ewma", "read_anthro", "as_matrix_delta",
      "identify_temp_sde",
      "get_dop", "calc_otl_evil_twins",
      "calc_and_recenter_z_scores",
      ".child_exc", ".cf_rescue_lookup",
      "ewma_cache_init", "ewma_cache_update",
      # Adult algorithm functions
      "cleanadult",
      "permissiveness_presets", "resolve_permissiveness",
      "check_between", "round_pt",
      "compute_et_limit", "compute_perc_limit", "compute_wtallow",
      "as.matrix.delta_dn", "ewma_dn",
      "adult_ewma_cache_init", "adult_ewma_cache_update",
      "remove_biv", "remove_biv_high", "remove_biv_low",
      "identify_rv", "temp_sde", "redo_identify_rv",
      "ht_allow", "ht_change_groups", "ht_3d_growth_compare",
      "detect_runs", "compute_trajectory_fails",
      "evil_twins", "propagate_to_rv",
      "remove_ewma_wt", "remove_mod_ewma_wt",
      "eval_2d_nonord", "eval_1d", "eval_error_load"
    )
    cl <- makeCluster(num.batches)
    clusterExport(cl = cl, varlist = var_for_par, envir = environment())
    registerDoParallel(cl)
  } else {
    if (is.na(num.batches)) num.batches <- 1
  }

  ### Loop Begin ###

  for (id_batch in unique(batches$batch)) {
    if (!quietly) message(sprintf("Batch %d start.", i))
    ids <- batches$subjid[batches$batch == id_batch]

    data.all <- copy(
      data.all.ages[
        subjid %in% ids & agedays < adult_cutpoint * 365.25
      ]
    )

    data.adult <- copy(
      data.all.ages[
        subjid %in% ids & agedays >= adult_cutpoint * 365.25
      ]
    )
  # Batch-level reference for result assembly (all ages for this batch's subjects)
  data.batch <- copy(data.all.ages[subjid %in% ids])

  # if there's no pediatric data, no need to go through this rigamarole
  if (nrow(data.all) > 0){

    # pediatric: height velocity calculations and preprocessing ----

    # for pediatric data, convert in and lbs to cm and kg (adult is done within algo)
    data.all[param == "HEIGHTIN", v := v*2.54]
    data.all[param == "HEIGHTIN", param := "HEIGHTCM"]
    data.all[param == "WEIGHTLBS", v := v/2.2046226]
    data.all[param == "WEIGHTLBS", param := "WEIGHTKG"]

    setkey(data.all, subjid)
    subjid.unique <- data.all[j = unique(subjid)]
    batches.all <- data.table(
      subjid = subjid.unique,
      batch = (seq_along(subjid.unique) - 1L) %% num.batches + 1L,
      key = 'subjid'
    )
    data.all <- batches.all[data.all]

    if (!quietly){
      message(sprintf("[%s] Begin processing pediatric data...", Sys.time()))
    }

    # Velocity reference tables loaded once before the batch loop (batch-invariant)

    # Algorithm design principles:
    # - All steps operate separately by parameter unless otherwise noted.
    # - Data sorted by subjid, param, agedays, internal_id; only non-excluded,
    #   non-missing rows participate (controlled by .child_valid()).
    # - "Next" and "previous" refer to adjacent included values for the same
    #   subject and parameter by age.
    # - Exclusions are recorded in the `exclude` factor column. Once a row is
    #   excluded it stays excluded unless explicitly rescued (e.g., CF rescue).
    # - CSD z-scores (not LMS z-scores) form the basis for the algorithm.
    #   CSD is more sensitive to extreme high values than standard LMS z-scores.

    # recategorize linear parameters as 'HEIGHTCM'
    # NOTE: this will be changed in future to consider this difference
    data.all[param == 'LENGTHCM', param := 'HEIGHTCM']

    # calculate z/sd scores
    if (!quietly)
      message(sprintf("[%s] Calculating z-scores...", Sys.time()))
      # removing z calculations, as they are not used
      # for infants, use z and who
      measurement.to.z <- if (!is.null(ref_tables)) ref_tables$mtz_cdc_prelim else
        read_anthro(ref.data.path, cdc.only = TRUE)
      measurement.to.z_who <- if (!is.null(ref_tables)) ref_tables$mtz_who_prelim else
        read_anthro(ref.data.path, cdc.only = FALSE)

      # calculate "standard deviation" scores
      if (!quietly)
        message(sprintf("[%s] Calculating SD-scores...", Sys.time()))
      data.all[, sd.orig_cdc := measurement.to.z(param, agedays, sex, v, TRUE)]
      # Note: sd.orig_who contains CDC values for ages >= 5y (the WHO closure
      # falls back to CDC when WHO data is unavailable). The name is historical.
      # The blending logic below gives zero weight to WHO at age 5+, so these
      # CDC values in sd.orig_who are never incorrectly used as WHO.
      data.all[, sd.orig_who := measurement.to.z_who(param, agedays, sex, v, TRUE)]

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


      # keep the original, uncorrected, unrecentered zscores
      data.all[,sd.orig_uncorr := sd.orig]

      # 2b: corrected z scores ----

      # keep the original column names -- we're adding a ton of columns that we
      # want to filter out after correction
      orig_colnames <- copy(colnames(data.all))

      # start by reading in fenton data only if potcorr subjects exist (see below)

      # add age in months and years (ageyears_2b used throughout Step 2b)
      data.all[, agemonths := agedays / 30.4375]
      data.all[, ageyears_2b := agedays / 365.25]

      # Sort by id, age-dependent direction
      # Age-dependent ID selection for SDEs:
      # At age 0: prefer LOWEST id (earliest measurement, before fluid/interventions)
      # At age > 0: prefer HIGHEST id (consistent with other SDE handling)
      # Create sort key: for age 0 use id ascending, for age > 0 use id descending
      data.all[, id_sort := ifelse(agedays == 0, internal_id, -internal_id)]
      setorder(data.all, subjid, param, agedays, id_sort)

      # initialize the column so it's always present
      # Create row-level potcorr_wt flag
      data.all[, potcorr := FALSE]
      data.all[, potcorr_wt := FALSE]

      # Mark infants needing Fenton correction
      # potcorr_wt: row-level flag for FIRST weight per subject (earliest agedays, not necessarily birth)
      #             that qualifies: Z<-2 AND age<10 months
      # After sort by agedays + id_sort, seq_len(.N)==1L selects first weight with age-dependent ID
      # NA-safe: sd.orig is NA when measurement is NA (Missing rows)
      data.all[param == "WEIGHTKG",
               potcorr_wt := (seq_len(.N) == 1L &
                              !is.na(sd.orig) & sd.orig < -2 &
                              agemonths < 10),
               by = subjid]

      # Clean up sort key
      data.all[, id_sort := NULL]

      # propagate that flag to all rows for that subject (subject-level)
      data.all[, potcorr := any(potcorr_wt, na.rm = TRUE), by = subjid]

      # ensure deterministic order
      setkey(data.all, subjid, param, agedays, internal_id)

      # --- Potcorr optimization: skip all reference merges if no potcorr subjects ---
      has_potcorr <- any(data.all$potcorr, na.rm = TRUE)

      if (has_potcorr) {
        # Read Fenton reference data (only when potcorr subjects exist)
        fent_foraga <- fread(
          system.file(file.path("extdata", "fent_foraga.csv.gz"),
                      package = "growthcleanr"))
        # Fenton 2025 reference: M, S_upper, S_lower from extracted curves (CSD method)
        fenton_ref <- fread(
          system.file(file.path("extdata", "fenton2025_ms_lookup_smoothed.csv"),
                      package = "growthcleanr"))

        # Subset to potcorr subject rows only (typically <1% of data)
        pc_ids <- unique(data.all[potcorr == TRUE, subjid])
        pc <- copy(data.all[subjid %in% pc_ids])

        # integer weight is in grams, rounded to the nearest 10
        pc[potcorr == TRUE, intwt := trunc(v*100)*10]
        # Floor weights in [100, 500) g up to the Fenton table minimum of 500 g,
        # so users who configure lower BIV thresholds can still get
        # Fenton-corrected z-scores for very low weights.
        pc[intwt >= 100 & intwt < 500, intwt := 500]

        # Fenton merge 1: weight -> estimated gestational age
        pc <- merge(
          pc, fent_foraga, by = c("sex", "intwt"),
          all.x = TRUE)

        # Propagate fengadays ONLY from potcorr_wt rows
        # Only use fengadays from the qualifying first weight (potcorr_wt==TRUE), not from any row
        pc[, fengadays_subj := ifelse(any(!is.na(fengadays) & potcorr_wt),
                                             min(fengadays[potcorr_wt], na.rm = TRUE),
                                             NA_real_),
                 by = subjid]
        # Use subject-level fengadays for all calculations
        pc[!is.na(fengadays_subj), fengadays := fengadays_subj]
        pc[, fengadays_subj := NULL]  # Clean up temporary variable

        pc[fengadays < 259, pmagedays := agedays + fengadays]
        pc[fengadays >= 259 | is.na(fengadays), pmagedays := NA_real_]

        pc[fengadays < 259, cagedays := pmagedays - 280]
        # replace fengadays with pmagedays to facilitate merging
        pc[, fengadays := pmagedays]

        # Fenton merge 2: gestational age -> Fenton 2025 reference (M, S_upper, S_lower)
        # Map gc param names to Fenton reference param names for long-format merge
        pc[, fen_param := fcase(
          param == "WEIGHTKG", "weight",
          param == "HEIGHTCM", "length",
          param == "HEADCM",   "headcirc"
        )]
        pc <- merge(
          pc, fenton_ref[, .(sex, ga_days, param, M, S_upper, S_lower)],
          by.x = c("sex", "fengadays", "fen_param"),
          by.y = c("sex", "ga_days", "param"),
          all.x = TRUE)

        # Reset potcorr_wt when Fenton merge fails
        # Only check weight rows (potcorr_wt is only TRUE on WEIGHTKG rows)
        pc[potcorr_wt == TRUE & is.na(M), potcorr_wt := FALSE]
        # Recalculate subject-level potcorr
        pc[, potcorr := any(potcorr_wt), by = subjid]

        # --- Fenton 2025 z-score calculation (CSD method) ---
        # Uses S_upper and S_lower extracted from plotted Fenton 2025 curves.
        # CSD: z = (v - M) / (S * M), with S_upper when v >= M, S_lower when v < M.
        # Units: weight M is in grams, length/HC M is in cm.
        pc[, v_fenton := fifelse(param == "WEIGHTKG", v * 1000, v)]

        # CSD z-score: S_upper/S_lower are proportional SDs (CV-like),
        # so absolute SD = S * M
        pc[, unmod_zscore := fifelse(
          v_fenton >= M,
          (v_fenton - M) / (S_upper * M),
          (v_fenton - M) / (S_lower * M)
        )]

        # Clean up temporary variables
        pc[, c("v_fenton", "fen_param", "M", "S_upper", "S_lower") := NULL]

        # Assign Fenton z-score as initial sd.corr for potcorr subjects
        pc[potcorr == TRUE & !is.na(unmod_zscore),
                 sd.corr := unmod_zscore]

        # For everyone else (non-potcorr rows of potcorr subjects), fall back to original z
        pc[is.na(sd.corr), sd.corr := sd.orig]

        # --- Corrected WHO/CDC z-scores via closures (no merge needed) ---
        # The measurement.to.z and measurement.to.z_who closures contain their own
        # internal reference data (loaded in read_anthro). The previous WHO/CDC merges
        # at this point added 40+ columns to data.all that were never referenced —
        # the closures compute z-scores independently. Merges removed for efficiency.
        pc[, cwho_cv := v]
        pc[, ccdc_cv := v]
        pc[param == "HEIGHTCM" & agedays > 730 & cagedays <= 730,
                 cwho_cv := cwho_cv + 0.8]
        pc[param == "HEIGHTCM" & agedays > 730 & cagedays <= 730,
                 ccdc_cv := ccdc_cv + 0.7]

        # create the corrected z scores
        pc[, sd.c_cdc :=
                   measurement.to.z(param, cagedays, sex, ccdc_cv, TRUE)]
        pc[, sd.c_who :=
                   measurement.to.z_who(param, cagedays, sex, cwho_cv, TRUE)]

        # Use WHO corrected z-score if pmagedays >= 350 and chronological age <= 730 days (~2 years)
        pc[potcorr == TRUE & pmagedays >= 350 & agedays <= 730 & !is.na(sd.c_who),
                 sd.corr := sd.c_who]

        # Assign sd.c (corrected WHO/CDC z-score) in one pass
        # HC uses WHO at all ages; WT/HT blend WHO->CDC between ages 2-5
        pc[, sd.c := fcase(
          param == "HEADCM" & !is.na(sd.c_who),
            sd.c_who,
          param != "HEADCM" & ageyears_2b <= 2 & !is.na(sd.c_who),
            sd.c_who,
          param != "HEADCM" & ageyears_2b > 2 & ageyears_2b < 5 &
            !is.na(sd.c_who) & !is.na(sd.c_cdc),
            (sd.c_who * (5 - ageyears_2b) + sd.c_cdc * (ageyears_2b - 2)) / 3,
          param != "HEADCM" & ageyears_2b >= 5 & !is.na(sd.c_cdc),
            sd.c_cdc
        )]

        # Prefer Fenton over corrected WHO for ages <=2. Only overwrite sd.c
        # when Fenton (unmod_zscore) is available; the := is a no-op otherwise,
        # which preserves the existing corrected WHO/CDC value as the fallback.
        pc[potcorr == TRUE & ageyears_2b <= 2 & !is.na(unmod_zscore),
                 sd.c := unmod_zscore]

        # Assign final sd.corr in one pass
        # Priority: <2y potcorr -> corrected; 2-4y -> smoothed blend; else -> original
        pc[, sd.corr := fcase(
          ageyears_2b <= 2 & potcorr & !is.na(sd.c),
            sd.c,
          ageyears_2b > 2 & ageyears_2b <= 4 & !is.na(sd.c) & !is.na(sd.orig),
            (sd.orig * (4 - ageyears_2b) + sd.c * (ageyears_2b - 2)) / 2,
          default = sd.orig
        )]

        # --- Evaluate whether correction makes z-scores more extreme (uncorr) ---
        examine_only <- pc$param == "WEIGHTKG" & pc$potcorr
        tmp <- copy(pc[examine_only,])

        # Filter SDE-Identicals BEFORE uncorr calculation (pre-Step-13 cleanup).
        # At birth (age 0): Keep LOWEST id (earliest, before fluid/interventions)
        # At age > 0: Keep HIGHEST id (consistent with other SDE handling)

        # Identify same-day identical values
        tmp[, `:=`(
          n_on_day = .N,
          n_unique_vals = uniqueN(v)
        ), by = .(subjid, agedays)]

        tmp[, all_identical := n_on_day > 1 & n_unique_vals == 1]

        # Age-dependent ID selection for identicals
        tmp[, keep_id := {
          ids <- internal_id
          if (agedays[1L] == 0L) min(ids) else max(ids)
        }, by = .(subjid, agedays)]

        # Filter out non-selected identicals
        tmp <- tmp[!(all_identical & internal_id != keep_id)]

        # Clean up temporary columns
        tmp[, `:=`(n_on_day = NULL, n_unique_vals = NULL, all_identical = NULL, keep_id = NULL)]

        # Now create seq_win on the filtered data (without SDE-Identicals)
        # Age-dependent ID sorting for consistent sequence numbering:
        # At age 0: sort by id ascending (lowest first)
        # At age > 0: sort by id descending (highest first)
        tmp[, id_sort := ifelse(agedays == 0, internal_id, -internal_id)]
        tmp <- tmp[order(subjid, agedays, id_sort),]
        tmp[, id_sort := NULL]

        # Create seq_win (sequence number per subject, after identicals removed)
        tmp[, seq_win := sequence(.N), by = subjid]
        # we're only looking at the first 4 values, and they need to be < 2 years
        tmp <- tmp[seq_win <= 4 & ageyears_2b < 2,]
        # don't look at subjects where there is only 1 score and there is no
        # value for either
        tmp <- tmp[!(is.na(sd.corr) & is.na(sd.orig)),]
        tmp <- tmp[subjid %in% names(which(table(subjid) > 1)),]

        # create differences and take the absolute value of their sum:
        #   abs(sum(sd.corr[1] - sd.corr)) == abs(sum of diffs from the first value)
        tmp[, sd.corr_abssumdiff := abs(sum(sd.corr[1] - sd.corr)), by = subjid]
        tmp[, sd.orig_abssumdiff := abs(sum(sd.orig[1] - sd.orig)), by = subjid]
        # find subjects where corrected value needs to be replaced

        # Only evaluate uncorr for the first weight per subject (earliest agedays)
        tmp[, is_first := seq_win == 1]

        sub_replace <- unique(tmp[sd.corr_abssumdiff > sd.orig_abssumdiff &
                                  !is.na(sd.orig_abssumdiff) &
                                  is_first == TRUE,
                                  subjid])

        if (!quietly)
          message(sprintf("  GA correction: %d subjects reverted to uncorrected z-scores", length(sub_replace)))

        # Revert correction for flagged subjects
        pc[subjid %in% sub_replace, sd.corr := sd.orig]
        pc[, uncorr := as.integer(subjid %in% sub_replace)]

        # --- Write results back to data.all ---
        data.all[, sd.corr := sd.orig]
        data.all[, uncorr := 0L]

        # Update potcorr subjects from subset results
        pc_result <- pc[, .(id, sd.corr, uncorr, potcorr)]
        data.all[pc_result, on = "id", `:=`(
          sd.corr = i.sd.corr,
          uncorr = i.uncorr,
          potcorr = i.potcorr
        )]

      } else {
        # No potcorr subjects -- fast path: skip all reference file I/O and merges
        data.all[, sd.corr := sd.orig]
        data.all[, uncorr := 0L]
      }


      # Add sd.corr, potcorr, uncorr to orig_colnames so they survive column cleanup
      orig_colnames <- c(orig_colnames, "sd.corr", "potcorr", "uncorr")

      # remove temporary columns (agemonths, ageyears_2b, potcorr_wt, etc.)
      data.all <- data.all[, colnames(data.all) %in% c(orig_colnames, "id"),
                           with = FALSE]


    # sort by subjid, param, agedays, internal_id for deterministic SDE order
    setkey(data.all, subjid, param, agedays, internal_id)

    # add a new convenience index for bookkeeping
    data.all[, index := 1:.N]

    # Mark missing values for exclusion
    data.all[, exclude := factor(
      ifelse(is.na(v) | agedays < 0, 'Exclude-Missing', 'Include'),
    levels = exclude.levels,
    ordered = TRUE)]
    # also mark certain measurements to not consider
    data.all[param == "HEADCM" & agedays > (3*365.25), exclude := "Exclude-Not-Cleaned"]

    # Initialize cf_rescued column to track CF rescue status
    # Populated in Step 6; rescued CFs get set back to "Include" but this column
    # preserves which rescue category they matched
    data.all[, cf_rescued := ""]

    # define field names needed by helper functions
    ewma.fields <- c('ewma.all', 'ewma.before', 'ewma.after')

    # SD-score recentering
    #
    # The code below reads a precomputed recentering file
    # (rcfile-2023-08-15_format.csv.gz) and subtracts its per-
    # (param, sex, agedays) median SD-score from sd.orig to produce
    # tbc.sd (and from sd.corr to produce ctbc.sd). Recentering is
    # needed because the population mean SD-score changes with age.
    #
    # The rcfile itself was built once using the procedure implemented
    # in sd_median() (see that function for the current-day code):
    #   a. Determine the median sd.orig for each param by year of age
    #      (sexes combined).
    #   b. Treat each year's median as applying to midyear-age (day
    #      floor(365.25 * year + 365.25 / 2)).
    #   c. Linearly interpolate between midyear medians by day of age;
    #      clamp to the earliest/latest year-median outside the covered
    #      range.
    # Users can override with a custom sd.recenter data.table if desired.

    if (!quietly)
      message(sprintf("[%s] Re-centering data...", Sys.time()))

    # Recentering: use built-in reference file or user-supplied data.table
    if (!is.data.table(sd.recenter)) {
      rc_path <- ifelse(
        ref.data.path == "",
        system.file(file.path("extdata",
                              "rcfile-2023-08-15_format.csv.gz"),
                    package = "growthcleanr"),
        file.path(ref.data.path, "rcfile-2023-08-15_format.csv.gz")
      )
      sd.recenter <- fread(rc_path)
      if (!quietly)
        message(sprintf("[%s] Using built-in reference medians...", Sys.time()))
    } else {
      if (!quietly)
        message(sprintf("[%s] Using user-supplied re-centering medians...", Sys.time()))
    }

    # ensure recentering medians are sorted correctly
    setkey(sd.recenter, param, sex, agedays)

    # add sd.recenter to data, and recenter
    setkey(data.all, param, sex, agedays)
    data.all <- sd.recenter[data.all]

    setkey(data.all, subjid, param, agedays, internal_id)
    data.all[, tbc.sd := sd.orig - sd.median]
    data.all[, ctbc.sd := sd.corr - sd.median]

    # notification: ensure awareness of small subsets in data
    if (!quietly) {
      year.counts <- data.all[, .N, floor(agedays / 365.25)]
      if (year.counts[N < 100, .N] > 0) {
        message(sprintf("[%s] Note: input data has at least one age-year with < 100 subjects...", Sys.time()))
      }
    }

    # safety check: treat observations where tbc.sd cannot be calculated as missing
    data.all[is.na(tbc.sd), exclude := 'Exclude-Missing']

    # HC >= 5 years: WHO HC reference only goes to 5 years, so no z-score data
    # exists. Combined with the pre-recentering assignment of
    # `Exclude-Not-Cleaned` for HC > 3*365.25, this keeps a single consistent
    # code across both "we don't clean HC >3y" and ">=5y has no reference".
    data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Exclude-Not-Cleaned']

    # pediatric: cleanchild (most of steps) ----

    # NOTE: the rest of cleangrowth's steps are done through cleanchild().

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
      message(sprintf("[%s] Cleaning growth data in %d batch(es)...", Sys.time(), num.batches))
    if (num.batches == 1) {
      ret.df <- cleanchild(
        data.all,
        log.path = log.path,
        quietly = quietly,
        parallel = parallel,
        measurement.to.z = measurement.to.z,
        ewma.fields = ewma.fields,
        cf_rescue = cf_rescue,
        cf_detail = cf_detail,
        biv.z.wt.low.young = biv.z.wt.low.young,
        biv.z.wt.low.old = biv.z.wt.low.old,
        biv.z.wt.high = biv.z.wt.high,
        biv.z.ht.low.young = biv.z.ht.low.young,
        biv.z.ht.low.old = biv.z.ht.low.old,
        biv.z.ht.high = biv.z.ht.high,
        biv.z.hc.low = biv.z.hc.low,
        biv.z.hc.high = biv.z.hc.high,
        exclude.levels = exclude.levels,
        tanner.ht.vel = tanner.ht.vel,
        error.load.threshold = error.load.threshold,
        error.load.mincount = error.load.mincount,
        ref.data.path = ref.data.path,
        ewma_window = ewma_window,
        ref_tables = ref_tables)
    } else {
      # create log directory if necessary
      if (!is.na(log.path)) {
        message(sprintf("[%s] Writing batch logs to '%s'...", Sys.time(), log.path))
        ifelse(!dir.exists(log.path), dir.create(log.path, recursive = TRUE), FALSE)
      }

      ret.df <- ddply(
        data.all,
        .(batch),
        cleanchild,
        .parallel = parallel,
        .paropts = list(.packages = c("data.table", "growthcleanr")),
        log.path = log.path,
        quietly = quietly,
        parallel = parallel,
        measurement.to.z = measurement.to.z,
        ewma.fields = ewma.fields,
        cf_rescue = cf_rescue,
        cf_detail = cf_detail,
        biv.z.wt.low.young = biv.z.wt.low.young,
        biv.z.wt.low.old = biv.z.wt.low.old,
        biv.z.wt.high = biv.z.wt.high,
        biv.z.ht.low.young = biv.z.ht.low.young,
        biv.z.ht.low.old = biv.z.ht.low.old,
        biv.z.ht.high = biv.z.ht.high,
        biv.z.hc.low = biv.z.hc.low,
        biv.z.hc.high = biv.z.hc.high,
        exclude.levels = exclude.levels,
        tanner.ht.vel = tanner.ht.vel,
        error.load.threshold = error.load.threshold,
        error.load.mincount = error.load.mincount,
        ref.data.path = ref.data.path,
        ewma_window = ewma_window,
        ref_tables = ref_tables
      )
    }


    if (!quietly)
      message(sprintf("[%s] Done with pediatric data!", Sys.time()))
  } else {
    ret.df <- data.table()

    if (!quietly)
      message(sprintf("[%s] No pediatric data. Moving to adult data...", Sys.time()))
  }

  # adult: send to cleanadult to do most of the work ----

  # no need to do this if there's no data
  if (nrow(data.adult) > 0){
    if (!quietly){
      message(sprintf("[%s] Begin processing adult data...", Sys.time()))
    }

    # this is where we do most of the adult work
    subjid.unique <- data.adult[j = unique(subjid)]
    batches.adult <- data.table(
      subjid = subjid.unique,
      newbatch = (seq_along(subjid.unique) - 1L) %% num.batches + 1L
    )
    data.adult <- merge(data.adult, batches.adult, by = "subjid")

    # add age in years
    data.adult[, age_years := agedays/365.25]
    # Provide measurement column expected by cleanadult(). Uses v_adult (original
    # units, pre-conversion) because cleanadult() handles imperial conversion
    # internally. The child path uses v (already converted to metric above).
    data.adult[, measurement := v_adult]

    if (num.batches == 1) {
      # do the cleaning
      res <- cleanadult(data.adult,
                        permissiveness = adult_permissiveness,
                        scale_max_lbs = adult_scale_max_lbs,
                        quietly = quietly)
    } else {
      res <- ddply(
        data.adult,
        .(newbatch),
        cleanadult,
        .parallel = parallel,
        .paropts = list(.packages = c("data.table", "growthcleanr")),
        permissiveness = adult_permissiveness,
        scale_max_lbs = adult_scale_max_lbs,
        quietly = quietly
      )

      res <- as.data.table(res)
    }

    # Apply Exclude-Missing AFTER cleanadult() returns (adult BIV handles NAs).
    # This differs from the child path, where Exclude-Missing is set BEFORE
    # dispatch (lines ~1333, ~1337). Both are correct — the child algorithm
    # requires pre-set Missing codes for .child_valid(), while the adult
    # algorithm handles NAs internally via BIV step 1.
    res[is.na(measurement) | agedays < 0, result := "Exclude-Missing"]

    if (!quietly)
      message(sprintf("[%s] Done with adult data!", Sys.time()))
  } else {
    res <- data.table()

    if (!quietly){
      message(sprintf("[%s] No adult data. Moving to postprocessing...", Sys.time()))
    }
  }


  # --- FINAL RETURN (modified to return all columns) ---

  # combine pediatric and adult results if present
  if (any(nrow(data.all) > 0, nrow(data.adult) > 0)) {

    # join pediatric and adult outputs together
    # cf_rescued only exists in the child algorithm; use empty strings for legacy path
    cf_rescued_peds <- if (!is.null(ret.df$cf_rescued)) {
      as.character(ret.df$cf_rescued)
    } else {
      rep("", nrow(ret.df))
    }

    # adult-specific columns: mean_ht and bin_result (NA for child rows)
    adult_mean_ht <- if (nrow(res) > 0 && "mean_ht" %in% names(res)) {
      res$mean_ht
    } else {
      rep(NA_real_, nrow(res))
    }
    adult_bin_result <- if (nrow(res) > 0 && "bin_result" %in% names(res)) {
      res$bin_result
    } else {
      rep(NA_character_, nrow(res))
    }

    full_out <- data.table(
      line = c(ret.df$line, res$line),
      exclude = c(as.character(ret.df$exclude), res$result),
      cf_rescued = c(cf_rescued_peds, rep("", nrow(res))),
      mean_ht = c(rep(NA_real_, nrow(ret.df)), adult_mean_ht),
      bin_result = c(rep(NA_character_, nrow(ret.df)), adult_bin_result)
    )

    # Add cf_detail columns if present in child output
    if (cf_detail && "cf_status" %in% names(ret.df)) {
      full_out[, cf_status := c(as.character(ret.df$cf_status), rep(NA_character_, nrow(res)))]
      full_out[, cf_deltaZ := c(ret.df$cf_deltaZ, rep(NA_real_, nrow(res)))]
    }

    # preserve original exclude levels
    full_out[, exclude := factor(exclude, levels = exclude.levels)]

    # join back to batch-level reference to preserve all columns
    all_results <- merge(
      data.batch,              # this batch's subjects (all ages)
      full_out,
      by = "line",
      all.x = TRUE,
      sort = FALSE
    )

    # Merge checkpoint diagnostics by ID (child algorithm only).
    # checkpoint_data uses the user's original id. Adult rows get NA for
    # checkpoint columns via all.x = TRUE, which is correct.
    if (exists("checkpoint_data") && nrow(checkpoint_data) > 0) {
      drop_cols <- intersect(c("subjid", "param", "agedays"), names(checkpoint_data))
      cp_merge <- checkpoint_data[, !..drop_cols]
      all_results <- merge(all_results, cp_merge, by = "id", all.x = TRUE)
    }

    # Reorder to match data.batch row order (id is guaranteed unique by input contract)
    all_results <- all_results[match(data.batch$id, all_results$id)]

    # restore original row order
    setorder(all_results, line)

    # if you have an 'id' field in data.work, it will now be preserved automatically
    # (since it existed in the original input)

    results_list[[i]] <- all_results
    if (!quietly) message(sprintf("Batch %d complete.", i))
    i <- i + 1


  } else {
    # if there was no data at all, return an empty table with consistent structure
    results_list[[i]] <- data.table::data.table()
    if (!quietly) message(sprintf("Batch %d complete.", i))
    i <- i + 1
  }

  }

  # --- Parallel teardown (once, after outer batch loop) ---
  if (parallel) stopCluster(cl)

  all_results <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
  setorder(all_results, line)

  # Partial-run mode: merge new results for changed subjects with cached results
  # for unchanged subjects, then sort by user id for a stable, semantically-
  # meaningful output order. (internal_id is renumbered 1..K for the partial
  # subset and would collide with the cached 1..N; user id is preserved
  # untouched through both paths and is contract-guaranteed unique.)
  if (do_partial) {
    all_results <- rbind(
      cached_results[!(as.character(cached_results$subjid) %in%
                         as.character(changed_subjids))],
      all_results,
      fill = TRUE
    )
    all_results <- all_results[order(id)]
  }

  # Derive bin_exclude and tri_exclude from final exclude codes.
  # These are added after all processing (including partial-run merge)
  # so they are consistent regardless of code path.

  # Check for Exclude-C-Temp-Same-Day surviving to output —
  # this is an internal working code that should be resolved before output.
  if (any(as.character(all_results$exclude) == "Exclude-C-Temp-Same-Day",
          na.rm = TRUE)) {
    stop("Exclude-C-Temp-Same-Day found in final output. ",
         "This is an internal code that should have been resolved during processing. ",
         "This indicates a bug in the algorithm.")
  }

  excl_char <- as.character(all_results$exclude)

  # bin_exclude: Include vs Exclude (always present)
  all_results[, bin_exclude := fifelse(excl_char == "Include", "Include", "Exclude")]

  # tri_exclude: Include vs Same-Day vs Exclude (opt-in via tri_exclude parameter)
  if (tri_exclude) {
    # Child SDE codes: Identical and Extraneous for each param
    # Adult SDE codes: same pattern with -A- prefix
    child_sde_codes <- c("Exclude-C-Identical", "Exclude-C-Extraneous")
    adult_sde_codes <- c("Exclude-A-Identical", "Exclude-A-Extraneous")
    sde_codes <- c(child_sde_codes, adult_sde_codes)

    all_results[, tri_exclude := fifelse(
      excl_char == "Include", "Include",
      fifelse(excl_char %in% sde_codes, "Same-Day", "Exclude")
    )]
  }

  return(all_results)
}

#' Function to calculate z-scores and csd-scores based on anthro tables.
#'
#' @param path Path to supplied reference anthro data. Defaults to package anthro tables.
#' @param cdc.only Whether or not only CDC data should be used. Defaults to false.
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
read_anthro <- function(path = "", cdc.only = FALSE) {
  # avoid "no visible bindings" warning
  src <- param <- sex <- age <- ret <- m <- NULL
  csdneg <- csdpos <- s <- NULL

  # Child algorithm reference tables (consolidated WHO + infants CDC)
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
  growth_cdc_ext <- read.csv(gzfile(growth_cdc_ext_path))

  l <- list(
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

#' Pre-load growthcleanr reference closures
#'
#' \code{gc_preload_refs} loads all reference table closures used by
#' \code{\link{cleangrowth}} once, avoiding repeated disk reads across
#' repeated calls (e.g., simulation loops). Pass the result as the
#' \code{ref_tables} argument to \code{\link{cleangrowth}}.
#'
#' @param path Path to the reference data directory. Defaults to the package
#'   extdata directory. Pass an explicit path only if using custom reference
#'   files.
#'
#' @return A named list with two \code{read_anthro} closures:
#'   \describe{
#'     \item{mtz_cdc_prelim}{CDC-only (used by child algorithm)}
#'     \item{mtz_who_prelim}{WHO+CDC (used by child algorithm)}
#'   }
#'
#' @export
#' @examples
#' \donttest{
#' refs <- gc_preload_refs()
#' result <- cleangrowth(
#'   subjid = syngrowth$subjid,
#'   param  = syngrowth$param,
#'   agedays = syngrowth$agedays,
#'   sex = syngrowth$sex,
#'   measurement = syngrowth$measurement,
#'   ref_tables = refs
#' )
#' }
gc_preload_refs <- function(path = "") {
  list(
    mtz_cdc_prelim = read_anthro(path, cdc.only = TRUE),
    mtz_who_prelim = read_anthro(path, cdc.only = FALSE)
  )
}

# Convert vector of ages (days) into pairwise absolute-difference matrix.
# Used by ewma() for age-gap weighting.
as_matrix_delta <- function(agedays) {
  n <- length(agedays)
  delta <- abs(matrix(rep(agedays, n), n, byrow = TRUE) - agedays)
  return(delta)
}

#' Exponentially Weighted Moving Average (EWMA)
#'
#' \code{ewma} calculates the exponentially weighted moving average (EWMA) for a set of numeric observations over time.
#'
#' @param agedays Vector of age in days for each z score (potentially transformed to adjust weighting).
#'
#' @param z Input vector of numeric z-score data.
#'
#' @param ewma.exp Exponent used for age-gap weighting. May be a single
#'   scalar (recycled across all observations) or a per-observation vector
#'   matched to \code{agedays} / \code{z}; the algorithm's main internal
#'   callers pass a per-observation vector that varies the exponent by the
#'   widest neighbor age gap.
#'
#' @param ewma.adjacent Specify whether EWMA values excluding adjacent measurements should be calculated.  Defaults to TRUE.
#'
#' @param window Maximum number of observations on each side for EWMA weighting.
#'   Default 15. Set to Inf to disable windowing.
#'
#' @param cache_env Optional environment for caching EWMA intermediate results.
#'   Used internally for performance. Defaults to NULL.
#'
#' @return A named list of numeric vectors, each the same length as
#'   \code{agedays}.
#' * \code{ewma.all} — EWMA at each observation's age, excluding only the
#'   observation itself.
#' * \code{ewma.before} — also excludes the immediate prior observation
#'   (sorted by \code{agedays}); equals \code{ewma.all} for the first
#'   observation and when \code{n <= 2}.
#' * \code{ewma.after} — also excludes the immediate subsequent observation;
#'   equals \code{ewma.all} for the last observation and when \code{n <= 2}.
#'
#' When \code{ewma.adjacent = FALSE}, only \code{ewma.all} is returned (still
#'   wrapped in a one-element named list).
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
ewma <- function(agedays, z, ewma.exp, ewma.adjacent = TRUE, window = 15, cache_env = NULL) {
  # For each observation i in (agedays, z), compute a weighted average of the
  # other observations' z-values using weights w_ij = (5 + |agedays_i - agedays_j|) ^ ewma.exp_i,
  # with w_ii = 0 so observation i is excluded from its own EWMA. The "+5"
  # prevents nearly-coincident agedays from blowing up the weight; the negative
  # exponent (typically -1.5 to -3.5) damps the contribution of distant ages.
  #
  # Three EWMA variants are returned (when ewma.adjacent = TRUE) so that
  # callers can decide which neighbors to suppress when checking whether an
  # observation is consistent with its trajectory:
  #   ewma.all     EWMA of all other observations
  #   ewma.before  also drops obs i-1 (immediate predecessor by agedays)
  #   ewma.after   also drops obs i+1 (immediate successor by agedays)
  # For the first/last observation, ewma.before/ewma.after equal ewma.all
  # (no predecessor/successor exists). For n <= 2, both equal ewma.all.
  #
  # Caller supplies one ewma.exp value per observation; sweep() applies the
  # exponent row-wise so each row of the weight matrix uses that row's
  # observation's exponent. A scalar ewma.exp is recycled by sweep() and
  # works for the simple single-call form shown in the example.
  #
  # window limits the EWMA at each observation to the ±window nearest
  # positions (default 15); set window = Inf to disable.
  #
  # cache_env (optional environment) lets callers reuse the windowed weight
  # matrix across paired calls on the same agedays + exponents (e.g.,
  # tbc.sd then ctbc.sd) to avoid recomputing the matrix.

  n <- length(agedays)
  # initialize response variables
  ewma.all <- ewma.before <- ewma.after <- vector('numeric', 0)
  if (n > 0) {
    # organize into data frame and sort into order of increasing age,
    # but retain original sort order information in index
    if (!all(agedays == cummax(agedays)))
      warning("EWMA ordering is not sorted; double check") #add in a check to make sure the inputs are already sorted (they should be)
    index <- order(agedays)

    # Use cached delta matrix if available, otherwise compute and optionally cache
    if (!is.null(cache_env) && !is.null(cache_env$delta)) {
      delta <- cache_env$delta
    } else {
      # Build weight matrix: row-wise exponent application means each
      # observation's EWMA uses that observation's own exponent for all
      # weights in its row. sweep(..., 1, ...) targets margin 1 (rows).
      delta <- as_matrix_delta(agedays)
      delta <- sweep(delta + 5, 1, ewma.exp, FUN = "^")
      diag(delta) <- 0  # self-weight is zero so obs i is excluded from its own EWMA

      # Apply windowing: zero out weights beyond window positions.
      if (!is.null(window) && is.finite(window)) {
        pos_diff <- abs(outer(1:n, 1:n, "-"))  # position difference matrix
        delta[pos_diff > window] <- 0
      }

      # Cache delta if environment provided
      if (!is.null(cache_env)) cache_env$delta <- delta
    }

    # Compute main EWMA using matrix multiply
    weighted_sums <- as.vector(delta %*% z)
    row_sums <- rowSums(delta)
    ewma.all[index] <- weighted_sums / row_sums

    if (ewma.adjacent) {
      if (n > 2) {
        # Before / After: subtract predecessor's / successor's contribution
        # from each observation's weighted sum and row sum directly. This
        # is O(n) per variant — no need to rebuild a separate weight matrix
        # for each adjacent-neighbor exclusion.
        pred_weights <- c(0, delta[cbind(2:n, 1:(n-1))])  # subdiagonal: delta[i, i-1]
        pred_z <- c(0, z[1:(n-1)])
        ewma.before[index] <- (weighted_sums - pred_weights * pred_z) / (row_sums - pred_weights)

        succ_weights <- c(delta[cbind(1:(n-1), 2:n)], 0)  # superdiagonal: delta[i, i+1]
        succ_z <- c(z[2:n], 0)
        ewma.after[index] <- (weighted_sums - succ_weights * succ_z) / (row_sums - succ_weights)
      } else {
        ewma.before <- ewma.after <- ewma.all
      }
    }
  }
  # Return as a named list (not a data.frame): list is data.table :=
  # compatible and avoids data.frame construction overhead in the per-group
  # callers in Steps 11, 13, 15/16, and 17.
  return(if (ewma.adjacent)
    list(ewma.all = ewma.all, ewma.before = ewma.before, ewma.after = ewma.after)
    else
      list(ewma.all = ewma.all))
}

# --- Incremental EWMA cache functions ---
# These allow O(n) updates when observations are excluded from EWMA groups,
# instead of O(n^2) full rebuild each iteration.

#' Create initial EWMA cache for a group (O(n^2) — first iteration only)
#' @keywords internal
#' @noRd
ewma_cache_init <- function(agedays, z_tbc, z_ctbc, exp_vals, ids, window = 15) {
  n <- length(agedays)

  # Check if ctbc == tbc (skip redundant ctbc computation for most groups)
  skip_ctbc <- identical(z_tbc, z_ctbc)

  # Build delta (weight) matrix
  delta <- as_matrix_delta(agedays)
  delta <- sweep(delta + 5, 1, exp_vals, FUN = "^")
  diag(delta) <- 0

  # Apply windowing
  if (!is.null(window) && is.finite(window)) {
    pos_diff <- abs(outer(1:n, 1:n, "-"))
    delta[pos_diff > window] <- 0
  }

  # Compute weighted sums and row sums
  ws_tbc <- as.vector(delta %*% z_tbc)
  rs <- rowSums(delta)

  # Only compute ctbc weighted sums if ctbc differs from tbc
  if (skip_ctbc) {
    ws_ctbc <- ws_tbc
  } else {
    ws_ctbc <- as.vector(delta %*% z_ctbc)
  }

  # Compute EWMA values
  ewma.all <- ws_tbc / rs
  c.ewma.all <- if (skip_ctbc) ewma.all else ws_ctbc / rs

  if (n > 2) {
    pred_w <- c(0, delta[cbind(2:n, 1:(n - 1))])
    pred_z <- c(0, z_tbc[1:(n - 1)])
    ewma.before <- (ws_tbc - pred_w * pred_z) / (rs - pred_w)
    succ_w <- c(delta[cbind(1:(n - 1), 2:n)], 0)
    succ_z <- c(z_tbc[2:n], 0)
    ewma.after <- (ws_tbc - succ_w * succ_z) / (rs - succ_w)
  } else {
    ewma.before <- ewma.after <- ewma.all
  }

  list(
    delta = delta, ws_tbc = ws_tbc, ws_ctbc = ws_ctbc, rs = rs,
    ids = ids, agedays = agedays, exp_vals = exp_vals,
    z_tbc = z_tbc, z_ctbc = z_ctbc, window = window,
    skip_ctbc = skip_ctbc,
    ewma.all = ewma.all, ewma.before = ewma.before,
    ewma.after = ewma.after, c.ewma.all = c.ewma.all
  )
}

#' Update EWMA cache after excluding one observation (O(n) instead of O(n^2))
#' @keywords internal
#' @noRd
ewma_cache_update <- function(cache, excluded_id) {
  pos_j <- which(cache$ids == excluded_id)
  if (length(pos_j) != 1L) {
    # Cache mismatch — return NULL to signal caller should do full rebuild
    return(NULL)
  }

  n <- length(cache$ids)
  new_n <- n - 1L
  keep <- seq_len(n) != pos_j
  skip_ctbc <- isTRUE(cache$skip_ctbc)

  # Remove j's contribution from weighted sums and row sums — O(n)
  col_j <- cache$delta[keep, pos_j]
  ws_tbc <- cache$ws_tbc[keep] - col_j * cache$z_tbc[pos_j]
  rs <- cache$rs[keep] - col_j

  # Only update ctbc sums if ctbc differs from tbc
  if (skip_ctbc) {
    ws_ctbc <- ws_tbc
  } else {
    ws_ctbc <- cache$ws_ctbc[keep] - col_j * cache$z_ctbc[pos_j]
  }

  # Remove row/col from delta matrix
  delta <- cache$delta[keep, keep, drop = FALSE]
  agedays <- cache$agedays[keep]
  z_tbc <- cache$z_tbc[keep]
  z_ctbc <- cache$z_ctbc[keep]
  ids <- cache$ids[keep]
  exp_vals <- cache$exp_vals[keep]
  window <- cache$window

  # Check if neighbors' exponents changed after removal.
  # Must check up to 2 positions on each side in the new index:
  # removing pos_j changes the max-gap for pos_j-1 (now adjacent
  # to pos_j-2 and pos_j) AND for pos_j-2 (whose gap-after is
  # now to pos_j instead of pos_j-1). Same logic applies on the
  # other side.
  neighbors <- integer(0)
  if (pos_j > 2L)     neighbors <- c(neighbors, pos_j - 2L)
  if (pos_j > 1L)     neighbors <- c(neighbors, pos_j - 1L)
  if (pos_j <= new_n) neighbors <- c(neighbors, pos_j)
  if (pos_j < new_n)  neighbors <- c(neighbors, pos_j + 1L)
  neighbors <- unique(
    neighbors[neighbors >= 1L & neighbors <= new_n])

  for (nb in neighbors) {
    # Recompute max gap for this neighbor
    diff_before <- if (nb > 1L) abs(agedays[nb] - agedays[nb - 1L]) else NA
    diff_after <- if (nb < new_n) abs(agedays[nb + 1L] - agedays[nb]) else NA
    maxdiff <- max(diff_before, diff_after, na.rm = TRUE)
    ageyears <- maxdiff / 365.25
    new_exp <- if (ageyears <= 1) -1.5 else if (ageyears >= 3) -3.5 else -1.5 - (ageyears - 1)

    if (abs(new_exp - exp_vals[nb]) > 1e-10) {
      # Exponent changed — rebuild this row of delta and recompute its sums
      exp_vals[nb] <- new_exp
      new_row <- (abs(agedays - agedays[nb]) + 5) ^ new_exp
      new_row[nb] <- 0  # self-weight

      # Apply windowing
      if (!is.null(window) && is.finite(window)) {
        pos_diffs <- abs(seq_len(new_n) - nb)
        new_row[pos_diffs > window] <- 0
      }

      delta[nb, ] <- new_row
      ws_tbc[nb] <- sum(new_row * z_tbc)
      if (!skip_ctbc) ws_ctbc[nb] <- sum(new_row * z_ctbc)
      rs[nb] <- sum(new_row)
    }
  }

  # If skip_ctbc, keep ws_ctbc in sync after neighbor rebuilds
  if (skip_ctbc) ws_ctbc <- ws_tbc

  # Compute EWMA values from cached sums — O(n)
  ewma.all <- ws_tbc / rs
  c.ewma.all <- if (skip_ctbc) ewma.all else ws_ctbc / rs

  if (new_n > 2L) {
    pred_w <- c(0, delta[cbind(2:new_n, 1:(new_n - 1L))])
    pred_z <- c(0, z_tbc[1:(new_n - 1L)])
    ewma.before <- (ws_tbc - pred_w * pred_z) / (rs - pred_w)
    succ_w <- c(delta[cbind(1:(new_n - 1L), 2:new_n)], 0)
    succ_z <- c(z_tbc[2:new_n], 0)
    ewma.after <- (ws_tbc - succ_w * succ_z) / (rs - succ_w)
  } else {
    ewma.before <- ewma.after <- ewma.all
  }

  list(
    delta = delta, ws_tbc = ws_tbc, ws_ctbc = ws_ctbc, rs = rs,
    ids = ids, agedays = agedays, exp_vals = exp_vals,
    z_tbc = z_tbc, z_ctbc = z_ctbc, window = window,
    skip_ctbc = skip_ctbc,
    ewma.all = ewma.all, ewma.before = ewma.before,
    ewma.after = ewma.after, c.ewma.all = c.ewma.all
  )
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
  # group all measurements above age 19 together
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

#' Look up the designated other parameter (DOP) for a growth parameter.
#'
#' Maps a single parameter name to its paired "other parameter" used in
#' cross-parameter plausibility checks (Child Steps 5, 6, 11, 13, 15/16,
#' 19). The DOP assignments are: WEIGHTKG -> HEIGHTCM (height is weight's
#' anchor), HEIGHTCM -> WEIGHTKG (weight is height's anchor), and
#' HEADCM -> HEIGHTCM (height is HC's anchor; there is no reverse
#' mapping).
#'
#' Scalar-only: the internal `if / else if / else` chain requires a
#' scalar logical test, so callers must pass a scalar param string. In
#' practice both call sites pass `df$param[1]` from a single-param
#' working subset.
#'
#' @param param_name character scalar: one of `"WEIGHTKG"`,
#'   `"HEIGHTCM"`, or `"HEADCM"`. No validation — anything other than
#'   `"WEIGHTKG"` or `"HEIGHTCM"` falls through to the HEADCM branch
#'   and returns `"HEIGHTCM"`.
#'
#' @return Character scalar: the DOP for the given param.
#'
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
#' Identify temporary SDE (same-day extraneous) "losers" on each
#' (subjid, param, agedays) group.
#'
#' Identical same-day values are removed earlier (Early Step 13
#' SDE-Identicals), so by the time identify_temp_sde() runs, all
#' same-day values in a given (subjid, param, agedays) group are
#' dissimilar. For each such group with more than one valid
#' measurement, selects one value to keep temporarily and flags the
#' remainder for exclusion. Selection uses the recentered z-score
#' (tbc.sd) and two per-subject reference points:
#'   median.spz  - median of tbc.sd for the same (subjid, param),
#'                 taken over valid rows (all non-excluded
#'                 measurements).
#'   median.dopz - median of tbc.sd for the subject's designated other
#'                 parameter (DOP).
#'
#' For each SDE group, the kept row is the one whose tbc.sd is closest
#' to median.spz, provided the subject has at least one non-SDE day.
#' Otherwise the row closest to median.dopz is kept. If no DOP data
#' exists for the subject, one row is chosen at random.
#'
#' Called from Steps 5, 6, 7, 9, 11 (mid-loop and end-of-step), and 13
#' of \code{cleanchild()}. All callers except Step 13 final SDE pass
#' \code{exclude_from_dop_ids = NULL}; the Step 13 caller passes the
#' ids of current temp SDEs so those rows do not contribute to the DOP
#' median anchor used to resolve remaining same-day groups.
#'
#' @param df data.table column-subset containing id, internal_id,
#'   subjid, param, agedays, tbc.sd, and exclude. Not mutated.
#' @param exclude_from_dop_ids optional vector of \code{id} values to
#'   remove from the DOP median calculation only (they still contribute
#'   to the SP median). NULL by default; non-NULL only at Step 13.
#' @return A logical vector aligned to the caller's input row order.
#'   TRUE marks rows the caller should flag
#'   \code{'Exclude-C-Temp-Same-Day'}; FALSE means leave the exclude
#'   column unchanged.
#'
#' @keywords internal
#' @noRd

identify_temp_sde <- function(df, exclude_from_dop_ids = NULL) {
  # exclude_from_dop_ids: when non-NULL, the passed ids are removed from the
  # DOP median calculation (but still contribute to the SP median). Only the
  # Step 13 caller passes ids; all earlier callers (Step 5 initial pass and
  # the recalc callers in Steps 6, 7, 9, and 11) leave it NULL.

  # avoid "no visible binding" warnings
  agedays <- absdmedian.spz <- absdmedian.dopz <- extraneous <- NULL
  extraneous.this.day <- index <- median.spz <- median.dopz <- NULL
  param <- subjid <- tbc.sd <- internal_id <- NULL

  # Make copy before modifying to avoid data.table alloccol error
  df <- copy(df)

  # Record the caller's original row positions before any reordering, so the
  # final result can be mapped back. Row ordering is set deterministically
  # below by explicit order(subjid, param, agedays, internal_id).
  df[, orig_row := .I]

  # Narrow df to the columns we need, grouped by (subjid, param,
  # agedays) and then sorted with internal_id as the tiebreaker so
  # same-day rows have a deterministic order for the winner selection
  # below.
  df <- df[, .(tbc.sd, exclude, id, internal_id, orig_row), by = .(subjid, param, agedays)]
  df <- df[order(subjid, param, agedays, internal_id)]

  # Compute valid.rows on the sorted data.
  valid.rows <- .child_valid(df, include.temporary.extraneous = TRUE)

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

  # Assign each valid row its SP median: median of tbc.sd within
  # (subjid, param), computed over all valid rows (not just SDE days).
  df[valid.rows, median.spz := median(tbc.sd, na.rm = TRUE), by = .(subjid, param)]

  # Compact (subjid, param) -> median.spz lookup table. Used below to
  # anchor each row's DOP median to the median of its other parameter.
  sp_medians <- df[valid.rows & !is.na(median.spz),
                   .(median.spz = median.spz[1]),
                   by = .(subjid, param)]

  # Assign each row its DOP median. When exclude_from_dop_ids is NULL
  # (all callers except Step 13) the DOP median is just the SP median
  # of the other parameter. When non-NULL (Step 13), recompute medians
  # for each param excluding the passed ids first, so the Step 13 final
  # SDE resolution is not biased by rows that are themselves under
  # consideration as temp SDEs.
  df[, median.dopz := as.double(NA)]

  if (is.null(exclude_from_dop_ids)) {
    # Default path (all callers except Step 13): reuse sp_medians
    # directly as each parameter's DOP median lookup.
    for (p in c("WEIGHTKG", "HEIGHTCM", "HEADCM")) {
      dop <- dop_map[p]
      dop_medians <- sp_medians[param == dop, .(subjid, dop_median = median.spz)]
      dop_lookup <- setNames(dop_medians$dop_median, dop_medians$subjid)
      df[param == p, median.dopz := dop_lookup[as.character(subjid)]]
    }
  } else {
    # Step 13 path: recompute per-(subjid, param) medians excluding the
    # passed ids, then use those as the DOP median lookup.
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

  # For rows on an SDE day, compute the absolute distance of tbc.sd
  # from the two median anchors. These distances drive the selection
  # below.
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

  # Within each SDE group, pick the winner by sorting on absdmedian.spz
  # (asc), then absdmedian.dopz (asc), then age-dependent internal_id.
  # The first row after sorting is the keeper; all others are marked
  # extraneous.
  df[valid.rows & extraneous.this.day, extraneous := {
    # Age-dependent internal_id tiebreaker:
    #   At agedays == 0: keep lowest internal_id (earliest, before fluid shifts)
    #   At agedays  > 0: keep highest internal_id (later, more careful measurement)
    # internal_id is used (not user id) because user id may be character.
    tiebreaker <- if (agedays[1] == 0) internal_id else -internal_id
    ord <- order(absdmedian.spz,
                 absdmedian.dopz, tiebreaker)
    keep_id <- id[ord[1]]
    id != keep_id
  }, by = .(subjid, param, agedays)]

  # Map result back to the caller's original row order. df was
  # re-ordered by the by-group work above; orig_row preserves each
  # row's caller-side position. df$extraneous is only ever set TRUE by
  # the guarded `extraneous := ...` assignment above, which restricts
  # to valid rows on SDE days, so we can use it directly without
  # re-masking on valid.rows.
  result <- logical(nrow(df))
  result[df$orig_row] <- df$extraneous
  return(result)
}

# evil twins ----

#' Function to calculate over-the-limit (OTL) measurements for Evil Twins step
#'
#' @param df data table with all parameters
#'
#' @keywords internal
#' @noRd
calc_otl_evil_twins <- function(df){
  # # start by determining if a measurement is over the limit (otl)

  # Handle minimal datasets with < 2 rows
  # When there are 0 or 1 rows, there are no adjacent pairs to compare for evil twins
  if (nrow(df) < 2) {
    df[, "otl" := FALSE]
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

  otl <- (tbc_next_diff > 5 & ctbc_next_diff > 5 & same_sp_next) |
         (tbc_prev_diff > 5 & ctbc_prev_diff > 5 & same_sp_prev)

  df[, "otl" := otl]

  return(df)
}

# moderate ewma ----

#' Recalculate recentered z-scores for a modified measurement column.
#'
#' Used by Child Step 15/16 (EWMA2 / moderate trajectory outliers) to
#' compute `tbc.p_plus` and `tbc.p_minus` from pre-filled `p_plus` /
#' `p_minus` columns — the ±5%-weight / ±1 cm perturbations of each
#' included value. The same WHO/CDC blending window and NA-fallback
#' logic as the main `cleangrowth()` z-score calculation are applied
#' here so that the perturbed z-scores are on the same footing as
#' `tbc.sd`.
#'
#' The `sd.median` column (populated by the main recentering in
#' `cleangrowth()`) must already be present on `df`. The helper looks
#' it up per row rather than recomputing or re-merging the recentering
#' file.
#'
#' Called twice per batch from the Step 15/16 pre-loop inside
#' `cleanchild()`, once with `cn = "p_plus"` and once with
#' `cn = "p_minus"`. Both callers pass pre-built `measurement.to.z` /
#' `measurement.to.z_who` closures so the expensive reference-table
#' reads happen once per batch rather than once per call.
#'
#' @param df data.table with columns `param`, `agedays`, `sex`,
#'   `sd.median`, and the measurement column named by `cn`. Modified
#'   in place: intermediate columns `cn.orig_cdc`, `cn.orig_who`,
#'   `cn.orig`, and a new `tbc.<cn>` column are added.
#' @param cn character scalar naming the measurement column to
#'   z-score and recenter (e.g. `"p_plus"` or `"p_minus"`). The output
#'   column is named `paste0("tbc.", cn)`.
#' @param ref.data.path path to reference-table directory; `""` uses
#'   the installed package's `inst/extdata/`. Only consulted when
#'   `measurement.to.z` / `measurement.to.z_who` are NULL.
#' @param measurement.to.z optional pre-built CDC-only closure from
#'   `read_anthro(..., cdc.only = TRUE)`. If NULL, is constructed here
#'   (one disk read).
#' @param measurement.to.z_who optional pre-built WHO-blended closure
#'   from `read_anthro(..., cdc.only = FALSE)`. If NULL, is
#'   constructed here (one disk read).
#'
#' @return The input `df` with an additional `tbc.<cn>` column
#'   containing the recentered, age-blended CSD z-score of the
#'   measurement column. The intermediate columns `cn.orig_cdc`,
#'   `cn.orig_who`, and `cn.orig` are also left on `df` but are not
#'   consumed downstream.
#'
#' @keywords internal
#' @noRd
calc_and_recenter_z_scores <- function(df, cn, ref.data.path,
                                       measurement.to.z = NULL,
                                       measurement.to.z_who = NULL){
  # avoid no visible warning errors
  cn.orig_cdc <- param <- agedays <- sex <- cn.orig_who <- cn.orig <- subjid <-
    tbc.cn <- sd.median <- NULL

  # Use pre-built closures if provided (the Step 15/16 callers always do).
  # The NULL fallback builds them on the fly — one disk read per closure —
  # and exists for direct / debugging use only.
  if (is.null(measurement.to.z))
    measurement.to.z <- read_anthro(ref.data.path, cdc.only = TRUE)
  if (is.null(measurement.to.z_who))
    measurement.to.z_who <- read_anthro(ref.data.path, cdc.only = FALSE)

  # CSD z-scores against each reference, one closure call per row of df.
  df[, cn.orig_cdc := measurement.to.z(param, agedays, sex, get(cn), TRUE)]
  df[, cn.orig_who := measurement.to.z_who(param, agedays, sex, get(cn), TRUE)]

  # Age-blending weights. WHO gets more weight at younger ages, CDC at
  # older ages; window boundaries and the /3 divisor match the main
  # z-score blend in cleangrowth() exactly.
  who_weight <- 5 - (df$agedays/365.25)
  cdc_weight <- (df$agedays/365.25) - 2

  smooth_val <- df$agedays/365.25 >= 2 & df$agedays/365.25 <= 5 & df$param != "HEADCM"

  # Inside the 2-5 year smoothing window (HT/WT only), blend the two
  # z-scores. Only cdc_weight and who_weight need to be re-subset with
  # [smooth_val]; the column references cn.orig_cdc / cn.orig_who on
  # the RHS are already restricted to smooth_val rows by the outer
  # df[smooth_val, ...] subset.
  df[smooth_val,
     cn.orig := (cn.orig_cdc * cdc_weight[smooth_val] +
                   cn.orig_who * who_weight[smooth_val])/3]

  # Under 2y (and HEADCM at any age): pure WHO. who_val and smooth_val
  # are mutually exclusive (HEADCM is excluded from smooth; <2y is
  # excluded from smooth), so this cannot overwrite the blended value.
  who_val <- df$param == "HEADCM" | df$agedays/365.25 < 2
  df[who_val, cn.orig := df$cn.orig_who[who_val]]

  # Over 5y HT/WT: pure CDC. cdc_val and smooth_val are similarly
  # mutually exclusive.
  cdc_val <- df$param != "HEADCM" & df$agedays/365.25 > 5
  df[cdc_val, cn.orig := df$cn.orig_cdc[cdc_val]]

  # NA fallbacks inside the smooth zone: if one reference is unavailable
  # for a row, use the other. Matches the main z-score path.
  df[smooth_val & is.na(cn.orig_cdc), cn.orig := cn.orig_who]
  df[smooth_val & is.na(cn.orig_who), cn.orig := cn.orig_cdc]

  # Recenter against the sd.median column populated by the main
  # recentering step; df inherits sd.median from its source data.df.
  df[, tbc.cn := cn.orig - sd.median]

  # rename ending column
  setnames(df, "tbc.cn", paste0("tbc.", cn))

  return(df)
}

# CF rescue threshold lookup table
# Returns the deltaZ threshold for CF rescue given age, interval, param, and rounding type.
# Rescue if |deltaZ| < threshold; exclude if |deltaZ| >= threshold.
# Returns 0 for NR (no rescue) cells and NA for impossible cells (--).
# See cf-rescue-thresholds.md for full documentation.
.cf_rescue_lookup <- function() {
  # Build lookup table once; keyed for fast joins
  # Age bins: agedays boundaries
  age_breaks <- c(0, 91, 183, 366, 731, 1827, 3653, 5479, Inf)
  age_labels <- c("0-3mo", "3-6mo", "6-12mo", "1-2y", "2-5y", "5-10y", "10-15y", "15-20y")
  # Interval bins: interval_days boundaries
  int_breaks <- c(1, 7, 30, 183, 365, Inf)
  int_labels <- c("<1wk", "1wk-1mo", "1-6mo", "6mo-1y", ">1y")

  # Threshold values: 0 = NR (no rescue), NA = impossible cell (--)
  # Format: list of param -> type -> matrix[age_row, int_col]
  #                           <1wk  1wk-1mo  1-6mo  6mo-1y  >1y
  ht_other <- matrix(c(
    0.40, 0.40,   NA,   NA,   NA,  # 0-3mo
    0.20, 0.40,   NA,   NA,   NA,  # 3-6mo
    0.05, 0.40, 0.00,   NA,   NA,  # 6-12mo
    0.05, 0.40, 0.40,   NA,   NA,  # 1-2y
    0.05, 0.20, 0.40,   NA,   NA,  # 2-5y
    0.05, 0.05, 0.40, 0.40,   NA,  # 5-10y
    0.05, 0.05, 0.20, 0.40, 0.00,  # 10-15y
    0.05, 0.05, 0.05, 0.20, 0.20   # 15-20y
  ), nrow = 8, ncol = 5, byrow = TRUE)

  ht_imperial <- matrix(c(
    NA,   NA,   NA,   NA,   NA,  # 0-3mo (not applicable)
    NA,   NA,   NA,   NA,   NA,  # 3-6mo
    NA,   NA,   NA,   NA,   NA,  # 6-12mo
    NA,   NA,   NA,   NA,   NA,  # 1-2y
    0.05, 0.20, 0.40, 0.40,   NA,  # 2-5y
    0.05, 0.05, 0.40, 0.40,   NA,  # 5-10y
    0.05, 0.05, 0.20, 0.40, 0.00,  # 10-15y
    0.05, 0.05, 0.05, 0.20, 0.20   # 15-20y
  ), nrow = 8, ncol = 5, byrow = TRUE)

  wt_other <- matrix(c(
    0.40, 0.40,   NA,   NA,   NA,  # 0-3mo
    0.20, 0.40, 0.00,   NA,   NA,  # 3-6mo
    0.05, 0.20, 0.00,   NA,   NA,  # 6-12mo
    0.05, 0.20, 0.40,   NA,   NA,  # 1-2y
    0.05, 0.05, 0.40, 0.40,   NA,  # 2-5y
    0.05, 0.05, 0.20, 0.40,   NA,  # 5-10y
    0.05, 0.05, 0.20, 0.40, 0.00,  # 10-15y
    0.05, 0.05, 0.20, 0.20, 0.20   # 15-20y
  ), nrow = 8, ncol = 5, byrow = TRUE)

  wt_imperial <- matrix(c(
    NA,   NA,   NA,   NA,   NA,  # 0-3mo
    NA,   NA,   NA,   NA,   NA,  # 3-6mo
    NA,   NA,   NA,   NA,   NA,  # 6-12mo
    NA,   NA,   NA,   NA,   NA,  # 1-2y
    0.05, 0.05, 0.40, 0.40,   NA,  # 2-5y
    0.05, 0.05, 0.20, 0.40,   NA,  # 5-10y
    0.05, 0.05, 0.20, 0.40, 0.00,  # 10-15y
    0.05, 0.05, 0.20, 0.20, 0.20   # 15-20y
  ), nrow = 8, ncol = 5, byrow = TRUE)

  list(
    age_breaks = age_breaks,
    age_labels = age_labels,
    int_breaks = int_breaks,
    int_labels = int_labels,
    tables = list(
      HEIGHTCM = list(other = ht_other, imperial = ht_imperial),
      HEADCM   = list(other = ht_other, imperial = ht_imperial),  # placeholder: use HT thresholds
      WEIGHTKG = list(other = wt_other, imperial = wt_imperial)
    )
  )
}

# Vectorized CF threshold lookup: given vectors of agedays, interval_days, param, wholehalfimp,
# returns a numeric vector of thresholds (0 = NR, NA = impossible cell)
.cf_get_thresholds <- function(agedays, interval_days, param, wholehalfimp, lookup = NULL) {
  if (is.null(lookup)) lookup <- .cf_rescue_lookup()

  n <- length(agedays)
  thresholds <- rep(NA_real_, n)

  # Bin agedays
  age_idx <- findInterval(agedays, lookup$age_breaks, rightmost.closed = TRUE)
  # Bin interval_days
  int_idx <- findInterval(interval_days, lookup$int_breaks, rightmost.closed = TRUE)

  # Clamp to valid range

  age_idx <- pmin(pmax(age_idx, 1L), length(lookup$age_labels))
  int_idx <- pmin(pmax(int_idx, 1L), length(lookup$int_labels))

  # Determine rounding type: imperial only applies at age > 2y (731+ days)
  rtype <- ifelse(wholehalfimp & agedays >= 731, "imperial", "other")

  # Look up each row
  for (i in seq_len(n)) {
    tbl <- lookup$tables[[param[i]]]
    if (is.null(tbl)) next
    mat <- tbl[[rtype[i]]]
    if (is.null(mat)) next
    thresholds[i] <- mat[age_idx[i], int_idx[i]]
  }

  thresholds
}

# Main child growthcleanr function (cleanchild)

#' Function to clean data (optionally in batches):
#'
#' NOTE: to use multiple processes we process patients in batches using the plyr
#' ddply function. Subjects are split into batches at the wrapper level; all
#' observations for a given subject stay together within one batch. This function
#' processes one batch.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats median embed
#' @noRd
cleanchild <- function(data.df,
                               log.path,
                               quietly,
                               parallel,
                               measurement.to.z,
                               ewma.fields,
                               cf_rescue = "standard",
                               cf_detail = FALSE,
                               biv.z.wt.low.young,
                               biv.z.wt.low.old,
                               biv.z.wt.high,
                               biv.z.ht.low.young,
                               biv.z.ht.low.old,
                               biv.z.ht.high,
                               biv.z.hc.low,
                               biv.z.hc.high,
                               exclude.levels,
                               tanner.ht.vel,
                               error.load.threshold,
                               error.load.mincount,
                               ref.data.path,
                               ewma_window = 15,
                               ref_tables = NULL) {
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

  cf_rescued <- cf_status <- cf_deltaZ <- cf_threshold <- NULL
  otl <- sd_med <- med_diff <- max_diff <- sum_otl <- i.exclude <- NULL


  # avoid no visible warning errors
  sum_sde <- no_sde <- cf <- wholehalfimp <- seq_win <- cs <- absdiff <-
    sd.orig_uncorr <- ageyears <- ctbc.sd <-
    originator_row <- cf_string_id <- ageday_has_include <-
    ..col_replace <- c.ewma.all <- pot_excl <- c.dewma.all <- p_plus <-
    p_minus <- tbc_diff_next <- tbc_diff_prior <-
    tbc_diff_plus_next <- tbc.p_plus <- tbc_diff_plus_prior <-
    tbc_diff_minus_next <- tbc.p_minus <- tbc_diff_minus_prior <- addcrithigh <-
    addcritlow <- tbc_dop <- i.tbc.sd <- rowind <- abssum <-
    whoagegrp_ht <- d_agedays <- mindiff <- maxdiff <- who_mindiff_ht <-
    who_maxdiff_ht <- mindiff_prior <- maxdiff_prior <- whoinc.age.hc <-
    who_maxdiff_hc <- who_mindiff_hc <- diff_prev <-
    diff_next <- aft.g.aftm1 <- val_excl <-
    absval <- comp_diff <- err_ratio <-
    NULL

  # Use internal_id for deterministic SDE order. internal_id is an integer
  # assigned in id-sorted order by cleangrowth(), so results are deterministic
  # regardless of input row order.
  data.df <- data.table(data.df, key = c('subjid', 'param', 'agedays', 'internal_id'))
  # Recreate index for batch processing
  # index was created on the full dataset in cleangrowth() before batching, so
  # batches have non-contiguous indices. This breaks merge/subset operations
  # that use index, so we reassign 1:.N here within the batch.
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
    message(sprintf(
      "[%s] Processing Batch #%s...",
      Sys.time(),
      data.df$batch[1]
    ))

  # NOTE: in each step, we redo temp SDEs

  # EARLY STEP 13: SDE-Identicals (before Steps 5/6) ----
  # Same-day identical values are excluded before CF detection so that CF logic
  # only sees distinct repeated values, not trivial same-day duplicates.
  # Tiebreaking uses age-dependent id preference:
  # At birth (age 0): Keep LOWEST internal_id (earliest, before fluid/interventions)
  # At age > 0: Keep HIGHEST internal_id (consistent with other SDE handling)

  if (!quietly)
    message(sprintf(
      "[%s] Exclude same-day identical measurements (early SDE-Identical)...",
      Sys.time()
    ))

  # SDE-Identical handles partial identicals - duplicates are marked even when
  # mixed with non-duplicates. Groups by subjid/param/agedays/v to find
  # duplicates of each specific value.
  data.df <- data.df[order(subjid, param, agedays, internal_id)]

  # Count included values with same measurement on same day
  data.df[, n_same_value := sum(exclude == "Include"), by = .(subjid, param, agedays, v)]
  data.df[, has_dup := (n_same_value > 1)]
  # Age-dependent: keep lowest internal_id at birth, highest internal_id otherwise
  # Guard against groups with no included rows (empty vector → -Inf/Inf warning)
  data.df[, keep_id := {
    incl_ids <- internal_id[exclude == "Include"]
    if (length(incl_ids) == 0L) NA_integer_
    else if (agedays[1L] == 0L) min(incl_ids)
    else max(incl_ids)
  }, by = .(subjid, param, agedays, v)]
  data.df[has_dup & exclude == "Include" & internal_id != keep_id,
          exclude := .child_exc(param, "Identical")]
  data.df[, c("n_same_value", "has_dup", "keep_id") := NULL]

  # Include internal_id for deterministic SDE order
  setkey(data.df, subjid, param, agedays, internal_id)

  # 5: temporary SDEs ----

  # save a copy of all original measurement values before any transformation
  data.df[, v.orig := v]

  if (!quietly)
    message(sprintf(
      "[%s] Preliminarily identify potential extraneous...",
      Sys.time()
    ))
  data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)]), exclude := 'Exclude-C-Temp-Same-Day']

  # capture a list of subjects with possible extraneous for efficiency later
  subj.dup <- data.df[exclude == 'Exclude-C-Temp-Same-Day', unique(subjid)]

  # Safety check after Step 5 - warn if duplicate Includes found
  step5_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(step5_dupe_check) > 0) {
    warning("Step 5: Found ", nrow(step5_dupe_check),
            " Include duplicates after temp SDE resolution. ",
            "This may indicate a bug in identify_temp_sde().")
  }

  # 6:  carried forwards ----

  # Initialize sde_identical_rows before CF block so it exists when referenced later
  sde_identical_rows <- data.df[0]  # Empty data.table with same structure

  if (!quietly)
    message(sprintf(
    "[%s] Exclude measurements carried forward...",
      Sys.time()
  ))

  # Only run CF detection on valid values with >1 measurement per subject-param
  sp_multi <- data.df[, .(is_multi = .N > 1), by = .(subjid, param)]
  not_single <- sp_multi[data.df, on = .(subjid, param), is_multi]
  valid_set <- .child_valid(data.df, include.temporary.extraneous = TRUE) &
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
    message(sprintf("  CF pre-filter: %d/%d subject-params have potential CFs (%.1f%%)",
                n_with_cf, n_total_sp, 100*n_with_cf/n_total_sp))

  # Initialize cf=FALSE for all, only process those with potential CFs
  data.sub[, cf := FALSE]

  # Only process subject-params with duplicate values
  if (length(sp_with_potential_cf) > 0) {
    cf_subset <- data.sub[sp_key %in% sp_with_potential_cf]
    # Add id for consistent SDE order
    cf_subset <- cf_subset[order(subjid, param, agedays, internal_id)]

    # CF detection: compare each value to the prior ageday's single value.
    # Only compare if the prior ageday has exactly ONE value — if there are
    # multiple values (SDEs) on the prior day, skip the CF comparison.

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

    # Step 6: CF = TRUE if prior day had exactly 1 value AND current matches it.
    # Uses exact equality (no numeric tolerance).
    cf_subset[, cf := !is.na(prior_single_val) & v.orig == prior_single_val]

    # Merge cf results back to data.sub by index
    data.sub[cf_subset, cf := i.cf, on = "index"]
  }
  # End of pre-filter if block - subject-params without potential CFs already have cf=FALSE

  # Cleanup temp columns
  data.sub[, sp_key := NULL]

  # Merge CF flags back to data.df
  cf_idx <- data.sub$index[data.sub$cf]
  data.df[index %in% cf_idx, exclude := .child_exc(param, "CF")]

  # Check if CFs exist before rescue processing.
  # This optimization skips rescue logic if no CFs are present.
  any_cf <- any(data.df$exclude == "Exclude-C-CF")
  if (!quietly)
    message(sprintf("  CF rescue pre-filter: CFs exist = %s", any_cf))

  # Only process CF rescue if CFs exist
  if (any_cf) {

  # Redo temporary SDEs after CF identification: temp SDEs should be
  # re-evaluated now that CFs are excluded.
  any_sde <- any(data.df$exclude == "Exclude-C-Temp-Same-Day")

  if (any_sde) {
    # Temporarily convert Temp SDEs back to Include for re-evaluation.
    data.df[exclude == "Exclude-C-Temp-Same-Day", exclude := "Include"]

    # Re-run temp SDE logic (now CFs are excluded, so SDE evaluation will differ)
    data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)]), exclude := 'Exclude-C-Temp-Same-Day']
  }

  # Determine if measurements are in whole or half imperial units (row-level flag)
  # WEIGHTKG: whole pounds; HEIGHTCM/HEADCM: whole or half inches
  # Tolerance: 0.01 of the imperial unit
  data.df[, wholehalfimp := FALSE]
  data.df[param == "WEIGHTKG",
          wholehalfimp := abs((v.orig * 2.20462262) %% 1) < 0.01]
  data.df[param == "HEIGHTCM",
          wholehalfimp := abs((v.orig / 2.54) %% 0.5) < 0.01]
  data.df[param == "HEADCM",
          wholehalfimp := abs((v.orig / 2.54) %% 0.5) < 0.01]

  # SDE-Identicals break multiple parts of the CF logic (rle, originator
  # detection, cs assignment) so temporarily remove them, do CF calculations,
  # then add them back after CF rescue below.
  sde_identical_rows <- data.df[exclude == "Exclude-C-Identical"]
  data.df <- data.df[!exclude == "Exclude-C-Identical"]

  # CFs on days that also have at least one Include are NOT eligible for
  # rescue and are excluded from string detection. ageday_has_include only
  # checks for Includes (temp SDEs do not block CF rescue eligibility).
  data.df[, ageday_has_include := any(as.character(exclude) == "Include"),
          by = c("subjid", "param", "agedays")]

  # POSITIONAL STRING DETECTION
  # Originators are Includes where the NEXT value is a CF. Strings propagate
  # forward (originator, then consecutive CFs) until interrupted by a non-CF.
  # The ageday_has_include check restricts CF rescue eligibility only; it does
  # not disqualify originators.

  # Initialize variables
  data.df[, cf_binary := exclude == "Exclude-C-CF"]

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

  # Store originator z-scores (uncorrected, pre-GA-correction).
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

  # Create seq_win variable (position in string: 0 for originator, 1/2/3... for CFs)
  data.df[, seq_win := NA_integer_]
  data.df[originator == TRUE, seq_win := 0]
  data.df[!is.na(cf_string_num) & cf_binary == TRUE,
          seq_win := seq_len(.N),
          by = c("subjid", "param", "cf_string_num")]

  # Create cs variable for compatibility with rescue code logic
  data.df[, cs := cf_string_num]

  # Calculate absdiff (absolute z-score difference from originator)
  data.df[!is.na(seq_win), absdiff := abs(sd.orig_uncorr - originator_z)]

  # Clean up temporary variables
  data.df[, c("cf_binary", "nextcf", "priorcf",
              "originator", "originator_seq", "cf_string_num", "originator_z") := NULL]

  # CF rescue: use age/interval/param-specific lookup thresholds
  # CFs with ageday_has_include are already excluded from cs/seq_win assignment
  # (they stay excluded regardless of rescue mode)

  # Compute interval_days for each CF: agedays - prior_ageday
  # We need the originator's ageday for each CF string
  # Reconstruct: originator is seq_win == 0, CFs are seq_win > 0
  # For each CF, the originator's ageday is the ageday of the row with seq_win == 0 in the same cs group
  data.df[, orig_ageday := NA_integer_]
  data.df[seq_win == 0, orig_ageday := agedays]
  # Propagate originator ageday forward within each string (same approach as originator_z)
  if (max_iterations > 0) {
    for (i in 1:max_iterations) {
      data.df[, orig_ageday := ifelse(
        !is.na(cs) & seq_win > 0 & is.na(orig_ageday),
        shift(orig_ageday, type = "lag"),
        orig_ageday
      ), by = c("subjid", "param")]
    }
  }

  # For CFs, interval = agedays - originator agedays
  # (This is the interval from the originator, which is the measurement being compared to)
  data.df[!is.na(seq_win) & seq_win > 0, cf_interval := agedays - orig_ageday]

  if (cf_rescue == "none") {
    # No rescue: all CFs stay excluded
    if (!quietly) message("  CF rescue mode: none (all CFs excluded)")

  } else if (cf_rescue == "all") {
    # Rescue every detected CF, including CFs on a SPA that also has a
    # non-CF Include. This is the "ignore CFs" / legacy-compatible mode:
    # a CF that is consistent with the trajectory should not be
    # preferentially excluded just because another value landed on the
    # same day. Step 13 final SDE resolution handles any resulting
    # multi-Include SPAs on downstream.
    cf_mask <- data.df$exclude == "Exclude-C-CF"
    if (any(cf_mask)) {
      data.df[cf_mask, cf_rescued := "Rescued-All"]
      data.df[cf_mask, exclude := "Include"]
      if (!quietly)
        message(sprintf("  CF rescue mode: all (%d measurements re-included)", sum(cf_mask)))
    }

  } else {
    # Standard rescue: use lookup table thresholds
    cf_lookup <- .cf_rescue_lookup()

    # Get thresholds for all CF rows
    cf_rows <- !is.na(data.df$seq_win) & data.df$seq_win > 0 &
               (data.df$ageday_has_include == FALSE | is.na(data.df$ageday_has_include))

    if (any(cf_rows)) {
      cf_dt <- data.df[cf_rows, .(index, agedays, cf_interval, param, wholehalfimp, absdiff)]
      cf_dt[, cf_threshold := .cf_get_thresholds(
        agedays = agedays,
        interval_days = cf_interval,
        param = as.character(param),
        wholehalfimp = wholehalfimp,
        lookup = cf_lookup
      )]

      # Rescue if absdiff < threshold (threshold > 0)
      # NR cells have threshold = 0: no rescue (absdiff >= 0 is always true)
      # NA threshold (impossible cell): treat as NR (no rescue)
      cf_dt[, rescued := !is.na(cf_threshold) & cf_threshold > 0 & absdiff < cf_threshold]

      # Apply rescues back to data.df
      rescued_idx <- cf_dt[rescued == TRUE, index]
      if (length(rescued_idx) > 0) {
        data.df[index %in% rescued_idx, cf_rescued := "Rescued"]
        data.df[index %in% rescued_idx, exclude := "Include"]
        if (!quietly)
          message(sprintf("  CF rescue: %d measurements re-included (lookup thresholds)",
                      length(rescued_idx)))
      }
    }
  }

  } # End if (any_cf)

  # Add SDE-Identical rows back after CF rescue
  if (nrow(sde_identical_rows) > 0) {
    data.df <- rbind(data.df, sde_identical_rows, fill = TRUE)
    data.df <- data.df[order(subjid, param, agedays, internal_id)]
  }

  # Populate cf_status and cf_deltaZ before dropping temp columns
  if (cf_detail) {
    # Initialize as NA (not a CF candidate)
    data.df[, cf_status := NA_character_]
    data.df[, cf_deltaZ := NA_real_]

    # CF candidates are rows that were flagged as CF (seq_win > 0)
    # They are either still excluded (CF-NR) or rescued (CF-Resc)
    if ("seq_win" %in% names(data.df)) {
      cf_candidate <- !is.na(data.df$seq_win) & data.df$seq_win > 0
      # CF-Resc: rows that were rescued (cf_rescued is non-empty)
      is_rescued <- cf_candidate & !is.na(data.df$cf_rescued) & data.df$cf_rescued != ""
      # CF-NR: CF candidates that were NOT rescued
      is_nr <- cf_candidate & !is_rescued
      data.df[is_rescued, cf_status := "CF-Resc"]
      data.df[is_nr, cf_status := "CF-NR"]
      # deltaZ for all CF candidates
      if ("absdiff" %in% names(data.df)) {
        data.df[cf_candidate, cf_deltaZ := absdiff]
      }
    }
  }

  # Drop columns no longer needed after Step 6
  # v.orig: only used for SDE tiebreaking and CF detection (Steps 5-6)
  # wholehalfimp: only used in CF exclusion logic (Step 6)
  cols_to_drop_6 <- intersect(
    c("v.orig", "wholehalfimp", "seq_win", "cs", "absdiff", "ageday_has_include",
      "orig_ageday", "cf_interval"),
    names(data.df)
  )
  if (length(cols_to_drop_6) > 0L) data.df[, (cols_to_drop_6) := NULL]

  # Step 7: BIV ----


  if (!quietly)
    message(sprintf(
      "[%s] Exclude BIVs...",
      Sys.time()
    ))

  # add age in years
  data.df[, ageyears := agedays/365.25]

  valid_set <- .child_valid(data.df, include.temporary.extraneous = TRUE)

  # identify absolute cutoffs — all use Exclude-C-BIV (not param-specific; param is in the data)
  # Min weight: <0.2 kg for first year, <1 kg after
  data.df[valid_set & param == "WEIGHTKG" & v < 0.2 & agedays <= 365,
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "WEIGHTKG" & v < 1 & agedays > 365,
          exclude := "Exclude-C-BIV"]
  # Max weight at birth
  data.df[valid_set & param == "WEIGHTKG" & v > 10.5 &
            agedays == 0,
          exclude := "Exclude-C-BIV"]
  # Max weight for <2y
  data.df[valid_set & param == "WEIGHTKG" & v > 35 &
            ageyears < 2,
          exclude := "Exclude-C-BIV"]
  # Max weight for all based on published data
  data.df[valid_set & param == "WEIGHTKG" & v > 600,
          exclude := "Exclude-C-BIV"]

  # Min/max HT: 18 is z=-6 for 22 0/7 in Fenton and 65 is z=6 for 40 0/7.
  data.df[valid_set & param == "HEIGHTCM" & v < 18,
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEIGHTCM" & v > 244,
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEIGHTCM" & v > 65 &
            agedays == 0,
          exclude := "Exclude-C-BIV"]

  # Min/max HC: 13 is z=-6 for 22 0/7 in Fenton.
  data.df[valid_set & param == "HEADCM" & v < 13,
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEADCM" & v > 75,
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEADCM" & v > 50 &
            agedays == 0,
          exclude := "Exclude-C-BIV"]

  # Standardized BIV — same exclusion code as absolute BIV (Exclude-C-BIV).
  # The !grepl(biv_pattern, exclude) guard skips rows just assigned
  # Exclude-C-BIV by the absolute block above, since valid_set was
  # computed before absolute BIV ran and so still treats those rows as valid.
  biv_pattern <- "^Exclude-C-BIV$"

  # z cutoffs use unrecentered CSD z-scores (sd.orig_uncorr).
  # Cutoffs are set from cleangrowth() parameters (biv.z.* family).
  # Weight
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < biv.z.wt.low.young &
            ageyears < 1 & !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr < biv.z.wt.low.old &
            ageyears >= 1 & !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "WEIGHTKG" & sd.orig_uncorr > biv.z.wt.high &
            !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]

  # Height. Default upper HT cutoff (biv.z.ht.high = 8) is tighter than the
  # weight upper cutoff based on analysis of CHOP data showing the
  # +/-15 / +/-25 range was too loose for heights.
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < biv.z.ht.low.young &
            ageyears < 1 & !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr < biv.z.ht.low.old &
            ageyears >= 1 & !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEIGHTCM" & sd.orig_uncorr > biv.z.ht.high &
            !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]

  # head circumference
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr < biv.z.hc.low &
            !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]
  data.df[valid_set & param == "HEADCM" & sd.orig_uncorr > biv.z.hc.high &
            !grepl(biv_pattern, exclude),
          exclude := "Exclude-C-BIV"]

  # 7d. Re-evaluate temp SDEs after BIV exclusions: reset all Exclude-C-Temp-Same-Day
  # rows to Include, then rerun identify_temp_sde(). Rationale: absolute or
  # standardized BIV may have excluded the prior temp-SDE "keeper" on an SPA;
  # another value in that SPA should now be flagged instead.
  data.df[exclude == 'Exclude-C-Temp-Same-Day', exclude := 'Include']
  data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)]), exclude := 'Exclude-C-Temp-Same-Day']

  # Drop column no longer needed after Step 7
  # ageyears: only used for BIV age-threshold checks (Step 7)
  data.df[, ageyears := NULL]

  # Step 9: Evil Twins ----
  # Evil Twins identifies pairs (or longer runs) of consecutive measurements
  # that are both extreme relative to the subject's trajectory. Such pairs
  # distort each other's EWMA enough that neither looks out of line on a
  # neighbor-weighted basis, so EWMA-based exclusions alone can miss them.
  # This step runs before EWMA processing to catch them.
  #
  # Per-(subjid, param) group processing: groups are small and independent
  # (typically 3-30 rows), so processing each group in a local copy avoids
  # repeatedly copying the full valid dataset on every while-loop iteration.
  # Same pattern is used in Steps 11, 15, and 17.

  # MUST reorder BEFORE computing valid_set.
  # valid_set is a boolean vector aligned to row positions; reordering
  # after would misalign it, causing temp SDEs to be included and
  # non-temp-SDEs to be excluded.
  # Sort key must include agedays (so calc_otl_evil_twins compares
  # temporally adjacent rows) and internal_id (so same-ageday rows have a
  # deterministic order across sequential and parallel runs).
  data.df <- data.df[order(subjid, param, agedays, internal_id),]

  # Evil Twins requires 3+ measurements per subject-param
  data.df[, sp_count_9 := .N, by = .(subjid, param)]
  not_single_pairs <- data.df$sp_count_9 > 2L
  valid_set <- .child_valid(data.df, include.temporary.extraneous = FALSE) &
    not_single_pairs

  # 9a. Global early-exit check: find out if any possible evil twins exist
  # at all (cheap vectorized check) before entering per-group loop.
  start_df <- calc_otl_evil_twins(data.df[valid_set,])

  if (any(start_df$otl, na.rm = TRUE)) {
    if (!quietly)
      message(sprintf(
      "[%s] Exclude evil twins...",
        Sys.time()
    ))

    # Identify which subject-params have OTL values — only process those groups
    sp_with_otl <- unique(start_df[otl == TRUE, .(subjid, param)])

    if (!quietly)
      message(sprintf("  Evil twins pre-filter: %d subject-params have OTL values",
                  nrow(sp_with_otl)))

    # 9b. Process each subject-param group independently via for loop
    # (Avoids data.table by+:= closure mechanics; groups are small, typically 3-30 rows)
    et_excl_lines <- integer(0)
    for (sp_i in seq_len(nrow(sp_with_otl))) {
      s <- sp_with_otl$subjid[sp_i]
      p <- sp_with_otl$param[sp_i]

      # Extract this group's valid rows
      grp_idx <- which(data.df$subjid == s & data.df$param == p & valid_set)
      if (length(grp_idx) < 2L) next

      df <- data.df[grp_idx, .(line, id, internal_id, subjid, param, agedays, tbc.sd, ctbc.sd, exclude)]
      df <- copy(df)

      # calc_otl for this group
      df <- calc_otl_evil_twins(df)

      while (any(df$otl, na.rm = TRUE)) {
        # 9c. median and distance from median (for Include rows only)
        incl_tbc <- df$tbc.sd[df$exclude == "Include"]
        sd_med <- median(incl_tbc, na.rm = TRUE)
        df[exclude == "Include", med_diff := abs(tbc.sd - sd_med)]

        # Find worst OTL value using tiebreaker hierarchy:
        #   1. Furthest from median (highest med_diff)
        #   2. Most extreme overall (highest abs(tbc.sd))
        #   3. Lowest internal_id (deterministic)
        otl_rows <- df[otl == TRUE]
        if (nrow(otl_rows) == 0L) break
        ord <- order(-otl_rows$med_diff, -abs(otl_rows$tbc.sd), otl_rows$internal_id)
        worst_line <- otl_rows$line[ord[1L]]

        # Mark exactly one exclusion per iteration
        df[line == worst_line, exclude := .child_exc(p, "Evil-Twins")]

        # Recalculate OTL on remaining Include values
        incl <- df[exclude == "Include"]
        if (nrow(incl) < 2L) break
        # Must remove otl column before calling calc_otl_evil_twins, because
        # df[, "otl" := otl] inside that function self-assigns from the existing
        # column (data.table scope resolution) instead of using the local variable.
        if ("otl" %in% names(incl)) incl[, otl := NULL]
        incl <- calc_otl_evil_twins(incl)
        # Map otl back to df (reset all, then set TRUE for OTL rows)
        df[, otl := FALSE]
        if (any(incl$otl)) {
          df[incl[otl == TRUE], otl := TRUE, on = .(line)]
        }
      }

      # Collect lines marked for exclusion
      excl <- df$line[df$exclude == "Exclude-C-Evil-Twins"]
      if (length(excl) > 0L) et_excl_lines <- c(et_excl_lines, excl)
    }

    # Apply all exclusions back to data.df
    if (length(et_excl_lines) > 0L) {
      data.df[line %in% et_excl_lines, exclude := .child_exc(param, "Evil-Twins")]
    }
  }
  data.df[, sp_count_9 := NULL]

  # 9d. Re-evaluate temp SDEs after Evil Twins exclusions: reset all
  # Exclude-C-Temp-Same-Day rows to Include, then rerun identify_temp_sde().
  # Rationale: an Evil Twins exclusion may have removed the prior temp-SDE
  # "keeper" on an SPA, so another value in that SPA should now be flagged.
  data.df[exclude == 'Exclude-C-Temp-Same-Day', exclude := 'Include']
  data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)]), exclude := 'Exclude-C-Temp-Same-Day']

  # Step 11: Extreme EWMA ----
  # Flag extreme outliers using an exponentially weighted moving average
  # (EWMA). For each measurement, the EWMA predicts what the value should
  # be based on the subject's other measurements for the same parameter;
  # values that deviate extremely from the prediction are flagged
  # Exclude-C-Traj-Extreme. Only Include rows participate in the EWMA and
  # as exclusion candidates — temp SDEs do not participate. Because an
  # erroneous measurement can distort the EWMA of its neighbours, at most
  # one value is flagged per subject-param per pass; the global while
  # loop then re-processes any subject-param that produced a new exclusion
  # until no further exclusions occur.

  if (!quietly)
    message(sprintf("[%s] Exclude extreme measurements based on EWMA...", Sys.time()))

  # 11a. Pre-filter / setup -----------------------------------------------
  # Two filters narrow the work; a third captures subjects with existing
  # temp SDEs for later targeted recalculation.

  # Count filter: EWMA needs at least 3 Include values per subject-param.
  data.df[, sp_key := paste0(subjid, "_", param)]
  include_mask <- .child_valid(data.df, include.temporary.extraneous = FALSE)
  include_counts <- data.df[include_mask, .(n_include = .N), by = sp_key]
  sp_to_process <- include_counts[n_include > 2, sp_key]

  # Pre-filter: EWMA1 exclusion requires |tbc.sd| > 3.5.
  # Skip groups where no Include observation has |tbc.sd| > 3.5.
  extreme_sp <- data.df[include_mask & sp_key %in% sp_to_process,
                        .(max_abs_tbc = if (any(!is.na(tbc.sd)))
                            max(abs(tbc.sd), na.rm = TRUE) else 0),
                        by = sp_key]
  sp_to_process <- extreme_sp[max_abs_tbc > 3.5, sp_key]
  if (!quietly)
    message(sprintf("  EWMA1 pre-filter: %d/%d groups have |tbc.sd| > 3.5",
                length(sp_to_process), nrow(extreme_sp)))

  # Track subjects with SDEs for targeted temp SDE recalculation.
  subj_with_sde <- unique(data.df[exclude == "Exclude-C-Temp-Same-Day", subjid])

  # 11b. Iteration loop ---------------------------------------------------
  # Each pass processes one subject-param at a time via a per-group
  # closure operating on a copy of the group's rows. At most one value is
  # flagged per subject-param per iteration. After each pass, temp SDEs
  # are recalculated only for subjects that had a new exclusion and
  # already carried temp SDEs. The loop terminates when no subject-param
  # produces a new exclusion.
  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      message(sprintf("  EWMA1 iteration %d: %d subject-params", iteration, length(sp_to_process)))

    # Track which rows had EWMA1 exclusions before this iteration
    data.df[sp_key %in% sp_to_process, had_ewma1_before := exclude == "Exclude-C-Traj-Extreme"]

    # Process each subject-param - ONE exclusion per subject-param per iteration
    data.df[sp_key %in% sp_to_process,
            exclude := (function(df) {
              # Only Include values participate (no temp SDEs)
              include_set <- .child_valid(df, include.temporary.extraneous = FALSE)

              if (sum(include_set) > 2) {
                # Initialize EWMA fields
                df[, (ewma.fields) := as.double(NaN)]

                # Calculate exp_vals based on age gaps (Include values only).
                # Uses linear interpolation for the exponent.
                inc_ages <- df$agedays[include_set]
                diff_before <- c(NA, diff(inc_ages))
                diff_after <- c(diff(inc_ages), NA)
                maxdiff <- pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)
                ageyears <- maxdiff / 365.25
                exp_vals <- fcase(ageyears <= 1, -1.5, ageyears >= 3, -3.5,
                                  default = -1.5 - (ageyears - 1))
                df[include_set, exp_vals := exp_vals]

                # Calculate EWMA for Include values only (no temp SDEs)
                # Share delta matrix between tbc.sd and ctbc.sd calls (same agedays/exponents)
                ewma_cache <- new.env(parent = emptyenv())
                df[include_set, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE, window = ewma_window, cache_env = ewma_cache)]
                # Skip ctbc EWMA if ctbc.sd == tbc.sd for this group (vast majority)
                if (all(df$ctbc.sd[include_set] == df$tbc.sd[include_set])) {
                  df[include_set, paste0("c.",ewma.fields) := .SD, .SDcols = ewma.fields]
                } else {
                  df[include_set, paste0("c.",ewma.fields) := ewma(agedays, ctbc.sd, exp_vals, TRUE, window = ewma_window, cache_env = ewma_cache)]
                }

                # Calculate dewma for Include values only
                df[include_set, `:=`(
                  dewma.all = tbc.sd - ewma.all,
                  dewma.before = tbc.sd - ewma.before,
                  dewma.after = tbc.sd - ewma.after,
                  c.dewma.all = ctbc.sd - c.ewma.all
                )]

                # Identify potential exclusions - Include values only.
                # c.dewma.all is tested in the same direction as dewma.all
                # (> 3.5 for positive outliers, < -3.5 for negative).
                df[, pot_excl := FALSE]
                df[include_set, pot_excl :=
                     (dewma.all > 3.5 & dewma.before > 3 & dewma.after > 3 & tbc.sd > 3.5 &
                        ((!is.na(ctbc.sd) & c.dewma.all > 3.5) | is.na(ctbc.sd))
                     ) |
                     (dewma.all < -3.5 & dewma.before < -3 & dewma.after < -3 & tbc.sd < -3.5 &
                        ((!is.na(ctbc.sd) & c.dewma.all < -3.5) | is.na(ctbc.sd))
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
                  df[pot_excl == TRUE, exclude := .child_exc(param, "Traj-Extreme")]
                } else if (num.exclude > 1) {
                  # Select worst: highest abs(tbc.sd + dewma.all), lowest internal_id as tiebreaker
                  worst.row <- with(df, order(pot_excl, abs(tbc.sd + dewma.all), -internal_id, decreasing = TRUE))[1]
                  df[worst.row, exclude := .child_exc(param, "Traj-Extreme")]
                }
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('internal_id', 'param', 'agedays', 'tbc.sd', 'ctbc.sd', 'exclude')]

    # Find subject-params with NEW exclusions this iteration
    data.df[sp_key %in% sp_to_process, has_new_excl :=
              exclude == "Exclude-C-Traj-Extreme" & !had_ewma1_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])
    subjects_with_new_excl <- unique(data.df[has_new_excl == TRUE, subjid])

    # Recalculate temp SDEs only for subjects that had exclusions AND have SDEs
    affected_sde_subj <- intersect(subjects_with_new_excl, subj_with_sde)
    if (length(affected_sde_subj) > 0) {
      # Reset temp SDEs to Include for affected subjects
      data.df[subjid %in% affected_sde_subj & exclude == 'Exclude-C-Temp-Same-Day',
              exclude := 'Include']
      # Recalculate temp SDEs for subset
      sde_subset <- data.df[subjid %in% affected_sde_subj]
      sde_result <- identify_temp_sde(sde_subset[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)])
      # Apply results back using index
      sde_indices_to_mark <- sde_subset$index[sde_result]
      data.df[index %in% sde_indices_to_mark, exclude := 'Exclude-C-Temp-Same-Day']
    }

    # Next iteration: only subject-params with new exclusions
    sp_to_process <- sp_with_new_excl

    # Cleanup iteration tracking columns
    data.df[, c("had_ewma1_before", "has_new_excl") := NULL]
  }

  data.df[, sp_key := NULL]

  # 11c. End-of-step temp SDE refresh -------------------------------------
  # Reset all Exclude-C-Temp-Same-Day rows to Include, then rerun
  # identify_temp_sde() across the full dataset. Rationale: the per-
  # iteration targeted recalc only covers subjects that had both an
  # existing temp SDE and a new EWMA1 exclusion in that pass; this final
  # pass catches any residual drift before Step 13. Does NOT use
  # exclude_from_dop_ids — that biasing is specific to Step 13.
  data.df[exclude == 'Exclude-C-Temp-Same-Day', exclude := 'Include']
  data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)]), exclude := 'Exclude-C-Temp-Same-Day']

  # 13: SDEs ----

  if (!quietly)
    message(sprintf(
      "[%s] Exclude same day extraneous...",
      Sys.time()
    ))

  # Pass temp SDE ids to exclude from DOP median.
  # In Step 13 the DOP median uses only fully-included values — temp SDEs are
  # still under consideration and should not bias the cross-parameter anchor.
  temp_sde_ids_step13 <- data.df[exclude == 'Exclude-C-Temp-Same-Day', id]

  # Safety check: no Include duplicates should exist before restoring temp SDEs
  pre_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(pre_dupe_check) > 0) {
    warning("Step 13: Found ", nrow(pre_dupe_check),
            " Include duplicates before restoring temp SDEs.")
  }

  data.df[exclude == 'Exclude-C-Temp-Same-Day', exclude := 'Include']
  data.df[identify_temp_sde(data.df[, .(id, internal_id, subjid, param, agedays, tbc.sd, exclude)],
                                        exclude_from_dop_ids = temp_sde_ids_step13),
          exclude := 'Exclude-C-Temp-Same-Day']

  # Safety check: no Include duplicates after temp SDE marking
  post_dupe_check <- data.df[exclude == "Include", .N, by = .(subjid, param, agedays)][N > 1]
  if (nrow(post_dupe_check) > 0) {
    warning("Step 13: Found ", nrow(post_dupe_check),
            " Include duplicates after temp SDE marking.")
  }
  keep_cols_sde <- names(data.df)

  # Pre-filter to subjects with same-day measurements (potential SDEs)
  # If a subject has no same-day measurements, SDE steps can do nothing - skip entirely
  # Exclude SDE-Identicals from count - they're already handled and shouldn't trigger SDE processing
  sde_day_counts <- data.df[exclude %in% c("Include", "Exclude-C-Temp-Same-Day"),
                            .(n_on_day = .N), by = .(subjid, param, agedays)]
  subj_with_sde_days <- unique(sde_day_counts[n_on_day > 1, subjid])
  n_total_subj <- uniqueN(data.df$subjid)
  n_with_sde <- length(subj_with_sde_days)
  if (!quietly)
    message(sprintf("  SDE pre-filter: %d/%d subjects have same-day measurements (%.1f%%)",
                n_with_sde, n_total_subj, 100*n_with_sde/n_total_subj))

  # Only process subjects with same-day measurements
  if (n_with_sde == 0) {
    # No SDEs to process - skip entire SDE section
    if (!quietly) message("  No SDEs to process, skipping SDE step")
    data.sde <- data.df[0, ]  # Empty data.table with same structure
  } else {
    # Filter to subjects with potential SDEs before dplyr chain
    data.df_sde_subset <- data.df[subjid %in% subj_with_sde_days]

    # Mark groups that have a temp SDE
    data.df_sde_subset[, sde_this := any(exclude == "Exclude-C-Temp-Same-Day"),
                        by = .(subjid, param, agedays)]

    # Filter to valid rows (Include or Temp SDE), then keep only subjects with any SDE
    data.sde <- data.df_sde_subset[.child_valid(data.df_sde_subset, include.temporary.extraneous = TRUE)]
    data.sde[, has_sde_subj := any(sde_this == TRUE), by = .(subjid)]
    data.sde <- data.sde[has_sde_subj == TRUE]
    data.sde[, has_sde_subj := NULL]

    data.sde <- data.sde[order(subjid, param, agedays, internal_id)]

    # Mark SDE-Identical: all values on same day are identical
    # Age 0: keep lowest internal_id, age > 0: keep highest internal_id
    data.sde[, keep_id := {
      ids <- internal_id
      if (agedays[1L] == 0L) min(ids) else max(ids)
    }, by = .(subjid, param, agedays)]
    data.sde[uniqueN(v) == 1L & internal_id != keep_id,
             exclude := .child_exc(param, "Identical"),
             by = .(subjid, param, agedays)]

    # Mark SDE-Identical for duplicate values within same day
    # Age-dependent internal_id selection for duplicate values
    data.sde[, keep_id_dup := {
      ids <- internal_id
      if (agedays[1L] == 0L) min(ids) else max(ids)
    }, by = .(subjid, param, agedays, v)]
    data.sde[, dup_count := .N, by = .(subjid, param, agedays, v)]
    # Mark currently-Included duplicates as SDE-Identical (only Include rows
    # are eligible — other exclusion codes are already handled).
    data.sde[dup_count > 1L & internal_id != keep_id_dup & exclude == "Include",
             exclude := .child_exc(param, "Identical")]
    data.sde[, c("keep_id", "keep_id_dup", "dup_count") := NULL]

    # Include id for deterministic SDE order
    setkey(data.sde, subjid, param, agedays, internal_id)


  # Track temp SDE status before restoration
  data.sde[, was_temp_sde := exclude == 'Exclude-C-Temp-Same-Day']
  data.sde[exclude == 'Exclude-C-Temp-Same-Day', exclude := 'Include']
  data.sde <- data.sde[order(subjid, param, agedays, internal_id)]

  # -------------------------------------------
  # Phase B2: One-Day SDE identification + SDE-All-Extreme
  # -------------------------------------------
  # Flag one-day SDE subject-param pairs
  data.sde[, n_days_with_data := uniqueN(agedays[exclude == "Include"]),
           by = .(subjid, param)]
  data.sde[, has_sde_day := any(was_temp_sde == TRUE), by = .(subjid, param)]
  data.sde[, one_day_sde_flag := (n_days_with_data == 1L & has_sde_day)]

  # Median tbc and distance from median (by subjid, param, agedays)
  data.sde[, median_tbc := as.numeric(median(tbc.sd[exclude == "Include"], na.rm = TRUE)),
           by = .(subjid, param, agedays)]
  data.sde[, absdiff_rel_to_median := fifelse(
    exclude == "Include", abs(tbc.sd - median_tbc), NA_real_)]
  data.sde[, min_absdiff_rel_to_median := suppressWarnings(
    min(absdiff_rel_to_median, na.rm = TRUE)),
    by = .(subjid, param, agedays)]
  data.sde[is.infinite(min_absdiff_rel_to_median),
           min_absdiff_rel_to_median := NA_real_]

  # SDE-All-Extreme: if closest value to median is still > 2 SD away (one-day only)
  data.sde[one_day_sde_flag & !is.na(absdiff_rel_to_median) &
             !is.na(min_absdiff_rel_to_median) &

               min_absdiff_rel_to_median > 2,
           exclude := .child_exc(param, "Extraneous")]

  # DOP medians (cross-parameter grouping by subjid, agedays)
  # DOP median uses FULLY included values only
  data.sde[, HT_dop_med := fifelse(
    one_day_sde_flag & param == "HEIGHTCM",
    as.numeric(median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "WEIGHTKG"],
                      na.rm = TRUE)),
    NA_real_),
    by = .(subjid, agedays)]
  data.sde[, WT_dop_med := fifelse(
    one_day_sde_flag & param == "WEIGHTKG",
    as.numeric(median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "HEIGHTCM"],
                      na.rm = TRUE)),
    NA_real_),
    by = .(subjid, agedays)]
  data.sde[, HC_dop_med := fifelse(
    one_day_sde_flag & param == "HEADCM",
    as.numeric(median(tbc.sd[exclude == "Include" & !was_temp_sde & param == "HEIGHTCM"],
                      na.rm = TRUE)),
    NA_real_),
    by = .(subjid, agedays)]
  data.sde[, absdiff_dop_med := fcase(
    one_day_sde_flag & exclude == "Include" & param == "HEIGHTCM", abs(tbc.sd - HT_dop_med),
    one_day_sde_flag & exclude == "Include" & param == "WEIGHTKG", abs(tbc.sd - WT_dop_med),
    one_day_sde_flag & exclude == "Include" & param == "HEADCM", abs(tbc.sd - HC_dop_med),
    default = NA_real_
  )]

  data.sde[, absdiff_dop_for_sort := fifelse(
    is.na(absdiff_dop_med), Inf, as.numeric(absdiff_dop_med))]

  # Age-dependent id tiebreaker: agedays==0 picks lowest internal_id, otherwise highest
  data.sde[, tiebreaker_oneday := fifelse(
    agedays == 0L, internal_id, -internal_id)]

  # Select one value to keep per SDE group (sort by absdiff_median, absdiff_dop, tiebreaker)
  data.sde[, keep_id_oneday := {
    eligible_mask <- one_day_sde_flag & exclude == "Include"
    if (sum(eligible_mask, na.rm = TRUE) < 2L) NA_integer_
    else {
      eligible_ids <- id[eligible_mask]
      eligible_median <- absdiff_rel_to_median[eligible_mask]
      eligible_dop <- absdiff_dop_for_sort[eligible_mask]
      eligible_tiebreaker <- tiebreaker_oneday[eligible_mask]
      ord <- order(eligible_median, eligible_dop, eligible_tiebreaker)
      eligible_ids[ord[1L]]
    }
  }, by = .(subjid, param, agedays)]

  # Assign SDE-One-Day exclusions
  data.sde[one_day_sde_flag & exclude == "Include" & !is.na(keep_id_oneday) &
             id != keep_id_oneday,
           exclude := .child_exc(param, "Extraneous")]

  # -------------------------------------------
  # Phase B3: SDE-EWMA resolution
  # -------------------------------------------

  # a. Calculate EWMA for fully Included values only (exclude temp SDEs). This
  # ensures (1) only one value per same-day contributes to EWMA, and (2) the
  # diff calculation for the exponent is correct (no same-day 0 diffs).
  ewma_df <- data.sde[exclude == "Include" & !(was_temp_sde)]
  setkey(ewma_df, subjid, param, agedays, internal_id)

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
  ewma_df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_val, TRUE, window = ewma_window), by = .(subjid, param)]

  # Merge EWMA columns back to data.sde
  ewma_merge_cols <- c("id", grep("ewma", names(ewma_df), value = TRUE))
  data.sde <- merge(data.sde, ewma_df[, ..ewma_merge_cols], by = "id", all.x = TRUE)

  # b. Assign EWMAs to temp SDEs (they need EWMA for SDE-EWMA resolution)
  # Compute max ewma.all from non-temp-SDE rows per group
  data.sde[, ewma_fill := suppressWarnings(
    max(ewma.all[!was_temp_sde], na.rm = TRUE)),
    by = .(subjid, param, agedays)]
  data.sde[is.infinite(ewma_fill), ewma_fill := NA_real_]
  data.sde[is.na(ewma.all) & was_temp_sde, ewma.all := ewma_fill]
  data.sde[, ewma_fill := NULL]

  # Use max(ewma) as the SPA (subject-parameter-ageday) reference point.
  # suppressWarnings: max/min on empty subsets produce -Inf/Inf with a warning;
  # converted to NA on the next line. Could use explicit length guards instead
  # (see commit 64c186b for the pattern), but these are in by-group expressions
  # where empty groups are rare and the Inf→NA conversion is correct.
  data.sde[, spa_ewma := suppressWarnings(max(ewma.all, na.rm = TRUE)),
           by = .(subjid, param, agedays)]
  data.sde[is.infinite(spa_ewma), spa_ewma := NA_real_]
  data.sde[, absdewma := abs(tbc.sd - spa_ewma)]

  # SDE-EWMA selection using sort logic
  data.sde[, n_available := sum(
    exclude %in% c("Exclude-C-Temp-Same-Day", "Include")),
    by = .(subjid, param, agedays)]
  data.sde[, min_absdewma := {
    vals <- absdewma[exclude %in% c("Exclude-C-Temp-Same-Day", "Include")]
    m <- suppressWarnings(min(vals, na.rm = TRUE))
    if (is.infinite(m)) NA_real_
    else m
  }, by = .(subjid, param, agedays)]

  # SDE-All-Extreme if min absdewma > 1
  data.sde[n_available >= 2L &
             exclude %in% c("Exclude-C-Temp-Same-Day", "Include") &
             !is.na(min_absdewma) & min_absdewma > 1,
           exclude := .child_exc(param, "Extraneous")]

  # Age-dependent tiebreaker for EWMA selection
  data.sde[, tiebreaker_ewma := fifelse(
    agedays == 0L, internal_id, -internal_id)]

  # Find the id to keep: sort eligible values by absdewma (asc), then tiebreaker (asc)
  data.sde[, keep_id_ewma := {
    eligible_mask <- exclude %in% c("Exclude-C-Temp-Same-Day", "Include")
    if (sum(eligible_mask) == 0L) NA_integer_
    else {
      eligible_ids <- id[eligible_mask]
      eligible_absdewma <- absdewma[eligible_mask]
      eligible_tiebreaker <- tiebreaker_ewma[eligible_mask]
      ord <- order(eligible_absdewma, eligible_tiebreaker)
      eligible_ids[ord[1L]]
    }
  }, by = .(subjid, param, agedays)]

  # Assign SDE-Extraneous exclusions: keep_id gets Include, others get Exclude-C-Extraneous
  data.sde[exclude %in% c("Exclude-C-Temp-Same-Day", "Include") &
             id == keep_id_ewma,
           exclude := "Include"]
  data.sde[exclude %in% c("Exclude-C-Temp-Same-Day", "Include") &
             id != keep_id_ewma,
           exclude := .child_exc(param, "Extraneous")]

  }  # End of SDE pre-filter else block

  # -------------------------------------------
  # Phase B4: Merge SDE results back to main data
  # -------------------------------------------
  sde_results <- data.sde[, .(id, sde_exclude = exclude)]
  data.df <- merge(data.df, sde_results, by = "id", all.x = TRUE)
  # Only overwrite rows that were in scope for SDE processing
  # (Include or Temp-SDE). Do not overwrite rows already carrying
  # permanent exclusion codes assigned before this block (e.g., BIV,
  # Evil Twins, early SDE-Identical).
  data.df[!is.na(sde_exclude) &
            exclude %in% c("Include",
                           "Exclude-C-Temp-Same-Day"),
          exclude := sde_exclude]
  data.df[, sde_exclude := NULL]
  # Keep only original columns (drop any extras from SDE processing)
  extra_cols <- setdiff(names(data.df), keep_cols_sde)
  if (length(extra_cols) > 0L) data.df[, (extra_cols) := NULL]
  setkey(data.df, subjid, param, agedays, internal_id)

  # 15-16: moderate EWMA ----
  # Iteration is global rather than per-(subjid, param): p_plus / p_minus
  # neighbor age-gaps and their z-scores are computed once for all valid rows,
  # one shared while loop runs across the whole batch, and after each
  # iteration only the subject-param groups that produced a new exclusion are
  # re-processed. This avoids re-walking groups that have already converged.

  if (!quietly)
    message(sprintf("[%s] Exclude moderate EWMA...", Sys.time()))

  # Order data for processing
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, internal_id),]

  # Pre-identify subject-params that need processing: Include values, >2 values
  # Include temp SDEs in count so they get p_plus/p_minus calculated
  data.df[, sp_key := paste0(subjid, "_", param)]
  include_counts <- data.df[.child_valid(data.df, include.temporary.extraneous = TRUE),
                            .(n_include = .N), by = sp_key]
  sp_to_process_15 <- include_counts[n_include > 2, sp_key]

  # 15A: Pre-calculate plus/minus values for ALL valid rows
  data.df[, p_plus := NA_real_]
  data.df[, p_minus := NA_real_]
  data.df[param == "WEIGHTKG" & sp_key %in% sp_to_process_15, `:=`(p_plus = 1.05*v, p_minus = 0.95*v)]
  data.df[param == "HEIGHTCM" & sp_key %in% sp_to_process_15, `:=`(p_plus = v+1, p_minus = v-1)]
  data.df[param == "HEADCM" & sp_key %in% sp_to_process_15, `:=`(p_plus = v+1, p_minus = v-1)]

  # 15B: Pre-calculate z-scores for p_plus/p_minus (once for all data)
  # This is expensive, so we do it once upfront
  # Build measurement.to.z_who once here (avoids 2 disk reads inside calc_and_recenter_z_scores)
  measurement.to.z_who_15 <- if (!is.null(ref_tables)) ref_tables$mtz_who_prelim else
    read_anthro(ref.data.path, cdc.only = FALSE)
  valid_for_zscore <- data.df$sp_key %in% sp_to_process_15 & !is.na(data.df$p_plus)
  if (sum(valid_for_zscore) > 0) {
    zscore_subset <- data.df[valid_for_zscore]
    zscore_subset <- calc_and_recenter_z_scores(zscore_subset, "p_plus", ref.data.path,
                                                measurement.to.z, measurement.to.z_who_15)
    zscore_subset <- calc_and_recenter_z_scores(zscore_subset, "p_minus", ref.data.path,
                                                measurement.to.z, measurement.to.z_who_15)
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

  # Get initial subject-params to process
  step15_filter <- data.df$sp_key %in% sp_to_process_15 &
    .child_valid(data.df, include.temporary.extraneous = FALSE) &
    !((data.df$param == "HEIGHTCM" | data.df$param == "HEADCM") & data.df$agedays == 0)
  sp_counts_15 <- data.df[step15_filter, .(n = .N), by = sp_key]
  sp_to_process <- sp_counts_15[n > 2, sp_key]

  # Pre-filter: All EWMA2 rules require addcrit (dewma.before > 1 AND dewma.after > 1).
  # dewma = tbc.sd - EWMA, where EWMA is a weighted average of other obs in the group.
  # |dewma| is bounded by range(tbc.sd). If range <= 1, all |dewma| <= 1,
  # so addcrit can never be TRUE and no exclusion rule can fire.
  tbc_range_15 <- data.df[step15_filter & sp_key %in% sp_to_process,
                           .(tbc_range = if (any(!is.na(tbc.sd)))
                               max(tbc.sd, na.rm = TRUE) - min(tbc.sd, na.rm = TRUE)
                             else 0),
                           by = sp_key]
  n_before_filter <- length(sp_to_process)
  sp_to_process <- tbc_range_15[tbc_range > 1, sp_key]
  if (!quietly)
    message(sprintf("  EWMA2 pre-filter: %d/%d groups have tbc.sd range > 1",
                length(sp_to_process), n_before_filter))

  # Incremental EWMA cache: persists across iterations, keyed by sp_key
  ewma2_caches <- new.env(parent = emptyenv())

  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      message(sprintf("  EWMA2 iteration %d: %d subject-params", iteration, length(sp_to_process)))

    # DOP snapshot: keyed table for fast O(log n) lookups inside group function
    # Refreshed each iteration because excludes change between iterations
    dop_snap <- data.df[exclude == "Include", .(subjid, param, agedays, tbc.sd)]
    setkey(dop_snap, subjid, param)

    # RECALCULATE filter each iteration - valid() depends on current exclusions
    # Must recalculate filter inside loop
    step15_filter <- data.df$sp_key %in% sp_to_process &
      .child_valid(data.df, include.temporary.extraneous = FALSE) &
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
    ewma2_codes <- "Exclude-C-Traj"
    data.df[sp_key %in% sp_to_process, had_ewma2_before := (exclude %in% ewma2_codes)]

    # Process each subject-param - ONE exclusion per subject-param per iteration
    data.df[step15_filter,
            exclude := (function(df) {
              if (nrow(df) < 3) return(df$exclude)

              sp <- paste0(df$subjid[1], "_", df$param[1])

              # Initialize EWMA fields
              df[, (ewma.fields) := as.double(NaN)]

              # Try incremental EWMA update from cache
              used_cache <- FALSE
              if (exists(sp, envir = ewma2_caches)) {
                prev_cache <- get(sp, envir = ewma2_caches)
                excluded_ids <- setdiff(prev_cache$ids, df$id)
                if (length(excluded_ids) == 1L) {
                  cache <- ewma_cache_update(prev_cache, excluded_ids[1L])
                  if (!is.null(cache)) {
                    # Map cache results to df (both in agedays order)
                    df[, `:=`(ewma.all = cache$ewma.all, ewma.before = cache$ewma.before,
                              ewma.after = cache$ewma.after, c.ewma.all = cache$c.ewma.all)]
                    df[, exp_vals := cache$exp_vals]
                    used_cache <- TRUE
                  }
                }
              }

              if (!used_cache) {
                # Full computation (first iteration or cache mismatch)
                df[, c("diff_before", "diff_after") := .(c(NA, diff(agedays)), c(diff(agedays), NA))]
                df[, maxdiff := pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)]
                df[, ageyears := maxdiff / 365.25]
                df[, exp_vals := fcase(ageyears <= 1, -1.5, ageyears >= 3, -3.5, default = -1.5 - (ageyears - 1))]

                cache <- ewma_cache_init(df$agedays, df$tbc.sd, df$ctbc.sd, df$exp_vals, df$id, window = ewma_window)
                df[, `:=`(ewma.all = cache$ewma.all, ewma.before = cache$ewma.before,
                          ewma.after = cache$ewma.after, c.ewma.all = cache$c.ewma.all)]
              }

              # Save cache for next iteration
              assign(sp, cache, envir = ewma2_caches)

              # Calculate dewma
              df[, `:=`(dewma.all = tbc.sd - ewma.all,
                        dewma.before = tbc.sd - ewma.before,
                        dewma.after = tbc.sd - ewma.after,
                        c.dewma.all = ctbc.sd - c.ewma.all)]

              # Calculate prior/next differences
              df[, `:=`(tbc_diff_next = tbc.sd - c(tbc.sd[2:.N], NA),
                        tbc_diff_prior = tbc.sd - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_plus_next = tbc.p_plus - c(tbc.sd[2:.N], NA),
                        tbc_diff_plus_prior = tbc.p_plus - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_minus_next = tbc.p_minus - c(tbc.sd[2:.N], NA),
                        tbc_diff_minus_prior = tbc.p_minus - c(NA, tbc.sd[1:(.N-1)]))]

              # Additional criteria
              df[, addcrithigh := dewma.before > 1 & dewma.after > 1 &
                   ((tbc_diff_next > 1 & tbc_diff_plus_next > 1 & tbc_diff_minus_next > 1) | is.na(tbc_diff_next)) &
                   ((tbc_diff_prior > 1 & tbc_diff_plus_prior > 1 & tbc_diff_minus_prior > 1) | is.na(tbc_diff_prior))]
              df[, addcritlow := dewma.before < -1 & dewma.after < -1 &
                   ((tbc_diff_next < -1 & tbc_diff_plus_next < -1 & tbc_diff_minus_next < -1) | is.na(tbc_diff_next)) &
                   ((tbc_diff_prior < -1 & tbc_diff_plus_prior < -1 & tbc_diff_minus_prior < -1) | is.na(tbc_diff_prior))]

              # DOP lookup (using pre-computed keyed snapshot — O(log n) vs O(n) scan)
              compare_df <- dop_snap[.(df$subjid[1], get_dop(df$param[1]))]
              compare_df <- compare_df[!is.na(tbc.sd)]
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
              df[dewma.all > 1 & (c.dewma.all > 1 | is.na(c.dewma.all)) & addcrithigh & !(rowind %in% c(1, n)),
                 pot_excl := .child_exc(param, "Traj")]
              df[dewma.all < -1 & (c.dewma.all < -1 | is.na(c.dewma.all)) & addcritlow & !(rowind %in% c(1, n)),
                 pot_excl := .child_exc(param, "Traj")]

              # Birth WT
              df[agedays == 0 & c(agedays[2:.N], NA) < 365.25 & dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & c(agedays[2:.N], NA) < 365.25 & dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & c(agedays[2:.N], NA) >= 365.25 & dewma.all > 4 & (c.dewma.all > 4 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & c(agedays[2:.N], NA) >= 365.25 & dewma.all < -4 & (c.dewma.all < -4 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]

              # First
              df[first_meas & (c(agedays[2:.N], NA) - agedays < 365.25) & dewma.all > 2 & (c.dewma.all > 2 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[first_meas & (c(agedays[2:.N], NA) - agedays < 365.25) & dewma.all < -2 & (c.dewma.all < -2 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]

              df[first_meas & (c(agedays[2:.N], NA) - agedays >= 365.25) & dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[first_meas & (c(agedays[2:.N], NA) - agedays >= 365.25) & dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]

              # Last
              if (n >= 2) {
                gap_last <- df$agedays[n] - df$agedays[n-1]
                tbc_prev <- df$tbc.sd[n-1]
                df[rowind == n & gap_last < 365.25*2 & abs(tbc_prev) < 2 & dewma.all > 2 & (c.dewma.all > 2 | is.na(c.dewma.all)) & addcrithigh,
                   pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last < 365.25*2 & abs(tbc_prev) < 2 & dewma.all < -2 & (c.dewma.all < -2 | is.na(c.dewma.all)) & addcritlow,
                   pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last < 365.25*2 & abs(tbc_prev) >= 2 & dewma.all > abs(tbc_prev) & (c.dewma.all > 3 | is.na(c.dewma.all)) & addcrithigh,
                   pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last < 365.25*2 & abs(tbc_prev) >= 2 & dewma.all < -abs(tbc_prev) & (c.dewma.all < -3 | is.na(c.dewma.all)) & addcritlow,
                   pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last >= 365.25*2 & abs(tbc_prev) < 2 & dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) &
                     (tbc.sd - tbc_dop > 4 | is.na(tbc_dop)) & addcrithigh, pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last >= 365.25*2 & abs(tbc_prev) < 2 & dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) &
                     (tbc.sd - tbc_dop < -4 | is.na(tbc_dop)) & addcritlow, pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last >= 365.25*2 & abs(tbc_prev) >= 2 & dewma.all > 1+abs(tbc_prev) & (c.dewma.all > 3 | is.na(c.dewma.all)) &
                     (tbc.sd - tbc_dop > 4 | is.na(tbc_dop)) & addcrithigh, pot_excl := .child_exc(param, "Traj")]
                df[rowind == n & gap_last >= 365.25*2 & abs(tbc_prev) >= 2 & dewma.all < -1-abs(tbc_prev) & (c.dewma.all < -3 | is.na(c.dewma.all)) &
                     (tbc.sd - tbc_dop < -4 | is.na(tbc_dop)) & addcritlow, pot_excl := .child_exc(param, "Traj")]
              }

              # Exclude the worst candidate (highest abssum)
              # All pot_excl candidates are eligible - birth HT/HC is already excluded
              # from Step 15 processing (temporarily excluded before iteration starts)
              candidates <- df[pot_excl != ""]
              if (nrow(candidates) > 0) {
                candidates[, abssum := abs(tbc.sd + dewma.all)]
                # Use internal_id tie-breaker for deterministic selection
                # which.max returns first tie, which depends on data order
                # Use order() with internal_id as tie-breaker for reproducibility
                ord <- order(-candidates$abssum, candidates$internal_id)
                worst_idx <- candidates$index[ord[1]]
                df[index == worst_idx, exclude := pot_excl]
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('index', 'id', 'internal_id', 'subjid', 'param', 'agedays', 'v', 'sex', 'tbc.sd', 'ctbc.sd',
                        'tbc.p_plus', 'tbc.p_minus', 'first_meas', 'exclude')]

    # Find subject-params with NEW exclusions this iteration
    data.df[sp_key %in% sp_to_process, has_new_excl := (exclude %in% ewma2_codes) & !had_ewma2_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])

    # Next iteration: only subject-params with new exclusions
    sp_to_process <- sp_with_new_excl

    # Cleanup iteration tracking
    data.df[, c("had_ewma2_before", "has_new_excl") := NULL]
  }
  rm(ewma2_caches)  # free cache memory

  # Step 16: Birth HT/HC ----

  if (!quietly)
    message(sprintf("[%s] Exclude moderate EWMA for birth height and head circumference...", Sys.time()))

  # Filter for Step 16: HT/HC only, subjects with birth measurement, >2 values
  # Get subjects with birth measurements (HT/HC at agedays=0)
  subj_with_birth <- unique(data.df[param %in% c("HEIGHTCM", "HEADCM") &
                                     .child_valid(data.df, include.temporary.extraneous = FALSE) &
                                     agedays == 0, subjid])

  # Initial filter for counting
  step16_filter <- data.df$param %in% c("HEIGHTCM", "HEADCM") &
    .child_valid(data.df, include.temporary.extraneous = FALSE) &
    data.df$subjid %in% subj_with_birth
  sp_counts_16 <- data.df[step16_filter, .(n = .N), by = sp_key]
  sp_to_process <- sp_counts_16[n > 2, sp_key]

  # Pre-filter: same logic as Step 15 — all rules require addcrit (dewma > 1),
  # and |dewma| is bounded by tbc.sd range. Skip groups with range <= 1.
  tbc_range_16 <- data.df[step16_filter & sp_key %in% sp_to_process,
                           .(tbc_range = if (any(!is.na(tbc.sd)))
                               max(tbc.sd, na.rm = TRUE) - min(tbc.sd, na.rm = TRUE)
                             else 0),
                           by = sp_key]
  n_before_filter_16 <- length(sp_to_process)
  sp_to_process <- tbc_range_16[tbc_range > 1, sp_key]
  if (!quietly)
    message(sprintf("  EWMA2-birth-HT-HC pre-filter: %d/%d groups have tbc.sd range > 1",
                length(sp_to_process), n_before_filter_16))

  # Incremental EWMA cache for Step 16
  ewma2b_caches <- new.env(parent = emptyenv())

  iteration <- 0
  while (length(sp_to_process) > 0) {
    iteration <- iteration + 1
    if (!quietly)
      message(sprintf("  EWMA2-birth-HT-HC iteration %d: %d subject-params", iteration, length(sp_to_process)))

    # RECALCULATE filter each iteration - valid() depends on current exclusions
    # Must recalculate filter inside loop
    step16_filter <- data.df$sp_key %in% sp_to_process &
      data.df$param %in% c("HEIGHTCM", "HEADCM") &
      .child_valid(data.df, include.temporary.extraneous = FALSE) &
      data.df$subjid %in% subj_with_birth

    ewma2_hthc_codes <- "Exclude-C-Traj"
    data.df[sp_key %in% sp_to_process, had_ewma2_before := (exclude %in% ewma2_hthc_codes)]

    data.df[step16_filter,
            exclude := (function(df) {
              if (nrow(df) < 3) return(df$exclude)

              sp <- paste0(df$subjid[1], "_", df$param[1])

              df[, (ewma.fields) := as.double(NaN)]

              # Try incremental EWMA update from cache
              used_cache <- FALSE
              if (exists(sp, envir = ewma2b_caches)) {
                prev_cache <- get(sp, envir = ewma2b_caches)
                excluded_ids <- setdiff(prev_cache$ids, df$id)
                if (length(excluded_ids) == 1L) {
                  cache <- ewma_cache_update(prev_cache, excluded_ids[1L])
                  if (!is.null(cache)) {
                    df[, `:=`(ewma.all = cache$ewma.all, ewma.before = cache$ewma.before,
                              ewma.after = cache$ewma.after, c.ewma.all = cache$c.ewma.all)]
                    df[, exp_vals := cache$exp_vals]
                    used_cache <- TRUE
                  }
                }
              }

              if (!used_cache) {
                df[, c("diff_before", "diff_after") := .(c(NA, diff(agedays)), c(diff(agedays), NA))]
                df[, maxdiff := pmax(abs(diff_before), abs(diff_after), na.rm = TRUE)]
                df[, ageyears := maxdiff / 365.25]
                df[, exp_vals := fcase(ageyears <= 1, -1.5, ageyears >= 3, -3.5, default = -1.5 - (ageyears - 1))]

                cache <- ewma_cache_init(df$agedays, df$tbc.sd, df$ctbc.sd, df$exp_vals, df$id, window = ewma_window)
                df[, `:=`(ewma.all = cache$ewma.all, ewma.before = cache$ewma.before,
                          ewma.after = cache$ewma.after, c.ewma.all = cache$c.ewma.all)]
              }

              # Save cache for next iteration
              assign(sp, cache, envir = ewma2b_caches)

              df[, `:=`(dewma.all = tbc.sd - ewma.all, dewma.before = tbc.sd - ewma.before,
                        dewma.after = tbc.sd - ewma.after, c.dewma.all = ctbc.sd - c.ewma.all)]

              df[, `:=`(tbc_diff_next = tbc.sd - c(tbc.sd[2:.N], NA),
                        tbc_diff_prior = tbc.sd - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_plus_next = tbc.p_plus - c(tbc.sd[2:.N], NA),
                        tbc_diff_plus_prior = tbc.p_plus - c(NA, tbc.sd[1:(.N-1)]))]
              df[, `:=`(tbc_diff_minus_next = tbc.p_minus - c(tbc.sd[2:.N], NA),
                        tbc_diff_minus_prior = tbc.p_minus - c(NA, tbc.sd[1:(.N-1)]))]

              df[, addcrithigh := dewma.before > 1 & dewma.after > 1 &
                   ((tbc_diff_next > 1 & tbc_diff_plus_next > 1 & tbc_diff_minus_next > 1) | is.na(tbc_diff_next)) &
                   ((tbc_diff_prior > 1 & tbc_diff_plus_prior > 1 & tbc_diff_minus_prior > 1) | is.na(tbc_diff_prior))]
              df[, addcritlow := dewma.before < -1 & dewma.after < -1 &
                   ((tbc_diff_next < -1 & tbc_diff_plus_next < -1 & tbc_diff_minus_next < -1) | is.na(tbc_diff_next)) &
                   ((tbc_diff_prior < -1 & tbc_diff_plus_prior < -1 & tbc_diff_minus_prior < -1) | is.na(tbc_diff_prior))]

              df[, pot_excl := ""]
              next_age <- if (nrow(df) > 1) df$agedays[2] else NA

              df[agedays == 0 & !is.na(next_age) & next_age < 365.25 & dewma.all > 3 & (c.dewma.all > 3 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & !is.na(next_age) & next_age < 365.25 & dewma.all < -3 & (c.dewma.all < -3 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & !is.na(next_age) & next_age >= 365.25 & dewma.all > 4 & (c.dewma.all > 4 | is.na(c.dewma.all)) & addcrithigh,
                 pot_excl := .child_exc(param, "Traj")]
              df[agedays == 0 & !is.na(next_age) & next_age >= 365.25 & dewma.all < -4 & (c.dewma.all < -4 | is.na(c.dewma.all)) & addcritlow,
                 pot_excl := .child_exc(param, "Traj")]

              candidates <- df[pot_excl != ""]
              if (nrow(candidates) > 0) {
                candidates[, abssum := abs(tbc.sd + dewma.all)]
                # Use id tie-breaker for deterministic selection
                ord <- order(-candidates$abssum, candidates$internal_id)
                worst_idx <- candidates$index[ord[1]]
                df[index == worst_idx, exclude := pot_excl]
              }

              return(df$exclude)
            })(copy(.SD)),
            by = .(subjid, param),
            .SDcols = c('index', 'id', 'internal_id', 'subjid', 'param', 'agedays', 'v', 'sex', 'tbc.sd', 'ctbc.sd',
                        'tbc.p_plus', 'tbc.p_minus', 'exclude')]

    data.df[sp_key %in% sp_to_process, has_new_excl := (exclude %in% ewma2_hthc_codes) & !had_ewma2_before]
    sp_with_new_excl <- unique(data.df[has_new_excl == TRUE, sp_key])
    sp_to_process <- sp_with_new_excl
    data.df[, c("had_ewma2_before", "has_new_excl") := NULL]
  }
  rm(ewma2b_caches)  # free cache memory

  # Cleanup
  data.df[, sp_key := NULL]

  # Drop columns no longer needed after Step 15
  # p_plus, p_minus, tbc.p_plus, tbc.p_minus: only used in EWMA2 +/-5% rule (Step 15)
  # first_meas: only used in EWMA2 first-measurement exclusion logic (Step 15)
  cols_to_drop_15 <- intersect(c("p_plus", "p_minus", "tbc.p_plus", "tbc.p_minus", "first_meas"),
                                names(data.df))
  if (length(cols_to_drop_15) > 0L) data.df[, (cols_to_drop_15) := NULL]

  # 17: raw differences ----

  if (!quietly)
    message(sprintf(
      "[%s] Exclude raw differences...",
      Sys.time()
    ))

  # tanner.ht.vel was already loaded (same file, same setnames/setkey)
  # in cleangrowth() and is passed as a parameter to this function.
  # Bug fix: was redundantly re-reading from disk each batch.
  tanner.ht.vel.rev <- tanner.ht.vel

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
  who.ht.vel <- merge(who.ht.vel, who.max.ht.vel, by = c('sex', 'whoagegrp_ht'), all = TRUE)

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
  who.hc.vel <- merge(who.hc.vel, who.max.hc.vel, by = c('sex', 'whoagegrp_ht'), all = TRUE)

  setnames(who.hc.vel, colnames(who.hc.vel), gsub('_', '.', colnames(who.hc.vel)))
  setkey(who.hc.vel, sex, whoagegrp.ht)

  # order just for ease later
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, internal_id),]

  # Compute valid_set AFTER sort, not before
  # valid_set is a logical vector indexed by row position, so it must be computed
  # on the sorted data, not before sorting (which would cause index misalignment)
  # create the valid set
  # we only run on valid values, non single values, and non weight
  tmp <- table(paste0(data.df$subjid, "_", data.df$param))
  not_single <- paste0(data.df$subjid, "_", data.df$param) %in% names(tmp)[tmp > 1]
  valid_set <- .child_valid(data.df, include.temporary.extraneous = FALSE) &
    not_single &
    data.df$param != "WEIGHTKG" # do not run on weight

  # ---- Step 17 Pre-filter: identify groups with raw-diff violations ----
  # Compute violation thresholds in one vectorized pass over all valid rows.
  # Only groups with at least one violation need the expensive per-group while loop.
  # Groups without violations are skipped entirely (saves copy(.SD), merges, EWMA, etc.)
  data.df[, sp_key := paste0(subjid, "_", param)]
  {
    pf <- data.df[valid_set, .(sp_key, subjid, param, agedays, id, internal_id, sex, v)]
    pf <- pf[order(subjid, param, agedays, internal_id)]

    # Compute gaps and raw value diffs within each group
    pf[, `:=`(
      d_agedays = shift(agedays, n = 1L, type = "lead") - agedays,
      diff_prev =
        v - shift(v, n = 1L, type = "lag"),
      diff_next =
        shift(v, n = 1L, type = "lead") - v
    ), by = .(subjid, param)]

    # Initialize threshold columns
    pf[, `:=`(mindiff = NA_real_, maxdiff = NA_real_)]

    # ---- HEIGHTCM thresholds ----
    ht_idx <- pf$param == "HEIGHTCM"
    if (any(ht_idx)) {
      # Tanner months (17B)
      pf[ht_idx, tanner.months := 6L + 12L * as.integer(round(
        .5 * (agedays + shift(agedays, n = 1L, type = "lead")) / 365.25)),
        by = .(subjid)]
      pf[ht_idx & (agedays / 30.4375) < 30, tanner.months := NA_integer_]

      # Merge tanner reference (update-on-join preserves row order)
      pf[tanner.ht.vel.rev, `:=`(min.ht.vel = i.min.ht.vel, max.ht.vel = i.max.ht.vel),
         on = .(sex, tanner.months)]

      # Apply max.ht.vel floors (17C-a)
      pf[ht_idx & !is.na(max.ht.vel) & max.ht.vel < 2.54, max.ht.vel := 2.54]
      pf[ht_idx & !is.na(d_agedays) & d_agedays > 2 * 30.4375 & !is.na(max.ht.vel) & max.ht.vel < 2 * 2.54,
         max.ht.vel := 2 * 2.54]
      pf[ht_idx & !is.na(d_agedays) & d_agedays > .5 * 365.25 & !is.na(max.ht.vel) & max.ht.vel < 4 * 2.54,
         max.ht.vel := 4 * 2.54]
      pf[ht_idx & !is.na(d_agedays) & d_agedays > 365.25 & !is.na(max.ht.vel) & max.ht.vel < 8 * 2.54,
         max.ht.vel := 8 * 2.54]

      # Tanner-based mindiff/maxdiff (17C-b)
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 365.25,
         mindiff := .5 * min.ht.vel * (d_agedays / 365.25)^2 - 3]
      pf[ht_idx & !is.na(d_agedays) & d_agedays > 365.25,
         mindiff := .5 * min.ht.vel - 3]
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 365.25,
         maxdiff := 2 * max.ht.vel * (d_agedays / 365.25)^0.33 + 5.5]
      pf[ht_idx & !is.na(d_agedays) & d_agedays > 365.25,
         maxdiff := 2 * max.ht.vel * (d_agedays / 365.25)^1.5 + 5.5]

      # WHO age groups (17D)
      pf[ht_idx, whoagegrp.ht := NA_integer_]
      pf[ht_idx & agedays / 30.4375 <= 24, whoagegrp.ht := as.integer(round(agedays / 30.4375))]
      pf[ht_idx, whoagegrp.ht.lead := shift(whoagegrp.ht, n = 1L, type = "lead"), by = .(subjid)]
      pf[ht_idx & ((!is.na(whoagegrp.ht) & whoagegrp.ht > 24) |
                    (!is.na(whoagegrp.ht.lead) & whoagegrp.ht.lead > 24)),
         whoagegrp.ht := NA_integer_]

      # WHO increment age groups (17E)
      pf[ht_idx, whoinc.age.ht := NA_integer_]
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 20, whoinc.age.ht := 1L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 20 & d_agedays < 46, whoinc.age.ht := 1L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 46 & d_agedays < 76, whoinc.age.ht := 2L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 76 & d_agedays < 107, whoinc.age.ht := 3L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 107 & d_agedays < 153, whoinc.age.ht := 4L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 153 & d_agedays < 199, whoinc.age.ht := 6L]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 200, whoinc.age.ht := 6L]

      # Merge WHO reference data (17F) — update-on-join
      who_ht_cols <- setdiff(colnames(who.ht.vel), c("sex", "whoagegrp.ht"))
      pf[who.ht.vel, (who_ht_cols) := mget(paste0("i.", who_ht_cols)), on = .(sex, whoagegrp.ht)]

      # Extract WHO mindiff/maxdiff based on whoinc.age.ht
      pf[ht_idx & !is.na(whoinc.age.ht) & !is.na(whoagegrp.ht),
         who_mindiff_ht := fcase(
           whoinc.age.ht == 1L, whoinc.1.ht,
           whoinc.age.ht == 2L, whoinc.2.ht,
           whoinc.age.ht == 3L, whoinc.3.ht,
           whoinc.age.ht == 4L, whoinc.4.ht,
           whoinc.age.ht == 6L, whoinc.6.ht,
           default = NA_real_)]
      pf[ht_idx & !is.na(whoinc.age.ht) & !is.na(whoagegrp.ht),
         who_maxdiff_ht := fcase(
           whoinc.age.ht == 1L, max.whoinc.1.ht,
           whoinc.age.ht == 2L, max.whoinc.2.ht,
           whoinc.age.ht == 3L, max.whoinc.3.ht,
           whoinc.age.ht == 4L, max.whoinc.4.ht,
           whoinc.age.ht == 6L, max.whoinc.6.ht,
           default = NA_real_)]

      # Scale WHO values based on d_agedays (17G)
      pf[ht_idx & !is.na(d_agedays) & !is.na(whoinc.age.ht) & d_agedays < whoinc.age.ht * 30.4375,
         who_mindiff_ht := who_mindiff_ht * d_agedays / (whoinc.age.ht * 30.4375)]
      pf[ht_idx & !is.na(d_agedays) & !is.na(whoinc.age.ht) & d_agedays > whoinc.age.ht * 30.4375,
         who_maxdiff_ht := who_maxdiff_ht * d_agedays / (whoinc.age.ht * 30.4375)]
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 9 * 30.4375,
         `:=`(who_mindiff_ht = who_mindiff_ht * .5 - 3,
              who_maxdiff_ht = who_maxdiff_ht * 2 + 3)]

      # Choose WHO vs Tanner (17H)
      # For gap < 9 months: use WHO (already transformed)
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 9 * 30.4375 & !is.na(who_mindiff_ht),
         mindiff := who_mindiff_ht]
      pf[ht_idx & !is.na(d_agedays) & d_agedays < 9 * 30.4375 & !is.na(who_maxdiff_ht),
         maxdiff := who_maxdiff_ht]
      # For gap >= 9 months: Tanner preferred, WHO fallback with transformation
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 9 * 30.4375 & is.na(min.ht.vel) & !is.na(who_mindiff_ht),
         mindiff := who_mindiff_ht * .5 - 3]
      pf[ht_idx & !is.na(d_agedays) & d_agedays >= 9 * 30.4375 & is.na(min.ht.vel) & !is.na(who_maxdiff_ht),
         maxdiff := who_maxdiff_ht * 2 + 3]

      # Fill defaults
      pf[ht_idx & is.na(mindiff), mindiff := -3]
      # Birth adjustments
      pf[ht_idx & agedays == 0, `:=`(mindiff = mindiff - 1.5, maxdiff = maxdiff + 1.5)]
    }

    # ---- HEADCM thresholds ----
    # HC uses WHO velocity tables only (no Tanner); starts at 2-month intervals.
    # Tolerance: ±1.5 cm (vs ±3 cm for HT).
    hc_idx <- pf$param == "HEADCM"
    if (any(hc_idx)) {
      # WHO age group (same floor-month formula as HT)
      pf[, whoagegrp.hc := NA_integer_]
      pf[hc_idx & agedays / 30.4375 <= 24,
         whoagegrp.hc := as.integer(round(agedays / 30.4375))]

      # Interval selection: HC has no 1-month interval; smallest available is 2-month
      pf[, whoinc.age.hc := NA_integer_]
      pf[hc_idx & !is.na(d_agedays) & d_agedays >= 46  & d_agedays < 76,  whoinc.age.hc := 2L]
      pf[hc_idx & !is.na(d_agedays) & d_agedays >= 76  & d_agedays < 107, whoinc.age.hc := 3L]
      pf[hc_idx & !is.na(d_agedays) & d_agedays >= 107 & d_agedays < 153, whoinc.age.hc := 4L]
      pf[hc_idx & !is.na(d_agedays) & d_agedays >= 153 & d_agedays < 200, whoinc.age.hc := 6L]

      # Rename HC velocity columns to avoid collision with already-merged HT columns
      # (both tables use "whoinc.N.ht" naming after gsub; rename HC copies to "whoinc.N.hc")
      who.hc.vel.r <- copy(who.hc.vel)
      ht_suffix_cols <- grep("whoinc", names(who.hc.vel.r), value = TRUE)
      setnames(who.hc.vel.r,
        old = ht_suffix_cols,
        new = gsub("\\.ht$", ".hc", ht_suffix_cols))
      setnames(who.hc.vel.r, "whoagegrp.ht", "whoagegrp.hc")
      setkey(who.hc.vel.r, sex, whoagegrp.hc)

      # Merge HC reference — non-HC rows have whoagegrp.hc = NA so they are not matched
      hc_vel_cols <- setdiff(names(who.hc.vel.r), c("sex", "whoagegrp.hc"))
      pf[who.hc.vel.r, (hc_vel_cols) := mget(paste0("i.", hc_vel_cols)),
         on = .(sex, whoagegrp.hc)]

      # Extract interval-specific WHO mindiff/maxdiff for HC
      pf[hc_idx & !is.na(whoinc.age.hc) & !is.na(whoagegrp.hc),
         who_mindiff_hc := fcase(
           whoinc.age.hc == 2L, whoinc.2.hc,
           whoinc.age.hc == 3L, whoinc.3.hc,
           whoinc.age.hc == 4L, whoinc.4.hc,
           whoinc.age.hc == 6L, whoinc.6.hc,
           default = NA_real_)]
      pf[hc_idx & !is.na(whoinc.age.hc) & !is.na(whoagegrp.hc),
         who_maxdiff_hc := fcase(
           whoinc.age.hc == 2L, max.whoinc.2.hc,
           whoinc.age.hc == 3L, max.whoinc.3.hc,
           whoinc.age.hc == 4L, max.whoinc.4.hc,
           whoinc.age.hc == 6L, max.whoinc.6.hc,
           default = NA_real_)]

      # Scale by d_agedays vs interval length (same logic as HT)
      pf[hc_idx & !is.na(d_agedays) & !is.na(whoinc.age.hc) &
           d_agedays < whoinc.age.hc * 30.4375,
         who_mindiff_hc := who_mindiff_hc * d_agedays / (whoinc.age.hc * 30.4375)]
      pf[hc_idx & !is.na(d_agedays) & !is.na(whoinc.age.hc) &
           d_agedays > whoinc.age.hc * 30.4375,
         who_maxdiff_hc := who_maxdiff_hc * d_agedays / (whoinc.age.hc * 30.4375)]

      # Apply HC tolerance: ±1.5 cm (vs ±3 cm for HT), whenever WHO data
      # exists (not restricted by age).
      pf[hc_idx & !is.na(who_mindiff_hc),
         `:=`(who_mindiff_hc = who_mindiff_hc * .5 - 1.5,
              who_maxdiff_hc = who_maxdiff_hc * 2 + 1.5)]

      # Assign to mindiff/maxdiff (WHO only — no Tanner for HC)
      pf[hc_idx & !is.na(d_agedays) & !is.na(who_mindiff_hc), mindiff := who_mindiff_hc]
      pf[hc_idx & !is.na(d_agedays) & !is.na(who_maxdiff_hc), maxdiff := who_maxdiff_hc]

      # Fallback: HC > 24 months or no WHO data → fixed -1.5.
      pf[hc_idx & is.na(mindiff), mindiff := -1.5]

      # Birth adjustments: ±0.5 cm for HC.
      pf[hc_idx & agedays == 0, `:=`(mindiff = mindiff - .5, maxdiff = maxdiff + .5)]

      # Clean up temporary HC columns
      hc_tmp <- c("whoagegrp.hc", "whoinc.age.hc", "who_mindiff_hc", "who_maxdiff_hc",
                  hc_vel_cols)
      pf[, (intersect(hc_tmp, names(pf))) := NULL]
    }

    # Compute lagged thresholds (previous row's mindiff/maxdiff)
    pf[, `:=`(
      mindiff_prior = shift(mindiff, n = 1L, type = "lag"),
      maxdiff_prior = shift(maxdiff, n = 1L, type = "lag")
    ), by = .(subjid, param)]

    # Check for any raw-diff violations per group
    pf[, has_violation :=
      (!is.na(diff_prev) & !is.na(mindiff_prior) & diff_prev < mindiff_prior) |
      (!is.na(diff_next) & !is.na(mindiff) & diff_next < mindiff) |
      (!is.na(diff_prev) & !is.na(maxdiff_prior) & diff_prev > maxdiff_prior) |
      (!is.na(diff_next) & !is.na(maxdiff) & diff_next > maxdiff)]

    sp17_to_process <- pf[has_violation == TRUE, unique(sp_key)]
    if (!quietly)
      message(sprintf("  Step 17 pre-filter: %d/%d groups have violations",
                  length(sp17_to_process), length(unique(pf$sp_key))))
    rm(pf)
  }

  if (length(sp17_to_process) > 0)
  data.df[valid_set & sp_key %in% sp17_to_process, exclude := (function(df) {
    # work inside a closure to drop added column values

    # save initial exclusions to keep track
    ind_all <- copy(df$index)
    id_all <- copy(df$id)

    exclude_all <- copy(df$exclude)

    # Pre-merge WHO velocity columns BEFORE the while loop (once per group, not per iteration).
    # This matches HT's pattern and avoids column-duplication on iteration 2+ when merge
    # encounters columns already added in iteration 1 and creates .x/.y suffixes.
    if (df$param[1] == "HEIGHTCM") {
      # For HEIGHTCM: pre-compute whoagegrp.ht (static across iterations; agedays don't change)
      df[agedays/30.4375 <= 24, whoagegrp.ht := round(agedays/30.4375)]
      df <- merge(df, who.ht.vel, by = c("sex", "whoagegrp.ht"), all.x = TRUE, sort = FALSE)
    } else if (df$param[1] == "HEADCM") {
      # For HEADCM: pre-compute whoagegrp.hc and pre-merge WHO HC velocity columns.
      # whoinc.2.hc etc. are observation-level constants (based on agedays, not d_agedays)
      # so they only need to be merged once before the while loop.
      who.hc.vel.r <- copy(who.hc.vel)
      hc_suffix_cols <- grep("whoinc|whoagegrp", names(who.hc.vel.r), value = TRUE)
      setnames(who.hc.vel.r, old = hc_suffix_cols,
               new = gsub("\\.ht$", ".hc", hc_suffix_cols))
      df[, whoagegrp.hc := NA_real_]
      df[agedays / 30.4375 <= 24, whoagegrp.hc := round(agedays / 30.4375)]
      df <- merge(df, who.hc.vel.r, by = c("sex", "whoagegrp.hc"), all.x = TRUE, sort = FALSE)
    }

    testing <- TRUE
    iter_count <- 0

    while (testing & nrow(df) > 1){
      iter_count <- iter_count + 1

      # sort df since it got reordered with keys
      # Include id for deterministic SDE order
      df <- df[order(agedays, internal_id),]

      # 17A
      df[, d_agedays := shift(agedays, n = 1L, type = "lead") - agedays]

      # 17E -- only applies to height
      if (df$param[1] == "HEIGHTCM"){
        # 17B
        df[, tanner.months := 6+12*(round(.5*(agedays + shift(agedays, n = 1L, type = "lead"))/365.25))]
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
        df[d_agedays < 365.25, mindiff := .5*min.ht.vel*(d_agedays/365.25)^2-3 ]
        df[d_agedays > 365.25, mindiff := .5*min.ht.vel-3 ]
        # maxdiff exponents: ^0.33 (cube root) for < 365.25 days, ^1.5 for > 365.25.
        df[d_agedays < 365.25,
           maxdiff := 2*max.ht.vel*(d_agedays/365.25)^0.33 + 5.5 ]
        df[d_agedays > 365.25,
           maxdiff := 2*max.ht.vel*(d_agedays/365.25)^1.5 + 5.5 ]

        # 17D: whoagegrp.ht pre-computed outside loop (static for HEIGHTCM)

        # 17E
        df[d_agedays >= 20 & d_agedays < 46, whoinc.age.ht := 1]
        df[d_agedays >= 46 & d_agedays < 76, whoinc.age.ht := 2]
        df[d_agedays >= 76 & d_agedays < 107, whoinc.age.ht := 3]
        df[d_agedays >= 107 & d_agedays < 153, whoinc.age.ht := 4]
        df[d_agedays >= 153 & d_agedays < 199, whoinc.age.ht := 6]


        # Chris updated this to >= from ==
        # update the edge intervals
        df[d_agedays < 20, whoinc.age.ht := 1]
        # Note: d_agedays == 200 is covered by the >= 200 line below
        df[d_agedays >= 200, whoinc.age.ht := 6]
        # Chris updated this to >= from ==
        # 17F: WHO velocity lookup (columns pre-merged outside loop; fcase selects by whoinc.age.ht)
        df[, who_mindiff_ht := fcase(
          !is.na(whoagegrp.ht) & whoinc.age.ht == 1, whoinc.1.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 2, whoinc.2.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 3, whoinc.3.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 4, whoinc.4.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 6, whoinc.6.ht,
          default = NA_real_
        )]
        df[, who_maxdiff_ht := fcase(
          !is.na(whoagegrp.ht) & whoinc.age.ht == 1, max.whoinc.1.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 2, max.whoinc.2.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 3, max.whoinc.3.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 4, max.whoinc.4.ht,
          !is.na(whoagegrp.ht) & whoinc.age.ht == 6, max.whoinc.6.ht,
          default = NA_real_
        )]

        # 17G: Scale WHO velocity limits by actual gap relative to reference interval
        df[d_agedays < whoinc.age.ht*30.4375, who_mindiff_ht :=
             who_mindiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        df[d_agedays > whoinc.age.ht*30.4375, who_maxdiff_ht :=
             who_maxdiff_ht * d_agedays/(whoinc.age.ht*30.4375)]
        # For gaps < 9 months: widen limits (allow more loss, allow more gain)
        df[d_agedays < 9*30.4375, who_mindiff_ht := who_mindiff_ht*.5-3]
        df[d_agedays < 9*30.4375, who_maxdiff_ht := who_maxdiff_ht*2+3]

        # 17H: Choose between Tanner and WHO velocity limits
        # For gap < 9 months: use WHO (already transformed at lines 4030-4031)
        df[d_agedays < 9*30.4375 & !is.na(who_mindiff_ht), mindiff := who_mindiff_ht]
        df[d_agedays < 9*30.4375 & !is.na(who_maxdiff_ht), maxdiff := who_maxdiff_ht]
        # Apply transformation for >= 9 month gaps.
        # When gap >= 9 months, Tanner is preferred. But if Tanner is not
        # available (min.ht.vel is NA), use WHO as fallback WITH the
        # transformation (same scaling as the < 9 month path).
        df[d_agedays >= 9*30.4375 & is.na(min.ht.vel) & !is.na(who_mindiff_ht),
           mindiff := who_mindiff_ht*.5-3]
        df[d_agedays >= 9*30.4375 & is.na(min.ht.vel) & !is.na(who_maxdiff_ht),
           maxdiff := who_maxdiff_ht*2+3]

        # Fallback mindiff default
        df[is.na(mindiff), mindiff := -3]
        # for birth measurements, add allowance of 1.5cm
        df[agedays == 0, mindiff := mindiff - 1.5]
        df[agedays == 0, maxdiff := maxdiff + 1.5]

        # 17I
        # sort df since it got reordered with keys
        # Include id for deterministic SDE order
        df <- df[order(agedays, internal_id),]
        df[, mindiff_prior := shift(mindiff, n = 1L, type = "lag")]
        df[, maxdiff_prior := shift(maxdiff, n = 1L, type = "lag")]
      } else { # head circumference
        # WHO HC velocity columns (whoinc.2.hc, max.whoinc.2.hc, etc.) were pre-merged
        # before the while loop. Here we compute only per-iteration logic (depends on d_agedays).

        # 17J: Interval selection — HC starts at 2-month (no 1-month interval).
        df[, whoinc.age.hc := NA_integer_]
        df[!is.na(d_agedays) & d_agedays >= 46  & d_agedays < 76,  whoinc.age.hc := 2L]
        df[!is.na(d_agedays) & d_agedays >= 76  & d_agedays < 107, whoinc.age.hc := 3L]
        df[!is.na(d_agedays) & d_agedays >= 107 & d_agedays < 153, whoinc.age.hc := 4L]
        df[!is.na(d_agedays) & d_agedays >= 153 & d_agedays < 200, whoinc.age.hc := 6L]

        # 17K: Extract interval-specific WHO mindiff/maxdiff using renamed columns
        df[, who_mindiff_hc := fcase(
          whoinc.age.hc == 2L, whoinc.2.hc,
          whoinc.age.hc == 3L, whoinc.3.hc,
          whoinc.age.hc == 4L, whoinc.4.hc,
          whoinc.age.hc == 6L, whoinc.6.hc,
          default = NA_real_)]
        df[, who_maxdiff_hc := fcase(
          whoinc.age.hc == 2L, max.whoinc.2.hc,
          whoinc.age.hc == 3L, max.whoinc.3.hc,
          whoinc.age.hc == 4L, max.whoinc.4.hc,
          whoinc.age.hc == 6L, max.whoinc.6.hc,
          default = NA_real_)]

        # 17L: Scale by actual gap vs. reference interval length (same logic as HT)
        df[!is.na(d_agedays) & !is.na(whoinc.age.hc) &
             d_agedays < whoinc.age.hc * 30.4375,
           who_mindiff_hc := who_mindiff_hc * d_agedays / (whoinc.age.hc * 30.4375)]
        df[!is.na(d_agedays) & !is.na(whoinc.age.hc) &
             d_agedays > whoinc.age.hc * 30.4375,
           who_maxdiff_hc := who_maxdiff_hc * d_agedays / (whoinc.age.hc * 30.4375)]

        # 17M: Apply HC tolerance ±1.5 cm whenever WHO data exists.
        df[!is.na(who_mindiff_hc),
           `:=`(who_mindiff_hc = who_mindiff_hc * .5 - 1.5,
                who_maxdiff_hc = who_maxdiff_hc * 2 + 1.5)]

        df[, mindiff := who_mindiff_hc]
        df[, maxdiff := who_maxdiff_hc]
        # Fallback: HC > 24 months or no WHO data → fixed -1.5.
        df[is.na(mindiff), mindiff := -1.5]

        # Birth adjustments: ±0.5 cm for HC.
        df[agedays == 0, `:=`(mindiff = mindiff - .5, maxdiff = maxdiff + .5)]

        # 17N: Sort and lag thresholds for pairwise violation check
        # Include id for deterministic SDE order
        df <- df[order(agedays, internal_id),]
        df[, mindiff_prior := shift(mindiff, n = 1L, type = "lag")]
        df[, maxdiff_prior := shift(maxdiff, n = 1L, type = "lag")]
      }

      # 17O: generate ewma
      df[, (ewma.fields) := as.double(NaN)]

      # Calculate exponent based on max age gap to nearest neighbor
      tmp <- data.frame(
        "before" = abs(df$agedays - c(NA, df$agedays[1:(nrow(df)-1)])),
        "after" = abs(df$agedays - c(df$agedays[2:(nrow(df))], NA))
      )
      maxdiff_e <- pmax(tmp$before, tmp$after, na.rm = TRUE)
      exp_vals <- rep(-1.5, nrow(tmp))
      maxdiff_years <- maxdiff_e / 365.25
      exp_vals[maxdiff_years > 1 & maxdiff_years < 3] <-
        -1.5 - (maxdiff_years[maxdiff_years > 1 & maxdiff_years < 3] - 1)
      exp_vals[maxdiff_years >= 3] <- -3.5
      df[, exp_vals := exp_vals]

      # calculate ewma
      df[, (ewma.fields) := ewma(agedays, tbc.sd, exp_vals, TRUE, window = ewma_window)]

      # calculate dewma
      df[, `:=`(
        dewma.all = tbc.sd - ewma.all,
        dewma.before = tbc.sd - ewma.before,
        dewma.after = tbc.sd - ewma.after
      )]

      # add differences for convenience
      df[, diff_prev := v-shift(v, n = 1L, type = "lag")]
      df[, diff_next := shift(v, n = 1L, type = "lead")-v]



      if (nrow(df) > 2){
        # 17P/R: identify pairs and calculate exclusions
        df[, pair := diff_prev < mindiff_prior |
             diff_next < mindiff |
             diff_prev > maxdiff_prior |
             diff_next > maxdiff
        ]
        df[is.na(pair), pair := FALSE]
        # Reset tie-breaker flags each iteration (they persist from prior)
        df[, bef.g.aftm1 := NA]
        df[, aft.g.aftm1 := NA]
        df[(pair & shift(pair, n = 1L, type = "lag")) & abs(dewma.before) > shift(abs(dewma.after), n = 1L, type = "lag"),
           bef.g.aftm1 := TRUE]
        df[(pair & shift(pair, n = 1L, type = "lead")) & abs(dewma.after) > shift(abs(dewma.before), n = 1L, type = "lead"),
           aft.g.aftm1 := TRUE]

        # Q
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prior & bef.g.aftm1, val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_next < mindiff & aft.g.aftm1, val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_prev > maxdiff_prior & bef.g.aftm1, val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_next > maxdiff & aft.g.aftm1, val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_prev < mindiff_prior & bef.g.aftm1, val_excl_code := "1"]
        df[diff_next < mindiff & aft.g.aftm1, val_excl_code := "2"]
        df[diff_prev > maxdiff_prior & bef.g.aftm1, val_excl_code := "3"]
        df[diff_next > maxdiff & aft.g.aftm1, val_excl_code := "4"]
      }else { # only 2 values
        # 17Q/R -- exclusions for pairs
        df[, val_excl := exclude]
        df[diff_prev < mindiff_prior & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lag"),
           val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_next < mindiff & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lead"),
           val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_prev > maxdiff_prior & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lag"),
           val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_next > maxdiff & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lead"),
           val_excl := .child_exc(param, "Abs-Diff")]
        df[diff_prev < mindiff_prior & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lag"),
           val_excl_code := "5"]
        df[diff_next < mindiff & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lead"),
           val_excl_code := "6"]
        df[diff_prev > maxdiff_prior & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lag"),
           val_excl_code := "7"]
        df[diff_next > maxdiff & abs(tbc.sd) > shift(abs(tbc.sd), n = 1L, type = "lead"),
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
        ord <- order(-candidates$absval, candidates$internal_id)
        idx <- candidates$index[ord[1]]


        exclude_all[ind_all == idx] <- df[index == idx, val_excl]

        # Continue iterating — after excluding 1 candidate, new violations may emerge.
        # count_exclude > 0 outer guard makes count_exclude >= 1 always true here.
        testing <- TRUE
        df <- df[index != idx, ]
      } else {
        testing <- FALSE
      }
    }

    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param),
  .SDcols = c('index', 'id', 'internal_id', 'agedays', 'sex', 'param', 'v', 'tbc.sd', 'exclude')]

  data.df[, sp_key := NULL]

  # 19: 1 or 2 measurements ----
  if (!quietly)
    message(sprintf(
      "[%s] Exclude 1 or 2 measurements...",
      Sys.time()
    ))

  # order just for ease later
  # Add id for consistent SDE order
  data.df <- data.df[order(subjid, param, agedays, internal_id),]

  # Compute valid_set AFTER sort, not before
  # Create the valid set: only Include rows contribute to the singles/pairs
  # count (so subjects whose excluded measurements leave them with 1 or 2
  # remaining Includes still flow into this step).
  include_df <- data.df[exclude == "Include"]
  sp_counts <- include_df[, .(sp_n = .N), by = .(subjid, param)]
  only_single_pairs <- sp_counts[data.df, on = .(subjid, param), fifelse(is.na(sp_n), FALSE, sp_n <= 2L)]
  valid_set <- .child_valid(data.df, include.temporary.extraneous = FALSE) &
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
    dop <- dop_snapshot[.(df$subjid[1], get_dop(df$param[1]))]

    # 19D: calculate the voi comparison
    if (nrow(dop) > 0){
      for (i in seq_len(nrow(df))){
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

      # 19E: which is larger
      # Use id tie-breaker for deterministic selection
      max_ind <- if (!all(is.na(df$comp_diff))){
        ord <- order(-df$comp_diff, df$internal_id)
        ord[1]
      } else {
        ord <- order(-abs(df$tbc.sd), df$internal_id)
        ord[1]
      }

      # 19F/G: which exclusion. Use absolute differences for the z-score
      # threshold comparisons.
      # NA-safe: diff_tbc.sd can be NA if z-scores are NA (edge cases).
      if (isTRUE(abs(diff_tbc.sd) > 4 &
          (abs(diff_ctbc.sd) > 4 | is.na(diff_ctbc.sd)) &
          diff_agedays >=365.25)){
        df[max_ind, exclude := .child_exc(param, "Pair")]
      } else if (isTRUE(abs(diff_tbc.sd) > 2.5 &
                 (abs(diff_ctbc.sd) > 2.5 | is.na(diff_ctbc.sd)) &
                 diff_agedays < 365.25)){
        df[max_ind, exclude := .child_exc(param, "Pair")]
      }

      # save the results
      exclude_all <- df$exclude

      # 19H
      # if one needs to get removed, we want to reevaluate as a single
      df <- df[exclude == "Include",]
    }

    if (nrow(df) == 1){
      # Check if 1-meas exclusion applies
      # NA-safe: tbc.sd can be NA in edge cases
      one_meas_cond <- isTRUE(
        (abs(df$tbc.sd) > 3 & !is.na(df$comp_diff) & df$comp_diff > 5) |
        (abs(df$tbc.sd) > 5 & is.na(df$comp_diff))
      )
      if (one_meas_cond) {
        df[, exclude := .child_exc(param, "Single")]
        # Only update exclude_all if 1-meas exclusion happens
        # Previous code at line 4529 unconditionally overwrote exclude_all, losing 2-meas results
        exclude_all[ind_all == df$index] <- .child_exc(df$param, "Single")
      }
    }

    return(exclude_all)
  })(copy(.SD)), by = .(subjid, param),
  .SDcols = c('index', 'id', 'internal_id', 'agedays', 'subjid', 'param', 'tbc.sd', 'ctbc.sd', 'exclude')]

  # 21: error load ----

  if (!quietly)
    message(sprintf(
      "[%s] Exclude error load...",
      Sys.time()
    ))

  # Error-load ratio: errors / (errors + includes). SDE codes, CF codes,
  # Missing, and Not-Cleaned rows are excluded from both the numerator and
  # the denominator (see non_error_codes below).
  valid_set <- rep(TRUE, nrow(data.df))

  # Non-error codes that should be excluded from both numerator AND denominator
  # CF rescue codes removed — rescued CFs are now "Include" (stored in cf_rescued column)
  non_error_codes <- c("Exclude-C-Identical",
                       "Exclude-C-Extraneous",
                       "Exclude-C-CF",
                       "Exclude-Missing",
                       "Exclude-Not-Cleaned")

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

  # error.load.mincount guards against tripping the step on a subject with
  # a small number of errors (default: need at least 2 errors before the
  # ratio is considered).
  data.df[valid_set & err_ratio > error.load.threshold & n_errors >= error.load.mincount & exclude == "Include",
          exclude := .child_exc(param, "Too-Many-Errors")]

  # Cleanup
  data.df[, c("err_ratio", "n_errors") := NULL]


  # end ----

  if (!quietly)
    message(sprintf("[%s] Completed Batch #%s...", Sys.time(), data.df$batch[1]))
  if (!quietly & parallel)
    sink()

  # Return appropriate z-score based on potcorr status
  # For potcorr subjects, return ctbc.sd (corrected); for others, return tbc.sd (uncorrected)
  # First check if we have potcorr and ctbc.sd columns
  has_potcorr <- "potcorr" %in% colnames(data.df)
  has_ctbc <- "ctbc.sd" %in% colnames(data.df)

  if (has_potcorr && has_ctbc) {
    # Create final_tbc column: use ctbc.sd for potcorr subjects, tbc.sd for others
    data.df[, final_tbc := ifelse(!is.na(potcorr) & potcorr == TRUE & !is.na(ctbc.sd), ctbc.sd, tbc.sd)]
  } else {
    # Fallback to tbc.sd if columns not available
    data.df[, final_tbc := tbc.sd]
  }

  # Assemble return columns: essentials + z-scores + final_tbc + optional
  # cf_detail columns.
  return_cols <- c("id", "line", "exclude", "param", "cf_rescued")

  # Add z-score columns if they exist
  zscore_cols <- c("sd.orig_who", "sd.orig_cdc", "sd.orig", "tbc.sd", "ctbc.sd")
  for (col in zscore_cols) {
    if (col %in% names(data.df)) {
      return_cols <- c(return_cols, col)
    }
  }

  # Also add final_tbc for reference
  return_cols <- c(return_cols, "final_tbc")

  # Add cf_detail columns if requested
  if (cf_detail) {
    for (col in c("cf_status", "cf_deltaZ")) {
      if (col %in% names(data.df)) {
        return_cols <- c(return_cols, col)
      }
    }
  }

  return(data.df[, ..return_cols])
}

# Supporting pediatric growthcleanr functions
# Supporting functions for pediatric piece of algorithm

#' Identify rows currently eligible for child-algorithm processing.
#'
#' Returns a logical vector the same length as the input, marking which
#' rows are currently "valid" for a given step of `cleanchild()`. A row
#' is valid if its exclusion code does not start with `"Exclude-"`, so
#' `"Include"` rows pass unconditionally and every `Exclude-*` code is
#' filtered out by default (including `Exclude-Missing`,
#' `Exclude-Not-Cleaned`, and all algorithm-assigned `Exclude-C-*`
#' codes).
#'
#' Three optional flags add back specific soft-excluded categories so a
#' step can opt them in without changing the default for the rest of
#' the algorithm. `include.temporary.extraneous` also admits rows
#' flagged `"Exclude-C-Temp-Same-Day"` by Child Step 5;
#' `include.extraneous` also admits rows flagged `"Exclude-C-Extraneous"`
#' by Child Step 13; `include.carryforward` also admits rows flagged
#' `"Exclude-C-CF"` by Child Step 6. The flags are additive.
#'
#' Accepts either a data.frame / data.table (in which case `df$exclude`
#' is used) or a plain character / factor vector of exclusion codes.
#' The exclusion column is coerced to character so factor inputs and
#' character inputs match identically under `^Exclude` pattern
#' matching.
#'
#' @param df data.frame / data.table containing an `exclude` column, or
#'   a character or factor vector of exclusion codes.
#' @param include.temporary.extraneous if TRUE, also keep rows with
#'   `exclude == "Exclude-C-Temp-Same-Day"`. Defaults to FALSE.
#' @param include.extraneous if TRUE, also keep rows with
#'   `exclude == "Exclude-C-Extraneous"`. Defaults to FALSE.
#' @param include.carryforward if TRUE, also keep rows with
#'   `exclude == "Exclude-C-CF"`. Defaults to FALSE.
#'
#' @return Logical vector the same length as the input, TRUE for rows
#'   that the calling step should process.
#'
#' @keywords internal
#' @noRd
.child_valid <- function(df,
                  include.temporary.extraneous = FALSE,
                  include.extraneous = FALSE,
                  include.carryforward = FALSE) {
  exclude <- if (is.data.frame(df)) df$exclude else df
  exclude <- as.character(exclude)

  # Base set: non-excluded rows. All exclusion codes start with "Exclude-"
  # (including Exclude-Missing and Exclude-Not-Cleaned), so this single
  # grepl catches everything. "Include" does not match, so included rows pass.
  keep <- !grepl("^Exclude", exclude)

  if (include.temporary.extraneous)
    keep <- keep | exclude == "Exclude-C-Temp-Same-Day"
  if (include.extraneous)
    keep <- keep | exclude == "Exclude-C-Extraneous"
  if (include.carryforward)
    keep <- keep | exclude == "Exclude-C-CF"

  return(keep)
}

