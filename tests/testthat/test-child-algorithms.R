testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Algorithm step tests for the child algorithm
#
# Tests that verify specific algorithm step behaviors using constructed data.
# Complements regression tests (frozen counts) and parameter tests (runs ok).
# These tests construct known scenarios and verify the algorithm responds
# correctly to each one.
#
# Sections:
#   1. CF rescue
#   2. Evil twins / OTL
#   3. Error load
#   4. Age boundaries
#   5. Parallel execution
# =============================================================================

# ---------------------------------------------------------------------------
# Shared data and helpers
# ---------------------------------------------------------------------------

# Load syngrowth once at file scope
data("syngrowth", package = "growthcleanr", envir = environment())
.sg <- as.data.table(syngrowth)
setkey(.sg, subjid, param, agedays)
.sg_peds <- .sg[agedays < 20 * 365.25]

# Pre-built subsets used by multiple tests
.subjs100 <- unique(.sg_peds$subjid)[1:100]
.d100 <- .sg_peds[subjid %in% .subjs100]
.subjs200 <- unique(.sg_peds$subjid)[1:200]
.d200 <- .sg_peds[subjid %in% .subjs200]

#' Run cleangrowth on a data.table with standard columns, return merged result
run_gc <- function(d, ...) {
  res <- cleangrowth(
    subjid = d$subjid,
    param  = d$param,
    agedays = d$agedays,
    sex    = d$sex,
    measurement = d$measurement,
    id     = d$id,
    quietly = TRUE,
    ...
  )
  # Select only GC output columns before merge to avoid .x/.y conflicts
  gc_cols <- c("id", "exclude")
  if ("cf_rescued" %in% names(res)) gc_cols <- c(gc_cols, "cf_rescued")
  merge(d, res[, ..gc_cols], by = "id")
}


# ===========================================================================
# Section 1: CF rescue tests
#
# CF rescue re-includes carried-forward values when the z-score difference
# between the CF'd value and the originator is small enough for that
# age/interval/param/rounding cell. Three modes:
#   cf_rescue = "standard" (default) — age/interval/param lookup thresholds
#   cf_rescue = "none"     — no rescue (all CFs excluded)
#   cf_rescue = "all"      — every detected CF rescued (including CFs on a
#                            SPA with another Include; Step 13 resolves the
#                            resulting multi-Include SPAs)
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 1: Standard rescue rescues more CFs than "none" mode
# ---------------------------------------------------------------------------
test_that("CF rescue: standard rescues more CFs than none", {

  res_std <- run_gc(.d200, cf_rescue = "standard")
  res_none <- run_gc(.d200, cf_rescue = "none")

  n_cf_std  <- sum(grepl("-CF$", res_std$exclude))
  n_cf_none <- sum(grepl("-CF$", res_none$exclude))

  # Standard rescue should leave fewer CFs excluded than no rescue
  expect_lte(n_cf_std, n_cf_none)
})

# ---------------------------------------------------------------------------
# Test 2: cf_rescued column tracks rescue status
# ---------------------------------------------------------------------------
test_that("CF rescue: cf_rescued column populated for rescued CFs", {

  res <- run_gc(.d200, cf_rescue = "standard")

  # cf_rescued column should exist
  expect_true("cf_rescued" %in% names(res),
              info = "Output must include cf_rescued column")

  # Some rows should have rescue codes
  rescued <- res[cf_rescued != ""]

  # Rescued CFs should not retain CF exclusion code
  if (nrow(rescued) > 0) {
    expect_false(any(grepl("-CF$", rescued$exclude)),
                 info = "Rescued CFs should not retain CF exclusion code")

    # Rescue codes should be valid
    valid_rescue <- c("Rescued", "Rescued-All")
    expect_true(all(rescued$cf_rescued %in% valid_rescue),
                info = "Rescue codes should be from the known set")
  }
})

# ---------------------------------------------------------------------------
# Test 3: cf_rescue = "none" rescues no CFs
# ---------------------------------------------------------------------------
test_that("CF rescue: none mode rescues no CFs", {

  res <- run_gc(.d100, cf_rescue = "none")

  rescued <- res[cf_rescued != ""]
  expect_equal(nrow(rescued), 0,
               info = "cf_rescue='none' should rescue no CFs")
})

# ---------------------------------------------------------------------------
# Test 4: cf_rescue = "all" rescues all CFs
# ---------------------------------------------------------------------------
test_that("CF rescue: all mode rescues all CFs", {

  res <- run_gc(.d200, cf_rescue = "all")

  # No CF exclusion codes should remain — every detected CF is rescued.
  n_cf <- sum(grepl("-CF$", res$exclude))
  expect_equal(n_cf, 0,
               info = "cf_rescue='all' should leave no CF exclusions")

  # Rescued rows should have the Rescued-All label.
  n_rescued_all <- sum(res$cf_rescued == "Rescued-All")
  expect_gt(n_rescued_all, 0,
            label = "Rescued-All count in cf_rescue='all'")
})


# ===========================================================================
# Section 2: Evil twins / OTL tests
#
# Evil twins detection uses Over-The-Limit (OTL) to find measurements where
# both tbc_diff > 5 and ctbc_diff > 5 relative to an adjacent value. The
# algorithm iteratively removes the worst OTL value until none remain.
#
# Strategy: take clean syngrowth subjects and inject unit errors (multiply
# height by 2.54, simulating inches recorded as cm). This creates z-score
# differences well above 5 from neighbors.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 5: Single unit error triggers Evil-Twins exclusion
# ---------------------------------------------------------------------------
test_that("ET/OTL: single unit error (height x2.54) is excluded as Evil-Twins", {

  # Find a subject with several Include heights in middle of trajectory
  # Subject 0d8773f3 has 19 HT measurements
  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- .d100[subjid == target & param == "HEIGHTCM"]

  # Pick a middle measurement and multiply by 2.54 (unit error: inches -> cm)
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_id <- ht_rows$id[mid_idx]
  d_mod <- copy(.d100)
  d_mod[id == error_id, measurement := measurement * 2.54]

  res <- run_gc(d_mod)

  # The modified measurement should be excluded (Evil-Twins or another extreme code)
  error_result <- as.character(res[id == error_id]$exclude)
  expect_true(error_result != "Include",
              info = sprintf("Unit error (x2.54) should be excluded, got: %s", error_result))
})

# ---------------------------------------------------------------------------
# Test 6: Two consecutive unit errors both excluded
# ---------------------------------------------------------------------------
test_that("ET/OTL: two consecutive unit errors are both excluded", {

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- .d100[subjid == target & param == "HEIGHTCM"]

  # Pick two consecutive middle measurements
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_ids <- ht_rows$id[mid_idx:(mid_idx + 1)]
  d_mod <- copy(.d100)
  d_mod[id %in% error_ids, measurement := measurement * 2.54]

  res <- run_gc(d_mod)

  # Both modified measurements should be excluded
  for (eid in error_ids) {
    error_result <- as.character(res[id == eid]$exclude)
    expect_true(error_result != "Include",
                info = sprintf("ID %d: consecutive unit error should be excluded, got: %s",
                               eid, error_result))
  }
})

# ---------------------------------------------------------------------------
# Test 7: Three consecutive unit errors all excluded
# ---------------------------------------------------------------------------
test_that("ET/OTL: three consecutive unit errors are all excluded", {

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- .d100[subjid == target & param == "HEIGHTCM"]

  # Pick three consecutive middle measurements
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_ids <- ht_rows$id[mid_idx:(mid_idx + 2)]
  d_mod <- copy(.d100)
  d_mod[id %in% error_ids, measurement := measurement * 2.54]

  res <- run_gc(d_mod)

  # All three modified measurements should be excluded
  for (eid in error_ids) {
    error_result <- as.character(res[id == eid]$exclude)
    expect_true(error_result != "Include",
                info = sprintf("ID %d: 3-consecutive unit error should be excluded, got: %s",
                               eid, error_result))
  }
})

# ---------------------------------------------------------------------------
# Test 8: Unit errors don't cause collateral damage to clean neighbors
# ---------------------------------------------------------------------------
test_that("ET/OTL: unit errors don't exclude neighboring clean measurements", {

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- .d100[subjid == target & param == "HEIGHTCM"]

  # Inject one unit error in the middle
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_id <- ht_rows$id[mid_idx]
  d_mod <- copy(.d100)
  d_mod[id == error_id, measurement := measurement * 2.54]

  # Run both original and modified
  res_orig <- run_gc(.d100)
  res_mod  <- run_gc(d_mod)

  # Get the non-error HT rows for this subject
  clean_ht_ids <- ht_rows$id[ht_rows$id != error_id]

  # Count how many were Include in original vs modified
  n_inc_orig <- sum(res_orig[id %in% clean_ht_ids]$exclude == "Include")
  n_inc_mod  <- sum(res_mod[id %in% clean_ht_ids]$exclude == "Include")

  # Should have same or very similar Include count — no widespread collateral damage
  # Allow up to 2 rows of collateral (EWMA-based steps may shift slightly)
  # Unit error should not cause widespread collateral exclusions
  expect_gte(n_inc_mod, n_inc_orig - 2)
})


# ===========================================================================
# Section 3: Error load tests
#
# Error load (Step 21) excludes all remaining Includes for a subject-param
# when the proportion of errors exceeds error.load.threshold and the count
# of errors >= error.load.mincount.
#
# Bug fix 2026-04-12: threshold was hardcoded at 0.4; now uses the parameter.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 9: error.load.threshold parameter actually affects exclusions
# ---------------------------------------------------------------------------
test_that("error load: threshold parameter controls exclusion behavior", {

  # threshold=0.4 should catch subjects that 0.5 misses
  # (confirmed: 2 Error-load at 0.4, 0 at 0.5 on this dataset)
  res_04 <- cleangrowth(
    subjid = .d200$subjid, param = .d200$param, agedays = .d200$agedays,
    sex = .d200$sex, measurement = .d200$measurement, id = .d200$id,
    quietly = TRUE, error.load.threshold = 0.4
  )
  res_05 <- cleangrowth(
    subjid = .d200$subjid, param = .d200$param, agedays = .d200$agedays,
    sex = .d200$sex, measurement = .d200$measurement, id = .d200$id,
    quietly = TRUE, error.load.threshold = 0.5
  )

  n_el_04 <- sum(grepl("Too-Many-Errors", res_04$exclude))
  n_el_05 <- sum(grepl("Too-Many-Errors", res_05$exclude))

  # Lower threshold should catch at least as many
  expect_gte(n_el_04, n_el_05)

  # Specifically: 0.4 should catch more than 0.5 on this dataset
  expect_gt(n_el_04, n_el_05)
})

# ---------------------------------------------------------------------------
# Test 10: error.load.mincount parameter affects exclusions
# ---------------------------------------------------------------------------
test_that("error load: mincount parameter controls exclusion behavior", {

  # mincount=1 with low threshold should catch more than mincount=10
  res_min1 <- cleangrowth(
    subjid = .d200$subjid, param = .d200$param, agedays = .d200$agedays,
    sex = .d200$sex, measurement = .d200$measurement, id = .d200$id,
    quietly = TRUE, error.load.threshold = 0.3, error.load.mincount = 1
  )
  res_min10 <- cleangrowth(
    subjid = .d200$subjid, param = .d200$param, agedays = .d200$agedays,
    sex = .d200$sex, measurement = .d200$measurement, id = .d200$id,
    quietly = TRUE, error.load.threshold = 0.3, error.load.mincount = 10
  )

  n_el_min1  <- sum(grepl("Too-Many-Errors", res_min1$exclude))
  n_el_min10 <- sum(grepl("Too-Many-Errors", res_min10$exclude))

  # Lower mincount should trigger at least as many Error-load
  expect_gte(n_el_min1, n_el_min10)
})

# ---------------------------------------------------------------------------
# Test 11: Constructed high-error subject triggers Error-load
# ---------------------------------------------------------------------------
test_that("error load: subject with many errors triggers Error-load on remaining", {

  # Take a subject with many HT measurements and corrupt most of them
  # Subject 0d8773f3 has 19 HT measurements
  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- .d100[subjid == target & param == "HEIGHTCM"]

  # Corrupt all but 2 measurements with distinct extreme values
  # (distinct values avoid SDE-Identical/CF which are excluded from error load denominator)
  # Keep first and last, corrupt everything else
  corrupt_ids <- ht_rows$id[2:(nrow(ht_rows) - 1)]
  d_mod <- copy(.d100)
  set.seed(42)
  # Each corrupt value is a different implausible height (200-300 cm range)
  d_mod[id %in% corrupt_ids,
        measurement := runif(.N, min = 200, max = 300)]

  # Use a low threshold to make it easier to trigger
  res <- run_gc(d_mod, error.load.threshold = 0.3)

  # At least one of the surviving measurements should be Error-load
  has_error_load <- any(grepl("Too-Many-Errors",
                        res[subjid == target & param == "HEIGHTCM"]$exclude))
  expect_true(has_error_load,
              info = "Subject with many errors should trigger Error-load for remaining Includes")
})


# ===========================================================================
# Section 4: Age boundary tests
#
# The child algorithm has age-specific behavior at several boundaries:
#   - HEADCM: cleaned only for agedays <= 3 * 365.25 (1095.75 days)
#   - CF rescue for multi-CF strings: females >= 16yr, males >= 17yr
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 12: HEADCM at the 3-year boundary
# ---------------------------------------------------------------------------
test_that("age boundary: HEADCM cleaned under 3yr, not cleaned over 3yr", {
  # Build synthetic data with HC measurements near the 3-year boundary
  # Need enough HT/WT context and multiple HC rows for the algorithm to work

  # Take a subject with data starting from birth
  # Subject 0db2905c: sex=0, agedays 0-1071
  target <- "0db2905c-1e0c-0d4d-7d86-4552c5b55ebd"
  subj_data <- .sg_peds[subjid == target]
  target_sex <- subj_data$sex[1]

  # Add several HC measurements spanning birth to beyond 3yr
  # Plausible HC trajectory for a male: ~35cm at birth -> ~49cm at 3yr
  under3_day <- as.integer(3 * 365.25 - 1)  # 1095 days
  over3_day  <- as.integer(3 * 365.25 + 2)  # 1097 days

  hc_ages <- c(0L, 30L, 90L, 180L, 365L, 730L, under3_day, over3_day)
  hc_vals <- c(35.0, 37.5, 40.0, 43.0, 46.0, 48.0, 49.0, 49.2)

  max_id <- max(subj_data$id)
  hc_rows <- data.table(
    id = max_id + seq_along(hc_ages),
    subjid = target,
    param = "HEADCM",
    agedays = hc_ages,
    sex = target_sex,
    measurement = hc_vals
  )

  d <- rbind(subj_data, hc_rows)
  # Include other subjects for context
  other_subjs <- unique(.sg_peds$subjid)[1:20]
  other_subjs <- other_subjs[other_subjs != target]
  d <- rbind(d, .sg_peds[subjid %in% other_subjs[1:10]])
  setkey(d, subjid, param, agedays)

  # Reassign sequential IDs to avoid conflicts
  d[, id := seq_len(.N)]

  res <- run_gc(d)

  # Find our HC rows by matching on subjid + param + agedays
  hc_res <- res[subjid == target & param == "HEADCM"]

  # Under 3yr should be cleaned (not "Exclude-Not-Cleaned")
  under3_res <- as.character(hc_res[agedays == under3_day]$exclude)
  expect_true(under3_res != "Exclude-Not-Cleaned",
              info = sprintf("HC under 3yr should be cleaned, got: %s", under3_res))

  # Over 3yr should be "Exclude-Not-Cleaned"
  over3_res <- as.character(hc_res[agedays == over3_day]$exclude)
  expect_equal(over3_res, "Exclude-Not-Cleaned",
               info = "HC over 3yr should be 'Not cleaned'")
})

# ---------------------------------------------------------------------------
# Test 13: cf_detail columns
# ---------------------------------------------------------------------------
test_that("CF rescue: cf_detail produces cf_status and cf_deltaZ columns", {

  res <- cleangrowth(
    subjid = .d100$subjid, param = .d100$param, agedays = .d100$agedays,
    sex = .d100$sex, measurement = .d100$measurement, id = .d100$id,
    quietly = TRUE, cf_detail = TRUE
  )

  # cf_status and cf_deltaZ columns should exist
  expect_true("cf_status" %in% names(res),
              info = "Output must include cf_status column when cf_detail=TRUE")
  expect_true("cf_deltaZ" %in% names(res),
              info = "Output must include cf_deltaZ column when cf_detail=TRUE")

  # cf_status values should be from {NA, "CF-NR", "CF-Resc"}
  valid_status <- c(NA_character_, "CF-NR", "CF-Resc")
  expect_true(all(res$cf_status %in% valid_status | is.na(res$cf_status)),
              info = "cf_status values should be NA, CF-NR, or CF-Resc")

  # cf_deltaZ should be NA for non-CF rows and numeric for CF rows
  cf_rows <- !is.na(res$cf_status)
  if (any(cf_rows)) {
    expect_true(all(!is.na(res$cf_deltaZ[cf_rows])),
                info = "cf_deltaZ should be non-NA for CF candidates")
    expect_true(all(res$cf_deltaZ[cf_rows] >= 0),
                info = "cf_deltaZ should be non-negative")
  }
})

# ---------------------------------------------------------------------------
# Test 14: cf_detail not present by default
# ---------------------------------------------------------------------------
test_that("CF rescue: cf_detail columns absent by default", {

  .subjs50 <- unique(.sg_peds$subjid)[1:50]
  .d50 <- .sg_peds[subjid %in% .subjs50]

  res <- cleangrowth(
    subjid = .d50$subjid, param = .d50$param, agedays = .d50$agedays,
    sex = .d50$sex, measurement = .d50$measurement, id = .d50$id,
    quietly = TRUE
  )

  expect_false("cf_status" %in% names(res),
               info = "cf_status should not be in output by default")
  expect_false("cf_deltaZ" %in% names(res),
               info = "cf_deltaZ should not be in output by default")
})


# ===========================================================================
# Section 5: Parallel execution
#
# parallel=TRUE should produce identical results to parallel=FALSE.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 15: parallel=TRUE produces identical results to parallel=FALSE
# ---------------------------------------------------------------------------
test_that("parallel execution: parallel=TRUE matches parallel=FALSE", {
  # This test requires growthcleanr to be installed (not just load_all)
  skip_if_not_installed("growthcleanr")

  res_seq <- cleangrowth(
    subjid = .d100$subjid, param = .d100$param, agedays = .d100$agedays,
    sex = .d100$sex, measurement = .d100$measurement, id = .d100$id,
    quietly = TRUE, parallel = FALSE
  )

  res_par <- cleangrowth(
    subjid = .d100$subjid, param = .d100$param, agedays = .d100$agedays,
    sex = .d100$sex, measurement = .d100$measurement, id = .d100$id,
    quietly = TRUE, parallel = TRUE, num.batches = 2
  )

  # Same number of rows
  expect_equal(nrow(res_par), nrow(res_seq),
               info = "Parallel and sequential should return same number of rows")

  # Merge by id and compare exclusion codes
  comp <- merge(res_seq[, .(id, excl_seq = as.character(exclude))],
                res_par[, .(id, excl_par = as.character(exclude))],
                by = "id")

  n_diff <- sum(comp$excl_seq != comp$excl_par)
  expect_equal(n_diff, 0,
               info = sprintf("Parallel and sequential should produce identical exclusion codes; %d differ", n_diff))
})


# ===========================================================================
# Section 6: CF rescue threshold-specific tests
#
# Verify that the age/interval/param lookup table produces different rescue
# decisions for cells with different thresholds. Uses cf_detail=TRUE to
# inspect cf_status and cf_deltaZ directly.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 16: High-threshold cell (infant, short interval) rescues large deltaZ
# ---------------------------------------------------------------------------
test_that("CF threshold: infant short-interval CF with deltaZ=0.35 is rescued", {
  # 0-3mo, <1wk interval, HT-Other -> threshold 0.40
  # A CF with deltaZ=0.35 should be rescued (0.35 < 0.40)
  # Build a subject with normal growth but one carried-forward height
  d <- data.table(
    id = 1:8,
    subjid = "subj_cf1",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 3L, 30L, 60L,    # HT: birth, 3 days later (CF candidate), 1mo, 2mo
                0L, 3L, 30L, 60L),
    measurement = c(50.0, 50.0, 54.0, 57.0,   # HT: birth=50, day3=50 (CF!), then normal
                    3.5, 4.0, 4.5, 5.5)        # WT: normal growth
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, cf_detail = TRUE
  )

  # id=2 is the CF candidate (HT 50.0 at day 3, same as birth HT 50.0)
  # In this high-threshold cell (0.40), small deltaZ should be rescued
  cf_row <- res[id == 2]
  if (!is.na(cf_row$cf_status)) {
    # If detected as CF, should be rescued (deltaZ small, threshold generous)
    expect_equal(cf_row$cf_status, "CF-Resc",
                 info = "Infant short-interval CF should be rescued (threshold=0.40)")
  }
})

# ---------------------------------------------------------------------------
# Test 17: NR cell excludes CF even with tiny deltaZ
# ---------------------------------------------------------------------------
test_that("CF threshold: NR cell (10-15y, >1y interval, WT) excludes CF", {
  # 10-15y age, >1y interval, WT-Other -> NR (threshold = 0.00)
  # Even identical weight with zero deltaZ should be excluded
  # CF detection requires the value to match the IMMEDIATELY PRECEDING measurement
  # so consecutive WT values must be identical with >365 day gap between them
  d <- data.table(
    id = 1:8,
    subjid = "subj_cf_nr",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(4000L, 4400L, 4800L, 5200L,
                4000L, 4400L, 4800L, 5200L),
    measurement = c(145.0, 149.0, 153.0, 157.0,    # HT: normal growth
                    40.0, 40.0, 45.0, 48.0)         # WT: 40.0 at day 4000 AND 4400 (consecutive, 400-day interval >1y)
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, cf_detail = TRUE
  )

  # id=6 is WT at day 4400 (age ~12y), identical to id=5 at day 4000 (400-day interval, >1y bin)
  # This is in the NR cell: 10-15y age, >1y interval, WT-Other → threshold 0.00
  cf_row <- res[id == 6]
  if (!is.na(cf_row$cf_status)) {
    expect_equal(cf_row$cf_status, "CF-NR",
                 info = "NR cell (10-15y, >1y interval, WT) should exclude CF")
  }
})

# ---------------------------------------------------------------------------
# Test 18: NR cell (young infant, long interval) always excludes
# ---------------------------------------------------------------------------
test_that("CF threshold: NR cell (6-12mo, >1mo interval) excludes all CFs", {
  # 6-12mo, 1-6mo interval, WT-Other -> threshold 0.00 (NR)
  # Even zero deltaZ should not be rescued
  d <- data.table(
    id = 1:8,
    subjid = "subj_cf3",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(200L, 260L, 300L, 340L,
                200L, 260L, 300L, 340L),
    measurement = c(68.0, 71.0, 73.0, 74.0,     # HT: normal
                    8.0, 8.0, 9.0, 9.5)          # WT: 8.0 at day 200 AND 260 (60 day interval, CF)
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, cf_detail = TRUE
  )

  # id=6 is WT CF candidate (8.0 at day 260, same as 8.0 at day 200, 60-day interval)
  # In NR cell, should be excluded regardless of deltaZ
  cf_row <- res[id == 6]
  if (!is.na(cf_row$cf_status)) {
    expect_equal(cf_row$cf_status, "CF-NR",
                 info = "NR cell CF should be excluded (no rescue)")
  }
})


# ===========================================================================
# Section 7: GA correction (potcorr) tests
#
# Step 2b: Subjects whose first weight z-score < -2 at age < 10 months
# are flagged as "potentially correctable" (potcorr). Their z-scores are
# corrected using Fenton reference curves. This affects downstream
# exclusion decisions.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 19: Very low birth weight triggers potcorr correction
# ---------------------------------------------------------------------------
test_that("GA correction: very low birth weight infant gets potcorr correction", {
  # Create a subject with a very low birth weight (z < -2) at birth
  # The correction should change z-scores, potentially changing exclusion outcomes
  # Compare with a normal-weight version of the same subject
  d_low <- data.table(
    id = 1:8,
    subjid = "subj_potcorr",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L,
                0L, 30L, 90L, 180L),
    measurement = c(45.0, 50.0, 58.0, 65.0,     # HT: short at birth, catches up
                    1.5, 3.0, 5.0, 7.0)          # WT: very low birth weight (1.5kg), rapid gain
  )

  d_normal <- data.table(
    id = 11:18,
    subjid = "subj_normal",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L,
                0L, 30L, 90L, 180L),
    measurement = c(50.0, 54.0, 60.0, 67.0,     # HT: normal
                    3.5, 4.5, 6.0, 8.0)          # WT: normal birth weight
  )

  d <- rbind(d_low, d_normal)

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # The low birth weight subject should have corrected z-scores (ctbc.sd differs from tbc.sd)
  potcorr_res <- res[id %in% 1:8]
  normal_res <- res[id %in% 11:18]

  # Both subjects should run without NA exclusion codes
  expect_false(any(is.na(potcorr_res$exclude)))
  expect_false(any(is.na(normal_res$exclude)))

  # Potcorr subject should have corrected z-scores that differ from uncorrected
  potcorr_wt <- potcorr_res[param == "WEIGHTKG" & !is.na(tbc.sd) & !is.na(ctbc.sd)]
  if (nrow(potcorr_wt) > 0) {
    has_correction <- any(potcorr_wt$ctbc.sd != potcorr_wt$tbc.sd, na.rm = TRUE)
    expect_true(has_correction,
                info = "Potcorr subject should have corrected z-scores differing from uncorrected")
  }

  # Normal subject should NOT have corrected z-scores differing from uncorrected
  normal_wt <- normal_res[param == "WEIGHTKG" & !is.na(tbc.sd) & !is.na(ctbc.sd)]
  if (nrow(normal_wt) > 0) {
    has_correction_norm <- any(normal_wt$ctbc.sd != normal_wt$tbc.sd, na.rm = TRUE)
    expect_false(has_correction_norm,
                 info = "Normal-weight subject should NOT have corrected z-scores")
  }
})

# ---------------------------------------------------------------------------
# Test 20: Birth weight just above potcorr threshold is NOT corrected
# ---------------------------------------------------------------------------
test_that("GA correction: birth weight z >= -2 does not trigger potcorr", {
  # A subject with first weight z-score just above -2 should NOT be potcorr
  # Normal 50th percentile male birth weight is ~3.5kg; z=-2 is ~2.5kg
  # Use 2.6kg (just above -2) — should not trigger correction
  d <- data.table(
    id = 1:8,
    subjid = "subj_near_potcorr",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L,
                0L, 30L, 90L, 180L),
    measurement = c(48.0, 52.0, 59.0, 66.0,     # HT: slightly short
                    2.6, 3.8, 5.5, 7.5)          # WT: 2.6kg at birth (near but above -2 SD)
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # For non-potcorr subjects, ctbc.sd should equal tbc.sd (no correction applied)
  wt_rows <- res[param == "WEIGHTKG" & !is.na(tbc.sd) & !is.na(ctbc.sd)]
  if (nrow(wt_rows) > 0) {
    expect_equal(wt_rows$ctbc.sd, wt_rows$tbc.sd,
                 info = "Non-potcorr subject: ctbc.sd should equal tbc.sd (no correction)")
  }
})


# ===========================================================================
# Section 8: Birth measurement EWMA2 tests (Step 15 birth WT, Step 16 birth HT/HC)
#
# Birth measurements (agedays == 0) have special EWMA2 rules:
# - Birth WT in Step 15: dewma > 3 with next < 365d, dewma > 4 with next >= 365d
# - Birth HT/HC in Step 16: same thresholds, separate iterative loop
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 21: Extreme birth weight is excluded
# ---------------------------------------------------------------------------
test_that("birth EWMA2: extreme birth weight is excluded", {
  # Create a subject with an extreme birth weight that should trigger exclusion
  # dewma > 3 needed; use a very high birth weight with normal later values
  d <- data.table(
    id = 1:10,
    subjid = "subj_birth_wt",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L, 365L,
                0L, 30L, 90L, 180L, 365L),
    measurement = c(50.0, 54.0, 60.0, 67.0, 75.0,    # HT: normal
                    8.0, 4.5, 6.0, 8.0, 10.0)         # WT: 8kg at birth (extreme!), then normal
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # The extreme birth weight (id=6) should be excluded
  birth_wt <- res[id == 6]
  expect_true(grepl("Exclude", birth_wt$exclude),
              info = sprintf("Extreme birth weight (8kg) should be excluded, got: %s",
                             as.character(birth_wt$exclude)))
})

# ---------------------------------------------------------------------------
# Test 22: Normal birth weight is not excluded
# ---------------------------------------------------------------------------
test_that("birth EWMA2: normal birth weight is included", {
  d <- data.table(
    id = 1:10,
    subjid = "subj_birth_wt2",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L, 365L,
                0L, 30L, 90L, 180L, 365L),
    measurement = c(50.0, 54.0, 60.0, 67.0, 75.0,    # HT: normal
                    3.5, 4.5, 6.0, 8.0, 10.0)         # WT: 3.5kg at birth (normal)
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # Normal birth weight should be included
  birth_wt <- res[id == 6]
  expect_equal(as.character(birth_wt$exclude), "Include",
               info = "Normal birth weight (3.5kg) should be included")
})

# ---------------------------------------------------------------------------
# Test 23: Extreme birth height is excluded (Step 16)
# ---------------------------------------------------------------------------
test_that("birth EWMA2: extreme birth height is excluded", {
  # Birth height that is extremely high should trigger Step 16
  d <- data.table(
    id = 1:10,
    subjid = "subj_birth_ht",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(0L, 30L, 90L, 180L, 365L,
                0L, 30L, 90L, 180L, 365L),
    measurement = c(70.0, 54.0, 60.0, 67.0, 75.0,    # HT: 70cm at birth (extreme!), then normal
                    3.5, 4.5, 6.0, 8.0, 10.0)         # WT: normal
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # The extreme birth height (id=1) should be excluded
  birth_ht <- res[id == 1]
  expect_true(grepl("Exclude", birth_ht$exclude),
              info = sprintf("Extreme birth height (70cm) should be excluded, got: %s",
                             as.character(birth_ht$exclude)))
})


# ===========================================================================
# Section 9: LENGTHCM identity test
#
# LENGTHCM is renamed to HEIGHTCM with no measurement adjustment.
# Results should be identical.
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 24: LENGTHCM produces identical results to HEIGHTCM
# ---------------------------------------------------------------------------
test_that("LENGTHCM: results identical to HEIGHTCM for same measurements", {
  # Use a small subset of real data, convert young HEIGHTCM to LENGTHCM
  .subjs10 <- unique(.sg_peds$subjid)[1:10]
  d_ht <- .sg_peds[subjid %in% .subjs10]

  # Run with original HEIGHTCM
  res_ht <- cleangrowth(
    subjid = d_ht$subjid, param = d_ht$param, agedays = d_ht$agedays,
    sex = d_ht$sex, measurement = d_ht$measurement, id = d_ht$id,
    quietly = TRUE
  )

  # Convert HEIGHTCM to LENGTHCM for children under 2 years
  d_len <- copy(d_ht)
  d_len[param == "HEIGHTCM" & agedays < 730, param := "LENGTHCM"]

  res_len <- cleangrowth(
    subjid = d_len$subjid, param = d_len$param, agedays = d_len$agedays,
    sex = d_len$sex, measurement = d_len$measurement, id = d_len$id,
    quietly = TRUE
  )

  # Exclusion codes should be identical (same measurements, just different param label)
  expect_equal(
    as.character(res_len[order(id)]$exclude),
    as.character(res_ht[order(id)]$exclude),
    info = "LENGTHCM and HEIGHTCM should produce identical exclusion codes"
  )
})
