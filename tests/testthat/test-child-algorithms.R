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
# Shared helpers
# ---------------------------------------------------------------------------

#' Load syngrowth pediatric data
load_syngrowth_peds <- function() {
  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- as.data.table(syngrowth)
  setkey(dt, subjid, param, agedays)
  dt[agedays < 20 * 365.25]
}

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

#' Extract a single subject's data for one param from syngrowth peds
#' Returns only Include rows from a baseline GC run (first n_subj subjects)
get_clean_subject <- function(dt_peds, target_subjid, target_param, n_subj = 200) {
  subjs <- unique(dt_peds$subjid)[seq_len(n_subj)]
  d <- dt_peds[subjid %in% subjs]
  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )
  d2 <- merge(d, res[, .(id, exclude)], by = "id")
  d2[subjid == target_subjid & param == target_param & exclude == "Include"]
}


# ===========================================================================
# Section 1: CF rescue tests
#
# CF rescue re-includes carried-forward values when the z-score difference
# between the CF'd value and the originator is small enough for that
# age/interval/param/rounding cell. Three modes:
#   cf_rescue = "standard" (default) — age/interval/param lookup thresholds
#   cf_rescue = "none"     — no rescue (all CFs excluded)
#   cf_rescue = "all"      — all CFs rescued (no CF exclusions)
# ===========================================================================

# ---------------------------------------------------------------------------
# Test 1: Standard rescue rescues more CFs than "none" mode
# ---------------------------------------------------------------------------
test_that("CF rescue: standard rescues more CFs than none", {
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:200]
  d <- dt_peds[subjid %in% subjs]

  res_std <- run_gc(d, cf_rescue = "standard")
  res_none <- run_gc(d, cf_rescue = "none")

  n_cf_std  <- sum(grepl("-CF$", res_std$exclude))
  n_cf_none <- sum(grepl("-CF$", res_none$exclude))

  # Standard rescue should leave fewer CFs excluded than no rescue
  expect_lte(n_cf_std, n_cf_none)
})

# ---------------------------------------------------------------------------
# Test 2: cf_rescued column tracks rescue status
# ---------------------------------------------------------------------------
test_that("CF rescue: cf_rescued column populated for rescued CFs", {
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:200]
  d <- dt_peds[subjid %in% subjs]

  res <- run_gc(d, cf_rescue = "standard")

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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  res <- run_gc(d, cf_rescue = "none")

  rescued <- res[cf_rescued != ""]
  expect_equal(nrow(rescued), 0,
               info = "cf_rescue='none' should rescue no CFs")
})

# ---------------------------------------------------------------------------
# Test 4: cf_rescue = "all" rescues all CFs
# ---------------------------------------------------------------------------
test_that("CF rescue: all mode rescues all CFs", {
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:200]
  d <- dt_peds[subjid %in% subjs]

  res <- run_gc(d, cf_rescue = "all")

  # No CF exclusion codes should remain
  n_cf <- sum(grepl("-CF$", res$exclude))
  expect_equal(n_cf, 0,
               info = "cf_rescue='all' should leave no CF exclusions")
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  # Find a subject with several Include heights in middle of trajectory
  # Subject 0d8773f3 has 19 HT measurements
  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- d[subjid == target & param == "HEIGHTCM"]

  # Pick a middle measurement and multiply by 2.54 (unit error: inches -> cm)
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_id <- ht_rows$id[mid_idx]
  d_mod <- copy(d)
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- d[subjid == target & param == "HEIGHTCM"]

  # Pick two consecutive middle measurements
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_ids <- ht_rows$id[mid_idx:(mid_idx + 1)]
  d_mod <- copy(d)
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- d[subjid == target & param == "HEIGHTCM"]

  # Pick three consecutive middle measurements
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_ids <- ht_rows$id[mid_idx:(mid_idx + 2)]
  d_mod <- copy(d)
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- d[subjid == target & param == "HEIGHTCM"]

  # Inject one unit error in the middle
  mid_idx <- ceiling(nrow(ht_rows) / 2)
  error_id <- ht_rows$id[mid_idx]
  d_mod <- copy(d)
  d_mod[id == error_id, measurement := measurement * 2.54]

  # Run both original and modified
  res_orig <- run_gc(d)
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:200]
  d <- dt_peds[subjid %in% subjs]

  # threshold=0.4 should catch subjects that 0.5 misses
  # (confirmed: 2 Error-load at 0.4, 0 at 0.5 on this dataset)
  res_04 <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, error.load.threshold = 0.4
  )
  res_05 <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:200]
  d <- dt_peds[subjid %in% subjs]

  # mincount=1 with low threshold should catch more than mincount=10
  res_min1 <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, error.load.threshold = 0.3, error.load.mincount = 1
  )
  res_min10 <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  # Take a subject with many HT measurements and corrupt most of them
  # Subject 0d8773f3 has 19 HT measurements
  target <- "0d8773f3-c18e-9736-0a78-f1fda9b4fa0a"
  ht_rows <- d[subjid == target & param == "HEIGHTCM"]

  # Corrupt all but 2 measurements with distinct extreme values
  # (distinct values avoid SDE-Identical/CF which are excluded from error load denominator)
  # Keep first and last, corrupt everything else
  corrupt_ids <- ht_rows$id[2:(nrow(ht_rows) - 1)]
  d_mod <- copy(d)
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
  dt_peds <- load_syngrowth_peds()

  # Take a subject with data starting from birth
  # Subject 0db2905c: sex=0, agedays 0-1071
  target <- "0db2905c-1e0c-0d4d-7d86-4552c5b55ebd"
  subj_data <- dt_peds[subjid == target]
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
  other_subjs <- unique(dt_peds$subjid)[1:20]
  other_subjs <- other_subjs[other_subjs != target]
  d <- rbind(d, dt_peds[subjid %in% other_subjs[1:10]])
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
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
  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:50]
  d <- dt_peds[subjid %in% subjs]

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
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

  dt_peds <- load_syngrowth_peds()
  subjs <- unique(dt_peds$subjid)[1:100]
  d <- dt_peds[subjid %in% subjs]

  res_seq <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, parallel = FALSE
  )

  res_par <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
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
