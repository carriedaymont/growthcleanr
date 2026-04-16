testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Layer 2: Parameter behavior tests for child algorithm
#
# Tests that function parameters change behavior as documented.
# Uses syngrowth built-in data.
# =============================================================================

# ---------------------------------------------------------------------------
# Shared data: load syngrowth once at file scope
# ---------------------------------------------------------------------------
data("syngrowth", package = "growthcleanr", envir = environment())
.sg <- as.data.table(syngrowth)
setkey(.sg, subjid, param, agedays)
.sg_peds <- .sg[agedays < 20 * 365.25]

# Helper: run child algorithm on N-subject pediatric subset
run_child_subset <- function(n_subjects, ...) {
  subjs <- unique(.sg$subjid)[seq_len(n_subjects)]
  d <- .sg_peds[subjid %in% subjs]
  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE,
    ...
  )
  return(list(input = d, result = res))
}

# ---------------------------------------------------------------------------
# Test 1: use_legacy_algorithm routes to legacy pediatric code
# ---------------------------------------------------------------------------
test_that("use_legacy_algorithm = TRUE runs legacy code and produces different results", {

  out_child <- run_child_subset(50)
  out_legacy <- run_child_subset(50, use_legacy_algorithm = TRUE)

  # Both return the same number of rows
  expect_equal(nrow(out_child$result), nrow(out_legacy$result))

  # Results should differ (algorithms are different)
  child_excl <- as.character(out_child$result$exclude)
  legacy_excl <- as.character(out_legacy$result$exclude)
  # At least some differences expected
  expect_true(any(child_excl != legacy_excl),
              info = "Child and legacy algorithms should produce at least some different results")

  # Both should have valid exclusion codes
  expect_false(any(is.na(out_legacy$result$exclude)))
})

# ---------------------------------------------------------------------------
# Test 2: sd.recenter options
# ---------------------------------------------------------------------------
test_that("sd.recenter parameter is accepted and runs without error", {

  out_nhanes <- run_child_subset(50, sd.recenter = "NHANES")
  out_none <- run_child_subset(50, sd.recenter = "none")

  # Both produce same row count and valid results
  expect_equal(nrow(out_nhanes$result), nrow(out_none$result))
  expect_false(any(is.na(out_nhanes$result$exclude)))
  expect_false(any(is.na(out_none$result$exclude)))

  # Recentering changes z-scores but may not always change exclusion codes
  # on small subsets — just verify both run without error
})

# ---------------------------------------------------------------------------
# Test 3: include.carryforward controls whether CFs are kept or excluded
#
# include.carryforward = TRUE means "keep CFs in output" (don't exclude them)
# include.carryforward = FALSE (default) means CFs are excluded
# ---------------------------------------------------------------------------
test_that("include.carryforward = TRUE keeps CFs, FALSE excludes them", {

  out_keep_cf <- run_child_subset(50, include.carryforward = TRUE)
  out_excl_cf <- run_child_subset(50, include.carryforward = FALSE)

  n_cf_kept <- sum(grepl("-CF$|-CF-deltaZ", out_keep_cf$result$exclude))
  n_cf_excl <- sum(grepl("-CF$|-CF-deltaZ", out_excl_cf$result$exclude))

  # include.carryforward = TRUE: CFs are kept as Include, so 0 CF exclusion codes
  expect_equal(n_cf_kept, 0,
               info = "include.carryforward=TRUE should keep CFs (no CF exclusion codes)")

  # include.carryforward = FALSE (default): CFs should be excluded
  expect_gt(n_cf_excl, 0,
            label = "CF count with include.carryforward = FALSE")
})

# ---------------------------------------------------------------------------
# Test 4: sd.extreme parameter changes extreme EWMA exclusions
# ---------------------------------------------------------------------------
test_that("sd.extreme parameter changes EWMA extreme exclusion behavior", {

  # Very strict threshold (lower sd.extreme = more exclusions)
  out_strict <- run_child_subset(50, sd.extreme = 3)
  # Very lenient threshold
  out_lenient <- run_child_subset(50, sd.extreme = 8)

  n_ewma_strict <- sum(grepl("EWMA", out_strict$result$exclude))
  n_ewma_lenient <- sum(grepl("EWMA", out_lenient$result$exclude))

  # Stricter threshold should exclude at least as many
  expect_gte(n_ewma_strict, n_ewma_lenient)
})

# ---------------------------------------------------------------------------
# Test 5: ewma_window parameter
# ---------------------------------------------------------------------------
test_that("ewma_window parameter is accepted and doesn't crash", {

  # Default is 15; try a different value
  out_w10 <- run_child_subset(20, ewma_window = 10)
  out_w20 <- run_child_subset(20, ewma_window = 20)

  # Both should complete without error

  expect_equal(nrow(out_w10$result), nrow(out_w10$input))
  expect_equal(nrow(out_w20$result), nrow(out_w20$input))

  # No NA exclusion codes
  expect_false(any(is.na(out_w10$result$exclude)))
  expect_false(any(is.na(out_w20$result$exclude)))
})

# ---------------------------------------------------------------------------
# Test 6: error.load.threshold and error.load.mincount
# ---------------------------------------------------------------------------
test_that("error load parameters affect error-load exclusions", {

  # Very aggressive error load threshold (0.1 = 10%)
  out_aggressive <- run_child_subset(100, error.load.threshold = 0.1,
                                      error.load.mincount = 1)
  # Very lenient (0.99 = 99%)
  out_lenient <- run_child_subset(100, error.load.threshold = 0.99,
                                   error.load.mincount = 10)

  n_errload_agg <- sum(out_aggressive$result$exclude == "Exclude-C-Too-Many-Errors")
  n_errload_len <- sum(out_lenient$result$exclude == "Exclude-C-Too-Many-Errors")

  # Aggressive should find at least as many error-load as lenient
  expect_gte(n_errload_agg, n_errload_len)
})

# ---------------------------------------------------------------------------
# Test 7: lt3.exclude.mode controls handling of subjects with < 3 measurements
# ---------------------------------------------------------------------------
# lt3.exclude.mode test removed — parameter was removed from cleangrowth()
# (walkthrough-todo-2026-04-13, D6)

# ---------------------------------------------------------------------------
# Test 8: recover.unit.error parameter
# ---------------------------------------------------------------------------
test_that("recover.unit.error parameter is accepted and doesn't crash", {

  out_recover <- run_child_subset(50, recover.unit.error = TRUE)
  out_no_recover <- run_child_subset(50, recover.unit.error = FALSE)

  # Both complete without error
  expect_equal(nrow(out_recover$result), nrow(out_recover$input))
  expect_equal(nrow(out_no_recover$result), nrow(out_no_recover$input))
})

# ---------------------------------------------------------------------------
# Test 9: Imperial unit conversion (HEIGHTIN, WEIGHTLBS)
# ---------------------------------------------------------------------------
test_that("imperial units (HEIGHTIN, WEIGHTLBS) are converted and processed", {

  d <- .sg_peds[subjid %in% unique(.sg$subjid)[1:10]]

  # Convert metric to imperial
  d_imp <- copy(d)
  d_imp[param == "HEIGHTCM", `:=`(
    measurement = measurement / 2.54,  # cm -> inches
    param = "HEIGHTIN"
  )]
  d_imp[param == "WEIGHTKG", `:=`(
    measurement = measurement * 2.20462,  # kg -> lbs
    param = "WEIGHTLBS"
  )]

  res_imp <- cleangrowth(
    subjid = d_imp$subjid,
    param = d_imp$param,
    agedays = d_imp$agedays,
    sex = d_imp$sex,
    measurement = d_imp$measurement,
    id = d_imp$id,
    quietly = TRUE
  )

  # All rows processed (none dropped as Missing due to unrecognized param)
  expect_equal(nrow(res_imp), nrow(d_imp))
  expect_equal(sum(res_imp$exclude == "Exclude-Missing"), 0,
               info = "No Missing codes should appear for valid imperial measurements")

  # Some should be Include
  expect_gt(sum(res_imp$exclude == "Include"), 0)
})

# ---------------------------------------------------------------------------
# Test 10: LENGTHCM treated same as HEIGHTCM
# ---------------------------------------------------------------------------
test_that("LENGTHCM param is accepted and processed", {

  d <- .sg_peds[subjid %in% unique(.sg$subjid)[1:10]]

  # Change all HEIGHTCM to LENGTHCM for young children
  d_len <- copy(d)
  d_len[param == "HEIGHTCM" & agedays < 730, param := "LENGTHCM"]

  res <- cleangrowth(
    subjid = d_len$subjid,
    param = d_len$param,
    agedays = d_len$agedays,
    sex = d_len$sex,
    measurement = d_len$measurement,
    id = d_len$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), nrow(d_len))
  # No Missing from unrecognized params
  expect_equal(sum(res$exclude == "Exclude-Missing"), 0)
})
