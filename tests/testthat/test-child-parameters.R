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

# (Test 1 and Test 2 removed — legacy algorithm removed in v3.0.0)

# ---------------------------------------------------------------------------
# Test 3: include.carryforward controls CF rescue behavior (deprecated alias
# for cf_rescue).
#
# include.carryforward = TRUE  → cf_rescue = "all"
#   Every detected CF is rescued (no Exclude-C-CF in output). Step 13
#   final-SDE resolution handles any resulting multi-Include SPAs.
# include.carryforward = FALSE → cf_rescue = "standard" (default)
#   Excludes CFs unless the lookup threshold rescues them.
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

# (Test 8 removed — recover.unit.error parameter removed in v3.0.0)

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
