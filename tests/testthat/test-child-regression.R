testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Regression + structural invariant tests for child algorithm (Infants_Main)
#
# These tests freeze known-good results from the child algorithm so that
# accidental behavioral changes are caught. When you intentionally change the
# algorithm, update the expected values here.
#
# Uses syngrowth built-in data, first 100 pediatric subjects (832 rows).
# Runtime: ~5-10 sec depending on machine.
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: run child algorithm on N-subject pediatric subset of syngrowth
# ---------------------------------------------------------------------------
run_child_subset <- function(n_subjects, sd_recenter = "NHANES") {
  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  data.table::setkey(dt, subjid, param, agedays)
  dt_peds <- dt[agedays < 20 * 365.25]
  subjs <- unique(dt$subjid)[seq_len(n_subjects)]
  d <- dt_peds[subjid %in% subjs]
  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    sd.recenter = sd_recenter,
    quietly = TRUE
  )
  return(list(input = d, result = res))
}

# Convenience: look up a single record's exclusion code by id
gcr_result <- function(res, rowid) {
  as.character(res[id == rowid]$exclude)
}

# ---------------------------------------------------------------------------
# Test 1: Structural invariants — these must ALWAYS hold
# ---------------------------------------------------------------------------
test_that("child algorithm: structural invariants on 100-subject subset", {

  out <- run_child_subset(100)
  d <- out$input
  res <- out$result

  # 1a. Output has same number of rows as input

  expect_equal(nrow(res), nrow(d),
               info = "Output row count must match input")

  # 1b. All input IDs are present in output
  expect_true(all(d$id %in% res$id),
              info = "Every input id must appear in output")

  # 1c. No NA exclusion codes
  expect_false(any(is.na(res$exclude)),
               info = "No exclusion codes should be NA")

  # 1d. Every exclusion code is a valid level
  valid_levels <- levels(res$exclude)
  expect_true(all(as.character(res$exclude) %in% valid_levels),
              info = "All exclusion codes must be valid factor levels")

  # 1e. Include count + all exclusion counts = total rows
  counts <- res[, .N, by = exclude]
  expect_equal(sum(counts$N), nrow(d),
               info = "Sum of all category counts must equal total rows")

  # 1f. Return is a data.table with expected columns
  expect_true(data.table::is.data.table(res))
  expect_true(all(c("id", "exclude") %in% names(res)),
              info = "Result must contain id and exclude columns")

  # 1g. Exclusion codes are a factor
  expect_true(is.factor(res$exclude))
})

# ---------------------------------------------------------------------------
# Test 2: Regression counts — frozen from 100-subject run
# ---------------------------------------------------------------------------
test_that("child algorithm: exclusion counts match expected on 100 subjects", {

  out <- run_child_subset(100)
  res <- out$result

  # Helper to get count for a category
  catcount <- function(category) {
    n <- res[exclude == category, .N]
    if (length(n) == 0) return(0L)
    return(n)
  }

  # Frozen expected counts (100 subjects, 832 rows, sd.recenter = "NHANES")
  expect_equal(nrow(res), 832)
  expect_equal(catcount("Include"), 599)
  expect_equal(catcount("Exclude-SDE-EWMA"), 127)
  expect_equal(catcount("Exclude-Carried-Forward"), 75)
  expect_equal(catcount("Exclude-SDE-Identical"), 13)
  expect_equal(catcount("Exclude-SDE-All-Extreme"), 6)
  expect_equal(catcount("Exclude-EWMA2-middle"), 3)
  expect_equal(catcount("Exclude-Absolute-BIV"), 3)
  expect_equal(catcount("Exclude-Min-diff"), 3)
  expect_equal(catcount("Exclude-Standardized-BIV"), 2)
  expect_equal(catcount("Exclude-Evil-Twins"), 1)

  # Exactly 10 distinct categories present (Error-load no longer triggers
  # at default threshold=0.5; previously triggered at hardcoded 0.4)
  expect_equal(length(unique(as.character(res$exclude))), 10)
})

# ---------------------------------------------------------------------------
# Test 3: Spot checks — individual records with known results
#
# These are chosen because their results are determined by within-subject
# patterns (CF, SDE-Identical, Absolute-BIV) or are stable across sample
# sizes with NHANES recentering. Verified consistent on 50- and 100-subject
# subsets.
# ---------------------------------------------------------------------------
test_that("child algorithm: spot-check individual records", {

  out <- run_child_subset(100)
  res <- out$result

  # Include — normal measurement
  expect_equal(gcr_result(res, 23751), "Include")

  # Exclude-SDE-EWMA — same-day extraneous based on EWMA
  expect_equal(gcr_result(res, 23755), "Exclude-SDE-EWMA")

  # Exclude-Carried-Forward — repeated value from prior visit
  expect_equal(gcr_result(res, 3), "Exclude-Carried-Forward")
  expect_equal(gcr_result(res, 4), "Exclude-Carried-Forward")

  # Exclude-SDE-Identical — same-day identical value
  expect_equal(gcr_result(res, 67275), "Exclude-SDE-Identical")

  # Exclude-SDE-All-Extreme — all same-day values are extreme
  expect_equal(gcr_result(res, 40091), "Exclude-SDE-All-Extreme")

  # Exclude-Absolute-BIV — biologically implausible value (313.6 cm at 434 days)
  expect_equal(gcr_result(res, 40108), "Exclude-Absolute-BIV")

  # Exclude-EWMA2-middle — EWMA pass 2 middle exclusion
  expect_equal(gcr_result(res, 7), "Exclude-EWMA2-middle")

  # Exclude-Min-diff — minimum height difference violation
  expect_equal(gcr_result(res, 38718), "Exclude-Min-diff")

  # ID 25251 was Exclude-Error-load when threshold was hardcoded at 0.4;
  # now Include at default threshold=0.5
  expect_equal(gcr_result(res, 25251), "Include")

  # Exclude-Standardized-BIV — standardized biologically implausible
  expect_equal(gcr_result(res, 25257), "Exclude-Standardized-BIV")

  # Exclude-Evil-Twins — evil twin detection (OTL in algorithm)
  expect_equal(gcr_result(res, 25258), "Exclude-Evil-Twins")
})

# ---------------------------------------------------------------------------
# Test 4: Stability across sample sizes (within-subject exclusions)
#
# Certain exclusion types are fully determined by within-subject patterns and
# should not change regardless of how many other subjects are in the batch.
# ---------------------------------------------------------------------------
test_that("child algorithm: within-subject exclusions stable at 50 vs 100 subjects", {

  out50 <- run_child_subset(50)
  out100 <- run_child_subset(100)

  # Records that exist in both 50- and 100-subject subsets
  # These within-subject results should be identical
  stable_ids <- c(
    3,      # CF
    4,      # CF
    67275,  # SDE-Identical
    40108,  # Absolute-BIV
    40091,  # SDE-All-Extreme
    23751,  # Include
    23755   # SDE-EWMA
  )

  for (sid in stable_ids) {
    r50 <- gcr_result(out50$result, sid)
    r100 <- gcr_result(out100$result, sid)
    expect_equal(r50, r100,
                 info = sprintf("id=%d should be stable across sample sizes", sid))
  }
})

# ---------------------------------------------------------------------------
# Test 5: Missing data handling
# ---------------------------------------------------------------------------
test_that("child algorithm: missing measurements return Missing", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  dt_peds <- dt[agedays < 20 * 365.25]

  # Use 20 subjects for more robust test (avoids edge cases with
  # same-day duplicates both being NA on tiny subsets)
  d20 <- dt_peds[subjid %in% unique(dt$subjid)[1:20]]

  # Introduce 10 missing values, avoiding same-day-same-param pairs
  # by picking one row per (subjid, param, agedays) group
  set.seed(42)
  d20_miss <- data.table::copy(d20)
  # Pick rows that are unique on subjid/param/agedays (no SDE duplicates)
  unique_spa <- d20_miss[, .I[1], by = .(subjid, param, agedays)]$V1
  missing_rows <- sample(unique_spa, 10)
  missing_ids <- d20_miss[missing_rows]$id
  d20_miss[missing_rows, measurement := NA]

  res <- suppressWarnings(cleangrowth(
    subjid = d20_miss$subjid,
    param = d20_miss$param,
    agedays = d20_miss$agedays,
    sex = d20_miss$sex,
    measurement = d20_miss$measurement,
    id = d20_miss$id,
    quietly = TRUE
  ))

  # Count of "Missing" should equal number of NA measurements
  n_missing_result <- sum(res$exclude == "Missing")
  expect_equal(n_missing_result, length(missing_ids),
               info = "Count of Missing results should match count of NA measurements")

  # Rows with NA measurement should never be coded "Include"
  na_rows <- res[id %in% missing_ids]
  expect_false(any(na_rows$exclude == "Include"),
               info = "NA measurements must not be coded Include")

  # Non-missing rows should NOT be "Missing"
  non_missing <- res[!id %in% missing_ids]
  expect_false(any(as.character(non_missing$exclude) == "Missing"))
})

# ---------------------------------------------------------------------------
# Test 6: HC (HEADCM) support
# ---------------------------------------------------------------------------
test_that("child algorithm: HEADCM measurements are processed", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  data.table::setkey(dt, subjid, param, agedays)
  dt_peds <- dt[agedays < 20 * 365.25]

  # Use 50 subjects — first 10 don't have infants, but first 50 has 8
  d50 <- dt_peds[subjid %in% unique(dt$subjid)[1:50]]
  young <- d50[param == "HEIGHTCM" & agedays < 1095]
  expect_gt(nrow(young), 0)

  # Add synthetic HC rows for young subjects (under 3 years)
  # Use plausible HC values: ~35cm at birth, growing ~1cm/month
  hc_rows <- young[, .(
    id = id + 100000L,
    subjid = subjid,
    sex = sex,
    param = "HEADCM",
    agedays = agedays,
    measurement = 35 + agedays * (10 / 365.25)  # ~10cm growth per year
  )]

  combined <- rbind(d50, hc_rows, fill = TRUE)

  res <- cleangrowth(
    subjid = combined$subjid,
    param = combined$param,
    agedays = combined$agedays,
    sex = combined$sex,
    measurement = combined$measurement,
    id = combined$id,
    quietly = TRUE
  )

  # All rows accounted for
  expect_equal(nrow(res), nrow(combined))

  # HC rows got exclusion codes (not NA)
  hc_results <- res[id %in% hc_rows$id]
  # All HC rows should have valid exclusion codes
  expect_false(any(is.na(hc_results$exclude)))

  # At least some HC rows should be "Include" (plausible values)
  expect_gt(sum(hc_results$exclude == "Include"), 0)
})

# ---------------------------------------------------------------------------
# Test 7: Backward compatibility — prelim_infants deprecation
# ---------------------------------------------------------------------------
test_that("prelim_infants = TRUE produces deprecation warning and runs child algorithm", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  dt_peds <- dt[agedays < 20 * 365.25]
  d5 <- dt_peds[subjid %in% unique(dt$subjid)[1:5]]

  # prelim_infants = TRUE should warn and produce child algorithm results
  expect_warning(
    res_pi <- cleangrowth(
      subjid = d5$subjid,
      param = d5$param,
      agedays = d5$agedays,
      sex = d5$sex,
      measurement = d5$measurement,
      id = d5$id,
      prelim_infants = TRUE,
      quietly = TRUE
    ),
    regexp = "prelim_infants.*deprecated"
  )

  # Should produce same results as default (child algorithm)
  res_default <- cleangrowth(
    subjid = d5$subjid,
    param = d5$param,
    agedays = d5$agedays,
    sex = d5$sex,
    measurement = d5$measurement,
    id = d5$id,
    quietly = TRUE
  )

  expect_equal(as.character(res_pi$exclude), as.character(res_default$exclude))
})

# ---------------------------------------------------------------------------
# Test 8: Performance benchmark — catch major regressions
#
# Runs 100-subject subset and reports elapsed time. Fails if runtime
# exceeds a generous ceiling (60 sec) that should only trigger if
# something has gone seriously wrong. The reported time is the useful
# signal — check it when making algorithm or batching changes.
# ---------------------------------------------------------------------------
test_that("child algorithm: 100-subject benchmark completes in reasonable time", {

  t0 <- proc.time()
  out <- run_child_subset(100)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  # Report timing — this is the main value of this test
  cat(sprintf(
    "\n=== BENCHMARK: 100 subjects (%d rows) completed in %.1f sec ===\n",
    nrow(out$result), elapsed
  ), file = stderr())

  # Sanity check that it actually ran
  expect_equal(nrow(out$result), 832)

  # Generous ceiling — not a tight benchmark, just catches catastrophic
  # slowdowns (e.g., batching bug causing N^2 behavior). Normal runtime
  # on this machine is ~2-5 sec for 100 subjects.
  expect_lt(elapsed, 60,
            label = sprintf("Runtime %.1f sec exceeded 60 sec ceiling", elapsed))
})

# ---------------------------------------------------------------------------
# Test 8: gc_preload_refs() + ref_tables produce identical results to full run
# ---------------------------------------------------------------------------
test_that("gc_preload_refs: ref_tables produces identical results to standard run", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  dt_peds <- dt[agedays < 20 * 365.25]
  subjs <- unique(dt_peds$subjid)[seq_len(30)]
  d <- dt_peds[subjid %in% subjs]

  # Standard run (reads refs from disk)
  res_standard <- cleangrowth(
    subjid = d$subjid,
    param  = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  # Pre-loaded refs run
  refs <- gc_preload_refs()
  res_preloaded <- cleangrowth(
    subjid = d$subjid,
    param  = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    ref_tables = refs,
    quietly = TRUE
  )

  # Results must be identical
  expect_equal(nrow(res_preloaded), nrow(res_standard),
               info = "Row counts must match")
  expect_equal(
    res_preloaded[order(id)]$exclude,
    res_standard[order(id)]$exclude,
    info = "Exclusion codes must be identical with and without ref_tables"
  )
})

# ---------------------------------------------------------------------------
# Test 9: changed_subjids partial run produces identical results to full run
# ---------------------------------------------------------------------------
test_that("changed_subjids: partial run produces identical results to full run", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- data.table::as.data.table(syngrowth)
  dt_peds <- dt[agedays < 20 * 365.25]
  subjs <- unique(dt_peds$subjid)[seq_len(30)]
  d <- dt_peds[subjid %in% subjs]

  # Full run (baseline cached results)
  res_full <- cleangrowth(
    subjid = d$subjid,
    param  = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  # Partial run: all subjects listed as "changed" — must equal full run
  res_partial <- cleangrowth(
    subjid = d$subjid,
    param  = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    cached_results = res_full,
    changed_subjids = unique(d$subjid),
    quietly = TRUE
  )

  expect_equal(nrow(res_partial), nrow(res_full),
               info = "Row counts must match")
  expect_equal(
    res_partial[order(id)]$exclude,
    res_full[order(id)]$exclude,
    info = "Exclusion codes must be identical: full run vs partial run with all subjids changed"
  )
})
