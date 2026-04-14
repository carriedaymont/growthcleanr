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
# Shared data: load syngrowth once, run 100-subject GC once
# ---------------------------------------------------------------------------
data("syngrowth", package = "growthcleanr", envir = environment())
.sg <- data.table::as.data.table(syngrowth)
data.table::setkey(.sg, subjid, param, agedays)
.sg_peds <- .sg[agedays < 20 * 365.25]

run_child_subset <- function(n_subjects, sd_recenter = "NHANES", ...) {
  subjs <- unique(.sg$subjid)[seq_len(n_subjects)]
  d <- .sg_peds[subjid %in% subjs]
  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    sd.recenter = sd_recenter,
    quietly = TRUE,
    ...
  )
  return(list(input = d, result = res))
}

# Cache the 100-subject default run — used by Tests 1-4 and benchmark
.out100 <- run_child_subset(100)

# Convenience: look up a single record's exclusion code by id
gcr_result <- function(res, rowid) {
  as.character(res[id == rowid]$exclude)
}

# ---------------------------------------------------------------------------
# Test 1: Structural invariants — these must ALWAYS hold
# ---------------------------------------------------------------------------
test_that("child algorithm: structural invariants on 100-subject subset", {

  d <- .out100$input
  res <- .out100$result

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

  res <- .out100$result

  # Helper to get count for a category
  catcount <- function(category) {
    n <- res[exclude == category, .N]
    if (length(n) == 0) return(0L)
    return(n)
  }

  # Frozen expected counts (100 subjects, 832 rows, sd.recenter = "NHANES")
  # Updated for CF rescue lookup thresholds (2026-04-14): WT-CF 43→30, Include 599→612
  expect_equal(nrow(res), 832)
  expect_equal(catcount("Include"), 612)
  expect_equal(catcount("Exclude-C-HT-Extraneous"), 67)
  expect_equal(catcount("Exclude-C-WT-Extraneous"), 66)
  expect_equal(catcount("Exclude-C-WT-CF"), 30)
  expect_equal(catcount("Exclude-C-HT-CF"), 32)
  expect_equal(catcount("Exclude-C-WT-Identical"), 8)
  expect_equal(catcount("Exclude-C-HT-BIV"), 5)
  expect_equal(catcount("Exclude-C-HT-Identical"), 5)
  expect_equal(catcount("Exclude-C-HT-Abs-Diff"), 3)
  expect_equal(catcount("Exclude-C-WT-Traj"), 2)
  expect_equal(catcount("Exclude-C-HT-Traj"), 1)
  expect_equal(catcount("Exclude-C-WT-Evil-Twins"), 1)

  # Exactly 12 distinct categories present
  expect_equal(length(unique(as.character(res$exclude))), 12)
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

  res <- .out100$result

  # Include — normal measurement
  expect_equal(gcr_result(res, 23751), "Include")

  # Exclude-C-HT-Extraneous — same-day extraneous (HEIGHTCM)
  expect_equal(gcr_result(res, 23755), "Exclude-C-HT-Extraneous")

  # Exclude-C-{param}-CF — carried forward (repeated value from prior visit)
  expect_equal(gcr_result(res, 3), "Exclude-C-HT-CF")
  expect_equal(gcr_result(res, 4), "Exclude-C-WT-CF")

  # Exclude-C-WT-Identical — same-day identical value (WEIGHTKG)
  expect_equal(gcr_result(res, 67275), "Exclude-C-WT-Identical")

  # Exclude-C-HT-Extraneous — same-day extraneous (HEIGHTCM)
  expect_equal(gcr_result(res, 40091), "Exclude-C-HT-Extraneous")

  # Exclude-C-HT-BIV — biologically implausible value (313.6 cm at 434 days)
  expect_equal(gcr_result(res, 40108), "Exclude-C-HT-BIV")

  # Exclude-C-WT-Traj — EWMA pass 2 middle exclusion (WEIGHTKG)
  expect_equal(gcr_result(res, 7), "Exclude-C-WT-Traj")

  # Exclude-C-HT-Abs-Diff — height difference violation
  expect_equal(gcr_result(res, 38718), "Exclude-C-HT-Abs-Diff")

  # ID 25251 was Exclude-Error-load when threshold was hardcoded at 0.4;
  # now Include at default threshold=0.5
  expect_equal(gcr_result(res, 25251), "Include")

  # Exclude-C-HT-BIV — standardized biologically implausible (HEIGHTCM)
  expect_equal(gcr_result(res, 25257), "Exclude-C-HT-BIV")

  # Exclude-C-WT-Evil-Twins — evil twin detection (WEIGHTKG)
  expect_equal(gcr_result(res, 25258), "Exclude-C-WT-Evil-Twins")
})

# ---------------------------------------------------------------------------
# Test 4: Stability across sample sizes (within-subject exclusions)
#
# Certain exclusion types are fully determined by within-subject patterns and
# should not change regardless of how many other subjects are in the batch.
# ---------------------------------------------------------------------------
test_that("child algorithm: within-subject exclusions stable at 50 vs 100 subjects", {

  out50 <- run_child_subset(50)
  out100 <- .out100

  # Records that exist in both 50- and 100-subject subsets
  # These within-subject results should be identical
  stable_ids <- c(
    3,      # CF (Exclude-C-HT-CF)
    4,      # CF (Exclude-C-WT-CF)
    67275,  # Identical (Exclude-C-WT-Identical)
    40108,  # BIV (Exclude-C-HT-BIV)
    40091,  # Extraneous (Exclude-C-HT-Extraneous)
    23751,  # Include
    23755   # Extraneous (Exclude-C-HT-Extraneous)
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

  # Use 20 subjects for more robust test (avoids edge cases with
  # same-day duplicates both being NA on tiny subsets)
  d20 <- .sg_peds[subjid %in% unique(.sg$subjid)[1:20]]

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
  n_missing_result <- sum(res$exclude == "Exclude-Missing")
  expect_equal(n_missing_result, length(missing_ids),
               info = "Count of Missing results should match count of NA measurements")

  # Rows with NA measurement should never be coded "Include"
  na_rows <- res[id %in% missing_ids]
  expect_false(any(na_rows$exclude == "Include"),
               info = "NA measurements must not be coded Include")

  # Non-missing rows should NOT be "Missing"
  non_missing <- res[!id %in% missing_ids]
  expect_false(any(as.character(non_missing$exclude) == "Exclude-Missing"))
})

# ---------------------------------------------------------------------------
# Test 6: HC (HEADCM) support
# ---------------------------------------------------------------------------
test_that("child algorithm: HEADCM measurements are processed", {

  # Use 50 subjects — first 10 don't have infants, but first 50 has 8
  d50 <- .sg_peds[subjid %in% unique(.sg$subjid)[1:50]]
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

  d5 <- .sg_peds[subjid %in% unique(.sg$subjid)[1:5]]

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
# Test 8: gc_preload_refs() + ref_tables produce identical results to full run
# ---------------------------------------------------------------------------
test_that("gc_preload_refs: ref_tables produces identical results to standard run", {

  subjs <- unique(.sg_peds$subjid)[seq_len(30)]
  d <- .sg_peds[subjid %in% subjs]

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

  subjs <- unique(.sg_peds$subjid)[seq_len(30)]
  d <- .sg_peds[subjid %in% subjs]

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
