testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Tests for cleangrowth() API: adult algorithm, edge cases, integration.
#
# Updated for v3.0.0: cleangrowth() returns a data.table (not a vector).
# =============================================================================

# Helper: merge cleangrowth results back to input by id
run_and_merge <- function(d, ...) {
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
  d_out <- copy(d)
  d_out[, gcr_result := res[match(d$id, res$id)]$exclude]
  return(d_out)
}

test_that("growthcleanr works as expected on adult synthetic data", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  data <- as.data.table(syngrowth)
  expect_equal(77721, data[, .N])
  setkey(data, subjid, param, agedays)

  data_adult <- copy(data[agedays >= 18 * 365.25, ])

  d100 <- copy(data_adult[subjid %in% unique(data_adult[, subjid])[1:100]])
  expect_equal(2023, d100[, .N])

  # Clean sample (default cutpoint = 20)
  cd100 <- run_and_merge(d100)

  # Clean with lower cutpoint = 18
  cd100cp <- run_and_merge(copy(d100), adult_cutpoint = 18)

  gcr_result <- function(dt, rowid) {
    as.character(dt[id == rowid]$gcr_result)
  }

  # These results should not change with cutpoint
  expect_equal("Include", gcr_result(cd100, 18695))
  expect_equal("Include", gcr_result(cd100cp, 18695))

  expect_equal("Exclude-A-Identical", gcr_result(cd100, 167))
  expect_equal("Exclude-A-Identical", gcr_result(cd100cp, 167))

  expect_equal("Exclude-A-BIV", gcr_result(cd100, 22009))
  expect_equal("Exclude-A-BIV", gcr_result(cd100cp, 22009))

  # Results should change due to younger cutpoint
  # With default cutpoint=20, records 18-20 yrs go through child algorithm
  expect_equal("Exclude-C-Extraneous", gcr_result(cd100, 69740))
  expect_equal("Include", gcr_result(cd100cp, 69740))

  expect_equal("Include", gcr_result(cd100cp, 55171))

  expect_equal("Exclude-C-BIV", gcr_result(cd100, 25259))
  expect_equal("Exclude-A-BIV", gcr_result(cd100cp, 25259))

  # Check counts (codes are no longer param-specific; use exact matching)
  catcount <- function(dt, category) {
    dt[gcr_result == category, .N]
  }

  expect_equal(1595, catcount(cd100, "Include"))
  expect_equal(334, catcount(cd100, "Exclude-A-Extraneous"))
  expect_equal(5, catcount(cd100, "Exclude-C-CF"))
  expect_equal(16, catcount(cd100, "Exclude-A-BIV"))

  expect_equal(1595, catcount(cd100cp, "Include"))
  expect_equal(356, catcount(cd100cp, "Exclude-A-Extraneous"))
  expect_equal(0, catcount(cd100cp, "Exclude-C-CF"))
  expect_equal(19, catcount(cd100cp, "Exclude-A-BIV"))

  # Missing data test for adults
  d5_miss <- copy(data_adult[subjid %in% unique(data[, subjid])[1:5]])
  set.seed(10)
  d5_miss$measurement[sample(seq_len(nrow(d5_miss)), 5)] <- NA

  res_miss <- cleangrowth(
    subjid = d5_miss$subjid,
    param = d5_miss$param,
    agedays = d5_miss$agedays,
    sex = d5_miss$sex,
    measurement = d5_miss$measurement,
    id = d5_miss$id,
    quietly = TRUE
  )

  expect_equal(sum(res_miss$exclude == "Exclude-Missing"), 5)
})

test_that("growthcleanr works without either adult or pediatric data", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- as.data.table(syngrowth)

  only_peds <- dt[agedays < 20 * 365.25][1:20]
  only_adult <- dt[agedays >= 20 * 365.25][1:20]
  nobody <- dt[agedays > 120 * 365.25]

  # Testing cleangrowth works without adult data
  peds_res <- cleangrowth(
    subjid = only_peds$subjid,
    param = only_peds$param,
    agedays = only_peds$agedays,
    sex = only_peds$sex,
    measurement = only_peds$measurement,
    id = only_peds$id,
    quietly = TRUE
  )
  expect_equal(nrow(peds_res), nrow(only_peds))

  # Testing cleangrowth works without pediatric data
  adult_res <- cleangrowth(
    subjid = only_adult$subjid,
    param = only_adult$param,
    agedays = only_adult$agedays,
    sex = only_adult$sex,
    measurement = only_adult$measurement,
    id = only_adult$id,
    quietly = TRUE
  )
  expect_equal(nrow(adult_res), nrow(only_adult))

  # Testing cleangrowth works with no data
  if (nrow(nobody) > 0) {
    no_res <- cleangrowth(
      subjid = nobody$subjid,
      param = nobody$param,
      agedays = nobody$agedays,
      sex = nobody$sex,
      measurement = nobody$measurement,
      id = nobody$id,
      quietly = TRUE
    )
    expect_equal(nrow(no_res), nrow(nobody))
  } else {
    # no data at all — skip
    expect_equal(0, nrow(nobody))
  }
})

# ===========================================================================
# Section: Child-adult spanning subjects
#
# Subjects with measurements on both sides of the adult_cutpoint (default 20y)
# have their data split: pediatric rows go through cleanchild(), adult rows
# through cleanadult(), then results are merged. These tests verify:
#   - No rows lost in the split/merge
#   - Child rows get child codes, adult rows get adult codes
#   - Constructed spanning subjects work correctly
# ===========================================================================

test_that("spanning subject: all rows preserved through child/adult split", {

  # Find a subject in syngrowth that spans the 20-year cutpoint
  data("syngrowth", package = "growthcleanr", envir = environment())
  sg <- as.data.table(syngrowth)
  setkey(sg, subjid, param, agedays)
  cutday <- 20 * 365.25

  # Find subjects with both child and adult data
  span_info <- sg[, .(n_child = sum(agedays < cutday),
                       n_adult = sum(agedays >= cutday)), by = subjid]
  span_subjs <- span_info[n_child > 0 & n_adult > 0]$subjid

  expect_true(length(span_subjs) > 0,
              info = "syngrowth should contain spanning subjects")

  # Take a spanning subject with substantial data on both sides
  # Use a571b120 which has 30 child + 23 adult rows
  target <- "a571b120-4ca7-86be-2ee3-a59d52edae59"
  d <- sg[subjid == target]

  n_child <- sum(d$agedays < cutday)
  n_adult <- sum(d$agedays >= cutday)

  expect_true(n_child > 0 & n_adult > 0,
              info = "Target subject should have both child and adult measurements")

  # Add some context subjects for GC to work with
  context_subjs <- unique(sg$subjid)[1:20]
  context_subjs <- context_subjs[context_subjs != target]
  d <- rbind(d, sg[subjid %in% context_subjs[1:10]])

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # All rows should be present in output
  expect_equal(nrow(res), nrow(d),
               info = "Output should have same number of rows as input")

  # All of the target subject's rows should be present
  target_res <- res[subjid == target]
  target_input <- d[subjid == target]
  expect_equal(nrow(target_res), nrow(target_input),
               info = "Spanning subject should have all rows in output")

  # No NA exclusion codes
  expect_false(any(is.na(target_res$exclude)),
               info = "No NA exclusion codes for spanning subject")
})

test_that("spanning subject: child rows get child codes, adult rows get adult codes", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  sg <- as.data.table(syngrowth)
  setkey(sg, subjid, param, agedays)
  cutday <- 20 * 365.25

  # Use subject with substantial data on both sides
  target <- "a571b120-4ca7-86be-2ee3-a59d52edae59"
  d <- sg[subjid == target]

  # Add context subjects
  context_subjs <- unique(sg$subjid)[1:20]
  context_subjs <- context_subjs[context_subjs != target]
  d <- rbind(d, sg[subjid %in% context_subjs[1:10]])

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  target_res <- res[subjid == target]

  # Split by age
  child_res <- target_res[agedays < cutday]
  adult_res <- target_res[agedays >= cutday]

  expect_true(nrow(child_res) > 0 & nrow(adult_res) > 0,
              info = "Should have results on both sides of cutpoint")

  # Child rows: exclusion codes should be Include, Exclude-Missing,
  # Exclude-Not-Cleaned, or Exclude-C-*
  child_codes <- as.character(child_res$exclude)
  child_ok <- grepl("^(Include|Exclude-Missing|Exclude-Not-Cleaned|Exclude-C-)", child_codes)
  expect_true(all(child_ok),
              info = sprintf("Child rows should have child codes, got: %s",
                             paste(unique(child_codes[!child_ok]), collapse = ", ")))

  # Adult rows: exclusion codes should be Include, Exclude-Missing, or Exclude-A-*
  adult_codes <- as.character(adult_res$exclude)
  adult_ok <- grepl("^(Include|Exclude-Missing|Exclude-A-)", adult_codes)
  expect_true(all(adult_ok),
              info = sprintf("Adult rows should have adult codes, got: %s",
                             paste(unique(adult_codes[!adult_ok]), collapse = ", ")))
})

test_that("spanning subject: constructed subject processed correctly across boundary", {

  # Construct a subject with clear growth trajectory across child/adult boundary
  # Male: normal growth from age 18y to 22y (crossing 20y cutpoint at 7305 days)
  d_span <- data.table(
    id = 1:12,
    subjid = "subj_span",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(6570L, 6935L, 7200L,    # HT: child (~18y, ~19y, ~19.7y)
                7400L, 7670L, 8000L,    # HT: adult (~20.3y, ~21y, ~21.9y)
                6570L, 6935L, 7200L,    # WT: child
                7400L, 7670L, 8000L),   # WT: adult
    measurement = c(175.0, 176.0, 176.5,        # HT: adolescent boy near adult height
                    177.0, 177.0, 177.0,        # HT: adult (stable)
                    70.0, 72.0, 73.0,           # WT: normal
                    74.0, 75.0, 76.0)           # WT: normal
  )

  res <- cleangrowth(
    subjid = d_span$subjid, param = d_span$param, agedays = d_span$agedays,
    sex = d_span$sex, measurement = d_span$measurement, id = d_span$id,
    quietly = TRUE
  )

  # All 12 rows should be present
  expect_equal(nrow(res), 12)

  # No NA exclusion codes
  expect_false(any(is.na(res$exclude)))

  cutday <- 20 * 365.25

  # Verify child/adult code assignment
  child_ids <- d_span[agedays < cutday]$id
  adult_ids <- d_span[agedays >= cutday]$id

  child_codes <- as.character(res[id %in% child_ids]$exclude)
  adult_codes <- as.character(res[id %in% adult_ids]$exclude)

  # Child codes should not have adult prefix
  expect_false(any(grepl("^Exclude-A-", child_codes)),
               info = "Child-age rows should not get adult exclusion codes")

  # Adult codes should not have child prefix
  expect_false(any(grepl("^Exclude-C-", adult_codes)),
               info = "Adult-age rows should not get child exclusion codes")

  # Normal values on both sides should mostly be Include
  # (small dataset may cause some exclusions, but not all)
  n_include <- sum(res$exclude == "Include")
  expect_true(n_include >= 6,
              info = "Normal spanning subject should have at least half rows included")
})


# ===========================================================================
# Section: Adult integration tests (cleangrowth → cleanadult)
#
# These tests verify that cleangrowth() correctly preprocesses data,
# passes parameters, and assembles output for adult-age rows. The adult
# unit tests (test-adult-clean.R) call cleanadult() directly — these
# tests exercise the full pipeline path.
# ===========================================================================

test_that("adult integration: output columns present and correctly typed", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  sg <- as.data.table(syngrowth)
  d_adult <- sg[agedays >= 20 * 365.25]
  d10 <- d_adult[subjid %in% unique(d_adult$subjid)[1:10]]

  res <- cleangrowth(
    subjid = d10$subjid, param = d10$param, agedays = d10$agedays,
    sex = d10$sex, measurement = d10$measurement, id = d10$id,
    quietly = TRUE
  )

  # All rows returned
  expect_equal(nrow(res), nrow(d10))

  # Adult-specific columns present
  expect_true("mean_ht" %in% names(res),
              info = "mean_ht column should be present for adult data")
  expect_true("bin_result" %in% names(res),
              info = "bin_result column should be present for adult data")
  expect_true("bin_exclude" %in% names(res),
              info = "bin_exclude column should be present")

  # mean_ht should be numeric and non-NA for subjects with included heights
  expect_true(is.numeric(res$mean_ht))
  has_incl_ht <- res[param == "HEIGHTCM" & exclude == "Include", unique(subjid)]
  if (length(has_incl_ht) > 0) {
    expect_true(any(!is.na(res[subjid %in% has_incl_ht]$mean_ht)),
                info = "mean_ht should be populated for subjects with included heights")
  }

  # bin_result should be Include or Exclude (no NAs for adult-only data)
  expect_true(all(res$bin_result %in% c("Include", "Exclude")),
              info = "bin_result should be Include or Exclude for adult data")

  # exclude should be a factor with no NAs
  expect_true(is.factor(res$exclude))
  expect_false(any(is.na(res$exclude)))
})

test_that("adult integration: imperial units processed correctly", {

  # HEIGHTIN and WEIGHTLBS should be converted internally by cleanadult()
  d <- data.table(
    id = 1:6,
    subjid = "imp_adult",
    sex = 0L,
    param = c("HEIGHTIN", "HEIGHTIN", "HEIGHTIN",
              "WEIGHTLBS", "WEIGHTLBS", "WEIGHTLBS"),
    agedays = c(8000L, 8400L, 8800L,
                8000L, 8400L, 8800L),
    measurement = c(70.0, 70.0, 70.5,        # ~178cm
                    170.0, 175.0, 180.0)      # ~77-82kg
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # All rows returned
  expect_equal(nrow(res), 6)

  # No NA exclusion codes
  expect_false(any(is.na(res$exclude)))

  # Param should be preserved as imperial in output
  expect_true(all(res[id %in% 1:3]$param == "HEIGHTIN"),
              info = "HEIGHTIN param should be preserved in output")
  expect_true(all(res[id %in% 4:6]$param == "WEIGHTLBS"),
              info = "WEIGHTLBS param should be preserved in output")

  # mean_ht should reflect metric conversion (70 inches ≈ 177.8 cm)
  mht <- unique(res$mean_ht[!is.na(res$mean_ht)])
  expect_true(length(mht) == 1 && mht > 170 && mht < 185,
              info = sprintf("mean_ht should be in cm (~178), got: %s",
                             paste(round(mht, 1), collapse = ", ")))

  # Normal values should be included
  expect_true(all(res$exclude == "Include"),
              info = "Normal imperial adult values should all be included")
})

test_that("adult integration: adult_permissiveness parameter passes through", {

  # 115cm height: BIV at looser (min 120cm) but passes BIV at loosest (min 50cm)
  d <- data.table(
    id = 1:6,
    subjid = "perm_test",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(8000L, 8400L, 8800L,
                8000L, 8400L, 8800L),
    measurement = c(115.0, 170.0, 170.0,
                    70.0, 72.0, 74.0)
  )

  res_looser <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, adult_permissiveness = "looser"
  )
  res_loosest <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, adult_permissiveness = "loosest"
  )

  # 115cm should be BIV at looser (outside 120-230 range)
  expect_equal(as.character(res_looser[id == 1]$exclude), "Exclude-A-BIV",
               info = "115cm height should be BIV at looser permissiveness")

  # 115cm should NOT be BIV at loosest (inside 50-244 range)
  expect_true(as.character(res_loosest[id == 1]$exclude) != "Exclude-A-BIV",
              info = "115cm height should pass BIV at loosest permissiveness")
})

test_that("adult integration: adult_scale_max_lbs passes through", {

  # 136.078 kg = 300 lbs — should be flagged with scale_max_lbs=300
  d <- data.table(
    id = 1:6,
    subjid = "cap_test",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(8000L, 8400L, 8800L,
                8000L, 8400L, 8800L),
    measurement = c(170.0, 170.0, 170.5,
                    70.0, 72.0, 136.078)   # 300 lbs
  )

  res_cap <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE, adult_scale_max_lbs = 300
  )
  res_nocap <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # With cap: scale max weight should be flagged
  expect_equal(as.character(res_cap[id == 6]$exclude), "Exclude-A-Scale-Max",
               info = "136kg (=300lbs) should be flagged with scale_max_lbs=300")

  # Without cap: same weight should not be scale-max flagged
  expect_true(as.character(res_nocap[id == 6]$exclude) != "Exclude-A-Scale-Max",
              info = "136kg should not be scale-max flagged without cap")
})

test_that("adult integration: weight_cap deprecation warns and works", {

  d <- data.table(
    id = 1:6,
    subjid = "dep_test",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(8000L, 8400L, 8800L,
                8000L, 8400L, 8800L),
    measurement = c(170.0, 170.0, 170.5,
                    70.0, 72.0, 136.078)   # 300 lbs
  )

  # weight_cap should warn about deprecation
  expect_warning(
    res <- cleangrowth(
      subjid = d$subjid, param = d$param, agedays = d$agedays,
      sex = d$sex, measurement = d$measurement, id = d$id,
      quietly = TRUE, weight_cap = 300
    ),
    regexp = "weight_cap.*deprecated"
  )

  # But should still work (passed through as adult_scale_max_lbs)
  expect_equal(as.character(res[id == 6]$exclude), "Exclude-A-Scale-Max",
               info = "Deprecated weight_cap should still function")
})

test_that("adult integration: character/UUID ids preserved", {

  d <- data.table(
    id = c("uuid-a1", "uuid-a2", "uuid-a3", "uuid-b1", "uuid-b2", "uuid-b3"),
    subjid = "subj1",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(8000L, 8400L, 8800L,
                8000L, 8400L, 8800L),
    measurement = c(170.0, 171.0, 170.5,
                    70.0, 72.0, 74.0)
  )

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  # All original ids should be present in output
  expect_true(all(d$id %in% res$id),
              info = "All input ids should be present in output")
  expect_true(all(res$id %in% d$id),
              info = "No extra ids should appear in output")

  # Id type should be preserved as character
  expect_true(is.character(res$id),
              info = "Character ids should remain character in output")
})

test_that("adult integration: spanning subject has correct output structure", {

  # Subject with data on both sides of cutpoint
  # Child rows should have NA mean_ht/bin_result; adult rows should have values
  d <- data.table(
    id = 1:12,
    subjid = "span_out",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(6570L, 6935L, 7200L,    # HT: child (~18y, ~19y, ~19.7y)
                7400L, 7670L, 8000L,    # HT: adult (~20.3y, ~21y, ~21.9y)
                6570L, 6935L, 7200L,    # WT: child
                7400L, 7670L, 8000L),   # WT: adult
    measurement = c(175.0, 176.0, 176.5,
                    177.0, 177.0, 177.0,
                    70.0, 72.0, 73.0,
                    74.0, 75.0, 76.0)
  )

  cutday <- 20 * 365.25

  res <- cleangrowth(
    subjid = d$subjid, param = d$param, agedays = d$agedays,
    sex = d$sex, measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  child_res <- res[agedays < cutday]
  adult_res <- res[agedays >= cutday]

  # Child rows: mean_ht should be NA (it's an adult-only column)
  expect_true(all(is.na(child_res$mean_ht)),
              info = "mean_ht should be NA for child rows")

  # Adult rows: mean_ht should be populated (this subject has included heights)
  expect_true(any(!is.na(adult_res$mean_ht)),
              info = "mean_ht should be populated for adult rows with included heights")

  # bin_exclude should be present for all rows (derived in postprocessing)
  expect_true(all(!is.na(res$bin_exclude)),
              info = "bin_exclude should be non-NA for all rows")
})


