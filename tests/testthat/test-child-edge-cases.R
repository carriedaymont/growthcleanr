testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Layer 3: Edge case tests for child algorithm
#
# Tests unusual inputs that could cause crashes or incorrect results.
# =============================================================================

# ---------------------------------------------------------------------------
# Shared data: load syngrowth once at file scope
# ---------------------------------------------------------------------------
data("syngrowth", package = "growthcleanr", envir = environment())
.sg <- as.data.table(syngrowth)
setkey(.sg, subjid, param, agedays)
.sg_peds <- .sg[agedays < 20 * 365.25]

# ---------------------------------------------------------------------------
# Test 1: Single subject
# ---------------------------------------------------------------------------
test_that("child algorithm handles single subject", {

  d1 <- .sg_peds[subjid == unique(.sg_peds$subjid)[1]]

  res <- cleangrowth(
    subjid = d1$subjid,
    param = d1$param,
    agedays = d1$agedays,
    sex = d1$sex,
    measurement = d1$measurement,
    id = d1$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), nrow(d1))
  expect_false(any(is.na(res$exclude)))
})

# ---------------------------------------------------------------------------
# Test 2: Subject with only 1 measurement per parameter
# ---------------------------------------------------------------------------
test_that("child algorithm handles subject with single measurement per param", {

  d <- data.table(
    id = 1:2,
    subjid = "subj001",
    sex = 0L,
    param = c("HEIGHTCM", "WEIGHTKG"),
    agedays = c(365, 365),
    measurement = c(75.0, 10.0)
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), 2)
  expect_false(any(is.na(res$exclude)))
})

# ---------------------------------------------------------------------------
# Test 3: Subject with only 2 measurements (triggers Step 19 singles/pairs)
# ---------------------------------------------------------------------------
test_that("child algorithm handles subject with exactly 2 measurements per param", {

  d <- data.table(
    id = 1:4,
    subjid = "subj001",
    sex = 1L,
    param = c("HEIGHTCM", "HEIGHTCM", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(365, 730, 365, 730),
    measurement = c(75.0, 85.0, 10.0, 12.0)
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), 4)
  expect_false(any(is.na(res$exclude)))
  # Plausible values should be Include
  expect_true(all(res$exclude == "Include"))
})

# ---------------------------------------------------------------------------
# Test 4: All measurements are NA (all Missing)
# ---------------------------------------------------------------------------
test_that("child algorithm handles all-NA measurements", {

  d <- data.table(
    id = 1:4,
    subjid = "subj001",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(365, 730, 365, 730),
    measurement = NA_real_
  )

  res <- suppressWarnings(cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  ))

  expect_equal(nrow(res), 4)
  expect_true(all(res$exclude == "Exclude-Missing"))
})

# ---------------------------------------------------------------------------
# Test 5: Mixed NA and valid measurements (same subject)
# ---------------------------------------------------------------------------
test_that("child algorithm handles mix of NA and valid measurements", {

  d <- .sg_peds[subjid %in% unique(.sg_peds$subjid)[1:5]]

  # Set half the measurements to NA
  d_half <- copy(d)
  set.seed(123)
  na_rows <- sample(seq_len(nrow(d_half)), nrow(d_half) %/% 2)
  na_ids <- d_half$id[na_rows]
  d_half[na_rows, measurement := NA_real_]

  res <- suppressWarnings(cleangrowth(
    subjid = d_half$subjid,
    param = d_half$param,
    agedays = d_half$agedays,
    sex = d_half$sex,
    measurement = d_half$measurement,
    id = d_half$id,
    quietly = TRUE
  ))

  expect_equal(nrow(res), nrow(d_half))
  # All NA rows should be Missing
  expect_true(all(res[id %in% na_ids]$exclude == "Exclude-Missing"))
  # Non-NA rows should NOT be Missing
  non_na_ids <- d_half$id[!seq_len(nrow(d_half)) %in% na_rows]
  expect_false(any(res[id %in% non_na_ids]$exclude == "Exclude-Missing"))
})

# ---------------------------------------------------------------------------
# Test 6: Same-day identical values (SDE-Identical)
# ---------------------------------------------------------------------------
test_that("child algorithm marks same-day identical measurements", {

  d <- data.table(
    id = 1:6,
    subjid = "subj001",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(365, 365, 730, 365, 365, 730),
    measurement = c(75.0, 75.0, 85.0, 10.0, 10.0, 12.0)
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  # Should have SDE-Identical exclusions for the duplicates
  n_sde_identical <- sum(grepl("Identical", res$exclude))
  expect_gt(n_sde_identical, 0,
            label = "Same-day identical values should get SDE-Identical")
})

# ---------------------------------------------------------------------------
# Test 7: Negative agedays marked as Missing
# ---------------------------------------------------------------------------
test_that("child algorithm marks negative agedays as Missing", {

  d <- data.table(
    id = 1:4,
    subjid = "subj001",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(-10, 365, -10, 365),
    measurement = c(50.0, 75.0, 3.5, 10.0)
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), 4)
  # Negative agedays should be Missing
  neg_rows <- res[agedays < 0]
  expect_true(all(neg_rows$exclude == "Exclude-Missing"))
})

# ---------------------------------------------------------------------------
# Test 8: HEADCM > 3 years excluded from cleaning
# ---------------------------------------------------------------------------
test_that("child algorithm excludes HEADCM > 3 years from cleaning", {

  # Use syngrowth subjects plus synthetic HC rows for a more realistic dataset
  d <- .sg_peds[subjid %in% unique(.sg_peds$subjid)[1:10]]

  # Add HC rows: some under 3 years, some over 3 years
  young_ht <- d[param == "HEIGHTCM" & agedays < 1000]
  old_ht <- d[param == "HEIGHTCM" & agedays > 1200]

  hc_young <- young_ht[, .(
    id = id + 200000L, subjid, sex, param = "HEADCM",
    agedays, measurement = 35 + agedays * (10 / 365.25)
  )]
  hc_old <- old_ht[, .(
    id = id + 300000L, subjid, sex, param = "HEADCM",
    agedays, measurement = 50 + agedays * (2 / 365.25)
  )]

  combined <- rbind(d, hc_young, hc_old, fill = TRUE)

  res <- suppressWarnings(cleangrowth(
    subjid = combined$subjid,
    param = combined$param,
    agedays = combined$agedays,
    sex = combined$sex,
    measurement = combined$measurement,
    id = combined$id,
    quietly = TRUE
  ))

  expect_equal(nrow(res), nrow(combined))

  # HEADCM under 3 years should NOT be "Exclude-Not-Cleaned"
  hc_young_res <- res[id %in% hc_young$id]
  if (nrow(hc_young_res) > 0) {
    expect_false(any(hc_young_res$exclude == "Exclude-Not-Cleaned"),
                 info = "HC under 3 years should not be 'Not cleaned'")
  }

  # HEADCM over 3 years should be excluded from cleaning
  # (marked "Exclude-Not-Cleaned" or equivalent exclusion)
  hc_old_res <- res[id %in% hc_old$id]
  if (nrow(hc_old_res) > 0) {
    expect_true(all(hc_old_res$exclude != "Include"),
                info = "HC over 3 years should not be Include")
  }
})

# ---------------------------------------------------------------------------
# Test 9: Extreme values get excluded
# ---------------------------------------------------------------------------
test_that("child algorithm excludes biologically implausible values", {

  # Include some extreme values alongside plausible ones
  d <- data.table(
    id = 1:8,
    subjid = rep("subj001", 8),
    sex = rep(0L, 8),
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(365, 730, 1095, 1460,
                365, 730, 1095, 1460),
    measurement = c(75.0, 85.0, 300.0, 100.0,   # 300cm is extreme
                    10.0, 12.0, 14.0, 200.0)     # 200kg at 4yrs is extreme
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  # Extreme values should be excluded (not Include)
  extreme_ht <- res[id == 3]  # 300cm height
  extreme_wt <- res[id == 8]  # 200kg weight
  expect_true(grepl("Exclude", extreme_ht$exclude),
              info = "300cm height at 3 years should be excluded")
  expect_true(grepl("Exclude", extreme_wt$exclude),
              info = "200kg weight at 4 years should be excluded")
})

# ---------------------------------------------------------------------------
# Test 10: Multiple subjects with varying data density
# ---------------------------------------------------------------------------
test_that("child algorithm handles mix of data-rich and data-sparse subjects", {

  # Get one data-rich subject and create a data-sparse one
  rich_subj <- unique(.sg_peds$subjid)[1]
  d_rich <- .sg_peds[subjid == rich_subj]

  d_sparse <- data.table(
    id = max(.sg_peds$id) + 1:2,
    subjid = "sparse_subj",
    sex = 0L,
    param = c("HEIGHTCM", "WEIGHTKG"),
    agedays = c(365, 365),
    measurement = c(75.0, 10.0)
  )

  combined <- rbind(d_rich, d_sparse, fill = TRUE)

  res <- cleangrowth(
    subjid = combined$subjid,
    param = combined$param,
    agedays = combined$agedays,
    sex = combined$sex,
    measurement = combined$measurement,
    id = combined$id,
    quietly = TRUE
  )

  expect_equal(nrow(res), nrow(combined))
  expect_false(any(is.na(res$exclude)))
})

# ---------------------------------------------------------------------------
# Test 11: Carried-forward values detected
# ---------------------------------------------------------------------------
test_that("child algorithm detects carried-forward values", {

  # Create a subject with obvious carry-forward (exact repeated weight)
  d <- data.table(
    id = 1:8,
    subjid = "subj001",
    sex = 0L,
    param = c("HEIGHTCM", "HEIGHTCM", "HEIGHTCM", "HEIGHTCM",
              "WEIGHTKG", "WEIGHTKG", "WEIGHTKG", "WEIGHTKG"),
    agedays = c(365, 730, 1095, 1460,
                365, 730, 1095, 1460),
    measurement = c(75.0, 85.0, 93.0, 100.0,  # height grows normally
                    10.0, 10.0, 10.0, 10.0)     # weight carried forward
  )

  res <- cleangrowth(
    subjid = d$subjid,
    param = d$param,
    agedays = d$agedays,
    sex = d$sex,
    measurement = d$measurement,
    id = d$id,
    quietly = TRUE
  )

  # Some of the repeated weights should be flagged as CF
  n_cf <- sum(grepl("-CF$|-CF-deltaZ", res$exclude))
  expect_gt(n_cf, 0, label = "Obvious carry-forwards should be detected")
})

# ---------------------------------------------------------------------------
# Test 12: Deterministic results (same input -> same output)
# ---------------------------------------------------------------------------
test_that("child algorithm produces deterministic results", {

  d <- .sg_peds[subjid %in% unique(.sg_peds$subjid)[1:20]]

  res1 <- cleangrowth(
    subjid = d$subjid, param = d$param,
    agedays = d$agedays, sex = d$sex,
    measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  res2 <- cleangrowth(
    subjid = d$subjid, param = d$param,
    agedays = d$agedays, sex = d$sex,
    measurement = d$measurement, id = d$id,
    quietly = TRUE
  )

  expect_identical(as.character(res1$exclude), as.character(res2$exclude))
})
