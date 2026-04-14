testthat::skip_on_cran()
library(growthcleanr)
library(data.table)

# =============================================================================
# Tests for cleangrowth() API: legacy algorithm, adult algorithm, edge cases.
#
# Updated for v3.0.0: cleangrowth() returns a data.table (not a vector).
# Legacy pediatric tests use use_legacy_algorithm = TRUE.
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

test_that("growthcleanr legacy algorithm works on pediatric synthetic data", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  data <- as.data.table(syngrowth)
  expect_equal(77721, data[, .N])
  setkey(data, subjid, param, agedays)

  data_peds <- copy(data[agedays < 20 * 365.25, ])

  d100_nhanes <- copy(data_peds[subjid %in% unique(data[, subjid])[1:100]])
  expect_equal(832, d100_nhanes[, .N])

  d100_derive <- copy(data_peds[subjid %in% unique(data[, subjid])[1:100]])
  expect_equal(832, d100_derive[, .N])

  # Run legacy algorithm with NHANES recentering
  cd100_nhanes <- run_and_merge(d100_nhanes,
                                 sd.recenter = "NHANES",
                                 use_legacy_algorithm = TRUE)

  # Run legacy algorithm with derived recentering
  cd100_derived <- run_and_merge(d100_derive,
                                  sd.recenter = "derive",
                                  use_legacy_algorithm = TRUE)

  # Spot check individual results
  gcr_result <- function(dt, rowid) {
    as.character(dt[id == rowid]$gcr_result)
  }

  # Results for these records should not change w/sample size
  expect_equal("Exclude-EWMA-8", gcr_result(cd100_nhanes, 35119))
  expect_equal("Exclude-EWMA-8", gcr_result(cd100_derived, 35119))

  expect_equal("Exclude-Min-Height-Change", gcr_result(cd100_nhanes, 38718))
  expect_equal("Exclude-Min-Height-Change", gcr_result(cd100_derived, 38718))

  expect_equal("Include", gcr_result(cd100_nhanes, 23766))
  expect_equal("Include", gcr_result(cd100_derived, 23766))

  # Results for these records can change w/NHANES vs. derived
  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd100_nhanes, 25))
  expect_equal("Include", gcr_result(cd100_derived, 25))

  expect_equal("Exclude-Carried-Forward", gcr_result(cd100_nhanes, 40094))
  expect_equal("Exclude-Carried-Forward", gcr_result(cd100_derived, 40094))

  expect_equal("Include", gcr_result(cd100_nhanes, 62606))
  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd100_derived, 62606))

  # Check counts of exclusions by category
  catcount <- function(dt, category) {
    dt[gcr_result == category, .N]
  }

  expect_equal(563, catcount(cd100_nhanes, "Include"))
  expect_equal(113, catcount(cd100_nhanes, "Exclude-Carried-Forward"))
  expect_equal(3, catcount(cd100_nhanes, "Exclude-EWMA-8"))

  expect_equal(563, catcount(cd100_derived, "Include"))
  expect_equal(113, catcount(cd100_derived, "Exclude-Carried-Forward"))
  expect_equal(3, catcount(cd100_derived, "Exclude-EWMA-8"))

  # Missing data test
  d5_miss <- copy(data_peds[subjid %in% unique(data[, subjid])[1:5]])
  set.seed(10)
  d5_miss$measurement[sample(seq_len(nrow(d5_miss)), 5)] <- NA

  res_miss <- cleangrowth(
    subjid = d5_miss$subjid,
    param = d5_miss$param,
    agedays = d5_miss$agedays,
    sex = d5_miss$sex,
    measurement = d5_miss$measurement,
    id = d5_miss$id,
    use_legacy_algorithm = TRUE,
    quietly = TRUE
  )

  expect_equal(sum(res_miss$exclude == "Exclude-Missing"), 5)
})

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

  expect_equal("Exclude-A-HT-Identical", gcr_result(cd100, 167))
  expect_equal("Exclude-A-HT-Identical", gcr_result(cd100cp, 167))

  expect_equal("Exclude-A-HT-BIV", gcr_result(cd100, 22009))
  expect_equal("Exclude-A-HT-BIV", gcr_result(cd100cp, 22009))

  # Results should change due to younger cutpoint
  # With default cutpoint=20, records 18-20 yrs go through child algorithm
  expect_equal("Exclude-C-WT-Extraneous", gcr_result(cd100, 69740))
  expect_equal("Include", gcr_result(cd100cp, 69740))

  expect_equal("Include", gcr_result(cd100cp, 55171))

  expect_equal("Exclude-C-HT-BIV", gcr_result(cd100, 25259))
  expect_equal("Exclude-A-HT-BIV", gcr_result(cd100cp, 25259))

  # Check counts (codes are now param-specific; use grepl for category totals)
  catcount <- function(dt, pattern) {
    dt[grepl(pattern, gcr_result), .N]
  }
  catcount_exact <- function(dt, category) {
    dt[gcr_result == category, .N]
  }

  expect_equal(1592, catcount_exact(cd100, "Include"))
  expect_equal(334, catcount(cd100, "Exclude-A-(HT|WT)-Extraneous"))
  expect_equal(8, catcount(cd100, "-CF$"))
  expect_equal(16, catcount(cd100, "Exclude-A-(HT|WT)-BIV"))

  expect_equal(1595, catcount_exact(cd100cp, "Include"))
  expect_equal(356, catcount(cd100cp, "Exclude-A-(HT|WT)-Extraneous"))
  expect_equal(0, catcount(cd100cp, "-CF$"))
  expect_equal(19, catcount(cd100cp, "Exclude-A-(HT|WT)-BIV"))

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

test_that("prelim_infants = TRUE runs child algorithm with deprecation warning", {

  data("syngrowth", package = "growthcleanr", envir = environment())
  dt <- as.data.table(syngrowth)
  setkey(dt, subjid, param, agedays)
  data_peds <- dt[agedays < 20 * 365.25]
  d5 <- data_peds[subjid %in% unique(dt$subjid)[1:5]]
  expect_equal(42, nrow(d5))

  # prelim_infants = TRUE should warn and run child algorithm
  expect_warning(
    res <- cleangrowth(
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

  expect_equal(nrow(res), 42)
  # Result should have valid exclusion codes
  expect_true(is.factor(res$exclude))
  expect_false(any(is.na(res$exclude)))
})
