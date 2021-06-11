test_that("growthcleanr works as expected on synthetic data", {

  # Run cleangrowth() on syngrowth data
  data <- as.data.table(syngrowth)

  # syngrowth hasn't changed in length
  expect_equal(77721, data[, .N])
  setkey(data, subjid, param, agedays)

  # subset to pediatric data
  data_peds <- copy(data[agedays < 20 * 365.25, ])

  # Create small samples; one small enough to trigger NHANES recentering,
  # one large enough not to.
  #
  # Note that we're creating distinct data tables to avoid accidentally
  # reusing the same by reference.
  d500_nhanes <- as.data.table(data_peds)[subjid %in% unique(data[, subjid])[1:500], ]
  expect_equal(5347, d500_nhanes[, .N])

  # And for overriding NHANES/derive option
  d500_derive <- as.data.table(data_peds)[subjid %in% unique(data[, subjid])[1:500], ]
  expect_equal(5347, d500_derive[, .N])

  # Clean samples: specify sd.recenter should use NHANES
  cd500_nhanes <-
    d500_nhanes[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      sd.recenter = "NHANES"
    )]

  # Specifying "derive" from data instead of NHANES
  cd500_derived <-
    d500_derive[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      sd.recenter = "derive"
    )]


  # Spot check individual results
  gcr_result <- function (dt, rowid) {
    return(as.character(dt[id == rowid]$gcr_result))
  }

  # Results for these records should not change w/sample size (NHANES vs. derived)
  expect_equal("Exclude-EWMA-8", gcr_result(cd500_nhanes, 10828))
  expect_equal("Exclude-EWMA-8", gcr_result(cd500_derived, 10828))

  expect_equal("Exclude-Min-Height-Change", gcr_result(cd500_nhanes, 29423))
  expect_equal("Exclude-Min-Height-Change", gcr_result(cd500_derived, 29423))

  expect_equal("Include", gcr_result(cd500_nhanes, 57638))
  expect_equal("Include", gcr_result(cd500_derived, 57638))

  # Results for these records can change w/NHANES vs. derived
  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd500_nhanes, 8456))
  expect_equal("Include", gcr_result(cd500_derived, 8456))

  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd500_nhanes, 6364))
  expect_equal("Exclude-Carried-Forward", gcr_result(cd500_derived, 6364))

  expect_equal("Exclude-EWMA-9", gcr_result(cd500_nhanes, 19328))
  expect_equal("Exclude-EWMA-Extreme", gcr_result(cd500_derived, 19328))

  # Check counts of exclusions by category
  catcount <- function (df, category) {
    return(as.numeric(df %>% filter(gcr_result == category) %>% select(n)))
  }

  d500_exclusions <-
    cd500_nhanes %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(3535, catcount(d500_exclusions, "Include"))
  expect_equal(751, catcount(d500_exclusions, "Exclude-Carried-Forward"))
  expect_equal(13, catcount(d500_exclusions, "Exclude-EWMA-8"))

  d500_derived_exclusions <-
    cd500_derived %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(3586, catcount(d500_derived_exclusions, "Include"))
  expect_equal(772, catcount(d500_derived_exclusions, "Exclude-Carried-Forward"))
  expect_equal(21, catcount(d500_derived_exclusions, "Exclude-EWMA-8"))

})

test_that("growthcleanr works without either adult or pediatric data", {
  # creating small only adult and only pediatric data
  # using default cutpoint -- 20
  only_peds <- syngrowth[syngrowth$agedays < 20*365.25,][1:50,]
  only_adult <- syngrowth[syngrowth$agedays >= 20*365.25,][1:50,]

  # testing cleangrowth works without adult data
  peds_res <- cleangrowth(
    only_peds$subjid,
    only_peds$param,
    only_peds$agedays,
    only_peds$sex,
    only_peds$measurement,
    quietly = T
  )

  expect_equal(length(peds_res), nrow(only_peds))

  # testing cleangrowth works without pediatric data
  adult_res <- cleangrowth(
    only_adult$subjid,
    only_adult$param,
    only_adult$agedays,
    only_adult$sex,
    only_adult$measurement,
    quietly = T
  )

  expect_equal(length(adult_res), nrow(only_adult))

})