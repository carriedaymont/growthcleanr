test_that("growthcleanr works as expected on pediatric synthetic data", {

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

test_that("growthcleanr works as expected on adult synthetic data", {

  # Run cleangrowth() on syngrowth data
  data <- as.data.table(syngrowth)

  # syngrowth hasn't changed in length
  expect_equal(77721, data[, .N])
  setkey(data, subjid, param, agedays)

  # subset to adult data
  data_adult <- copy(data[agedays >= 18 * 365.25, ])

  # Create small sample
  d500 <- as.data.table(data_adult)[subjid %in% unique(data_adult[, subjid])[1:500], ]
  expect_equal(12447, d500[, .N])

  # Clean sample
  cd500 <-
    d500[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  # Clean again with lower cutpoint
  cd500cp <-
    copy(d500)[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      adult_cutpoint = 18
    )]


  # Spot check individual results
  gcr_result <- function (dt, rowid) {
    return(as.character(dt[id == rowid]$gcr_result))
  }

  # These results should not change with cutpoint
  expect_equal("Include", gcr_result(cd500, 27166))
  expect_equal("Include", gcr_result(cd500cp, 27166))

  expect_equal("Exclude-Adult-Identical-Same-Day", gcr_result(cd500, 47596))
  expect_equal("Exclude-Adult-Identical-Same-Day", gcr_result(cd500cp, 47596))

  expect_equal("Exclude-Adult-BIV", gcr_result(cd500, 41872))
  expect_equal("Exclude-Adult-BIV", gcr_result(cd500cp, 41872))

  # Results for these records should change due to younger cutpoint
  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd500, 38722))
  expect_equal("Exclude-Adult-Extraneous-Same-Day", gcr_result(cd500cp, 38722))

  expect_equal("Exclude-Carried-Forward", gcr_result(cd500, 12923))
  expect_equal("Include", gcr_result(cd500cp, 12923))

  expect_equal("Exclude-Extraneous-Same-Day", gcr_result(cd500, 25259))
  expect_equal("Exclude-Adult-Distinct-3-Or-More", gcr_result(cd500cp, 25259))

  # Check counts of exclusions by category
  catcount <- function (df, category) {
    return(as.numeric(df %>% filter(gcr_result == category) %>% select(n)))
  }

  d500_exclusions <-
    cd500 %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(9745, catcount(d500_exclusions, "Include"))
  expect_equal(2090, catcount(d500_exclusions, "Exclude-Adult-Extraneous-Same-Day"))
  expect_equal(59, catcount(d500_exclusions, "Exclude-Adult-Distinct-3-Or-More"))
  expect_equal(43, catcount(d500_exclusions, "Exclude-Carried-Forward"))
  expect_equal(2, catcount(d500_exclusions, "Exclude-Adult-Transpositions"))

  d500cp_exclusions <-
    cd500cp %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(9774, catcount(d500cp_exclusions, "Include"))
  expect_equal(2200, catcount(d500cp_exclusions, "Exclude-Adult-Extraneous-Same-Day"))
  expect_equal(62, catcount(d500cp_exclusions, "Exclude-Adult-Distinct-3-Or-More"))
  expect_true(is.na(catcount(d500cp_exclusions, "Exclude-Carried-Forward")))
  expect_equal(2, catcount(d500cp_exclusions, "Exclude-Adult-Transpositions"))

})

test_that("growthcleanr works without either adult or pediatric data", {
  # creating small only adult and only pediatric data
  # using default cutpoint -- 20
  only_peds <- syngrowth[syngrowth$agedays < 20*365.25,][1:50,]
  only_adult <- syngrowth[syngrowth$agedays >= 20*365.25,][1:50,]
  nobody <- syngrowth[syngrowth$agedays > 120*365.25,]

  # testing cleangrowth works without adult data
  peds_res <- cleangrowth(
    only_peds$subjid,
    only_peds$param,
    only_peds$agedays,
    only_peds$sex,
    only_peds$measurement,
    quietly = TRUE
  )

  expect_equal(length(peds_res), nrow(only_peds))

  # testing cleangrowth works without pediatric data
  adult_res <- cleangrowth(
    only_adult$subjid,
    only_adult$param,
    only_adult$agedays,
    only_adult$sex,
    only_adult$measurement,
    quietly = TRUE
  )

  expect_equal(length(adult_res), nrow(only_adult))

  # testing cleangrowth works with no data
  no_res <- cleangrowth(
    nobody$subjid,
    nobody$param,
    nobody$agedays,
    nobody$sex,
    nobody$measurement,
    quietly = TRUE
  )

  expect_equal(length(no_res), nrow(nobody))

})
