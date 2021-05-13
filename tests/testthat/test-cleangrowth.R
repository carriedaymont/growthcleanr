test_that("growthcleanr works as expected on synthetic data", {

  # Run cleangrowth() on syngrowth data
  data <- as.data.table(syngrowth)
  # subset to pediatric data
  data_orig <- data <-
    copy(data[suppressWarnings(!is.na(as.numeric(data$subjid))),])
  data[, subjid := as.numeric(subjid)]

  # syngrowth hasn't changed in length
  expect_equal(56703, data[, .N])
  setkey(data, subjid, param, agedays)

  # Create small samples; one small enough to trigger NHANES recentering,
  # one large enough not to.
  #
  # Note that we're creating distinct data tables to avoid accidentally
  # reusing the same by reference.
  d100 <- as.data.table(data_orig)[subjid %in% unique(data[, subjid])[1:100], ]
  expect_equal(2215, d100[, .N])
  d300 <- as.data.table(data_orig)[subjid %in% unique(data[, subjid])[1:300], ]
  expect_equal(6647, d300[, .N])

  # And for overriding NHANES/derive option
  d100_derive <- as.data.table(data_orig)[subjid %in% unique(data[, subjid])[1:100], ]
  expect_equal(2215, d100_derive[, .N])
  d300_nhanes <- as.data.table(data_orig)[subjid %in% unique(data[, subjid])[1:300], ]
  expect_equal(6647, d300_nhanes[, .N])

  # Clean samples: cd100 w/o specifying sd.recenter should use NHANES
  cd100 <-
    d100[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  # cd100 specifying "derive" should derive
  cd100_derived <-
    d100_derive[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      sd.recenter = "derive"
    )]

  # cd300 w/o specifying sd.recenter should derive
  cd300 <-
    d300[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  # cd300_nhanes specifying "nhanes" should use NHANES
  cd300_nhanes <-
    d300_nhanes[, gcr_result := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement,
      sd.recenter = "NHANES"
    )]

  # Spot check individual results
  gcr_result <- function (dt, rowid) {
    return(as.character(dt[id == rowid]$gcr_result))
  }

  # Results for these records should not change w/sample size (NHANES vs. derived)
  expect_equal("Exclude-EWMA-8", gcr_result(cd100, 9652))
  expect_equal("Exclude-EWMA-8", gcr_result(cd100_derived, 9652))
  expect_equal("Exclude-EWMA-8", gcr_result(cd300, 9652))
  expect_equal("Exclude-EWMA-8", gcr_result(cd300_nhanes, 9652))

  expect_equal("Exclude-Min-Height-Change", gcr_result(cd100, 31450))
  expect_equal("Exclude-Min-Height-Change", gcr_result(cd100_derived, 31450))
  expect_equal("Exclude-Min-Height-Change", gcr_result(cd300, 31450))
  expect_equal("Exclude-Min-Height-Change", gcr_result(cd300_nhanes, 31450))

  expect_equal("Include", gcr_result(cd100, 15102))
  expect_equal("Include", gcr_result(cd100_derived, 15102))
  expect_equal("Include", gcr_result(cd300, 15102))
  expect_equal("Include", gcr_result(cd300_nhanes, 15102))

  # Results for these records can change w/NHANES vs. derived
  expect_equal("Include", gcr_result(cd100, 862))
  expect_equal("Include", gcr_result(cd100_derived, 862))
  expect_equal("Exclude-EWMA-8", gcr_result(cd300, 862))
  expect_equal("Include", gcr_result(cd300_nhanes, 862))

  expect_equal("Exclude-EWMA-8", gcr_result(cd100, 3453))
  expect_equal("Exclude-EWMA-Extreme", gcr_result(cd100_derived, 3453))
  expect_equal("Exclude-EWMA-Extreme", gcr_result(cd300, 3453))
  expect_equal("Exclude-EWMA-8", gcr_result(cd300_nhanes, 3453))

  expect_equal("Exclude-EWMA-9", gcr_result(cd100, 38081))
  expect_equal("Include", gcr_result(cd100_derived, 38081))
  expect_equal("Include", gcr_result(cd300, 38081))
  expect_equal("Exclude-EWMA-9", gcr_result(cd300_nhanes, 38081))

  # Check counts of exclusions by category
  catcount <- function (df, category) {
    return(as.numeric(df %>% filter(gcr_result == category) %>% select(n)))
  }

  d100_exclusions <-
    cd100 %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(1503, catcount(d100_exclusions, "Include"))
  expect_equal(275, catcount(d100_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d100_exclusions, "Exclude-EWMA-11"))

  d100_derived_exclusions <-
    cd100_derived %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(1505, catcount(d100_derived_exclusions, "Include"))
  expect_equal(275, catcount(d100_derived_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d100_derived_exclusions, "Exclude-EWMA-11"))

  d300_exclusions <-
    cd300 %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(4540, catcount(d300_exclusions, "Include"))
  expect_equal(809, catcount(d300_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d300_exclusions, "Exclude-EWMA-11"))

  d300_nhanes_exclusions <-
    cd300_nhanes %>% group_by(gcr_result) %>% tally(sort = TRUE)
  expect_equal(4545, catcount(d300_nhanes_exclusions, "Include"))
  expect_equal(809, catcount(d300_nhanes_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d300_nhanes_exclusions, "Exclude-EWMA-11"))

  # Verify that the 100-subject set which defaulted to NHANES and the same 100
  # subjects from the 300-subject set which specified NHANES had the same results
  cd300_nhanes_100 <- cd300_nhanes[subjid %in% unique(data[, subjid])[1:100], ]
  nhanes_combined <- merge(cd100, cd300_nhanes_100,
                           by=c("id", "subjid", "sex", "agedays", "param", "measurement"),
                           suffixes=c(".c100", ".c300n100"))
  expect_equal(0, nhanes_combined[gcr_result.c100 != gcr_result.c300n100, .N])
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