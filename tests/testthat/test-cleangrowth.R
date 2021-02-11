test_that("growthcleanr works as expected on synthetic data", {

  # Run cleangrowth() on syngrowth data
  data <- as.data.table(syngrowth)
  # syngrowth hasn't changed in length
  expect_equal(56703, data[, .N])
  setkey(data, subjid, param, agedays)

  # Create two small samples; one small enough to trigger NHANES recentering,
  # one large enough not to.
  d100_sample <- data[subjid %in% unique(data[, subjid])[1:100], ]
  expect_equal(2215, d100_sample[, .N])
  d300_sample <- data[subjid %in% unique(data[, subjid])[1:300], ]
  expect_equal(6647, d300_sample[, .N])

  # Clean both samples
  cd100 <-
    d100_sample[, clean_value := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  cd300 <-
    d300_sample[, clean_value := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  # Spot check individual results
  clean_value <- function (dt, rowid) {
    return(as.character(dt[id == rowid]$clean_value))
  }

  # Results for these records should not change w/sample size (NHANES vs. derived)
  expect_equal("Exclude-EWMA-8", clean_value(cd100, 9652))
  expect_equal("Exclude-EWMA-8", clean_value(cd300, 9652))

  expect_equal("Exclude-Min-Height-Change", clean_value(cd100, 31450))
  expect_equal("Exclude-Min-Height-Change", clean_value(cd300, 31450))

  expect_equal("Include", clean_value(cd100, 15102))
  expect_equal("Include", clean_value(cd300, 15102))

  # Results for these records should change w/sample size (NHANES vs. derived)
  expect_equal("Include", clean_value(cd100, 862))
  expect_equal("Exclude-EWMA-8", clean_value(cd300, 862))

  expect_equal("Exclude-EWMA-8", clean_value(cd100, 3453))
  expect_equal("Exclude-EWMA-Extreme", clean_value(cd300, 3453))

  expect_equal("Exclude-EWMA-9", clean_value(cd100, 38081))
  expect_equal("Include", clean_value(cd300, 38081))

  # Check counts of exclusions by category
  catcount <- function (df, category) {
    return(as.numeric(df %>% filter(clean_value == category) %>% select(n)))
  }

  d100_exclusions <-
    cd100 %>% group_by(clean_value) %>% tally(sort = TRUE)
  expect_equal(1503, catcount(d100_exclusions, "Include"))
  expect_equal(275, catcount(d100_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d100_exclusions, "Exclude-EWMA-11"))

  d300_exclusions <-
    cd300 %>% group_by(clean_value) %>% tally(sort = TRUE)
  expect_equal(4540, catcount(d300_exclusions, "Include"))
  expect_equal(809, catcount(d300_exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(d300_exclusions, "Exclude-EWMA-11"))

})
