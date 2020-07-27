library(dplyr)

test_that("growthcleanr works as expected on synthetic data", {

  # Run cleangrowth() on syngrowth data
  data <- as.data.table(syngrowth)
  # syngrowth hasn't changed in length
  expect_equal(56703, data[, .N])
  setkey(data, subjid, param, agedays)
  data_sample <- data[subjid %in% unique(data[, subjid])[1:100], ]
  # Nor has the simplistic sample
  expect_equal(2215, data_sample[, .N])
  cleaned_data <-
    data_sample[, clean_value := cleangrowth(
      subjid,
      param,
      agedays,
      sex,
      measurement
    )]

  # Spot check individual results
  clean_value <- function (df, rowid) {
    return(as.character(df[id == rowid]$clean_value))
  }
  expect_equal("Exclude-EWMA-8", clean_value(data_sample, 9652))
  expect_equal("Exclude-Min-Height-Change", clean_value(data_sample, 31450))
  expect_equal("Include", clean_value(data_sample, 15102))

  # Check counts of exclusions by category
  catcount <- function (df, category) {
    return(as.numeric(df %>% filter(clean_value == category) %>% select(n)))
  }
  exclusions <-
    cleaned_data %>% group_by(clean_value) %>% tally(sort = TRUE)
  expect_equal(1505, catcount(exclusions, "Include"))
  expect_equal(275, catcount(exclusions, "Exclude-Carried-Forward"))
  expect_equal(1, catcount(exclusions, "Exclude-EWMA-11"))

})