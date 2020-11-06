test_that("adjustcarryforward correctly reinserts strings of carried forwards using specific examples", {
  # load Carrie's specific example data
  ex_df <- read.csv(file.path("inst", "testdata", "test_acf_specific_examples.csv"))

  # assign to vectors for ease
  subjid <- samp_a$subjid
  param <- samp_a$param
  agedays <- samp_a$agedays
  sex <- samp_a$sex
  measurement <- samp_a$measurement

  # run growthcleanr using all defaults
  orig.exclude <- cleangrowth(subjid,
                              param,
                              agedays,
                              sex,
                              measurement)

  # run adjustcarriedforward

})