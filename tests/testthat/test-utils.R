test_that("recode_sex works as expected with defaults", {

  # create a dataframe according to recode_sex defaults
  num_obs <- 20
  set.seed(7) # for replicability
  df <- data.frame(matrix(NA, nrow = num_obs, ncol = 1))
  colnames(df) <- "sex"
  df$sex <- sample(c("0","1"), num_obs, replace = T)

  # run recode sex with all defaults
  r_df <- recode_sex(df)

  # check that column names are correct
  expect(all(colnames(r_df) %in% c("sex","sex_recoded")),
         "column names incorrect")

  # all observations are accounted for
  expect_equal(df$sex, r_df$sex)

  # check that sex was recoded according to specifications
  expect_equal(r_df$sex_recoded[r_df$sex == "0"], rep(1, sum(df$sex == "0")))
  expect_equal(r_df$sex_recoded[r_df$sex == "1"], rep(2, sum(df$sex == "1")))
})

test_that("recode_sex works as expected with custom inputs", {

  # create a dataframe with customization
  num_obs <- 31
  set.seed(7) # for replicability
  df <- data.frame(matrix(NA, nrow = num_obs, ncol = 1))
  colnames(df) <- "sex_type"
  df$sex_type <- sample(c("M","F"), num_obs, replace = T)

  # run recode sex with all defaults
  r_df <- recode_sex(df,
                     sourcecol = "sex_type",
                     sourcem = "M",
                     sourcef = "F",
                     targetcol = "sex_r",
                     targetm = "m",
                     targetf = "f")

  # check that column names are correct
  expect(all(colnames(r_df) %in% c("sex_type","sex_r")),
         "column names incorrect")

  # all observations are accounted for
  expect_equal(df$sex_type, r_df$sex_type)

  # check that sex was recoded according to specifications
  expect_equal(r_df$sex_r[r_df$sex_type == "M"],
               rep("m", sum(df$sex_type == "M")))
  expect_equal(r_df$sex_r[r_df$sex_type == "F"],
               rep("f", sum(df$sex_type == "F")))
})