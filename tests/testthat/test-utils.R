test_that("splitinput splits files correctly with default values", {
  create_df <- function(num_ids, num_obs) {
    df <- data.frame(matrix(NA, nrow = num_ids * num_obs, ncol = 5))
    colnames(df)[1] <- "subjid"
    df$subjid <- rep(paste0("subj", 1:num_ids), each = num_obs)

    return(df)
  }

  # create input dataframe for default value test
  num_ids <- 1000 * 6
  num_obs <- 7
  df <- create_df(num_ids, num_obs)
  # where to put output
  fcount <- splitinput(df,
                       fname = "df",
                       fdir = tempdir())

  # check that it yielded the correct number of files, with the correct name
  f_log <- grepl("df.*.csv", list.files(tempdir()))
  expect_equal(sum(f_log),
               ceiling(num_ids * num_obs / 10000))

  # check file contents
  sp_list <- lapply(list.files(tempdir())[f_log],
                    function(x) {
                      read.csv(paste0(tempdir(), "/", x))
                    })

  # check that each, except the last, is above the default limit
  # also check that all columns are accounted for
  for (x in sp_list[1:(length(sp_list) - 1)]) {
    expect_gte(nrow(x), 10000)
    expect_equal(ncol(x), 5)
  }

  # check that subjects are not split, between files
  all_subj <- lapply(sp_list, function(x) {
    unique(x$subjid)
  })
  expect_equal(length(Reduce(intersect, all_subj)), 0)

  # remove created csvs
  file.remove(list.files(tempdir(), full.names = T)[f_log])
})

test_that("splitinput splits files correctly with custom values", {
  create_df <- function(num_ids, num_obs) {
    df <- data.frame(matrix(NA, nrow = num_ids * num_obs, ncol = 5))
    colnames(df)[1] <- "subjid"
    df$subjid <- rep(paste0("subj", 1:num_ids), each = num_obs)

    return(df)
  }

  remove_files <- function(f_log) {
    file.remove(list.files(tempdir(), full.names = T)[f_log])
  }

  # create input dataframe
  num_ids <- 1
  num_obs <- 20
  df <- create_df(num_ids, num_obs)
  # run splitinput with new name, less than the default observations
  fcount <- splitinput(df,
                       fname = "onesub",
                       fdir = tempdir())

  # check that it yielded the correct number of files, with the correct name
  f_log <- grepl("onesub.*.csv", list.files(tempdir()))
  expect_equal(sum(f_log),
               ceiling(num_ids * num_obs / 10000))

  # check file contents
  sp_list <- lapply(list.files(tempdir())[f_log],
                    function(x) {
                      read.csv(paste0(tempdir(), "/", x))
                    })

  # check that the only file has all the observations
  expect_equal(nrow(sp_list[[1]]), num_ids * num_obs)
  expect_equal(ncol(sp_list[[1]]), 5)

  # remove created csvs
  remove_files(f_log)

  # try reducing the amount of minimum rows
  fcount <- splitinput(df,
                       min_nrow = 2,
                       fname = "lessrows",
                       fdir = tempdir())

  # check that it did not split the file
  f_log <- grepl("lessrows.*.csv", list.files(tempdir()))
  expect_equal(sum(f_log), 1)

  remove_files(f_log)

  # check that it splits correctly for a given minimum amount of rows
  df <- create_df(2, 10)
  df$X2 <- c(1:(2 * 10)) # creating fake observations
  df <- df[sample(1:nrow(df), nrow(df)), ] # reorder
  fcount <- splitinput(df,
                       min_nrow = 5,
                       fname = "multless",
                       fdir = tempdir())

  # check that it yielded the correct number of files, with the correct name
  f_log <- grepl("multless.*.csv", list.files(tempdir()))
  expect_equal(sum(f_log), 2)

  # check file contents
  sp_list <- lapply(list.files(tempdir())[f_log],
                    function(x) {
                      read.csv(paste0(tempdir(), "/", x))
                    })

  # check that each file is above the limit
  # also check that all columns are accounted for
  for (x in sp_list[1:(length(sp_list))]) {
    expect_gte(nrow(x), 5)
    expect_equal(ncol(x), 5)
  }

  # check that subjects are not split, between files
  all_subj <- lapply(sp_list, function(x) {
    unique(x$subjid)
  })
  expect_equal(length(Reduce(intersect, all_subj)), 0)

  # remove created csvs
  remove_files(f_log)
})

test_that("splitinput throws errors when expected", {
  # run splitinput with several wrong directory names
  expect_error(splitinput(data.frame(), fdir = "hello"))
  expect_error(splitinput(data.frame(), fdir = T))
  expect_error(splitinput(data.frame(), fdir = data.frame()))
})

test_that("recode_sex works as expected with defaults", {
  # create a dataframe according to recode_sex defaults
  num_obs <- 20
  set.seed(7) # for replicability
  df <- data.frame(matrix(NA, nrow = num_obs, ncol = 1))
  colnames(df) <- "sex"
  df$sex <- sample(c("0", "1"), num_obs, replace = T)

  # run recode sex with all defaults
  r_df <- recode_sex(df)

  # check that column names are correct
  expect(all(colnames(r_df) %in% c("sex", "sex_recoded")),
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
  df$sex_type <- sample(c("M", "F"), num_obs, replace = T)

  # run recode sex with all defaults
  r_df <- recode_sex(
    df,
    sourcecol = "sex_type",
    sourcem = "M",
    sourcef = "F",
    targetcol = "sex_r",
    targetm = "m",
    targetf = "f"
  )

  # check that column names are correct
  expect(all(colnames(r_df) %in% c("sex_type", "sex_r")),
         "column names incorrect")

  # all observations are accounted for
  expect_equal(df$sex_type, r_df$sex_type)

  # check that sex was recoded according to specifications
  expect_equal(r_df$sex_r[r_df$sex_type == "M"],
               rep("m", sum(df$sex_type == "M")))
  expect_equal(r_df$sex_r[r_df$sex_type == "F"],
               rep("f", sum(df$sex_type == "F")))
})

test_that("longwide works as expected with default values", {
  # use synthetic data, running cleaning on a subset
  data("syngrowth")
  sub_syn <-
    syngrowth[syngrowth$subjid %in% unique(syngrowth$subjid)[1:100], ]
  sub_syn <- cbind(
    sub_syn,
    "gcr_result" = cleangrowth(
      subjid = sub_syn$subjid,
      param = sub_syn$param,
      agedays = sub_syn$agedays,
      sex = sub_syn$sex,
      measurement = sub_syn$measurement
    )
  )

  # run longwide on changed data
  wide_syn <- longwide(sub_syn)

  # check that it has the correct amount of columns
  expect_equal(ncol(wide_syn), 9)

  # check that all but one subject are accounted for
  ss <- unique(sub_syn$subjid)
  ws <- unique(wide_syn$subjid)
  expect_false(all(ss %in% ws),
               "not all subjects appear in wide format")
  expect_equal(setdiff(ss, ws), c("542abc54-c79f-9895-0350-ded2bf04af6e"))

  obs_ids <- c(wide_syn$wt_id, wide_syn$ht_id)

  # check that all subjects' measurements with at least two occurrences appear
  all_obs <- sapply(unique(sub_syn$subjid), function(i) {
    sub_group <- sub_syn[sub_syn$gcr_result == "Include" &
                           sub_syn$agedays >= 730, ]
    sum(table(sub_group$agedays[sub_group$subjid == i]) >= 2)
  })
  total_obs <- sum(all_obs)
  expect(total_obs <= nrow(wide_syn),
         "observations are dropped in wide format")

  # check that it includes specified inclusion types
  expect(
    all(sub_syn$gcr_result[sub_syn$id %in% obs_ids] == "Include"),
    "longwide() includes inclusion values that were not specified"
  )

  # check that all sexes have been correctly recoded, using values from ws
  # because one will be missing in wide_syn otherwise
  orig_sex <- sub_syn$sex[ws]
  names(orig_sex) <- ws
  aft_sex <- wide_syn$sex[ws]
  names(aft_sex) <- ws

  expect_equal(aft_sex[names(orig_sex)], orig_sex + 1)

  # spot check that data is correct
  set.seed(7)
  # check height ids
  ht_sub <-
    sub_syn[sub_syn$param == "HEIGHTCM" & sub_syn$id %in% obs_ids, ]
  for (x in ht_sub$id[sample(1:nrow(ht_sub), 3)]) {
    w_idx <- wide_syn$ht_id == x
    ht_idx <- ht_sub$id == x

    # check ages
    expect_equal(wide_syn$agey[w_idx], round(ht_sub$agedays[ht_idx] / 365.25), 4)
    expect_equal(wide_syn$agem[w_idx],
                 round(round(ht_sub$agedays[ht_idx] / 365.25), 4) * 12, 4)
    expect_equal(wide_syn$agedays[w_idx], ht_sub$agedays[ht_idx])

    # check height
    expect_equal(wide_syn$ht[w_idx], ht_sub$measurement[ht_idx])

  }

  # check weight ids
  wt_sub <-
    sub_syn[sub_syn$param == "WEIGHTKG" & sub_syn$id %in% obs_ids, ]
  for (x in wt_sub$id[sample(1:nrow(wt_sub), 3)]) {
    w_idx <- wide_syn$wt_id == x
    wt_idx <- wt_sub$id == x

    # check ages
    expect_equal(wide_syn$agey[w_idx], round(wt_sub$agedays[wt_idx] / 365.25), 4)
    expect_equal(wide_syn$agem[w_idx],
                 round(round(wt_sub$agedays[wt_idx] / 365.25), 4) * 12, 4)
    expect_equal(wide_syn$agedays[w_idx], wt_sub$agedays[wt_idx])

    # check height
    expect_equal(wide_syn$wt[w_idx], wt_sub$measurement[wt_idx])
  }

})

test_that("longwide works as expected with custom values", {
  # just checking the custom-ness, so use a smaller subset for speed
  # use synthetic data, running cleaning on a subset
  data("syngrowth")
  sub_syn <-
    syngrowth[syngrowth$subjid %in% unique(syngrowth$subjid)[1:20], ]
  sub_syn <- cbind(
    sub_syn,
    "cv" = cleangrowth(
      subjid = sub_syn$subjid,
      param = sub_syn$param,
      agedays = sub_syn$agedays,
      sex = sub_syn$sex,
      measurement = sub_syn$measurement
    )
  )

  # run longwide on changed data with all exclusion types included
  wide_syn <- longwide(sub_syn,
                       gcr_result = "cv",
                       include_all = T)

  # check that it has the correct amount of columns
  expect_equal(ncol(wide_syn), 9)

  # check that all subjects are accounted for
  expect(all(unique(sub_syn$subjid) %in% unique(wide_syn$subjid)),
         "not all subjects appear in wide format")

  # check that all subjects' measurements with at least two occurrences appear
  all_obs <- sapply(unique(sub_syn$subjid), function(i) {
    sum(table(sub_syn$agedays[sub_syn$subjid == i]) >= 2)
  })
  total_obs <- sum(all_obs)

  expect(total_obs <= nrow(wide_syn),
         "observations are dropped in wide format")


  # run longwide on changed data with some exclusion types included
  inc_types <- c("Include",
                 "Exclude-Carried-Forward",
                 "Exclude-Extraneous-Same-Day")
  wide_syn <- longwide(sub_syn,
                       gcr_result = "cv",
                       inclusion_types = inc_types)

  # check that it has the correct amount of columns
  expect_equal(ncol(wide_syn), 9)

  # check that all subjects are accounted for
  expect(all(unique(sub_syn$subjid) %in% unique(wide_syn$subjid)),
         "not all subjects appear in wide format with custom")

  # check that all subjects' measurements with at least two occurrences appear
  all_obs <- sapply(unique(sub_syn$subjid), function(i) {
    sub_group <- sub_syn[sub_syn$gcr_result %in% inc_types &
                           sub_syn$agedays >= 730, ]
    sum(table(sub_group$agedays[sub_group$subjid == i]) >= 2)
  })
  total_obs <- sum(all_obs)
  expect(total_obs <= nrow(wide_syn),
         "observations are dropped in wide format")

  # check that it includes specified inclusion types
  obs_ids <- c(wide_syn$wt_id, wide_syn$ht_id)
  expect(
    all(sub_syn$gcr_result[sub_syn$id %in% obs_ids] == inc_types),
    "longwide() includes inclusion values that were not specified"
  )

})

test_that("longwide throws errors correctly", {
  # use synthetic data, running cleaning on a very small subset
  data("syngrowth")
  sub_syn <-
    syngrowth[syngrowth$subjid %in% unique(syngrowth$subjid)[1:5], ]
  sub_syn <- cbind(
    sub_syn,
    "gcr_result" = cleangrowth(
      subjid = sub_syn$subjid,
      param = sub_syn$param,
      agedays = sub_syn$agedays,
      sex = sub_syn$sex,
      measurement = sub_syn$measurement
    )
  )

  # test with deleting a necessary column
  expect_error(longwide(sub_syn[, -2]))

  # test include_all not being correct
  expect_error(longwide(sub_syn, include_all = "hello"))

  # test duplicated ids
  sub_syn$id <- 1
  expect_error(longwide(sub_syn))
})

test_that("simple_bmi works as expected", {
  data("syngrowth")

  # Similar strategy as for longwide, create subset for speed
  sub_syn <- syngrowth[syngrowth$subjid %in% unique(syngrowth$subjid)[101:200], ]
  sub_syn <- cbind(
    sub_syn,
    "cv" = cleangrowth(
      subjid = sub_syn$subjid,
      param = sub_syn$param,
      agedays = sub_syn$agedays,
      sex = sub_syn$sex,
      measurement = sub_syn$measurement
    )
  )

  wide_syn <- longwide(sub_syn, gcr_result = "cv", include_all = TRUE)
  bmi_syn <- simple_bmi(wide_syn)
  expect_equal(TRUE, "wt" %in% names(bmi_syn))
  expect_equal(bmi_syn$bmi,
               bmi_syn$wt / ((bmi_syn$ht * .01) ^ 2))

  # Verify that invalid column names throw an error
  expect_error(simple_bmi(wide_syn, ht="invalid_column"))
  expect_error(simple_bmi(wide_syn, wt="invalid_wt_col", ht="invalid_ht_col"))
})
