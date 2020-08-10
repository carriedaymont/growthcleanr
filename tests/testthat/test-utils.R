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
                       fname = paste0(tempdir(), "/df"))

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
                       fname = paste0(tempdir(), "/onesub"))

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
                       fname = paste0(tempdir(), "/lessrows"))

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
                       fname = paste0(tempdir(), "/multless"))

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
