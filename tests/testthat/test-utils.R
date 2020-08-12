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