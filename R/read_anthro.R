#' Function to calculate z-scores and csd-scores based on anthro tables.
#'
#' @param path Path to supplied reference anthro data. Defaults to package anthro tables.
#' @param cdc.only Whether or not only CDC data should be used. Defaults to false.
#'
#' @return Function for calculating BMI based on measurement, age in days, sex, and measurement value.
#'
#' @export
#' @import data.table
#'
#' @examples
#' # Return calculating function with all defaults
#' afunc <- read_anthro()
#'
#' # Return calculating function while specifying a path and using only CDC data
#' afunc <- read_anthro(
#'   path = system.file("extdata", package = "growthcleanr"),
#'   cdc.only = TRUE
#' )
read_anthro <- function(path = NULL, cdc.only = FALSE) {
  # set correct path based on input reference table path (if any)
  anthro_list <- lapply(
    X = c("weianthro.txt", "lenanthro.txt", "bmianthro.txt"),
    dir_path = path,
    FUN = function(ifile, dir_path) {
      ifile_path <- if (is.null(dir_path)) {
        system.file(file.path("extdata", ifile), package = "growthcleanr")
      } else {
        file.path(dir_path, ifile)
      }

      fread(ifile_path)[, .(
        src = "WHO",
        sex = sex - 1,
        age,
        param = c(
          "weianthro.txt" = "WEIGHTKG",
          "lenanthro.txt" = "HEIGHTCM",
          "bmianthro.txt" = "BMI"
        )[ifile],
        l,
        m,
        s,
        csdpos = NA_real_,
        csdneg = NA_real_
      )]

    }
  )

  cdc_wide <-  fread(
    file = ifelse(
      is.null(path),
      system.file(file.path("extdata", "growthfile_cdc_ext.csv"), package = "growthcleanr"),
      file.path(path, "growthfile_cdc_ext.csv")
    )
  )[,
    c("src", paste("cdc", c("wt", "ht", "hc", "bmi"), "param", sep = "_")) := list(
      "CDC", "WEIGHTKG", "HEIGHTCM", "hc", "BMI"
    )
  ]

  cdc_long <- melt(
    data = cdc_wide,
    id.vars = c("src", "sex", "agedays"),
    measure.vars = patterns(
      param = "_param$",
      l = "_l$",
      m = "_m$",
      s = "_s$",
      csdpos = "_csd_pos$",
      csdneg = "_csd_neg$"
    )
  )[, variable := NULL][param != "hc"]
  setnames(cdc_long, "agedays", "age")

  anthro_list[[length(anthro_list) + 1]] <- cdc_long
  anthro <- rbindlist(anthro_list)

  setkeyv(anthro, c("src", "param", "sex", "age"))

  return(function(param, agedays, sex, measurement, csd = FALSE) {
    # For now, we will only use CDC growth reference data, note that the cubically interpolated file
    # we are using has linear measurments derived from length data for children < 731 days, and height thereafter

    # keep column sequence the same fo efficient join
    dt <- data.table(
      src = fifelse(agedays < 731 & !cdc.only, "WHO", "CDC"),
      param, sex, agedays, measurement
    )
    dt <- anthro[dt]

    if (csd) {
      dt[, fifelse(measurement < m, (measurement - m) / csdneg, (measurement - m) / csdpos)]
    } else {
      dt[, fifelse(l == 0, log(measurement / m) / s, (((measurement / m)^l) - 1) / (l * s))]
    }
  })
}
