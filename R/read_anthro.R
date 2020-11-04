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
  weianthro_path <- ifelse(
    is.null(path),
    system.file(file.path("extdata", "weianthro.txt"), package = "growthcleanr"),
    file.path(path, "weianthro.txt")
  )
  lenanthro_path <- ifelse(
    is.null(path),
    system.file(file.path("extdata", "lenanthro.txt"), package = "growthcleanr"),
    file.path(path, "lenanthro.txt")
  )
  bmianthro_path <- ifelse(
    is.null(path),
    system.file(file.path("extdata", "bmianthro.txt"), package = "growthcleanr"),
    file.path(path, "bmianthro.txt")
  )
  growth_cdc_ext_path <- ifelse(
    is.null(path),
    system.file(file.path("extdata", "growthfile_cdc_ext.csv"), package = "growthcleanr"),
    file.path(path, "growthfile_cdc_ext.csv")
  )


  growth_cdc_ext <- fread(growth_cdc_ext_path)

  l <- list(
    with(
      fread(weianthro_path, header = TRUE),
      data.frame(
        src = "WHO",
        param = "WEIGHTKG",
        sex = sex - 1,
        age,
        l,
        m,
        s,
        csdpos = as.double(NA),
        csdneg = as.double(NA)
      )
    ),
    with(
      fread(lenanthro_path, header = TRUE),
      data.frame(
        src = "WHO",
        param = "HEIGHTCM",
        sex = sex - 1,
        age,
        l,
        m,
        s,
        csdpos = as.double(NA),
        csdneg = as.double(NA)
      )
    ),
    with(
      fread(bmianthro_path, header = TRUE),
      data.frame(
        src = "WHO",
        param = "BMI",
        sex = sex - 1,
        age,
        l,
        m,
        s,
        csdpos = as.double(NA),
        csdneg = as.double(NA)
      )
    ),
    with(
      growth_cdc_ext,
      data.frame(
        src = "CDC",
        param = "WEIGHTKG",
        sex,
        age = agedays,
        l = cdc_wt_l,
        m = cdc_wt_m,
        s = cdc_wt_s,
        csdpos = cdc_wt_csd_pos,
        csdneg = cdc_wt_csd_neg
      )
    ),
    with(
      growth_cdc_ext,
      data.frame(
        src = "CDC",
        param = "HEIGHTCM",
        sex,
        age = agedays,
        l = cdc_ht_l,
        m = cdc_ht_m,
        s = cdc_ht_s,
        csdpos = cdc_ht_csd_pos,
        csdneg = cdc_ht_csd_neg
      )
    ),
    with(
      growth_cdc_ext,
      data.frame(
        src = "CDC",
        param = "BMI",
        sex,
        age = agedays,
        l = cdc_bmi_l,
        m = cdc_bmi_m,
        s = cdc_bmi_s,
        csdpos = cdc_bmi_csd_pos,
        csdneg = cdc_bmi_csd_neg
      )
    )
  )

  anthro <- rbindlist(l)


  setkeyv(anthro, c("src", "param", "sex", "age"))

  return(function(param, agedays, sex, measurement, csd = FALSE) {
    # For now, we will only use CDC growth reference data, note that the cubically interpolated file
    # we are using has linear measurments derived from length data for children < 731 days, and height thereafter
    src <- ifelse(agedays < 731 & !cdc.only, "WHO", "CDC")

    # keep column sequence the same fo efficient join
    dt <- data.table(src, param, sex, agedays, measurement)
    dt <- anthro[dt]

    dt[, ret := as.double(NA)]
    if (csd) {
      dt[measurement < m, ret := (measurement - m) / csdneg]
      dt[measurement >= m, ret := (measurement - m) / csdpos]
    } else {
      dt[l == 0, ret := log(measurement / m) / s]
      dt[l != 0, ret := (((measurement / m)^l) - 1) / (l * s)]
    }

    return(dt$ret)
  })
}