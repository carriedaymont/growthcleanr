library(growthcleanr)
library(dplyr)
library(data.table)
library(stringr)
library(gt)

source_new_code <- function() {
  files <- list.files("../code/modified_source_code/", pattern = "\\.[rR]$", full.names = TRUE)
  invisible(lapply(files, source))
  tmp <- knitr::purl(
    "/Documents and Settings/penne/Documents/Github/Daymont_Contract/growthcleanr/cp/code/modified_source_code/Contracted_Fix_Main.Rmd",
    output = tempfile()
  )
  source(tmp, chdir = TRUE)
}


data_in <- fread(
  "/Documents and Settings/penne/Documents/Github/Daymont_Contract/growthcleanr/cp/data/2025-11-15-gc-Stata-chopdall-results.csv",
  select = c(
    "subjid", "agedays", "sex", "obsid",
    "exc_ht", "exc_wt", "exc_hc", "measurement", "param"
  )
) %>%
  mutate(
    gcr_result_carrie = dplyr::case_when(
      param == "WEIGHTKG" ~ exc_wt,
      param == "HEIGHTCM" ~ exc_ht,
      param == "HEADCM"   ~ exc_hc
    )
  ) 



# Arrange.
data_in_subs <- data_in %>% select(subjid) %>% distinct() %>% slice_sample(n=50)

data.work <- data_in %>% arrange(subjid, param, agedays, obsid) %>% as.data.table() %>% filter(subjid %in% data_in_subs$subjid)
# Source our modified code.
source_new_code()

setkey(data.work, subjid, param, agedays)



# Run cleangrowth.a
cleaned_data_updated <- cleangrowth(
  subjid = data.work$subjid,
  param = data.work$param,
  agedays = data.work$agedays,
  sex = data.work$sex,
  measurement = data.work$measurement,
  prelim_infants = TRUE,
  id = data.work$obsid, 
  sd.recenter = "nhanes",
  quietly = FALSE
)

write.csv(cleaned_data_updated, paste("/Documents and Settings/penne/Documents/Github/Daymont_Contract/growthcleanr/cp/data/100pct_data_test", Sys.Date(), ".csv", sep = ""))