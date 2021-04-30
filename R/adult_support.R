# Supporting Adult growthcleanr functions
# Supporting functions for adult piece of algorithm, ordered by step first used (if
# not a convenience function, or EWMA)

# convenience functions ----

#' convenience function -- see if numeric vector falls between two numbers
#' returns boolean vector
#' @keywords internal
#' @noRd
check_between <- function(vect, num_low, num_high, incl = T){
  return(
    if (incl){
      vect <= num_high & vect >= num_low
    } else {
      vect < num_high & vect > num_low
    }
  )
}

#' convenience function -- round to the nearest .x
#' @keywords internal
#' @noRd
round_pt <- function(val, pt){
  return(round(val/pt)*pt)
}

#' convenience function to get remainders for floats (a/b)
#' @keywords internal
#' @noRd
get_float_rem <- function(a, b){
  return(abs(round(a/b) - (a/b)))
}

# EWMA functions ----

#' function to calculate as delta matrix for adults
#' @keywords internal
#' @noRd
as.matrix.delta_dn <- function(agedays) {
  n <- length(agedays)
  delta <- abs(matrix(rep(agedays, n), n, byrow = T) - agedays)

  return(delta)
}

#' Exponentially Weighted Moving Average (EWMA) (daymont implementation)
#'
#' \code{ewma} calculates the exponentially weighted moving average (EWMA) for a set of numeric observations over time.
#'
#' @param agedays Vector of age in days for each z score (potentially transformed to adjust weighting).
#'
#' @param z Input vector of numeric MEASUREMENT data.
#'
#' @param ewma.exp Exponent to use for weighting.
#'
#' @param ewma.adjacent Specify whether EWMA values excluding adjacent measurements should be calculated.  Defaults to TRUE.
#'
#' @return Data frame with 3 variables:
#' * The first variable (ewma.all) contains the EWMA at observation time
#'   excluding only the actual observation for that time point.
#' * The second variable (ewma.before) contains the EWMA for each observation excluding both the actual observation
#'   and the immediate prior observation.
#' * The third variable (ewma.after) contains the EWMA for each observation excluding both the actual observation
#'   and the subsequent observation.
#' @keywords internal
#' @noRd
ewma_dn <- function(agedays, z, ewma.exp = 5, ewma.adjacent = T) {
  # 6.  EWMA calculation description: Most of the next steps will involve calculating the exponentially weighted moving average for each subject and parameter. I will
  #     describe how to calculate EWMASDs, and will describe how it needs to be varied in subsequent steps.
  # a.	The overall goal of the EWMASD calculation is to identify the difference between the SD-score and what we might predict that DS-score should be, in order to
  #     determine whether it should be excluded.
  # b.	Only nonmissing SD-scores for a parameter that are not designated for exclusion are included in the following calculations.
  # c.	For each SD-score SDi and associated agedaysi calculate the following for every other z-score (SDj...SDn) and associated agedays (agedaysj...agedaysn)  for the
  #     same subject and parameter
  #   i.	(delta)Agej=agedaysj-agedaysi
  #   ii.	EWMAZ=SDi=[(sigma)j->n(SDj*((5+(delta)Agej)^-1.5))]/[ (sigma)j->n((5+(delta)Agej)^-1.5)]
  #   iii.	For most EWMASD calculations, there are 3 EWMASDs that need to be calculated. I will note if not all of these need to be done for a given step.
  #     1.	EWMASDall calculated as above
  #     2.	EWMAZbef calculated excluding the SD-score just before the SD-score of interest (sorted by agedays). For the first observation for a parameter for a
  #         subject, this should be identical to EWMASDall rather than missing.
  #     3.	EWMAZaft calculated excluding the z-score just after the SD-score of interest (sorted by agedays). For the lastobservation for a parameter for a subject,
  #         this should be identical to EWMASDall rather than missing.
  #   iv.	For each of the three EWMASDs, calculate the dewma_*=SD-EWMASD
  # d.	EWMASDs and (delta)EWMASDs will change if a value is excluded or manipulated using one of the methods below, therefore EWMASDs and (delta)EWMASDs be recalculated for each
  #     step where they are needed.
  # e.	For these calculations, use variables that allow for precise storage of numbers (in Stata this is called 'double') because otherwise rounding errors can cause
  #     problems in a few circumstances

  n <- length(agedays)
  # initialize response variables
  ewma.all <- ewma.before <- ewma.after <- vector('numeric', 0)
  if (n > 0) {
    # organize into data frame and sort into order of increasing age,
    # but retain original sort order information in index
    if (!all(agedays == cummax(agedays)))
      warning("EWMA ordering is not sorted; double check") #add in a check to make sure the inputs are already sorted (they should be)
    index <- order(agedays)

    # calculate matrix of differences in age, and add 5 to each delta per Daymont algorithm
    delta <- as.matrix.delta_dn(agedays)
    delta <- ifelse(delta == 0, 0, (delta) ^ ewma.exp)

    # calculate EWMAs, and return in order of original data
    ewma.all[index] <- delta %*% z / apply(delta, 1, sum)

    if (ewma.adjacent) {
      if (n > 2) {
        delta2 = delta
        delta2[col(delta2) == row(delta2) - 1] = 0
        ewma.before[index] = delta2 %*% z / apply(delta2, 1, sum)
        delta3 = delta
        delta3[col(delta3) == row(delta3) + 1] = 0
        ewma.after[index] = delta3 %*% z / apply(delta3, 1, sum)
      } else {
        ewma.before <- ewma.after <- ewma.all
      }
    }
  }
  # return all 3 EWMAs as a data frame
  return(if (ewma.adjacent)
    data.frame(ewma.all, ewma.before, ewma.after)
    else
      data.frame(ewma.all))
}

# step 1w, W BIV ----

#' function to remove BIVs, based on cutoffs for the given method
#' inputs:
#'   subj_df: data frame with measurements of a given type
#'   type: height, weight, or bmi
#'   biv_df: data frame with BIV cutoffs for the given type
#'   include: default F, whether or not to include the endpoints
#' outputs:
#'   logical, true if the given record should be removed due to being a BIV
#' @keywords internal
#' @noRd
remove_biv <- function(subj_df, type, biv_df, include = F){
  too_low <- remove_biv_low(subj_df, type, biv_df, include)
  too_high <- remove_biv_high(subj_df, type, biv_df, include)

  return(too_low | too_high)
}

#' function to remove only the low end of BIVs, based on cutoffs for the given
#' method (for intermediate processing only)
#' inputs:
#'   subj_df: data frame with measurements of a given type
#'   type: height, weight, or bmi
#'   biv_df: data frame with BIV cutoffs for the given type
#'   include: default F, whether or not to include the endpoints
#' outputs:
#'   logical, true if the given record should be removed due to being a BIV
#' @keywords internal
#' @noRd
remove_biv_low <- function(subj_df, type, biv_df, include = F){
  if (!include){
    too_low <- subj_df$measurement < biv_df[type, "low"]
  } else {
    too_low <- subj_df$measurement <= biv_df[type, "low"]
  }

  return(too_low)
}

#' function to remove only the high end of BIVs, based on cutoffs for the given
#' method (for intermediate processing only)
#' inputs:
#'   subj_df: data frame with measurements of a given type
#'   type: height, weight, or bmi
#'   biv_df: data frame with BIV cutoffs for the given type
#'   include: default F, whether or not to include the endpoints
#' outputs:
#'   logical, true if the given record should be removed due to being a BIV
#' @keywords internal
#' @noRd
remove_biv_high <- function(subj_df, type, biv_df, include = F){
  if (!include){
    too_high <- subj_df$measurement > biv_df[type, "high"]
  } else {
    too_high <- subj_df$measurement >= biv_df[type, "high"]
  }

  return(too_high)
}

# step 2w, W repeated values ----

#' Function to identify repeated values in WEIGHT values
#' adds two colummns to w_subj_df:
#'   is_first_rv: is it the first repeated value?
#'   is_rv: is it a repeated value that is not a first rv
#' @keywords internal
#' @noRd
identify_rv <- function(w_subj_df){
  if (nrow(w_subj_df) > 0){
    # follows a similar process to temp sde (step 3)
    # identify which of these have duplicate values
    tab_vals <- table(w_subj_df$meas_m)
    dup_vals <- names(tab_vals)[tab_vals > 1] # coerces to character

    # if there are any duplicate days, we want to identify the first
    if (length(dup_vals) > 0){
      # flag repeated values, as well as if it's a first value
      w_subj_df$is_first_rv <- w_subj_df$is_rv <- F
      for (dv in dup_vals){
        first_rv <- which(as.character(w_subj_df$meas_m) == dv)[1]
        w_subj_df$is_first_rv[first_rv] <- T
        w_subj_df$is_rv[which(as.character(w_subj_df$meas_m) == dv)[-1]] <- T
      }
    } else {
      # no repeated values
      w_subj_df$is_first_rv <- w_subj_df$is_rv <- F
    }
  }

  return(w_subj_df)
}

# step 3, temp extraneous ----

#' Function to identify temporary same days extraneous
#' adds a column to subj_df : "extraneous" designating whether or not the row is
#' temporarily extraneous (not to be considered in the future)
#' ptype: height or weight, weight excludes repeated values
#' @keywords internal
#' @noRd
temp_sde <- function(subj_df, ptype = "height"){
  # identify which of these have duplicate days
  tab_days <- table(subj_df$age_days)
  dup_days <- names(tab_days)[tab_days > 1] # coerces to character

  if (nrow(subj_df) >= 2){
    # can't do this without having no duplicate days
    med_wo <-
      if (sum(!as.character(subj_df$age_days) %in% dup_days) > 0){
        # get the median without duplicate days
        median(subj_df$measurement[
          !as.character(subj_df$age_days) %in% dup_days &
            if (ptype == "weight"){
              !subj_df$is_rv
            } else {
              rep(T, nrow(subj_df))
            }
        ])
        # if it's weight, we also don't want to include RV in median calculation

        # distribute that median out duplicate days
      } else {
        # for those with other nonduplicate parameters, make the median 0
        0
      }

    subj_df$diff <- NA
    subj_df$diff[as.character(subj_df$age_days) %in% dup_days] <-
      abs(subj_df$measurement[as.character(subj_df$age_days) %in% dup_days] -
            med_wo)

    # flag extraneous values on the same day that is not the minimum distance
    # from the median sd score
    subj_df$extraneous <- F
    for (dd in dup_days){
      minz <-which.min(subj_df$diff[as.character(subj_df$age_days) == dd])
      subj_df$extraneous[as.character(subj_df$age_days) == dd][-minz] <-
        T
    }

  } else if (nrow(subj_df) > 0){
    # nothing is duplicate
    subj_df$extraneous <- F
  }

  return(subj_df)
}

# function to redo repeated values after you've done a temporary same day extraneous
# w_subj_df: weight subject data.table
#' @keywords internal
#' @noRd
redo_identify_rv <- function(w_subj_df){
  # redo RVs just if any first RVs became extraneous
  if (nrow(w_subj_df) > 0 & any(w_subj_df$extraneous & w_subj_df$is_first_rv)){
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous,])
    inc_df <- identify_rv(inc_df)

    # reassign new rvs to weight df -- ordered the same way
    w_subj_df$is_first_rv <- w_subj_df$is_rv <- F
    w_subj_df$is_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_rv
    w_subj_df$is_first_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_first_rv
  }

  return(w_subj_df)
}

# step 5, hundreds ----

#' function to identify hundred exclusions
#' inc_df: subj_df with temp extraneous/first rv removed
#' dewma: delta ewma, for metric
#' meas_col: meas_im or meas_m
#' hundreds: 100/200/300, etc.
#' ptype: param type, "weight" or "height"
#' mtype: m or imp (metric or imperial)
#' returns criteria (true if implausible)
#' @keywords internal
#' @noRd
rem_hundreds <- function(inc_df, dewma, meas_col, hundreds, ptype = "weight"){
  # calculate difference between values -- ENDS ARE PROTECTED ON EITHER SIDE
  inc_df$diff_prev <- c(NA, diff(inc_df[,..meas_col]))
  inc_df$diff_next <- c(diff(inc_df[,..meas_col]), NA)

  # state upper and lower limits (hundreds +/- 2)
  # modifier for height vs weight
  modifier <- if (ptype == "height"){
    5.08
  } else {
    2
  }
  # conversion for imperial for weight (none for height)
  div_modifier <- if (grepl("_m", meas_col) | ptype == "height"){
    1
  } else {
    2.2046226
  }
  # these are metric limits
  llimit <- (hundreds / div_modifier) - modifier
  ulimit <- (hundreds / div_modifier) + modifier
  # these are imperial limits
  llimit_imp <- if (ptype == "height" | grepl("_m", meas_col)){
    llimit
  } else {
    hundreds - 2*2.2046226
  }
  ulimit_imp <- if (ptype == "height" | grepl("_m", meas_col)){
    ulimit
  } else {
    hundreds + 2*2.2046226
  }

  # identify hundred exclusions -- only hundred down
  exc_up <-
    (check_between(dewma$dewma.all, llimit, ulimit)) &
    (check_between(dewma$dewma.before, llimit, ulimit)) &
    (check_between(dewma$dewma.after, llimit, ulimit)) &
    (check_between(inc_df$diff_prev, llimit_imp, ulimit_imp) |
       check_between(inc_df$diff_next, llimit_imp, ulimit_imp))
  exc_down <-
    (check_between(dewma$dewma.all, -ulimit, -llimit)) &
    (check_between(dewma$dewma.before, -ulimit, -llimit)) &
    (check_between(dewma$dewma.after, -ulimit, -llimit)) &
    (check_between(inc_df$diff_prev, -ulimit_imp, -llimit_imp) |
       check_between(inc_df$diff_next, -ulimit_imp, -llimit_imp))
  exc_hundred <- if (ptype == "height"){
    exc_down
  } else {
    exc_down | exc_up
  }
  # end criteria depends on the number of distinct values
  criteria <-
    if ((ptype == "height" & length(unique(inc_df$meas_m)) > 2) &
        (ptype == "weight" & length(inc_df$meas_m) > 2)){
      exc_hundred
    } else {
      exc_hundred &
        if (ptype == "height"){
          inc_df$meas_m < 100
        } else {
          inc_df$meas_m < 40 | inc_df$meas_m > 182
        }
    }
  criteria[is.na(criteria)] <- F

  return(criteria)
}

# step 6, unit errors ----

#' function to remove unit errors
#' inc_df: subj_df with temp extraneous/first rv removed
#' ptype: height or weight
#' returns criteria (true if implausible)
#' @keywords internal
#' @noRd
rem_unit_errors <- function(inc_df, ptype = "height"){
  # add "unit error": metric encoded as imperial
  inc_df$ue <- inc_df$meas_m * (if (ptype == "height"){ 2.54 } else {1/2.2046226})

  # calculate ewma (using metric)
  ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)
  dewma <- (inc_df$meas_m- ewma_res)
  # delta ewma with unit error
  absdewma_ue <- abs(ewma_res-inc_df$ue)
  colnames(dewma) <- colnames(absdewma_ue) <-
    paste0("d",colnames(ewma_res))

  # calculate difference between values
  inc_df$abs_ue_prev <- c(NA, abs(inc_df$ue[2:nrow(inc_df)] -
                                    inc_df$meas_m[1:(nrow(inc_df)-1)]))
  inc_df$abs_ue_next <- c(abs(inc_df$ue[1:(nrow(inc_df)-1)] -
                                inc_df$meas_m[2:nrow(inc_df)]), NA)

  # identify unit error exclusions
  exc_ue <- if (ptype == "height"){
    dewma$dewma.all < -80 &
      dewma$dewma.before < -80 &
      dewma$dewma.after < -80 &
      (absdewma_ue$dewma.all <= 2.54 | inc_df$abs_ue_prev <= 2.54 |
         inc_df$abs_ue_next <= 2.54)
  } else {
    dewma$dewma.all > 40 &
      dewma$dewma.before > 40 &
      dewma$dewma.after > 40 &
      (absdewma_ue$dewma.all <= 2 | inc_df$abs_ue_prev <= 2 |
         inc_df$abs_ue_next <= 2)
  }

  # end criteria depends on the number of distinct values
  criteria <-
    if ((ptype == "height" & length(unique(inc_df$meas_m)) > 2) &
        (ptype == "weight" & length(inc_df$meas_m) > 2)){
      exc_ue
    } else {
      exc_ue &
        if (ptype == "height"){
          inc_df$meas_m < 100
        } else {
          inc_df$meas_m > 182
        }
    }
  criteria[is.na(criteria)] <- F

  return(criteria)
}

# step 7h, transpositions ----

# get ones/tens places of a given number
# returns vector of desired digits
#' @keywords internal
#' @noRd
get_num_places <- function(num, place){
  place_map <- c(
    "ones" = 1,
    "tens" = 2
  )

  # there should always be a tens place -- BIVs filter out single digits
  res <- sapply(strsplit(as.character(num), ""), function(x){
    if ("." %in% x){
      as.numeric(x[which(x == ".") - place_map[place]])
    } else {
      as.numeric(x[length(x) + 1 - place_map[place]])
    }
  })

  return(res)
}

#' function to switch the ones and tens place digit of a number
#' returns vector of switched numbers
#' @keywords internal
#' @noRd
switch_tens_ones <- function(num){
  # gets 10s and ones digits
  tens <- get_num_places(num, "tens")
  ones <- get_num_places(num, "ones")

  # get everything to the left of the tens place
  left_num <- sapply(strsplit(as.character(num), ""), function(x){
    if (length(x) == 2){
      ""
    } else if ("." %in% x){
      paste(x[1:(which(x == ".") - 3)], collapse = "")
    } else {
      paste(x[1:(length(x) + 1 - 3)], collapse = "")
    }
  })
  # get everything to right of ones place
  right_num <- sapply(strsplit(as.character(num), ""), function(x){
    if ("." %in% x){
      paste(x[(which(x == ".")):length(x)], collapse = "")
    } else {
      ""
    }
  })

  # return number with tens and ones place switched
  return(as.numeric(paste0(left_num, ones, tens, right_num)))
}

#' function to calculate transpositions
#' inc_df: subj_df with temp extraneous/first rv removed
#' ptype: height or weight
#' returns criteria (true if implausible)
#' @keywords internal
#' @noRd
rem_transpositions <- function(inc_df, ptype = "height"){
  # calculate ewma (using metric)
  ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)
  dewma <- (inc_df$meas_m- ewma_res)
  colnames(dewma) <- paste0("d",colnames(ewma_res))

  criteria <- rep(F, nrow(inc_df))
  for (mtype in c("m", "im")){
    inc_df$transpo <- switch_tens_ones(
      unlist(inc_df[, paste0("meas_", mtype), with = F])
      )
    # if imperial, we want to convert to metric
    if (mtype == "im"){
      inc_df$transpo <- inc_df$transpo *
        (if (ptype == "height"){ 2.54 } else {1/2.2046226})
    }

    inc_df$ones <- get_num_places(
      unlist(inc_df[, paste0("meas_", mtype), with = F]), "ones"
    )
    inc_df$tens <- get_num_places(
      unlist(inc_df[, paste0("meas_", mtype), with = F]), "tens"
    )
    absdewma_transpo <- abs(inc_df$transpo - ewma_res)
    colnames(absdewma_transpo) <- paste0("d",colnames(ewma_res))

    # calculate difference between values
    inc_df$abs_transpo_prev <- c(NA, abs(inc_df$transpo[2:nrow(inc_df)] -
                                           inc_df$meas_m[1:(nrow(inc_df)-1)]))
    inc_df$abs_transpo_next <- c(abs(inc_df$transpo[1:(nrow(inc_df)-1)] -
                                       inc_df$meas_m[2:nrow(inc_df)]), NA)

    # transposition cutoff
    tcut <- if (ptype == "height"){10} else {30}
    # allowance for change in height/weight
    allowance <- if (ptype == "height"){2.54} else {2}

    # identify transposition exclusions
    exc_transpo <-
      ((abs(dewma$dewma.all) > tcut &
          abs(dewma$dewma.before) > abs(.9*dewma$dewma.all) &
          abs(dewma$dewma.after) > abs(.9*dewma$dewma.all)) |
         (abs(dewma$dewma.all) < -1*tcut &
            abs(dewma$dewma.before) < abs(-.9*dewma$dewma.all) &
            abs(dewma$dewma.after) < abs(-.9*dewma$dewma.all))) &
      (absdewma_transpo$dewma.all <= allowance |
         inc_df$abs_transpo_prev <= allowance |
         inc_df$abs_transpo_next <= allowance) &
      abs(inc_df$tens - inc_df$ones) >= 3

    criteria <- criteria | exc_transpo
    criteria[is.na(criteria)] <- F
  }

  return(criteria)
}

# step 10 hab, H distinct values ----

#' function to calculate height growth allowance
#' @keywords internal
#' @noRd
ht_allow <- function(velocity, ageyears1, ageyears2){
  return(
    velocity*(log(ageyears2 - 16.9)) - (velocity*log(ageyears1 - 16.9))
  )
}

#' function to generate height growth/loss groups
#' @keywords internal
#' @noRd
ht_change_groups <- function(h_subj_df, cutoff){
  # already ordered by age
  glist <- galist <- list()
  # keep track of some current group variables
  cg <-  1 # current group
  glist[[cg]] <- setNames(h_subj_df$meas_m[1], h_subj_df$id[1])
  galist[[cg]] <- h_subj_df$age_years[1]
  for (m in 2:nrow(h_subj_df)){
    cm <- h_subj_df$meas_m[m] # current measurement

    # find range of group with current measurement added
    crng <- max(c(glist[[cg]], cm)) - min(c(glist[[cg]], cm))

    # if the range is below 2 inches with the added value, add and move on
    # we're also going to set the names on the measurements as ids for ease later
    if (crng < (5.08 + .001)){
      glist[[cg]] <- setNames(c(glist[[cg]], cm),
                              c(names(glist[[cg]]), h_subj_df$id[m]))
      galist[[cg]] <- c(galist[[cg]], h_subj_df$age_years[m])
    } else {
      # otherwise, we add a new group for the current measurement
      cg <- cg + 1
      glist[[cg]] <- setNames(cm, h_subj_df$id[m])
      galist[[cg]] <- h_subj_df$age_years[m]
    }

    # if there are too many groups, we're going to stop -- we're not doing
    # this test
    if (cg > cutoff){
      break
    }
  }

  return(list(
    "meas" = glist,
    "age" = galist
  ))
}

#' function to compare growth for 3D height groups
#' compare: before or first
#' returns whether or not to use original exclusions
#' @keywords internal
#' @noRd
ht_3d_growth_compare <- function(mean_ht, min_age, glist,
                                 compare = "before"){
  # preallocating on whether or not we want to go by original exclusion
  origexc <- F
  for (i in 2:6){
    # if there are no members of the group, we want to skip
    if (i > length(glist)){
      next
    }

    # for ease, creating reference variables
    check_num <- if (compare == "before"){ i - 1} else {1}
    ageyears1 <- min_age[check_num]
    ageyears2 <- min_age[i]
    mh1 <- mean_ht[check_num]
    mh2 <- mean_ht[i]

    # check based on growth
    htcompare <- ifelse(ageyears2 > 25, 25, ageyears2)

    # using short circuiting
    hta <-
      if ((htcompare - ageyears1) < 1){
        ht_allow(20, ageyears1, htcompare)
      } else if ((htcompare - ageyears1) <= 3){
        ht_allow(15, ageyears1, htcompare)
      } else if ((htcompare - ageyears1) > 3){
        ht_allow(12, ageyears1, htcompare)
      }

    origexc <- origexc |
      (((mh2 - mh1) < 0 |
      (mh2 - mh1) > hta) &
        ageyears1 < 25)
  }

  return(origexc)
}

# Step 9w, W extreme EWMA ----

#' Function to remove data based on exponentially-weighted moving average
#' (Daymont, et al.) for WEIGHT. Cutoff defaults adjusted for adults.
#' inputs:
#' subj_df: subject data frame, which has age in days and z-score
#' ewma_cutoff: EWMA past which considered invalid (center value). left and right
#'   are .5 less.
#' outputs:
#'  logical indicating whether to exclude a record
#' @keywords internal
#' @noRd
remove_ewma_wt <- function(subj_df, ewma_cutoff_low = 60,
                           ewma_cutoff_high = 100){
  orig_subj_df <- subj_df

  # all three need to be beyond a cutoff for exclusion
  # exclude the most extreme, then recalculate again and again
  rem_ids <- c()
  change <- T
  iter <- 1
  while (change){
    # figure out time difference between points
    agedays_bef <- c(Inf, diff(subj_df$age_days))
    agedays_aft <- c(diff(subj_df$age_days), Inf)
    # both year: different cutoffs if they're at least a year apart
    both_year <- agedays_bef > 365.25 & agedays_aft > 365.25

    # calculate ewma
    ewma_res <- ewma_dn(subj_df$age_days, subj_df$meas_m)

    dewma <- subj_df$meas_m - ewma_res
    colnames(dewma) <- paste0("d",colnames(ewma_res))

    criteria_low <-
      !both_year &
      ((dewma$dewma.all > ewma_cutoff_low &
          dewma$dewma.before > .9*dewma$dewma.all &
          dewma$dewma.after > .9*dewma$dewma.all) |
         (dewma$dewma.all < -1*ewma_cutoff_low &
            dewma$dewma.before < -.9*dewma$dewma.all &
            dewma$dewma.after < -.9*dewma$dewma.all))
    criteria_high <-
      both_year &
      ((dewma$dewma.all > ewma_cutoff_high &
          dewma$dewma.before > .9*dewma$dewma.all &
          dewma$dewma.after > .9*dewma$dewma.all) |
         (dewma$dewma.all < -1*ewma_cutoff_high &
            dewma$dewma.before < -.9*dewma$dewma.all &
            dewma$dewma.after < -.9*dewma$dewma.all))
    criteria_new <- criteria_low | criteria_high
    criteria_new[is.na(criteria_new)] <- F

    if (all(!criteria_new)){
      # if none of them are to be removed
      change <- F
      # if using intermediate values, we want to keep some
    } else {
      # figure out the most extreme value and remove it and rerun
      to_rem <- which.max(abs(dewma$dewma.all)[criteria_new])

      # keep the ids that failed and remove
      rem_ids[length(rem_ids)+1] <- unlist(subj_df[criteria_new,][to_rem, "id"])
      subj_df <- subj_df[subj_df$id != rem_ids[length(rem_ids)],]
      # update iteration
      iter <- iter + 1

      # check if this is viable -- you need at least three points, otherwise
      # we're done
      if (nrow(subj_df) < 3){
        change <- F
      }
    }
  }

  # form results into a logical vector
  criteria <- rep(F, nrow(orig_subj_df))
  criteria[orig_subj_df$id %in% rem_ids] <- T

  return(criteria)
}






