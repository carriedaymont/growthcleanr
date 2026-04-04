# Adult growthcleanr support functions
# Internal functions for the adult algorithm (adult_clean.R)
# Supports both independent and linked repval_handling modes

# permissiveness presets ----

#' Return preset parameter values for each permissiveness level
#' @keywords internal
permissiveness_presets <- function() {
  list(
    loosest = list(
      # BIV (overall) limits
      overall_ht_min = 50, overall_ht_max = 244,
      overall_wt_min = 20, overall_wt_max = 500,
      overall_bmi_min = 5, overall_bmi_max = 300,
      # 1D (single) limits
      single_ht_min_bmi = 60, single_ht_max_bmi = 245,
      single_wt_min_bmi = 12, single_wt_max_bmi = 350,
      single_ht_min_nobmi = 122, single_ht_max_nobmi = 245,
      single_wt_min_nobmi = 30, single_wt_max_nobmi = 350,
      single_bmi_min = 10, single_bmi_max = 250,
      # Algorithm parameters
      wtallow_formula = "piecewise",
      ewma_cap_short = 50, ewma_cap_long = 80,
      wt_scale_et = 0.70, wt_scale_wtallow = 0.50,
      wt_scale_threshold = 120,
      perclimit_low = 0.5, perclimit_mid = 0.4, perclimit_high = 0.0,
      error_load_threshold = 0.41,
      mod_ewma_f = 0.75,
      ht_band = 3,
      allow_ht_loss = TRUE, allow_ht_gain = TRUE,
      repval_handling = "independent"
    ),
    looser = list(
      overall_ht_min = 120, overall_ht_max = 230,
      overall_wt_min = 30, overall_wt_max = 270,
      overall_bmi_min = 12, overall_bmi_max = 65,
      single_ht_min_bmi = 120, single_ht_max_bmi = 230,
      single_wt_min_bmi = 30, single_wt_max_bmi = 270,
      single_ht_min_nobmi = 120, single_ht_max_nobmi = 230,
      single_wt_min_nobmi = 30, single_wt_max_nobmi = 270,
      single_bmi_min = 12, single_bmi_max = 65,
      wtallow_formula = "piecewise",
      ewma_cap_short = 50, ewma_cap_long = 80,
      wt_scale_et = 0.70, wt_scale_wtallow = 0.50,
      wt_scale_threshold = 120,
      perclimit_low = 0.5, perclimit_mid = 0.4, perclimit_high = 0.0,
      error_load_threshold = 0.41,
      mod_ewma_f = 0.75,
      ht_band = 3,
      allow_ht_loss = FALSE, allow_ht_gain = TRUE,
      repval_handling = "independent"
    ),
    tighter = list(
      overall_ht_min = 142, overall_ht_max = 213,
      overall_wt_min = 36, overall_wt_max = 159,
      overall_bmi_min = 16, overall_bmi_max = 45,
      single_ht_min_bmi = 142, single_ht_max_bmi = 213,
      single_wt_min_bmi = 36, single_wt_max_bmi = 159,
      single_ht_min_nobmi = 142, single_ht_max_nobmi = 213,
      single_wt_min_nobmi = 36, single_wt_max_nobmi = 159,
      single_bmi_min = 16, single_bmi_max = 45,
      wtallow_formula = "piecewise-lower",
      ewma_cap_short = 40, ewma_cap_long = 60,
      wt_scale_et = 0, wt_scale_wtallow = 0,
      wt_scale_threshold = 120,
      perclimit_low = 0.7, perclimit_mid = 0.4, perclimit_high = 0.4,
      error_load_threshold = 0.29,
      mod_ewma_f = 0.60,
      ht_band = 2,
      allow_ht_loss = FALSE, allow_ht_gain = TRUE,
      repval_handling = "linked"
    ),
    tightest = list(
      overall_ht_min = 147, overall_ht_max = 208,
      overall_wt_min = 39, overall_wt_max = 136,
      overall_bmi_min = 18, overall_bmi_max = 40,
      single_ht_min_bmi = 147, single_ht_max_bmi = 208,
      single_wt_min_bmi = 39, single_wt_max_bmi = 136,
      single_ht_min_nobmi = 147, single_ht_max_nobmi = 208,
      single_wt_min_nobmi = 39, single_wt_max_nobmi = 136,
      single_bmi_min = 18, single_bmi_max = 40,
      wtallow_formula = "allofus15",
      ewma_cap_short = 40, ewma_cap_long = 40,
      wt_scale_et = 0, wt_scale_wtallow = 0,
      wt_scale_threshold = 120,
      perclimit_low = 0.7, perclimit_mid = 0.4, perclimit_high = 0.4,
      error_load_threshold = 0.29,
      mod_ewma_f = 0.60,
      ht_band = 2,
      allow_ht_loss = FALSE, allow_ht_gain = FALSE,
      repval_handling = "linked"
    )
  )
}

#' Resolve permissiveness: fill NULLs from preset, keep explicit values
#' @param permissiveness Character: "loosest", "looser", "tighter", "tightest"
#' @param ... Named parameter values (NULL = use preset, non-NULL = override)
#' @return Named list of all resolved parameter values
#' @keywords internal
resolve_permissiveness <- function(permissiveness, ...) {
  valid <- c("loosest", "looser", "tighter", "tightest")
  if (!permissiveness %in% valid) {
    stop(paste0("permissiveness must be one of: ",
                paste(valid, collapse = ", "),
                ". Got: '", permissiveness, "'"))
  }
  preset <- permissiveness_presets()[[permissiveness]]
  user <- list(...)
  # For each preset param, use user value if non-NULL, else preset
  resolved <- preset
  for (nm in names(user)) {
    if (!is.null(user[[nm]])) {
      resolved[[nm]] <- user[[nm]]
    }
  }
  resolved
}

# convenience functions ----

#' convenience function -- see if numeric vector falls between two numbers
#' returns boolean vector
#' @keywords internal
#' @noRd
check_between <- function(vect, num_low, num_high, incl = TRUE) {
  return(
    if (incl) {
      vect <= num_high & vect >= num_low
    } else {
      vect < num_high & vect > num_low
    }
  )
}

#' convenience function -- round to the nearest .x
#' @keywords internal
#' @noRd
round_pt <- function(val, pt) {
  return(round(val / pt) * pt)
}

# Dynamic ET/EWMA cap and perclimit helpers ----

#' Compute the dynamic ET/EWMA cap for given intervals and weights.
#' Vectorized: interval_months and maxwt can be vectors.
#' @param interval_months Interval in months (minimum neighbor gap)
#' @param maxwt Maximum of (weight, EWMA) for each observation/pair
#' @param cap_short Baseline cap for intervals < 6 months
#' @param cap_long Baseline cap for intervals >= 6 months
#' @param wt_scale_et Scaling factor for weight (0 = no scaling)
#' @param wt_scale_threshold Weight threshold above which scaling applies
#' @return Numeric vector of ET limits
#' @keywords internal
compute_et_limit <- function(interval_months, maxwt,
                             cap_short, cap_long,
                             wt_scale_et = 0, wt_scale_threshold = 120) {
  baseline <- ifelse(interval_months < 6, cap_short, cap_long)
  wt_addition <- wt_scale_et * pmax(0, maxwt - wt_scale_threshold)
  baseline + wt_addition
}

#' Compute observation-level percentage criterion limit.
#' Vectorized over meas.
#' @param meas Weight measurements in kg
#' @param perclimit_low Limit for wt <= 45 kg
#' @param perclimit_mid Limit for 45 < wt <= 80 kg
#' @param perclimit_high Limit for wt > 80 kg (0 = disabled: percewma < 0 is never TRUE)
#' @return Numeric vector of perc_limits
#' @keywords internal
compute_perc_limit <- function(meas, perclimit_low, perclimit_mid, perclimit_high) {
  ifelse(meas <= 45, perclimit_low, ifelse(meas <= 80, perclimit_mid, perclimit_high))
}

# EWMA functions ----

#' function to calculate as delta matrix for adults
#' @keywords internal
#' @noRd
as.matrix.delta_dn <- function(agedays) {
  n <- length(agedays)
  delta <- abs(matrix(rep(agedays, n), n, byrow = TRUE) - agedays)
  return(delta)
}

#' Exponentially Weighted Moving Average (EWMA) (daymont implementation)
#' Adult version uses |delta|^(-5) weighting (not pediatric (5+delta)^(-1.5))
#' @keywords internal
#' @noRd
ewma_dn <- function(agedays, meas, ewma.exp = -5, ewma.adjacent = TRUE,
                    ewma_window = 15) {
  n <- length(agedays)
  ewma.all <- ewma.before <- ewma.after <- vector('numeric', 0)
  if (n > 0) {
    if (!all(agedays == cummax(agedays)))
      warning("EWMA ordering is not sorted; double check")
    index <- order(agedays)

    delta <- as.matrix.delta_dn(agedays)
    delta <- ifelse(delta == 0, 0, (delta) ^ ewma.exp)

    # Apply position-based window: zero out entries beyond ewma_window positions
    if (!is.null(ewma_window)) {
      pos_dist <- abs(row(delta) - col(delta))
      delta[pos_dist > ewma_window] <- 0
    }

    ewma.all[index] <- delta %*% meas / apply(delta, 1, sum)

    if (ewma.adjacent) {
      if (n > 2) {
        delta2 = delta
        delta2[col(delta2) == row(delta2) - 1] = 0
        ewma.before[index] = delta2 %*% meas / apply(delta2, 1, sum)
        delta3 = delta
        delta3[col(delta3) == row(delta3) + 1] = 0
        ewma.after[index] = delta3 %*% meas / apply(delta3, 1, sum)
      } else {
        ewma.before <- ewma.after <- ewma.all
      }
    }
  }
  return(if (ewma.adjacent)
    data.frame(ewma.all, ewma.before, ewma.after)
    else
      data.frame(ewma.all))
}

# EWMA Cache functions ----
# EWMA cache: O(n) iterative updates instead of O(n²) full rebuild

#' Initialize EWMA cache for a set of observations.
#' Builds the full weight matrix once (O(n²)), computes EWMA values.
#' @param agedays Sorted age in days
#' @param meas Measurements (same order as agedays)
#' @param ewma_exp Exponent for distance weighting (default -5)
#' @return Cache list with delta matrix, weighted sums, row sums, and EWMA values
#' @keywords internal
adult_ewma_cache_init <- function(agedays, meas, ewma_exp = -5, ewma_window = 15) {
  n <- length(agedays)
  if (n == 0) return(NULL)

  # Build weight matrix
  delta <- as.matrix.delta_dn(agedays)
  delta <- ifelse(delta == 0, 0, delta ^ ewma_exp)

  # Apply position-based window: zero out entries beyond ewma_window positions
  if (!is.null(ewma_window)) {
    pos_dist <- abs(row(delta) - col(delta))
    delta[pos_dist > ewma_window] <- 0
  }

  # Weighted sums and row sums
  ws <- as.vector(delta %*% meas)
  rs <- rowSums(delta)

  # EWMA values
  ewma_all <- ws / rs

  if (n > 2) {
    # Subdiagonal trick: ewma_before removes predecessor, ewma_after removes successor
    pred_w <- c(0, delta[cbind(2:n, 1:(n - 1))])
    pred_m <- c(0, meas[1:(n - 1)])
    ewma_before <- (ws - pred_w * pred_m) / (rs - pred_w)

    succ_w <- c(delta[cbind(1:(n - 1), 2:n)], 0)
    succ_m <- c(meas[2:n], 0)
    ewma_after <- (ws - succ_w * succ_m) / (rs - succ_w)
  } else {
    ewma_before <- ewma_after <- ewma_all
  }

  list(
    delta = delta,
    ws = ws,
    rs = rs,
    meas = meas,
    agedays = agedays,
    n = n,
    ewma_all = ewma_all,
    ewma_before = ewma_before,
    ewma_after = ewma_after
  )
}

#' Update EWMA cache by removing one observation. O(n) instead of O(n²).
#' @param cache Cache list from adult_ewma_cache_init or previous update
#' @param pos_j Position (index) of the observation to remove
#' @return Updated cache list with n-1 observations, or NULL if n would be 0
#' @keywords internal
adult_ewma_cache_update <- function(cache, pos_j) {
  n <- cache$n
  if (n <= 1) return(NULL)

  keep <- seq_len(n)[-pos_j]

  # Subtract obs j's contribution from all other weighted sums and row sums
  col_j <- cache$delta[keep, pos_j]
  ws <- cache$ws[keep] - col_j * cache$meas[pos_j]
  rs <- cache$rs[keep] - col_j

  # Trim matrix and vectors
  delta <- cache$delta[keep, keep, drop = FALSE]
  meas <- cache$meas[keep]
  agedays <- cache$agedays[keep]
  n_new <- n - 1L

  # Recompute EWMA values from updated sums
  ewma_all <- ws / rs

  if (n_new > 2) {
    pred_w <- c(0, delta[cbind(2:n_new, 1:(n_new - 1))])
    pred_m <- c(0, meas[1:(n_new - 1)])
    ewma_before <- (ws - pred_w * pred_m) / (rs - pred_w)

    succ_w <- c(delta[cbind(1:(n_new - 1), 2:n_new)], 0)
    succ_m <- c(meas[2:n_new], 0)
    ewma_after <- (ws - succ_w * succ_m) / (rs - succ_w)
  } else {
    ewma_before <- ewma_after <- ewma_all
  }

  list(
    delta = delta,
    ws = ws,
    rs = rs,
    meas = meas,
    agedays = agedays,
    n = n_new,
    ewma_all = ewma_all,
    ewma_before = ewma_before,
    ewma_after = ewma_after
  )
}

# step 1w, W BIV ----

#' function to remove BIVs, based on cutoffs for the given method
#' @keywords internal
#' @noRd
remove_biv <- function(subj_df, type, biv_df, include = FALSE) {
  too_low <- remove_biv_low(subj_df, type, biv_df, include)
  too_high <- remove_biv_high(subj_df, type, biv_df, include)
  return(too_low | too_high)
}

#' @keywords internal
#' @noRd
remove_biv_low <- function(subj_df, type, biv_df, include = FALSE) {
  # 0.12 tolerance for rounding (0.1 cm/kg rounding + float precision)
  if (!include) {
    too_low <- subj_df$meas_m < biv_df[type, "low"] - 0.12
  } else {
    too_low <- subj_df$meas_m <= biv_df[type, "low"] - 0.12
  }
  return(too_low)
}

#' @keywords internal
#' @noRd
remove_biv_high <- function(subj_df, type, biv_df, include = FALSE) {
  # 0.12 tolerance for rounding (0.1 cm/kg rounding + float precision)
  if (!include) {
    too_high <- subj_df$meas_m > biv_df[type, "high"] + 0.12
  } else {
    too_high <- subj_df$meas_m >= biv_df[type, "high"] + 0.12
  }
  return(too_high)
}

# step 2w, W repeated values ----

#' Identify repeated weight values within a subject.
#' Marks the first occurrence (earliest age, id tiebreaker) as is_first_rv=TRUE
#' and all subsequent identical values as is_rv=TRUE. Unique values get both FALSE.
#' Relies on w_subj_df being pre-sorted by age/id. Exact numeric match only
#' (78.1 != 78.101). Does not filter by extraneous status; caller is responsible
#' for sequencing (see identify_rv -> temp_sde -> redo_identify_rv cycle).
#' @keywords internal
#' @noRd
identify_rv <- function(w_subj_df) {
  if (nrow(w_subj_df) > 0) {
    # is_rv: TRUE for all duplicates except the first occurrence
    w_subj_df$is_rv <- duplicated(w_subj_df$meas_m)
    # is_first_rv: TRUE for the first occurrence of values that have duplicates
    w_subj_df$is_first_rv <- duplicated(w_subj_df$meas_m, fromLast = TRUE) &
      !w_subj_df$is_rv
  }
  return(w_subj_df)
}

# step 3, temp extraneous ----

#' Temporarily flag same-day values that deviate most from patient median.
#' For weight (ptype="weight"): median of non-RV values (all included, not
#' filtered by same-day status). For height (ptype="height"): median of all
#' included values. On each same-day group, the value closest to median
#' survives (id tiebreaker); others are marked extraneous=TRUE.
#' @keywords internal
#' @noRd
temp_sde <- function(subj_df, ptype = "height") {
  tab_days <- table(subj_df$age_days)
  dup_days <- names(tab_days)[tab_days > 1]

  if (nrow(subj_df) >= 2) {
    # Median of all included values (non-RV only for weight, all for height)
    if (ptype == "weight") {
      med_val <- median(subj_df$meas_m[!subj_df$is_rv])
    } else {
      med_val <- median(subj_df$meas_m)
    }

    subj_df$diff <- NA
    subj_df$diff[as.character(subj_df$age_days) %in% dup_days] <-
      abs(subj_df$meas_m[as.character(subj_df$age_days) %in% dup_days] -
            med_val)

    subj_df$extraneous <- FALSE
    for (dd in dup_days) {
      day_diffs <- subj_df$diff[as.character(subj_df$age_days) == dd]
      # Keep highest id among ties (last position with min diff, matching
      # final SDE resolution which keeps highest id)
      keeper <- max(which(day_diffs == min(day_diffs)))
      subj_df$extraneous[as.character(subj_df$age_days) == dd][-keeper] <- TRUE
    }

    subj_df$diff <- NULL
  } else if (nrow(subj_df) > 0) {
    subj_df$extraneous <- FALSE
  }

  return(subj_df)
}

#' Re-identify repeated values after temp SDE resolution.
#' Only runs if an is_first_rv value was marked extraneous by temp_sde() —
#' in that case, a different value in the RV group may need to become first_rv.
#' Subsets to non-extraneous values, re-runs identify_rv(), maps results back.
#' Note: after non-SDE exclusions (weight cap, evil twins), the calling code
#' uses identify_rv() directly on the remaining rows, then temp_sde(), then
#' this function for SDE cleanup.
#' @keywords internal
#' @noRd
redo_identify_rv <- function(w_subj_df) {
  if (nrow(w_subj_df) > 0 & any(w_subj_df$extraneous & w_subj_df$is_first_rv)) {
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous, ])
    inc_df <- identify_rv(inc_df)
    w_subj_df$is_first_rv <- w_subj_df$is_rv <- FALSE
    w_subj_df$is_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_rv
    w_subj_df$is_first_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_first_rv
  }
  return(w_subj_df)
}

# step 10 hab, H distinct values ----

#' function to calculate height growth allowance
#' @keywords internal
#' @noRd
ht_allow <- function(velocity, ageyears1, ageyears2) {
  return(
    velocity * (log(ageyears2 - 16.9)) - (velocity * log(ageyears1 - 16.9))
  )
}

#' function to generate height growth/loss groups
#' @keywords internal
#' @noRd
ht_change_groups <- function(h_subj_df, cutoff, type = "loss") {
  glist <- galist <- list()
  cg <- 1
  glist[[cg]] <- setNames(h_subj_df$meas_m[1], h_subj_df$id[1])
  galist[[cg]] <- h_subj_df$ageyears[1]
  for (m in 2:nrow(h_subj_df)) {
    cm <- h_subj_df$meas_m[m]
    crng <- max(c(glist[[cg]], cm)) - min(c(glist[[cg]], cm))
    temp_mindiff <- min(glist[[cg]]) - cm
    temp_maxdiff <- max(glist[[cg]]) - cm

    if (crng < (5.08 + 0.12)) {
      glist[[cg]] <- setNames(c(glist[[cg]], cm),
                              c(names(glist[[cg]]), h_subj_df$id[m]))
      galist[[cg]] <- c(galist[[cg]], h_subj_df$ageyears[m])
    } else {
      cg <- cg + 1
      glist[[cg]] <- setNames(cm, h_subj_df$id[m])
      galist[[cg]] <- h_subj_df$ageyears[m]
    }

    if (cg > cutoff) {
      break
    }

    if (type == "loss") {
      if (temp_mindiff < -(5.08 + 0.12)) {
        glist <- galist <- list()
        break
      }
    } else {
      if (temp_maxdiff > (5.08 + 0.12)) {
        glist <- galist <- list()
        break
      }
    }
  }

  return(list(
    "meas" = glist,
    "age" = galist
  ))
}

#' function to compare growth for 3D height groups
#' @keywords internal
#' @noRd
ht_3d_growth_compare <- function(mean_ht, min_age, glist,
                                 compare = "before") {
  origexc <- FALSE
  for (i in 2:6) {
    if (i > length(glist)) {
      next
    }
    check_num <- if (compare == "before") { i - 1 } else { 1 }
    ageyears1 <- min_age[check_num]
    ageyears2 <- min_age[i]
    mh1 <- mean_ht[check_num]
    mh2 <- mean_ht[i]

    htcompare <- ifelse(ageyears2 > 25, 25, ageyears2)

    hta <-
      if ((htcompare - ageyears1) < 1) {
        ht_allow(20, ageyears1, htcompare)
      } else if ((htcompare - ageyears1) <= 3) {
        ht_allow(15, ageyears1, htcompare)
      } else if ((htcompare - ageyears1) > 3) {
        ht_allow(12, ageyears1, htcompare)
      }

    origexc <- origexc |
      ((mh2 - mh1) < -0.12 |
       (mh2 - mh1) > hta + 0.12)
  }

  return(origexc)
}

# compute_wtallow ----

#' Compute weight allowance from interval in months (vectorized)
#'
#' Built-in formulas: "piecewise" (default), "piecewise-lower", "allofus15".
#' Custom: pass a file path to a CSV with columns 'months' and 'wtallow'.
#' Values are linearly interpolated; months beyond the table are clamped to the
#' nearest row.
#'
#' @param months Numeric vector of interval lengths in months
#' @param formula Character: formula name or path to custom CSV
#' @param cap_params List with cap_short, cap_long, wt_scale_et,
#'   wt_scale_wtallow, wt_scale_threshold
#' @param maxwt Optional numeric vector of max(wt, ewma) for weight scaling
#' @return Numeric vector of weight allowances
#' @keywords internal
compute_wtallow <- function(months, formula = "piecewise",
                            cap_params = list(cap_short = 50, cap_long = 80,
                                              wt_scale_et = 0.70,
                                              wt_scale_wtallow = 0.50,
                                              wt_scale_threshold = 120),
                            maxwt = NULL) {
  cap_short <- cap_params$cap_short
  cap_long  <- cap_params$cap_long

  if (formula == "piecewise") {
    # Breakpoints: 20 @ 1mo, cap_short @ 6mo, cap_long @ 12mo, flat after
    slope_1_6  <- (cap_short - 20) / 5
    slope_6_12 <- (cap_long - cap_short) / 6
    wta <- ifelse(months <= 1, 10 + 10 * log(1 + 5 * months) / log(6),
           ifelse(months <= 6, 20 + slope_1_6 * (months - 1),
           ifelse(months <= 12, cap_short + slope_6_12 * (months - 6),
                                cap_long)))
  } else if (formula == "piecewise-lower") {
    # Tighter variant: reaches 2/3*cap_short at 6mo, caps at cap_short for >12mo
    target_6mo  <- cap_short * 2 / 3
    target_12mo <- cap_short * 2 / 3 + (cap_long - cap_short) * 2 / 3
    slope_1_6   <- (target_6mo - 20) / 5
    slope_6_12  <- (target_12mo - target_6mo) / 6
    wta <- ifelse(months <= 1, 10 + 10 * log(1 + 5 * months) / log(6),
           ifelse(months <= 6, 20 + slope_1_6 * (months - 1),
           ifelse(months <= 12, target_6mo + slope_6_12 * (months - 6),
                                cap_short)))
  } else if (formula == "allofus15") {
    # Short-interval breakpoints: 5kg (1-2 days), 10kg (3-7 days),
    # 15kg (8-180 days), then linear 15→cap_short from 180d–12mo, flat after
    days_2  <- 2 / 30.4375    # ~0.0657 months
    days_7  <- 7 / 30.4375    # ~0.230 months
    days_180 <- 180 / 30.4375 # ~5.913 months
    slope_180d_12mo <- (cap_short - 15) / (12 - days_180)
    wta <- ifelse(months <= days_2, 5,
           ifelse(months <= days_7, 10,
           ifelse(months <= days_180, 15,
           ifelse(months <= 12,
                  15 + slope_180d_12mo * (months - days_180),
                  cap_short))))
  } else {
    # Custom CSV: columns 'months' and 'wtallow', linearly interpolated
    if (!exists(".wtallow_custom_cache", envir = .GlobalEnv) ||
        !identical(attr(get(".wtallow_custom_cache", envir = .GlobalEnv), "path"), formula)) {
      if (!file.exists(formula)) {
        stop(paste0("wtallow_formula '", formula, "' is not a built-in formula ",
                     "('piecewise', 'piecewise-lower', 'allofus15') ",
                     "and file not found."))
      }
      custom <- read.csv(formula, stringsAsFactors = FALSE)
      if (!all(c("months", "wtallow") %in% names(custom))) {
        stop("Custom wtallow CSV must have columns 'months' and 'wtallow'.")
      }
      custom <- custom[order(custom$months), ]
      attr(custom, "path") <- formula
      assign(".wtallow_custom_cache", custom, envir = .GlobalEnv)
    }
    custom <- get(".wtallow_custom_cache", envir = .GlobalEnv)
    wta <- approx(x = custom$months, y = custom$wtallow, xout = months, rule = 2)$y
  }

  # Weight-scaled addition (loosest/looser: wt_scale_wtallow > 0)
  if (!is.null(maxwt) && cap_params$wt_scale_wtallow > 0) {
    wt_addition <- cap_params$wt_scale_wtallow *
                   pmax(0, maxwt - cap_params$wt_scale_threshold)
    wta <- wta + wt_addition
  }

  # Ceiling: wtallow cannot exceed the ET limit for this interval/weight
  if (!is.null(maxwt)) {
    et_ceiling <- compute_et_limit(months, maxwt,
                                   cap_short, cap_long,
                                   cap_params$wt_scale_et,
                                   cap_params$wt_scale_threshold)
    wta <- pmin(wta, et_ceiling)
  } else {
    # Without maxwt, use the static cap_long as ceiling (conservative)
    wta <- pmin(wta, cap_long)
  }

  wta
}

# detect_runs ----

#' Detect consecutive runs of TRUE in a logical vector
#' Returns list with run_id, run_len, run_pos for each element
#' @keywords internal
detect_runs <- function(flagged) {
  n <- length(flagged)
  run_id  <- rep(NA_integer_, n)
  run_len <- rep(NA_integer_, n)
  run_pos <- rep(NA_integer_, n)
  if (n == 0 || !any(flagged)) return(list(run_id = run_id, run_len = run_len, run_pos = run_pos))

  is_start <- flagged & c(TRUE, !flagged[-n])
  cum_id <- cumsum(is_start)
  cum_id[!flagged] <- NA

  if (any(!is.na(cum_id))) {
    tbl <- table(cum_id)
    for (rid in names(tbl)) {
      idx <- which(cum_id == as.integer(rid))
      run_id[idx]  <- as.integer(rid)
      run_len[idx] <- as.integer(tbl[rid])
      run_pos[idx] <- seq_along(idx)
    }
  }
  list(run_id = run_id, run_len = run_len, run_pos = run_pos)
}

# compute_trajectory_fails ----

#' Pre-computes trajectory rescue for all observations.
#' Returns TRUE = fails all rescue (not rescued).
#' @keywords internal
compute_trajectory_fails <- function(meas, age_days, err = 5) {
  n <- length(meas)
  if (n < 3) return(rep(TRUE, n))

  p1 <- c(NA, meas[-n])
  n1 <- c(meas[-1], NA)
  p2 <- c(NA, NA, meas[1:(n - 2)])
  n2 <- c(meas[3:n], NA, NA)
  ap1 <- c(NA, age_days[-n])
  an1 <- c(age_days[-1], NA)
  ap2 <- c(NA, NA, age_days[1:(n - 2)])
  an2 <- c(age_days[3:n], NA, NA)

  # METHOD 1: Interpolation between p1 and n1 (±err)
  lo <- pmin(p1, n1, na.rm = FALSE) - err
  hi <- pmax(p1, n1, na.rm = FALSE) + err
  rescued_interp <- !is.na(lo) & !is.na(hi) & meas >= lo & meas <= hi

  # METHOD 2: Extrapolation from prior (p2 → p1)
  slope_p <- (p1 - p2) / (ap1 - ap2)
  lepolate_p <- p1 + slope_p * (age_days - ap1)
  lepolate_p <- round(lepolate_p / 0.2) * 0.2
  lo_p <- pmin(p2, lepolate_p, na.rm = FALSE) - err
  hi_p <- pmax(p2, lepolate_p, na.rm = FALSE) + err
  rescued_prior <- !is.na(lo_p) & !is.na(hi_p) & meas >= lo_p & meas <= hi_p
  # Distance guard: don't trust extrapolation > 2× source interval
  dist_extrap_p <- abs(age_days - ap1)
  dist_source_p <- abs(ap1 - ap2)
  rescued_prior[!is.na(dist_extrap_p) & !is.na(dist_source_p) &
                dist_extrap_p > 2 * dist_source_p] <- FALSE

  # METHOD 3: Extrapolation from next (n2 → n1)
  slope_n <- (n1 - n2) / (an1 - an2)
  lepolate_n <- n1 + slope_n * (age_days - an1)
  lepolate_n <- round(lepolate_n / 0.2) * 0.2
  lo_n <- pmin(n2, lepolate_n, na.rm = FALSE) - err
  hi_n <- pmax(n2, lepolate_n, na.rm = FALSE) + err
  rescued_next <- !is.na(lo_n) & !is.na(hi_n) & meas >= lo_n & meas <= hi_n
  dist_extrap_n <- abs(an1 - age_days)
  dist_source_n <- abs(an2 - an1)
  rescued_next[!is.na(dist_extrap_n) & !is.na(dist_source_n) &
               dist_extrap_n > 2 * dist_source_n] <- FALSE

  rescued_interp[is.na(rescued_interp)] <- FALSE
  rescued_prior[is.na(rescued_prior)]   <- FALSE
  rescued_next[is.na(rescued_next)]     <- FALSE

  fails <- !rescued_interp & !rescued_prior & !rescued_next
  fails
}

# Step 9Wa: Evil Twins ----

#' Evil twins detection for one subject (independent mode)
#' @param w_subj_df Weight data for one subject (data.table with meas_m, age_days, id, is_rv)
#'   Must include all Inc values (including RVs in independent mode)
#' @return Character vector of ids to exclude with "Evil twins" code
#' @keywords internal
evil_twins <- function(w_subj_df, cap_params) {
  # Need at least 3 included values for pairs guard
  inc_df <- w_subj_df[order(w_subj_df$age_days, w_subj_df$id), ]
  exc_ids <- character(0)

  while (TRUE) {
    # Work with currently non-excluded values
    working <- inc_df[!inc_df$id %in% exc_ids, ]
    n <- nrow(working)
    if (n < 3) break

    working <- working[order(working$age_days, working$id), ]

    # Compute interval-based caps (etcap) between adjacent values
    age_diff_days <- diff(working$age_days)
    months_next <- age_diff_days / 30.4375
    wt_diff <- abs(diff(working$meas_m))

    # maxwt for each pair: max of the two weights (no EWMA in evil twins)
    maxwt_pair <- pmax(working$meas_m[-n], working$meas_m[-1])

    # Dynamic cap by interval and weight, with 0.12 kg rounding tolerance
    etcap <- compute_et_limit(months_next, maxwt_pair,
                              cap_params$cap_short, cap_params$cap_long,
                              cap_params$wt_scale_et,
                              cap_params$wt_scale_threshold) + 0.12

    # Find out-of-bounds pairs
    oob_pairs <- which(wt_diff > etcap)
    if (length(oob_pairs) == 0) break

    # Mark which observations are in OOB pairs
    oob_obs <- unique(c(oob_pairs, oob_pairs + 1))
    working$oob <- FALSE
    working$oob[oob_obs] <- TRUE

    if (!any(working$oob)) break

    # Independent mode: all values (Inc + RV) participate in median and pairs guard
    all_vals <- working$meas_m
    if (length(all_vals) == 0) break
    subj_median <- median(all_vals)

    # Absolute deviation from median
    working$absd_med <- abs(working$meas_m - subj_median)

    # Plausibility guardrail: <38 or >180 → maximum deviation
    working$absd_med[working$oob & (working$meas_m < 38 | working$meas_m > 180)] <- 99999

    # Among OOB observations, find the most deviant
    oob_working <- working[working$oob, ]
    # Sort: most deviant first, then highest weight, then lowest id
    oob_working <- oob_working[order(-oob_working$absd_med, -oob_working$meas_m, oob_working$id), ]

    # Exclude the most deviant one
    exc_ids <- c(exc_ids, as.character(oob_working$id[1]))
  }

  exc_ids
}

# Step 9Wb/11Wb: EWMA weight outlier removal ----

#' Remove EWMA weight outliers (used for both Extreme and Moderate EWMA)
#'
#' Iteratively excludes values whose EWMA deviation exceeds an interval-specific
#' threshold. Each round removes the single worst outlier (largest |dewma|, with
#' age then id as tiebreakers). Iterates until no candidates remain or <3 values.
#'
#' Threshold: 2-tier (<6m: cap_short, ≥6m: cap_long) with optional weight
#' scaling via compute_et_limit(). Per-observation threshold based on
#' min neighbor gap and max(wt, ewma).
#'
#' 90% rule: directional dewma (before/after) must exceed 90% of threshold.
#' Missing neighbors (edge values) are treated as Inf (confirming exclusion).
#'
#' @param subj_df Data frame with age_days, meas_m, id columns (pre-sorted)
#' @param cap_params List with cap_short, cap_long, wt_scale_et, wt_scale_threshold.
#' @param exc_label Base label for exclusion codes. Round number is appended.
#' @return Named character vector: id -> "<exc_label>-N" for excluded values
#' @keywords internal
remove_ewma_wt <- function(subj_df, cap_params,
                           exc_label = "Exclude-A-WT-Traj-Ext",
                           ewma_window = 15) {
  orig_subj_df <- subj_df
  rem_ids <- character(0)
  round_codes <- character(0)
  change <- TRUE
  round_num <- 1
  cache <- NULL

  while (change) {
    n <- nrow(subj_df)
    if (n < 3) break

    # Minimum neighbor gap in months (NA → Inf for edge values)
    agedays_bef <- c(NA, subj_df$age_days[-n])
    agedays_aft <- c(subj_df$age_days[-1], NA)
    gap_bef <- subj_df$age_days - agedays_bef
    gap_aft <- agedays_aft - subj_df$age_days
    gap_bef_check <- ifelse(is.na(gap_bef), Inf, gap_bef)
    gap_aft_check <- ifelse(is.na(gap_aft), Inf, gap_aft)
    min_gap_months <- pmin(gap_bef_check, gap_aft_check) / 30.4375

    # EWMA — rebuild each round (window shifts when observations are removed)
    cache <- adult_ewma_cache_init(subj_df$age_days, subj_df$meas_m,
                             ewma_window = ewma_window)
    dewma_all <- subj_df$meas_m - cache$ewma_all
    dewma_bef <- subj_df$meas_m - cache$ewma_before
    dewma_aft <- subj_df$meas_m - cache$ewma_after

    # Dynamic per-observation threshold: 2-tier + weight scaling
    maxwt <- pmax(subj_df$meas_m, cache$ewma_all)
    threshold <- compute_et_limit(min_gap_months, maxwt,
                                  cap_params$cap_short, cap_params$cap_long,
                                  cap_params$wt_scale_et,
                                  cap_params$wt_scale_threshold)

    # 90% rule: directional dewma must exceed 90% of threshold
    # NA → Inf for edge values (Stata missing-as-infinity: . > x is TRUE)
    dewma_bef_safe <- ifelse(is.na(dewma_bef), Inf, dewma_bef)
    dewma_aft_safe <- ifelse(is.na(dewma_aft), Inf, dewma_aft)

    # Positive direction
    # 0.12 kg rounding tolerance on all threshold comparisons
    criteria_pos <- !is.na(dewma_all) & dewma_all > threshold + 0.12 &
                    dewma_bef_safe > 0.9 * threshold + 0.12 &
                    dewma_aft_safe > 0.9 * threshold + 0.12
    # Negative direction
    criteria_neg <- !is.na(dewma_all) & dewma_all < -(threshold + 0.12) &
                    dewma_bef_safe < -(0.9 * threshold + 0.12) &
                    dewma_aft_safe < -(0.9 * threshold + 0.12)

    criteria_new <- criteria_pos | criteria_neg

    if (all(!criteria_new)) {
      change <- FALSE
    } else {
      absdewma <- abs(dewma_all)
      cand_idx <- which(criteria_new)
      # id tiebreaker for sort determinism
      ord <- order(-absdewma[cand_idx], subj_df$age_days[cand_idx],
                   as.numeric(subj_df$id[cand_idx]))
      to_rem <- cand_idx[ord[1]]

      rem_ids <- c(rem_ids, as.character(subj_df$id[to_rem]))
      round_codes <- c(round_codes, paste0(exc_label, "-", round_num))
      subj_df <- subj_df[-to_rem, ]
      round_num <- round_num + 1

      if (nrow(subj_df) < 3) {
        change <- FALSE
      }
    }
  }

  # Return named character vector
  result <- setNames(round_codes, rem_ids)
  result
}

# Linked mode: RV propagation ----

#' Propagate firstRV exclusions forward to RV copies and extraneous values
#' (linked mode). For each excluded value, finds all other values with the
#' same measurement that are still "Include" — including extraneous values
#' whose is_rv flag may have been cleared by redo_identify_rv(). Marks them
#' with "<code>-RV-propagated". This is forward-only propagation (firstRV →
#' later RVs + extraneous), unlike Step 4W's bidirectional scale-max propagation.
#' @param exc_codes Named character vector: id -> exclusion code (from firstRV pass)
#' @param w_subj_df Current working weight data (must have meas_m, id columns)
#' @param w_subj_keep Named vector of all exclusion codes for this subject's weights
#' @return List with updated w_subj_keep and propagated_ids
#' @keywords internal
propagate_to_rv <- function(exc_codes, w_subj_df, w_subj_keep) {
  propagated_ids <- character(0)
  if (length(exc_codes) == 0) {
    return(list(w_subj_keep = w_subj_keep, propagated_ids = propagated_ids))
  }

  for (exc_id in names(exc_codes)) {
    idx <- which(as.character(w_subj_df$id) == exc_id)
    if (length(idx) == 0) next
    exc_meas <- w_subj_df$meas_m[idx[1]]

    # Find all values with same meas_m (RV copies + extraneous with same value)
    match_idx <- which(w_subj_df$meas_m == exc_meas &
                       as.character(w_subj_df$id) != exc_id)
    if (length(match_idx) == 0) next

    match_ids <- as.character(w_subj_df$id[match_idx])
    match_ids <- match_ids[w_subj_keep[match_ids] == "Include"]

    if (length(match_ids) > 0) {
      prop_code <- paste0(exc_codes[exc_id], "-RV-propagated")
      w_subj_keep[match_ids] <- prop_code
      propagated_ids <- c(propagated_ids, match_ids)
    }
  }

  list(w_subj_keep = w_subj_keep, propagated_ids = propagated_ids)
}

# Step 11Wb: 7-step Moderate EWMA ----

#' Remove implausible weight values using 7-step moderate EWMA flow.
#' @param full_inc_df Data frame for one subject with columns:
#'   id, age_days, ageyears, meas_m (weight in kg). Must have >= 3 rows.
#' @param exc_label Character prefix for exclusion codes
#' @param max_rounds Maximum number of rounds (default 5)
#' @param wtallow_formula Formula for wtallow ("piecewise", "piecewise-lower", "allofus15")
#' @return Named character vector: id -> exclusion code for excluded values.
#' @keywords internal
remove_mod_ewma_wt <- function(full_inc_df, exc_label = "Exclude-A-WT-Traj-Moderate",
                               max_rounds = 100, wtallow_formula = "piecewise",
                               cap_params = NULL, ewma_window = 15,
                               mod_ewma_f = 0.75,
                               perclimit_low = 0.5, perclimit_mid = 0.4,
                               perclimit_high = 0.0) {
  inc_df <- full_inc_df[order(full_inc_df$ageyears, as.numeric(full_inc_df$id)), ]
  exclusions <- character(0)

  for (round_num in seq_len(max_rounds)) {
    n <- nrow(inc_df)
    if (n < 3) break

    ids   <- as.character(inc_df$id)
    meas  <- inc_df$meas_m
    adays <- inc_df$age_days
    ayrs  <- inc_df$ageyears

    # Age differences to neighbors (in years)
    ageyrs_bef <- c(Inf, diff(ayrs))
    ageyrs_aft <- c(diff(ayrs), Inf)
    agedays_bef <- round(ageyrs_bef * 365.25)
    agedays_aft <- round(ageyrs_aft * 365.25)

    # Minimum adjacent age diff (years), floored at 0
    minagediff <- pmin(ageyrs_bef, ageyrs_aft)
    minagediff[minagediff < 0] <- 0
    months <- minagediff * 12

    # Weight differences to neighbors
    wt_bef <- c(NA, diff(meas))
    wt_aft <- c(diff(meas), NA)

    # EWMA — rebuild each round (window shifts when observations are removed)
    cache <- adult_ewma_cache_init(adays, meas, ewma_window = ewma_window)
    dewma_all <- meas - cache$ewma_all
    dewma_bef <- meas - cache$ewma_before
    dewma_aft <- meas - cache$ewma_after

    # wtallow (with weight scaling via maxwt)
    maxwt_obs <- pmax(meas, cache$ewma_all)
    wta <- compute_wtallow(months, formula = wtallow_formula,
                           cap_params = cap_params, maxwt = maxwt_obs)

    # Observation-level perclimit (each observation's perclimit is based on
    # its own weight). This differs from 11Wa which uses subject-level (max).
    perc_limit <- compute_perc_limit(meas, perclimit_low, perclimit_mid, perclimit_high)

    sn <- seq_len(n)

    # Trajectory rescue
    traj_fails <- compute_trajectory_fails(meas, adays)

    # === STEP 1: Standard pathway + trajectory rescue ===
    # 0.12 kg rounding tolerance on all threshold comparisons
    f <- mod_ewma_f
    std_exc <- (dewma_all > wta + 0.12 & dewma_bef > (f * wta) + 0.12 & dewma_aft > (f * wta) + 0.12) |
               (dewma_all < -(wta + 0.12) & dewma_bef < -((f * wta) + 0.12) & dewma_aft < -((f * wta) + 0.12))
    std_exc[is.na(std_exc)] <- FALSE

    exc_stand <- std_exc & traj_fails

    # Initialize accumulators
    exc_wt_i <- rep(FALSE, n)
    error_load <- rep(FALSE, n)

    # === STEP 2: Standard run detection + pair scoring ===
    runs <- detect_runs(exc_stand)

    # 4+ runs → error load
    is_4plus <- !is.na(runs$run_len) & runs$run_len >= 4 & exc_stand
    error_load[is_4plus] <- TRUE
    exc_wt_i[is_4plus] <- TRUE

    # Isolated (run_len == 1): directly exc_wt_i
    is_isolated <- !is.na(runs$run_len) & runs$run_len == 1 & !error_load
    exc_wt_i[is_isolated] <- TRUE

    # Pairs/trios (run_len 2-3): score
    is_pair_trio <- !is.na(runs$run_len) & runs$run_len >= 2 &
                    runs$run_len <= 3 & !error_load & exc_stand

    if (any(is_pair_trio)) {
      # wtallow_unrel: use gap to non-pair-member
      p2_age <- c(NA, NA, adays[1:(n - 2)])
      n2_age <- c(adays[3:n], NA, NA)

      # When before neighbor is a pair member, use gap to p2
      d_bef_rm <- (adays - p2_age) / 365.25
      d_bef_keep <- ageyrs_aft
      minadiff_bef_rm <- pmin(d_bef_rm, d_bef_keep, na.rm = TRUE)
      minadiff_bef_rm[is.na(d_bef_rm) & is.na(d_bef_keep)] <- NA
      minadiff_bef_rm[is.na(d_bef_rm)] <- d_bef_keep[is.na(d_bef_rm)]
      minadiff_bef_rm[is.na(d_bef_keep)] <- d_bef_rm[is.na(d_bef_keep)]
      wta_bef_unrel <- compute_wtallow(minadiff_bef_rm * 12, formula = wtallow_formula,
                                     cap_params = cap_params, maxwt = maxwt_obs)

      # When after neighbor is a pair member, use gap to n2
      d_aft_rm <- (n2_age - adays) / 365.25
      d_aft_keep <- ageyrs_bef
      minadiff_aft_rm <- pmin(d_aft_keep, d_aft_rm, na.rm = TRUE)
      minadiff_aft_rm[is.na(d_aft_keep) & is.na(d_aft_rm)] <- NA
      minadiff_aft_rm[is.na(d_aft_keep)] <- d_aft_rm[is.na(d_aft_keep)]
      minadiff_aft_rm[is.na(d_aft_rm)] <- d_aft_keep[is.na(d_aft_rm)]
      wta_aft_unrel <- compute_wtallow(minadiff_aft_rm * 12, formula = wtallow_formula,
                                     cap_params = cap_params, maxwt = maxwt_obs)

      score_std <- rep(NA_real_, n)
      is_first <- is_pair_trio & runs$run_pos == 1
      score_std[is_first] <- abs(dewma_aft[is_first] / wta_aft_unrel[is_first])
      is_last <- is_pair_trio & runs$run_pos == runs$run_len
      score_std[is_last] <- abs(dewma_bef[is_last] / wta_bef_unrel[is_last])
      is_mid <- is_pair_trio & runs$run_pos > 1 & runs$run_pos < runs$run_len
      score_std[is_mid] <- abs(dewma_bef[is_mid] / wta_bef_unrel[is_mid] +
                               dewma_aft[is_mid] / wta_aft_unrel[is_mid])

      unique_runs <- unique(runs$run_id[is_pair_trio])
      for (rid in unique_runs) {
        in_run <- which(is_pair_trio & runs$run_id == rid)
        if (length(in_run) > 0) {
          best <- in_run[which.max(score_std[in_run])]
          exc_wt_i[best] <- TRUE
        }
      }
    }

    # === STEP 3: Alternate pathway ===
    # 0.12 kg rounding tolerance on all threshold comparisons
    alt_exc <- rep(FALSE, n)
    prior_unrel <- agedays_bef <= 14 & !is.na(wt_bef) & abs(wt_bef) > wta + 0.12
    alt_exc[prior_unrel] <-
      (dewma_all[prior_unrel] > wta[prior_unrel] + 0.12 &
       dewma_aft[prior_unrel] > (f * wta[prior_unrel]) + 0.12) |
      (dewma_all[prior_unrel] < -(wta[prior_unrel] + 0.12) &
       dewma_aft[prior_unrel] < -((f * wta[prior_unrel]) + 0.12))
    next_unrel <- agedays_aft <= 14 & !is.na(wt_aft) & abs(wt_aft) > wta + 0.12
    alt_exc[next_unrel] <- alt_exc[next_unrel] |
      (dewma_all[next_unrel] > wta[next_unrel] + 0.12 &
       dewma_bef[next_unrel] > (f * wta[next_unrel]) + 0.12) |
      (dewma_all[next_unrel] < -(wta[next_unrel] + 0.12) &
       dewma_bef[next_unrel] < -((f * wta[next_unrel]) + 0.12))
    alt_exc[is.na(alt_exc)] <- FALSE

    exc_pair <- alt_exc & !exc_wt_i

    # === STEP 4: Alternate run detection + pair scoring ===
    alt_runs <- detect_runs(exc_pair)

    is_alt_isolated <- !is.na(alt_runs$run_len) & alt_runs$run_len == 1 &
                       exc_pair & !error_load
    exc_wt_i[is_alt_isolated] <- TRUE

    is_alt_pt <- !is.na(alt_runs$run_len) & alt_runs$run_len >= 2 &
                 alt_runs$run_len <= 3 & !error_load & exc_pair

    if (any(is_alt_pt)) {
      p2_age <- c(NA, NA, adays[1:(n - 2)])
      n2_age <- c(adays[3:n], NA, NA)

      d_bef_rm2 <- (adays - p2_age) / 365.25
      d_bef_keep2 <- ageyrs_aft
      minadiff_bef2 <- pmin(d_bef_rm2, d_bef_keep2, na.rm = TRUE)
      minadiff_bef2[is.na(d_bef_rm2) & is.na(d_bef_keep2)] <- NA
      minadiff_bef2[is.na(d_bef_rm2)] <- d_bef_keep2[is.na(d_bef_rm2)]
      minadiff_bef2[is.na(d_bef_keep2)] <- d_bef_rm2[is.na(d_bef_keep2)]
      wta_bef_unrel2 <- compute_wtallow(minadiff_bef2 * 12, formula = wtallow_formula,
                                      cap_params = cap_params, maxwt = maxwt_obs)

      d_aft_rm2 <- (n2_age - adays) / 365.25
      d_aft_keep2 <- ageyrs_bef
      minadiff_aft2 <- pmin(d_aft_keep2, d_aft_rm2, na.rm = TRUE)
      minadiff_aft2[is.na(d_aft_keep2) & is.na(d_aft_rm2)] <- NA
      minadiff_aft2[is.na(d_aft_keep2)] <- d_aft_rm2[is.na(d_aft_keep2)]
      minadiff_aft2[is.na(d_aft_rm2)] <- d_aft_keep2[is.na(d_aft_rm2)]
      wta_aft_unrel2 <- compute_wtallow(minadiff_aft2 * 12, formula = wtallow_formula,
                                      cap_params = cap_params, maxwt = maxwt_obs)

      score_alt <- rep(NA_real_, n)
      prior_committed <- c(FALSE, exc_wt_i[-n] & !exc_pair[-n])
      next_committed  <- c(exc_wt_i[-1] & !exc_pair[-1], FALSE)

      eff_first  <- is_alt_pt & alt_runs$run_pos == 1 & !prior_committed
      eff_last   <- is_alt_pt & alt_runs$run_pos == alt_runs$run_len & !next_committed
      eff_middle <- is_alt_pt & !eff_first & !eff_last

      score_alt[eff_first]  <- abs(dewma_aft[eff_first] / wta_aft_unrel2[eff_first])
      score_alt[eff_last]   <- abs(dewma_bef[eff_last] / wta_bef_unrel2[eff_last])
      score_alt[eff_middle] <- abs(dewma_bef[eff_middle] / wta_bef_unrel2[eff_middle] +
                                   dewma_aft[eff_middle] / wta_aft_unrel2[eff_middle])

      unique_alt_runs <- unique(alt_runs$run_id[is_alt_pt])
      for (rid in unique_alt_runs) {
        in_run <- which(is_alt_pt & alt_runs$run_id == rid)
        if (length(in_run) > 0) {
          best <- in_run[which.max(score_alt[in_run])]
          exc_wt_i[best] <- TRUE
        }
      }
    }

    # === STEP 5: Percentage criterion ===
    percewma     <- meas / cache$ewma_all
    percewma_bef <- meas / cache$ewma_before
    percewma_aft <- meas / cache$ewma_after

    perc_flag <- percewma < perc_limit & percewma_bef < perc_limit &
                 percewma_aft < perc_limit & !exc_wt_i
    perc_flag[is.na(perc_flag)] <- FALSE
    exc_wt_i[perc_flag] <- TRUE

    # === STEP 6: 4+ consecutive exc_wt_i → error load ===
    consec_runs <- detect_runs(exc_wt_i)
    new_error_load <- !is.na(consec_runs$run_len) & consec_runs$run_len >= 4 & exc_wt_i
    error_load[new_error_load] <- TRUE

    # === STEP 7: Prioritized exclusion ===

    # Error load: all excluded immediately
    el_ids <- ids[error_load]
    if (length(el_ids) > 0) {
      el_code <- paste0(exc_label, "-Error-Load-", round_num)
      exclusions[el_ids] <- el_code
    }

    # Non-error-load candidates
    candidates <- exc_wt_i & !error_load
    if (!any(candidates)) {
      if (any(error_load)) {
        inc_df <- inc_df[!error_load, ]
      } else {
        break
      }
      next
    }

    # Score each candidate using graduated multiplier
    score_final <- rep(NA_real_, n)

    # Edge values: abs(dewma/wta)
    is_edge_first <- candidates & sn == 1
    is_edge_last  <- candidates & sn == n
    score_final[is_edge_first] <- abs(dewma_all[is_edge_first] / wta[is_edge_first])
    score_final[is_edge_last]  <- abs(dewma_all[is_edge_last] / wta[is_edge_last])

    # Interior values: graduated multiplier
    is_middle <- candidates & sn > 1 & sn < n
    if (any(is_middle)) {
      ratio_bef <- abs(dewma_bef[is_middle] / wta[is_middle])
      ratio_aft <- abs(dewma_aft[is_middle] / wta[is_middle])
      min_ratio <- pmin(ratio_bef, ratio_aft)
      multiplier <- pmin(1.0, 0.6 + 0.4 * min_ratio)
      score_final[is_middle] <- multiplier * abs(
        dewma_bef[is_middle] / wta[is_middle] + dewma_aft[is_middle] / wta[is_middle]
      )
    }

    # Tiebreakers: smallest wta → closest to median age → earliest position
    median_age <- median(adays)
    absdiff_median <- abs(adays - median_age)

    cand_idx <- which(candidates)
    ord <- order(-score_final[cand_idx], wta[cand_idx],
                 absdiff_median[cand_idx], sn[cand_idx],
                 as.numeric(ids[cand_idx]))
    best_idx <- cand_idx[ord[1]]

    exc_code <- paste0(exc_label, "-", round_num)
    exclusions[ids[best_idx]] <- exc_code

    # Remove excluded + error load
    to_remove <- error_load
    to_remove[best_idx] <- TRUE
    inc_df <- inc_df[!to_remove, ]
  }

  exclusions
}

# Step 11Wa2: 2D Non-Ordered Pairs ----

#' Evaluate 2D non-ordered weight pairs for one subject
#' @param w_subj_df Weight data with meas_m, age_days, id. Already filtered to Inc.
#' @param w_subj_keep Named vector of current exclusion codes for all weight obs
#' @param wtallow_formula Formula for wtallow
#' @return Character vector of ids to exclude with "2D Non-Ord" code
#' @keywords internal
eval_2d_nonord <- function(w_subj_df, w_subj_keep, wtallow_formula = "piecewise",
                           cap_params = NULL) {
  exc_ids <- character(0)
  if (nrow(w_subj_df) < 2) return(exc_ids)

  vals <- w_subj_df$meas_m
  ids <- as.character(w_subj_df$id)

  # Check: exactly 2 distinct values
  uvals <- unique(vals)
  if (length(uvals) != 2) return(exc_ids)

  # Check: NOT time-ordered (interleaved)
  v1_ages <- w_subj_df$age_days[vals == uvals[1]]
  v2_ages <- w_subj_df$age_days[vals == uvals[2]]
  # Time-ordered if all of one value come before all of the other
  if (max(v1_ages) < min(v2_ages) || max(v2_ages) < min(v1_ages)) {
    return(exc_ids)  # This is ordered → handled by 2D Ord
  }

  # Compute wtallow for each adjacent pair
  w_sorted <- w_subj_df[order(w_subj_df$age_days, as.numeric(w_subj_df$id)), ]
  n <- nrow(w_sorted)
  any_outside <- FALSE

  for (i in 1:(n - 1)) {
    if (w_sorted$meas_m[i] != w_sorted$meas_m[i + 1]) {
      age_diff_months <- abs(w_sorted$age_days[i + 1] - w_sorted$age_days[i]) / 30.4375
      maxwt_pair <- max(w_sorted$meas_m[i], w_sorted$meas_m[i + 1])
      wta <- compute_wtallow(age_diff_months, formula = wtallow_formula,
                             cap_params = cap_params, maxwt = maxwt_pair)
      if (abs(w_sorted$meas_m[i + 1] - w_sorted$meas_m[i]) > wta + 0.12) {
        any_outside <- TRUE
        break
      }
    }
  }

  if (!any_outside) return(exc_ids)  # Rule 1: all within wtallow → keep all

  # Check for prior non-SDE exclusions
  # SDE codes are Identical and Extraneous (ht or wt); all others count as non-SDE.
  # Use exact matching to avoid grepl("Identical") accidentally matching
  # "Exclude-A-WT-Scale-Max-Identical", which is a non-SDE exclusion.
  all_wt_ids <- names(w_subj_keep)
  sde_codes <- c(
    "Exclude-A-HT-Identical", "Exclude-A-WT-Identical",
    "Exclude-A-HT-Extraneous", "Exclude-A-WT-Extraneous"
  )
  prior_nonSDE <- any(
    w_subj_keep[all_wt_ids] != "Include" &
    w_subj_keep[all_wt_ids] != "temp extraneous" &
    !w_subj_keep[all_wt_ids] %in% sde_codes &
    w_subj_keep[all_wt_ids] != "" &
    !is.na(w_subj_keep[all_wt_ids])
  )

  if (prior_nonSDE) {
    # Rule 2: Outside wtallow + prior non-SDE exclusions → all excluded
    return(ids)
  }

  # Check dominance
  count_v1 <- sum(vals == uvals[1])
  count_v2 <- sum(vals == uvals[2])
  total <- count_v1 + count_v2
  dominant_pct <- max(count_v1, count_v2) / total

  if (dominant_pct > 0.65) {
    # Rule 3: dominant > 65% → exclude minority
    minority_val <- if (count_v1 > count_v2) uvals[2] else uvals[1]
    exc_ids <- ids[vals == minority_val]
  } else {
    # Rule 4: each ≤ 65% → all excluded
    exc_ids <- ids
  }

  exc_ids
}

# Step 13: 1D Evaluation ----

#' Evaluate single distinct value observations
#' @param subj_results data.table with result, meas_m, param, ageyears, id columns
#'   for ONE subject — all params included
#' @param params List of limit parameters
#' @return Named character vector of ids to mark as "1D"
#' @keywords internal
eval_1d <- function(subj_results, params) {
  # Two-pass: (1) evaluate with BMI, (2) re-evaluate if BMI lost due to pass 1 exclusion

  exc_ids <- character(0)

  for (pass in 1:2) {
    # Get currently included heights and weights (excluding pass-1 exclusions)
    ht_inc <- subj_results[subj_results$param %in% c("HEIGHTCM", "HEIGHTIN") &
                           subj_results$result == "Include" &
                           !subj_results$id %in% exc_ids, ]
    wt_inc <- subj_results[subj_results$param %in% c("WEIGHTKG", "WEIGHTLBS") &
                           subj_results$result == "Include" &
                           !subj_results$id %in% exc_ids, ]

    ht_1d <- length(unique(ht_inc$meas_m)) == 1 & nrow(ht_inc) > 0
    wt_1d <- length(unique(wt_inc$meas_m)) == 1 & nrow(wt_inc) > 0

    if (!ht_1d && !wt_1d) next

    # Check BMI availability
    ht_days <- if (nrow(ht_inc) > 0) unique(ht_inc$age_days) else integer(0)
    wt_days <- if (nrow(wt_inc) > 0) unique(wt_inc$age_days) else integer(0)
    bmi_days <- intersect(ht_days, wt_days)
    has_bmi <- length(bmi_days) > 0

    if (pass == 1 && has_bmi) {
      # Pass 1: evaluate with BMI
      ht_val <- ht_inc$meas_m[ht_inc$age_days %in% bmi_days][1]
      wt_val <- wt_inc$meas_m[wt_inc$age_days %in% bmi_days][1]
      bmi <- wt_val / ((ht_val / 100)^2)

      bmi_extreme <- bmi < params$bmi_min | bmi > params$bmi_max

      if (bmi_extreme) {
        if (ht_1d) exc_ids <- c(exc_ids, as.character(ht_inc$id))
        if (wt_1d) exc_ids <- c(exc_ids, as.character(wt_inc$id))
      } else {
        # 0.12 tolerance for rounding (0.1 cm/kg rounding + float precision)
        if (ht_1d) {
          ht_val_check <- unique(ht_inc$meas_m)
          if (ht_val_check < params$ht_min_bmi - 0.12 | ht_val_check > params$ht_max_bmi + 0.12) {
            exc_ids <- c(exc_ids, as.character(ht_inc$id))
          }
        }
        if (wt_1d) {
          wt_val_check <- unique(wt_inc$meas_m)
          if (wt_val_check < params$wt_min_bmi - 0.12 | wt_val_check > params$wt_max_bmi + 0.12) {
            exc_ids <- c(exc_ids, as.character(wt_inc$id))
          }
        }
      }
    } else if ((pass == 1 && !has_bmi) || (pass == 2 && !has_bmi)) {
      # No BMI: use tighter limits (0.12 tolerance)
      if (ht_1d) {
        ht_val_check <- unique(ht_inc$meas_m)
        if (ht_val_check < params$ht_min_nobmi - 0.12 | ht_val_check > params$ht_max_nobmi + 0.12) {
          exc_ids <- c(exc_ids, as.character(ht_inc$id))
        }
      }
      if (wt_1d) {
        wt_val_check <- unique(wt_inc$meas_m)
        if (wt_val_check < params$wt_min_nobmi - 0.12 | wt_val_check > params$wt_max_nobmi + 0.12) {
          exc_ids <- c(exc_ids, as.character(wt_inc$id))
        }
      }
    }
    # Pass 2 with BMI: nothing new to check (already evaluated with BMI in pass 1)
  }

  exc_ids
}

# Step 14: Error Load ----

#' Evaluate error load for one subject
#' @param subj_results data.table with result, param columns for ONE subject
#' @param error_threshold Threshold for error ratio (default 0.41)
#' @return Named character vector of ids to mark as "Error load"
#' @keywords internal
eval_error_load <- function(subj_results, error_threshold = 0.41) {
  exc_ids <- character(0)

  # Process height and weight separately
  for (p in c("ht", "wt")) {
    if (p == "ht") {
      p_rows <- subj_results[subj_results$param %in% c("HEIGHTCM", "HEIGHTIN"), ]
    } else {
      p_rows <- subj_results[subj_results$param %in% c("WEIGHTKG", "WEIGHTLBS"), ]
    }
    if (nrow(p_rows) == 0) next

    # Exclude SDEs from denominator
    sde_codes <- c("Exclude-A-HT-Identical", "Exclude-A-HT-Extraneous",
                   "Exclude-A-WT-Identical", "Exclude-A-WT-Extraneous")
    non_sde <- p_rows[!p_rows$result %in% sde_codes, ]
    if (nrow(non_sde) == 0) next

    # Count errors
    error_codes_ht <- c("Exclude-A-HT-BIV",
                        "Exclude-A-HT-Single",
                        "Exclude-A-HT-Ord-Pair",
                        "Exclude-A-HT-Ord-Pair-All",
                        "Exclude-A-HT-Window",
                        "Exclude-A-HT-Window-All")
    error_codes_wt <- c("Exclude-A-WT-BIV",
                        "Exclude-A-WT-Single",
                        "Exclude-A-WT-2D-Ordered",
                        "Exclude-A-WT-2D-Non-Ordered",
                        "Exclude-A-WT-Scale-Max",
                        "Exclude-A-WT-Scale-Max-Identical",
                        "Exclude-A-WT-Scale-Max-RV-Propagated",
                        "Exclude-A-Evil-Twins")

    # EWMA-RV-propagated codes don't count as errors (same underlying error as source).
    # Scale-Max-RV-propagated DOES count (matches Stata "RV 400" behavior).
    is_ewma_propagated <- grepl("^Exclude-A-WT-Traj.*-RV-propagated$", non_sde$result, ignore.case = TRUE)
    is_ewma_error <- grepl("^Exclude-A-WT-Traj", non_sde$result) & !is_ewma_propagated
    if (p == "ht") {
      n_errors <- sum((non_sde$result %in% error_codes_ht | is_ewma_error) & !is_ewma_propagated)
    } else {
      n_errors <- sum((non_sde$result %in% error_codes_wt | is_ewma_error) & !is_ewma_propagated)
    }

    n_inc <- sum(non_sde$result == "Include")
    denom <- n_errors + n_inc
    if (denom < 3) next

    ratio <- n_errors / denom
    if (ratio > error_threshold) {
      # Exclude all remaining Include values
      inc_ids <- as.character(non_sde$id[non_sde$result == "Include"])
      exc_ids <- c(exc_ids, inc_ids)
    }
  }

  exc_ids
}
