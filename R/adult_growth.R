# Main adult growthcleanr function
# internal supporting functions for adults can be found in: adult_support.R

# function to clean height and weight data for adults
# inputs:
# df: data.table with 7 columns:
#   id: row id, must be unique
#   subjid: subject id
#   sex: sex of subject
#   age_years: age, in years
#   param: HEIGHTCM or WEIGHTKG
#   measurement: height or weight measurement
# outputs:
#   df, with additional columns:
#     result, which specifies whether the height measurement should be included,
#       or is implausible (designated with an error code).
#     mean_sde, mean of similar same day extraneous values
cleanadult <- function(df, weight_cap = Inf){
  # method specific constants ----
  # this includes specified cutoffs, etc.

  # BIVs
  biv_df <- data.frame(
    "low" = c(50, 20),
    "high" = c(244, 500)
  )
  rownames(biv_df) <- c("height", "weight")

  # begin implementation ----

  # preallocate final designation
  df[, result := "Include"]
  df[, mean_sde := NA]
  # rownames(df) <- df$id # NO ROWNAMES IN DATA.TABLE
  # go through each subject
  for (i in unique(df$subjid)){
    slog <- df$subjid == i

    # start with height (steps 1 - 7) ----

    h_df <- copy(df[param %in% c("HEIGHTCM", "HEIGHTIN") & slog,])

    # order by age
    h_df <- h_df[order(age_years, id),]

    h_subj_keep <- rep("Include", nrow(h_df))
    h_subj_mean_sde <- rep(NA, nrow(h_df))
    names(h_subj_keep) <- names(h_subj_mean_sde) <-
      h_df$id

    h_subj_df <- copy(h_df)

    # if there are no valid heights, skip
    if (nrow(h_df) > 0){
      # add metric (m) and imperial (im) measurements
      # also keep original measurements (will overwrite for ease in functions)
      h_subj_df$meas_m <- h_subj_df$meas_im <- h_subj_df$meas_orig <-
        h_subj_df$measurement
      h_subj_df$meas_m[h_subj_df$param == "HEIGHTIN"] <-
        h_subj_df$measurement[h_subj_df$param == "HEIGHTIN"]*2.54
      h_subj_df$meas_im[h_subj_df$param == "HEIGHTCM"] <-
        h_subj_df$measurement[h_subj_df$param == "HEIGHTCM"]/2.54
      # convert age years to days -- already in there
      h_subj_df$age_days <- h_subj_df$agedays
    }

    # 1h, H BIV ----
    # 1h. remove biologically impossible height records
    step <- "Exclude-BIV"

    if (nrow(h_subj_df) > 0){
      # overwrite measurement with metric (bivs are in metric)
      h_subj_df$measurement <- h_subj_df$meas_m

      criteria <- remove_biv(h_subj_df, "height", biv_df)
      h_subj_keep[criteria] <- step

      h_subj_df <- h_subj_df[!criteria,]
    }

    # 3h, H temp extraneous ----
    # 3h. Calculate temporary extraneous for exclusion in intermediate steps
    step <- "temp extraneous" # no exclusions here

    if (nrow(h_subj_df) > 0){
      # overwrite measurement with metric (step done in metric)
      h_subj_df$measurement <- h_subj_df$meas_m

      # adds the "extraneous" column, designating whether or not the row is
      # temporarily extraneous (not to be considered in the future)
      # also adds age days and "diff" (unused outside of this function)
      h_subj_df <- temp_sde(h_subj_df)
    }

    # 5h, H hundreds ----
    # 5h. when height goes down by 100 cm -- is it valid?
    step <- "Exclude-Hundreds"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(h_subj_df[!h_subj_df$extraneous,])

    # only do this if there are at least two values
    if (nrow(inc_df) > 1){
      # calculate ewma (using metric)
      ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)
      # delta ewma
      dewma <- (inc_df$meas_m- ewma_res)
      colnames(dewma) <- paste0("d",colnames(ewma_res))

      criteria <- rem_hundreds(inc_df, dewma, "meas_m", 100, "height")

      # update and remove
      h_subj_keep[as.character(inc_df$id)][criteria] <- step

      # don't get rid of extraneous just yet
      h_subj_df <- h_subj_df[h_subj_df$id %in% inc_df$id[!criteria] |
                               h_subj_df$extraneous,]

      # reevaluate temp same day
      h_subj_df <- temp_sde(h_subj_df)
    }

    # 6h, H unit errors ----
    # 6h. checking whether or not height should be a different type of value
    step <- "Exclude-Unit-Errors"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(h_subj_df[!h_subj_df$extraneous,])

    # only do this if there are at least two values
    if (nrow(inc_df) > 1){

      criteria <- rem_unit_errors(inc_df, ptype = "weight")

      # update and remove
      h_subj_keep[as.character(inc_df$id)][criteria] <- step

      # don't get rid of extraneous just yet
      h_subj_df <- h_subj_df[h_subj_df$id %in% inc_df$id[!criteria] |
                               h_subj_df$extraneous,]

      # reevaluate temp same day
      h_subj_df <- temp_sde(h_subj_df)
    }

    # 7h, H transpositions ----
    # 7h. checking whether or not 10s and 1s digit should be switched
    step <- "Exclude-Transpositions"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(h_subj_df[!h_subj_df$extraneous,])
    # only do this if there are at least two unique values
    if (nrow(inc_df) > 1 & length(unique(inc_df$meas_m)) > 2){

      criteria <- rem_transpositions(inc_df, ptype = "height")

      # update and remove
      h_subj_keep[as.character(inc_df$id)][criteria] <- "step"

      # don't get rid of extraneous just yet
      h_subj_df <- h_subj_df[h_subj_df$id %in% inc_df$id[!criteria] |
                               h_subj_df$extraneous,]

      # reevaluate temp same day
      h_subj_df <- temp_sde(h_subj_df)

    }

    # then do weight (steps 1 - 7) ----

    w_df <- copy(df[df$param == "WEIGHTKG" & slog,])

    # order by age
    w_df <- w_df[order(age_years, id),]

    w_subj_keep <- rep("Include", nrow(w_df))
    w_subj_mean_sde <- rep(NA, nrow(w_df))
    names(w_subj_keep) <- names(w_subj_mean_sde) <-
      w_df$id

    w_subj_df <- copy(w_df)

    if (nrow(w_df) > 0){
      # add metric (m) and imperial (im) measurements
      # also keep original measurements (will overwrite for ease in functions)
      w_subj_df$meas_m <- w_subj_df$meas_im <- w_subj_df$meas_orig <-
        w_subj_df$measurement
      w_subj_df$meas_m[w_subj_df$param == "WEIGHTLBS"] <-
        w_subj_df$measurement[w_subj_df$param == "WEIGHTLBS"]/2.2046226
      w_subj_df$meas_im[w_subj_df$param == "WEIGHTKG"] <-
        w_subj_df$measurement[w_subj_df$param == "WEIGHTKG"]*2.2046226

      # convert age years to days -- already exists
      w_subj_df$age_days <- w_subj_df$agedays
    }

    # if there are no valid heights, skip
    if (nrow(w_subj_df) > 0){
      # 1w, W BIV ----
      # 1w. remove biologically impossible weight records
      step <- "Exclude-BIV"

      # overwrite measurement with metric (bivs are in metric)
      w_subj_df$measurement <- w_subj_df$meas_m

      criteria <- remove_biv(w_subj_df, "weight", biv_df)
      w_subj_keep[criteria] <- step

      w_subj_df <- w_subj_df[!criteria,]
    }

    # 2w, W repeated values ----
    # 2w. Identify repeated values (RV) -- values that are the same over different days
    step <- "repeated values" # no exclusions here

    w_subj_df <- identify_rv(w_subj_df)

    # 3w, W temp extraneous ----
    # 3w. Calculate temporary extraneous for exclusion in intermediate steps
    step <- "temp extraneous" # no exclusions here

    if (nrow(w_subj_df) > 0){
      # overwrite measurement with metric (step done in metric)
      w_subj_df$measurement <- w_subj_df$meas_m

      # adds the "extraneous" column, designating whether or not the row is
      # temporarily extraneous (not to be considered in the future)
      # also adds age days and "diff" (unused outside of this function)
      w_subj_df <- temp_sde(w_subj_df, ptype = "weight")

      # redo RVs just if any first RVs became extraneous
      if (any(w_subj_df$extraneous & w_subj_df$is_first_rv)){
        inc_df <- copy(w_subj_df[!w_subj_df$extraneous,])
        inc_df <- identify_rv(inc_df)

        # reassign new rvs to weight df -- ordered the same way
        w_subj_df$is_first_rv <- w_subj_df$is_rv <- F
        w_subj_df$is_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_rv
        w_subj_df$is_first_rv[w_subj_df$id %in% inc_df$id] <- inc_df$is_first_rv
      }
    }

    # 4w, W weight cap ----
    # 4w. excluding specified "weight caps" if a user specifies
    step <- "Exclude-Weight-Cap"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous,])

    # only do this if there are at least two values, and if there's a weight cap
    # to evaluate
    if (nrow(inc_df) > 1 & weight_cap < Inf){
      # weight cap is evaluated with +/ .1 (for precision)
      wc_low <- round(weight_cap, 1) - .1
      wc_high <- round(weight_cap, 1) + .1

      # if all are weight cap, implausible
      is_wc <- check_between(inc_df$meas_m, wc_low, wc_high)
      criteria <- if (all(is_wc)){
        rep(T, nrow(inc_df))
      } else {
        rep(F, nrow(inc_df))
      }

      # if we have somewhere between 1 and < all weight caps, we can evaluate
      if (any(is_wc) & !all(is_wc) & nrow(inc_df) > 2){
        # calculate ewma (using metric)
        ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)
        # delta ewma
        dewma <- (inc_df$meas_m- ewma_res)
        colnames(dewma) <- paste0("d",colnames(ewma_res))

        # figure out weight difference between points
        wt_next <- c(NA, diff(inc_df$meas_m))
        wt_prev <- c(diff(inc_df$meas_m), NA)

        exc_wc <-
          is_wc &
          ((dewma$dewma.all > 50 &
              dewma$dewma.before > .9*dewma$dewma.all &
              dewma$dewma.all > .9*dewma$dewma.all) |
             (dewma$dewma.all < -50 &
                dewma$dewma.before < -.9*dewma$dewma.all &
                dewma$dewma.all < -.9*dewma$dewma.all)) &
          ((wt_prev > 50 & wt_next > 50) |
             (wt_prev < -50 & wt_next < -50))

        criteria <- exc_wc
      }

      # implausible ids from the step
      impl_ids <- as.character(inc_df$id)[criteria]

      # update and remove
      w_subj_keep[c(impl_ids)] <- step

      # don't get rid of extraneous just yet
      w_subj_df <- w_subj_df[!w_subj_df$id %in% impl_ids,]

      # reevaluate first RV
      w_subj_df <- identify_rv(w_subj_df)
      # reevaluate temp same day
      w_subj_df <- temp_sde(w_subj_df, ptype = "weight")
    }

    # 5w, W hundreds ----
    # 5w. when weight goes up/down by 100/200 kg/100-300 lbs -- is it valid?
    step <- "Exclude-Hundreds"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous & !w_subj_df$is_rv,])

    # only do this if there are at least two values
    if (nrow(inc_df) > 1){
      # calculate ewma (using metric)
      ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)
      # delta ewma
      dewma <- (inc_df$meas_m- ewma_res)
      colnames(dewma) <- paste0("d",colnames(ewma_res))

      criteria <- rep(F, nrow(inc_df))
      for (mtype in c("m", "im")){
        # only test certain 100s for metric v imperial
        test_hundreds <- if(mtype == "m"){
          c(100, 200)
        } else {
          c(100, 200, 300)
        }

        # go through and test each combination of hundreds
        for (th in test_hundreds){
          criteria <- criteria |
            rem_hundreds(inc_df, dewma, paste0("meas_", mtype), th, "weight")
        }
      }

      # implausible ids from the step
      impl_ids <- as.character(inc_df$id)[criteria]
      # if it's a repeated value, we want to get rid of it as well
      rv_impl_ids <- as.character(
        w_subj_df$id[w_subj_df$meas_m %in% inc_df$meas_m[criteria &
                                                           inc_df$is_first_rv]]
      )

      # update and remove
      w_subj_keep[impl_ids] <- step
      w_subj_keep[rv_impl_ids] <- paste0(step, "-RV")

      # don't get rid of extraneous just yet
      w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

      # reevaluate temp same day
      w_subj_df <- temp_sde(w_subj_df, ptype = "weight")
    }

    # 6w, W unit errors ----
    # 6w. if a record recorded as metric should be imperial for interior values
    step <- "Exclude-Unit-Errors"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous & !w_subj_df$is_rv,])

    # only do this if there are at least two values
    if (nrow(inc_df) > 1){

      criteria <- rem_unit_errors(inc_df, ptype = "weight")

      # implausible ids from the step
      impl_ids <- as.character(inc_df$id)[criteria]
      # if it's a repeated value, we want to get rid of it as well
      rv_impl_ids <- as.character(
        w_subj_df$id[w_subj_df$meas_m %in% inc_df$meas_m[criteria &
                                                           inc_df$is_first_rv]]
      )

      # update and remove
      w_subj_keep[impl_ids] <- step
      w_subj_keep[rv_impl_ids] <- paste0(step, "-RV")

      # don't get rid of extraneous just yet
      w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

      # reevaluate temp same day
      w_subj_df <- temp_sde(w_subj_df, ptype = "weight")
    }

    # 7w, W transpositions ----
    # 7w. if a record should have swapped the 10s and 1s digits
    step <- "Exclude-Transpositions"

    # we only want to consider subjects without temp extraneous
    inc_df <- copy(w_subj_df[!w_subj_df$extraneous & !w_subj_df$is_rv,])

    # only do this if there are at least two values
    if (nrow(inc_df) > 1 & length(inc_df$meas_m) > 2){

      criteria <- rem_transpositions(inc_df, ptype = "weight")

      # implausible ids from the step
      impl_ids <- as.character(inc_df$id)[criteria]
      # if it's a repeated value, we want to get rid of it as well
      rv_impl_ids <- as.character(
        w_subj_df$id[w_subj_df$meas_m %in% inc_df$meas_m[criteria &
                                                           inc_df$is_first_rv]]
      )

      # update and remove
      w_subj_keep[impl_ids] <- step
      w_subj_keep[rv_impl_ids] <- paste0(step, "-RV")

      # don't get rid of extraneous just yet
      w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

      # reevaluate temp same day
      w_subj_df <- temp_sde(w_subj_df, ptype = "weight")
    }

    # do step 8: swaps (both height and weight) ----
    # 8. checking whether heights and weights should have been swapped
    step <- "Exclude-Swaps"

    if (nrow(h_subj_df) > 0 & nrow(w_subj_df) > 0){
      # we only want to consider subjects without temp extraneous
      h_inc_df <- copy(h_subj_df[!h_subj_df$extraneous,])
      w_inc_df <- copy(w_subj_df[!w_subj_df$extraneous,])

      # calculate ewma
      h_ewma_res <- ewma_dn(h_inc_df$age_days, h_inc_df$meas_m)
      w_ewma_res <- ewma_dn(w_inc_df$age_days, w_inc_df$meas_m)
      # add ewma to dataframes
      h_inc_df <- cbind(h_inc_df, h_ewma_res)
      w_inc_df <- cbind(w_inc_df, w_ewma_res)

      # possible removal of height/weights by bmi
      # h = height, w = weight
      comb_df <- comb_df_orig <-
        merge(h_inc_df, w_inc_df, by = "age_days", all = T,
              suffixes = c(".h", ".w"))
      # remove ones that don't match
      comb_df <- comb_df[!(is.na(comb_df$id.h) | (is.na(comb_df$id.w))),]

      # you need at least two values for this -- can't evaluate next/previous for
      # first and last values
      if (nrow(comb_df) >= 3){
        # identify safe swaps
        # check exact conversion multiples of weight kg to lbs
        # .1 lbs for higher precision, .5 lbs for values with lower precision)
        noswap <- get_float_rem(comb_df$meas_m.w, 2.2046226 * .1) < .01
        noswap <- noswap |
          (get_float_rem(comb_df$meas_m.w, 2.2046226 * .5) < .1 &
             get_float_rem(comb_df$meas_m.w, .1) == 0)
        # check exact conversion multiples of height cm to in
        # .25 in for higher precision (.01), 1 in for lower precision (.1)
        noswap <- noswap |
          get_float_rem(comb_df$meas_m.h, 2.54 * .25) < .01
        noswap <- noswap |
          (get_float_rem(comb_df$meas_m.h, 2.54) < .1 &
             get_float_rem(comb_df$meas_m.h, .1) == 0)

        opposite_map <- c("h" = "w", "w" = "h")

        for (typ in c("h", "w")){
          # create "swaps"
          comb_df[,paste0("swap_",typ)] <-
            comb_df[, paste0("meas_m.", opposite_map[typ]), with = F]

          # calculate difference between values -- DOING ABS HERE
          comb_df[,paste0("swap_prev_",typ)] <-
            c(NA, abs(unlist(
              comb_df[, paste0("swap_",typ), with = F][2:nrow(comb_df)] -
                comb_df[, paste0("meas_m.",typ), with = F][1:(nrow(comb_df)-1)]
            )))
          comb_df[,paste0("swap_next_",typ)] <-
            c(abs(unlist(
              comb_df[, paste0("swap_",typ), with = F][1:(nrow(comb_df)-1)] -
                comb_df[, paste0("meas_m.",typ), with = F][2:nrow(comb_df)]
              )), NA)

          # create dewmas
          for (dtyp in c("all", "before", "after")){
            comb_df[,paste0("dewma.",dtyp,".", typ)] <-
              comb_df[, paste0("meas_m.",typ), with = F] -
              comb_df[,paste0("ewma.",dtyp,".", typ), with = F]

            # dewma for swaps
            comb_df[,paste0("dewma.",dtyp,".swap_", typ)] <-
              comb_df[, paste0("swap_",typ), with = F] -
              comb_df[,paste0("ewma.",dtyp,".", typ), with = F]
          }
        }

        # criteria check - 4 parts
        # check if dewma and one or more diff of recorded values is high
        wt_far <-
          (comb_df$dewma.all.w > 30 &
             comb_df$dewma.before.w > (.9*comb_df$dewma.all.w) &
             comb_df$dewma.after.w > (.9*comb_df$dewma.all.w)) |
          (comb_df$dewma.all.w < -30 &
             comb_df$dewma.before.w < (-.9*comb_df$dewma.all.w) &
             comb_df$dewma.after.w < (-.9*comb_df$dewma.all.w))
        ht_far <-
          (comb_df$dewma.all.h > 12.7 &
             comb_df$dewma.before.h > (.9*comb_df$dewma.all.h) &
             comb_df$dewma.after.h > (.9*comb_df$dewma.all.h)) |
          (comb_df$dewma.all.h < -12.7 &
             comb_df$dewma.before.h < (-.9*comb_df$dewma.all.h) &
             comb_df$dewma.after.h < (-.9*comb_df$dewma.all.h))

        # check if the swapped value is close
        swap_wt_close <-
          (abs(comb_df$dewma.all.swap_w) <= 2 |
             abs(comb_df$swap_prev_w) <= 2 |
             abs(comb_df$swap_next_w) <= 2) #|
        # (abs(comb_df$dewma.all.swap_im_w) <= 2 |
        #    abs(comb_df$swap_im_prev_w) <= 2 |
        #    abs(comb_df$swap_im_next_w) <= 2)
        swap_ht_close <-
          (abs(comb_df$dewma.all.swap_h) <= 2.54 |
             abs(comb_df$swap_prev_h) <= 2.54 |
             abs(comb_df$swap_next_h) <= 2.54) #|
        # (abs(comb_df$dewma.all.swap_im_h) <= 2.54 |
        #    abs(comb_df$swap_im_prev_h) <= 2.54 |
        #    abs(comb_df$swap_im_next_h) <= 2.54)

        exc_swap <- wt_far & ht_far & swap_wt_close & swap_ht_close
        # replace noswaps
        exc_swap[noswap] <- F
        criteria <- exc_swap
        criteria[is.na(criteria)] <- F

        # implausible ids from the step
        impl_ids <- as.character(comb_df$id.w)[criteria]
        # if it's a repeated value, we want to get rid of it as well
        rv_impl_ids <- as.character(
          w_subj_df$id[w_subj_df$meas_m %in% comb_df$meas_m.w[criteria &
                                                                comb_df$is_first_rv]]
        )

        # update and remove -- weight
        w_subj_keep[impl_ids] <- step
        w_subj_keep[rv_impl_ids] <- paste0(step, "-RV")

        # don't get rid of extraneous just yet -- shouldn't be in
        w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

        # update and remove -- height
        h_subj_keep[as.character(comb_df$id.h)][criteria] <- step

        # don't get rid of extraneous just yet
        h_subj_df <- h_subj_df[h_subj_df$id %in% comb_df$id.h[!criteria] |
                                 h_subj_df$extraneous,]

        # reevaluate temp same day -- don't need to reevaluate if nothing has
        # changed
        if (any(criteria)){
          h_subj_df <- temp_sde(h_subj_df)
          w_subj_df <- temp_sde(w_subj_df, ptype = "weight")
        }
      }
    }

    # finish height steps (steps 9-11) ----
    # 9h, H same day extraneous ----
    # 9h. reevaluate same day extraneous for exclusion
    # step specified inside

    # don't do this if there aren't any non extraneous for the subject
    if (nrow(h_subj_df) > 0 & any(h_subj_df$extraneous)){
      step <- "Exclude-Same-Day-Identical"

      # identify duplicate days
      dup_days <- unique(h_subj_df$age_days[h_subj_df$extraneous])

      # first, exclude identical measurements on the same day
      ide_ids <- c() # where we keep the identical ones
      for (dd in dup_days){
        # count amount of unique values
        s_df <- copy(h_subj_df[h_subj_df$age_days == dd,])
        ide_tab <- table(s_df$meas_m)
        if (any(ide_tab > 1)){
          # for each identical, keep only the first one, by id
          # (ordered initially)
          ide_ids <- c(ide_ids, s_df$id[duplicated(
            s_df$meas_m[
              as.character(s_df$meas_m) %in% names(ide_tab[ide_tab > 1])
            ]
          )])
        }
      }
      criteria <- h_subj_df$id %in% ide_ids

      # we're going to update h_subj_df before moving on to the rest of this
      # subject
      h_subj_keep[as.character(h_subj_df$id)][criteria] <- step

      h_subj_df <- h_subj_df[!criteria,]

      if (any(criteria)){
        # reevaluate temp same day
        h_subj_df <- temp_sde(h_subj_df)
        # identify duplicate days
        dup_days <- unique(h_subj_df$age_days[h_subj_df$extraneous])
      }

      step <- "Exclude-Same-Day-Extraneous"
      # now the rest!

      # TODO: CHECK THIS
      # check if heights on duplicate days are trivially the same, keep both,
      # use mean for all of those
      rem_ids <- c()
      for (dd in dup_days){
        s_df <- copy(h_subj_df[h_subj_df$age_days == dd,])
        sde_range <- max(s_df$meas_m) - min(s_df$meas_m)
        if (sde_range < 2.541){ # 1 inch
          h_subj_df$mean_sde[h_subj_df$age_days == dd] <-
            h_subj_df$meas_m[h_subj_df$age_days == dd] <-
            mean(s_df$meas_m)
          # update imperial
          h_subj_df$meas_im[h_subj_df$age_days == dd] <-
            round(mean(s_df$meas_m)/2.54, 2)

          # remove all except the first by id
          rem_ids <- c(rem_ids, h_subj_df$id[h_subj_df$age_days == dd][
            duplicated(h_subj_df$meas_m[h_subj_df$age_days == dd])])
        }
      }
      criteria <- h_subj_df$id %in% rem_ids

      # we're going to update h_subj_df before moving on to the rest of this
      # subject
      h_subj_keep[as.character(h_subj_df$id)][criteria] <- step
      # also update mean sde
      h_subj_mean_sde[as.character(h_subj_df$id)] <- h_subj_df$mean_sde

      h_subj_df <- h_subj_df[!criteria,]

      if (any(criteria)){
        # reevaluate temp same day
        h_subj_df <- temp_sde(h_subj_df)
        # identify duplicate days
        dup_days <- unique(h_subj_df$age_days[h_subj_df$extraneous])
      }

      # next, calculate the duplicate ratio -- what proportion of days are
      # duplicated
      dup_ratio <- mean(!is.na(h_subj_df$diff))
      # also check whether or not any same-days are adjacent -- need 4 at least
      # rolling windows of day differences -- we are looking for 0,x,0
      if (nrow(h_subj_df) > 3){
        roll <- embed(diff(h_subj_df$age_days), 3)
        adjacent <- any(sapply(1:nrow(roll), function(x){
          all(c(roll[x,1] == 0, roll[x,2] != 0, roll[x,3] == 0))
        }))
      } else {
        adjacent <- F
      }

      # if dup ratio is too high, or any adjacent same days, we exclude all
      criteria <-
        if ((dup_ratio > .25) | adjacent){
          rep(T, nrow(h_subj_df))
        } else {
          rep(F, nrow(h_subj_df))
        }

      # if criteria didn't catch it, we now compare with medians
      if (!all(criteria) & any(h_subj_df$extraneous)){
        med <- median(h_subj_df$measurement[
          !h_subj_df$age_days %in% dup_days
        ])
        h_subj_df$absdiff <- NA
        h_subj_df$absdiff[h_subj_df$age_days %in% dup_days] <-
          abs(med - h_subj_df$diff[h_subj_df$age_days %in% dup_days])

        # go through each duplicate day and figure out criteria
        # keep the first one that satisfies this criteria
        # TODO: CHECK THIS
        mm <- h_subj_df$meas_m # shorthand
        rem_ids <- c()
        for (dd in dup_days){
          # find prior and next index
          curr_ind <- which(h_subj_df$age_days == dd)[1]
          p_ind <- which(h_subj_df$age_days == dd)[1] - 1
          n_ind <- which(h_subj_df$age_days == dd)[
            sum(h_subj_df$age_days == dd)] + 1
          # modify for first and last
          if (p_ind < 0){
            p_ind <- 1
          }
          if (n_ind > nrow(h_subj_df)){
            n_ind <- nrow(h_subj_df)
          }

          compare <-
            abs(mm[curr_ind] - mm[p_ind]) < 2.541 &
            abs(mm[curr_ind] - mm[n_ind]) < 2.541 &
            h_subj_df$diff[h_subj_df$age_days == dd] < 2.541

          comp_id <-
            if (sum(compare) == 1){
              h_subj_df$id[which(compare)]
            } else {
              ""
            }
          # the rest are to be excluded
          rem_ids <- h_subj_df$id[h_subj_df$age_days == dd &
                                    h_subj_df$id != comp_id]
        }
        # update criteria
        criteria <- h_subj_df$id %in% rem_ids
      }

      h_subj_keep[as.character(h_subj_df$id)][criteria] <- step

      h_subj_df <- h_subj_df[!criteria,]
    }

    # 10hab, H distinct values ----
    # 10ha-b. evaluate based on numbers of distinct values
    # steps specified in specific caveat

    # get number of distinct values
    num_distinct <- length(unique(h_subj_df$meas_m))

    # preallocate criteria
    criteria <- rep(F, nrow(h_subj_df))

    # go through each type of exclusion
    if (num_distinct == 2){
      # 10ha, H distinct pairs ----
      step <- "Exclude-Distinct-Pairs"

      # identify "height 1 and 2" and their corresponding ages
      ht_1 <- unique(h_subj_df$meas_m[order(h_subj_df$age_days)])[1]
      ht_2 <- unique(h_subj_df$meas_m[order(h_subj_df$age_days)])[2]
      ht_1_log <- h_subj_df$meas_m == ht_1
      ht_2_log <- h_subj_df$meas_m == ht_2
      ht_1_imp <- unique(h_subj_df$meas_im[ht_1_log])
      ht_2_imp <- unique(h_subj_df$meas_im[ht_2_log])
      ht_1_ageyears <- h_subj_df$age_years[ht_1_log]
      ht_2_ageyears <- h_subj_df$age_years[ht_2_log]

      # check if pairs outside two inch range
      # imperial will also be unique
      exc_2d <- abs(ht_1_imp - ht_2_imp) > 2

      # only if outside the range
      if (exc_2d){
        # calculate first and second ageyears
        ageyears1 <- sort(ht_1_ageyears)[1]
        ageyears2 <- sort(ht_2_ageyears)[1]

        # calculate "height allow": allowed growth for ages <= 25
        if (ageyears1 <= 25){
          htcompare <- ifelse(ageyears2 > 25, 25, ageyears2)

          # using short circuiting
          hta <-
            if ((ageyears2 - ageyears1) < 2){
              ht_allow(15.5, ageyears1, htcompare)
            } else if ((ageyears2 - ageyears1) <= 3){
              ht_allow(13, ageyears1, htcompare)
            } else if ((ageyears2 - ageyears1) > 3){
              ht_allow(12, ageyears1, htcompare)
            }

          # potential reallow: height gain for youth gain in height
          pairhtgain <-
            (ht_2_imp - ht_1_imp) <= (hta + 2) &
            (ht_2_imp - ht_1_imp) > 0 &
            min(ht_2_ageyears) > max(ht_1_ageyears)
        } else {
          # this is not a check to apply if they're outside the age range
          pairhtgain <- T
        }

        # generate keep ratios for height 1 and 2
        keepht1 <- sum(ht_1_log) >= (sum(ht_2_log)*(4/3))
        keepht2 <- sum(ht_2_log) >= (sum(ht_1_log)*(4/3))

        # potential reallow: falls are <= 3 (+2) in or <= 5 (+2) for ageyears > 50
        pairhtloss <-
          (ht_1 > ht_2) &
          ((ht_1_imp - ht_2_imp) <= (3+2) |
             ((ht_1_imp - ht_2_imp) <= (5+2) & ageyears2 >= 50)) &
          min(ht_2_ageyears) > max(ht_1_ageyears)

        # exclude if there are no reallow for loss and gain
        exc_pairs <- !(pairhtloss & pairhtgain)
        if (exc_pairs){
          criteria <- rep(T, nrow(h_subj_df))
        }

        # reallow if the ratios are satisfied
        if (keepht1){
          criteria[ht_1_log] <- F
        }
        if (keepht2){
          criteria[ht_2_log] <- F
        }
      }

    } else if (num_distinct >= 3){
      # 10ha, H distinct 3 or more ----
      step <- "Exclude-Distinct-3-Or-More"

      h_subj_df <- h_subj_df[order(h_subj_df$age_years),]
      # create w2 (w/in 2) and o2  (outside 2) groups
      w2_groups <- lapply(unique(h_subj_df$meas_m), function(x){
        h_subj_df$meas_m[check_between(h_subj_df$meas_m, x, x+5.081)]
      })
      o2_groups <- lapply(unique(h_subj_df$meas_m), function(x){
        h_subj_df$meas_m[!check_between(h_subj_df$meas_m, x, x+5.081)]
      })
      # will coerce to character
      names(w2_groups) <- names(o2_groups) <- unique(h_subj_df$meas_m)

      # calculate keep ratio (if nothing's outside, it's infinite = good)
      ratio_w2o2 <- sapply(unique(h_subj_df$meas_m), function(x){
        length(w2_groups[[as.character(x)]])/length(o2_groups[[as.character(x)]])
      })
      # calculate which ones are ok (1.5 x more in w2 than o2)
      okratio <- ratio_w2o2 > 3/2

      # if there are multiple ok, we need some tiebreakers
      if (sum(okratio) > 1){
        # w2 groups for consideration
        consider_w2 <- unique(h_subj_df$meas_m)[okratio]

        # create the total and distinct score
        score <- sapply(consider_w2, function(x){
          length(w2_groups[[as.character(x)]]) +
            .5*length(unique(w2_groups[[as.character(x)]]))
        })

        best_scores <- which(score == max(score))
        if (length(best_scores) > 1){
          # if there are multiple best scores, we move to our next tiebreaker
          # calculating mean absolute distance between w2 and o2
          mean_abs_dist <- sapply(consider_w2[best_scores], function(x){
            mean(abs(mean(w2_groups[[as.character(x)]]) -
                       o2_groups[[as.character(x)]]))
          })

          # this will choose the first if there's a tiebreaker
          best_w2 <- consider_w2[best_scores][which.max(mean_abs_dist)]
        } else {
          best_w2 <- consider_w2[best_scores]
        }
      } else if (sum(okratio) == 1){
        best_w2 <- unique(h_subj_df$meas_m)[okratio]
      } else {
        best_w2 <- "none"
      }

      if (best_w2 != "none"){ # one of the w2s was selected -- o2s are out
        criteria[!h_subj_df$meas_m %in% w2_groups[[as.character(best_w2)]]] <- T
      } else { # none -- all are false
        criteria <- rep(T, nrow(h_subj_df))
      }

      # 3D loss ----

      # if anything was excluded, we want to possibly reinclude with gain/loss code
      if (any(criteria)){
        # create 3 groups:
        # g1: ht 1 and all within 2 inches
        # g2: first height outside g1 + 2 inches
        # g3: first height outside g2 + 2 inches
        # ^ if any heights out, we want to go with what we had originally

        gtotal <- ht_change_groups(h_subj_df, 3)
        glist <- gtotal$meas # groups with measurements
        galist <- gtotal$age # groups with ages, in years

        # if there are any outside these groups, we don't attempt to fix
        # we need all of them to be in some group
        # if (all(Reduce("|", glist))){
        if (length(glist) <= 3){
          # create a mean ht and a min age for each group
          mean_ht <- sapply(glist, function(g){
            suppressWarnings(mean(g))
          })
          min_age <- sapply(galist, function(ga){
            suppressWarnings(min(ga))
          })

          # check g2 v g1 -- true indicates reinclude all
          g2_g1_incl_check <-
            if (!is.na(mean_ht[2])){
              (mean_ht[2] - mean_ht[1]) < 0 &
                ((min_age[2] < 50 &
                    (mean_ht[2] - mean_ht[1]) > ((-5 * 2.54) +.001)) |
                   (min_age[2] >= 50 &
                      (mean_ht[2] - mean_ht[1]) > ((-7 * 2.54) +.001)))
            } else {
              F
            }

          # check g3 v g2 -- true indicates use the original exclusions
          g3_g2_check <-
            if (!is.na(mean_ht[2]) & !is.na(mean_ht[3])){
              (mean_ht[3] - mean_ht[2]) > 0 |
                (min_age[3] < 50 &
                   (mean_ht[3] - mean_ht[2]) < ((-5 * 2.54) +.001)) |
                (min_age[3] >= 50 &
                   (mean_ht[3] - mean_ht[2]) < ((-7 * 2.54) +.001))
            } else {
              F
            }

          # check g3 v g1 -- true indicates use the original exclusions
          g3_g1_check <-
            if (!is.na(mean_ht[2]) & !is.na(mean_ht[3])){
              (min_age[3] < 50 &
                 (mean_ht[3] - mean_ht[1]) < ((-6 * 2.54) +.001)) |
                (min_age[2] < 50 & min_age[3] >= 50 &
                   (mean_ht[3] - mean_ht[1]) < ((-8 * 2.54) +.001)) |
                (min_age[2] >= 50 &
                   (mean_ht[3] - mean_ht[1]) < ((-9 * 2.54) +.001))
            } else {
              F
            }

          # if g2 g1 check passess, reinclude all
          # TODO: CHECK IF SHOULD BE ANY OR ALL
          if (any(g2_g1_incl_check)){
            criteria <- rep(F, nrow(h_subj_df))
          }

          # if all are false, reinclude all?
          # TODO: CHECK THIS
          if (all(!c(g3_g2_check, g3_g1_check))){
            criteria <- rep(F, nrow(h_subj_df))
          }
        }

      }

      # 3D gain ----

      # account for gain -- only for subjects with < 25 years
      if (any(criteria) & min(h_subj_df$age_years) < 25){
        # create 6 groups:
        # g1: ht 1 and all within 2 inches
        # g2: first height outside g1 + 2 inches ....
        # ^ if any heights out, we want to go with what we had originally

        gtotal <- ht_change_groups(h_subj_df, 6)
        glist <- gtotal$meas # groups with measurements
        galist <- gtotal$age # groups with ages

        # if there are any outside these groups, we don't attempt to fix
        # we need all of them to be in some group
        if (length(glist) <= 6){
          # create a mean ht and a min age for each group
          mean_ht <- sapply(glist, function(g){
            suppressWarnings(mean(g))
          })
          min_age <- sapply(galist, function(ga){
            suppressWarnings(min(ga))
          })

          # preallocating on whether or not we want to go by original exclusion
          origexc <- F
          # first, compare to the age before
          origexc <- origexc |
            ht_3d_growth_compare(mean_ht, min_age, glist, compare = "before")
          # then, compare to the first age group
          origexc <- origexc |
            ht_3d_growth_compare(mean_ht, min_age, glist, compare = "first")

          # if all are false, reinclude all?
          if (!origexc){
            criteria <- rep(F, nrow(h_subj_df))
          }
        }
      }

    }

    # update and remove
    h_subj_keep[as.character(h_subj_df$id)][criteria] <- step

    h_subj_df <- h_subj_df[!criteria,]

    # finish weight steps (steps 9-12) ----

    # 9w, W extreme EWMA ----
    # 9w. mark extreme values using EWMA method
    step <- "Exclude-Extreme-EWMA"

    # first, remove ewma without temp extraneous and repeated values
    # we only want to consider subjects without temp extraneous
    inc_df_first <- copy(w_subj_df[!w_subj_df$extraneous & !w_subj_df$is_rv,])
    # only do this if there are at least two values
    criteria_first <-
      if (nrow(inc_df_first) > 1){
        remove_ewma_wt(inc_df_first)
      } else {
        rep(F, nrow(inc_df_first))
      }

    # then, remove ewma just without temp extraneous
    inc_df_rv <- copy(w_subj_df[
      !w_subj_df$extraneous & !w_subj_df$id %in% inc_df_first$id[criteria_first],
      ])
    # only do this if there are at least two UNIQUE values
    criteria_rv <-
      if (length(unique(inc_df_rv)) > 1){
        remove_ewma_wt(inc_df_rv)
      } else {
        rep(F, nrow(inc_df_rv))
      }

    # implausible ids from the step
    impl_ids <- c(
      as.character(inc_df_first$id)[criteria_first],
      as.character(inc_df_rv$id)[criteria_rv]
    )
    # if it's a repeated value, we want to get rid of it as well
    rv_impl_ids <- as.character(
      w_subj_df$id[w_subj_df$meas_m %in%
                     inc_df_first$meas_m[criteria_first &
                                           inc_df_first$is_first_rv]],
      w_subj_df$id[w_subj_df$meas_m %in%
                     inc_df_rv$meas_m[criteria_rv &
                                        inc_df_rv$is_first_rv]]
    )

    # update and remove
    w_subj_keep[impl_ids] <- step
    w_subj_keep[rv_impl_ids] <- paste0(step, "-RV")

    # don't get rid of extraneous just yet
    w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

    # don't need to do this if we don't find anything to remove
    if (length(c(impl_ids, rv_impl_ids)) > 0){
      # reevaluate first rvs
      w_subj_df <- identify_rv(w_subj_df)

      # don't need to reevaluate temp same day -- next step we get rid of them
    }

    # 10w, W same day extraneous ----
    # 10w. check and remove same day extraneous values

    # includes all RV

    # don't do this if there aren't any non extraneous for the subject
    if (nrow(w_subj_df) > 0 & any(w_subj_df$extraneous)){
      step <- "10w, W same day identical"

      # identify duplicate days
      dup_days <- unique(w_subj_df$age_days[w_subj_df$extraneous])

      # first, exclude identical measurements on the same day
      ide_ids <- c() # where we keep the identical ones
      for (dd in dup_days){
        # count amount of unique values
        s_df <- w_subj_df[w_subj_df$age_days == dd,]
        ide_tab <- table(s_df$meas_m)
        if (any(ide_tab > 1)){
          # for each identical, keep only the first one, by id
          # (ordered initially)
          ide_ids <- c(ide_ids, s_df$id[duplicated(
            s_df$meas_m[
              as.character(s_df$meas_m) %in% names(ide_tab[ide_tab > 1])
            ]
          )])
        }
      }
      criteria <- w_subj_df$id %in% ide_ids

      # we're going to update h_subj_df before moving on to the rest of this
      # subject
      w_subj_keep[as.character(w_subj_df$id)][criteria] <-
        "Implausible"
      w_subj_reason[as.character(w_subj_df$id)][criteria] <-
        paste0("Implausible, Step ",step)

      w_subj_df <- w_subj_df[!criteria,]

      if (any(criteria)){
        # reevaluate temp same day
        w_subj_df <- temp_sde(w_subj_df)
        # identify duplicate days
        dup_days <- unique(w_subj_df$age_days[w_subj_df$extraneous])
      }

      step <- "10w, W same day extraneous"
      # now the rest!

      # check if heights on duplicate days are trivially the same, keep both,
      # use mean for all of those
      rem_ids <- c()
      for (dd in dup_days){
        s_df <- w_subj_df[w_subj_df$age_days == dd,]
        sde_range <- max(s_df$meas_m) - min(s_df$meas_m)
        if (sde_range < 1.001){ # 1 pound
          w_subj_df$mean_sde[w_subj_df$age_days == dd] <-
            w_subj_df$meas_m[w_subj_df$age_days == dd] <-
            mean(s_df$meas_m)
          # update imperial
          w_subj_df$meas_im[w_subj_df$age_days == dd] <-
            round(mean(s_df$meas_m)*2.2046226, 1)

          # remove all except the first by id
          rem_ids <- c(rem_ids, w_subj_df$id[w_subj_df$age_days == dd][
            duplicated(w_subj_df$meas_m[w_subj_df$age_days == dd])])
        }
      }
      criteria <- w_subj_df$id %in% rem_ids

      # we're going to update w_subj_df before moving on to the rest of this
      # subject
      w_subj_keep[as.character(w_subj_df$id)][criteria] <-
        "Implausible"
      w_subj_reason[as.character(w_subj_df$id)][criteria] <-
        paste0("Implausible, Step ",step)
      # also update mean sde
      w_subj_mean_sde[as.character(w_subj_df$id)] <- w_subj_df$mean_sde

      w_subj_df <- w_subj_df[!criteria,]

      if (any(criteria)){
        # reevaluate temp same day
        w_subj_df <- temp_sde(w_subj_df)
        # identify duplicate days
        dup_days <- unique(w_subj_df$age_days[w_subj_df$extraneous])
      }

      # next, calculate the duplicate ratio -- what proportion of days are
      # duplicated
      dup_ratio <- mean(!is.na(w_subj_df$diff))
      # also check whether or not any same-days are adjacent -- need 4 at least
      # rolling windows of day differences -- we are looking for 0,x,0
      if (nrow(w_subj_df) > 3){
        roll <- embed(diff(w_subj_df$age_days), 3)
        adjacent <- any(sapply(1:nrow(roll), function(x){
          all(c(roll[x,1] == 0, roll[x,2] != 0, roll[x,3] == 0))
        }))
      } else {
        adjacent <- F
      }

      # if dup ratio is too high, or any adjacent same days, we exclude all
      criteria <-
        if ((dup_ratio > .25) | adjacent){
          rep(T, nrow(w_subj_df))
        } else {
          rep(F, nrow(w_subj_df))
        }

      # if criteria didn't catch it, we now compare with medians
      if (!all(criteria)){
        # calculate ewma
        # calculate ewma (using metric)
        ewma_res <- ewma_dn(w_subj_df$age_days, w_subj_df$meas_m,
                            ewma.adjacent = F)
        # delta ewma
        dewma <- (w_subj_df$meas_m- ewma_res)
        colnames(dewma) <- paste0("d",colnames(ewma_res))

        # go through each duplicate day and figure out criteria
        # keep the first one that satisfies this criteria
        rem_ids <- c()
        for (dd in dup_days){
          # TODO: DOES COMPARE ALWAYS IMPLY LESS 5?
          compare <- abs(dewma$dewma.all[w_subj_df$age_days == dd]) <= 2.001
          # TODO: CHECK IF THIS IS ONLY FOR THE SPECIFIC DAY AND SUBJECT
          less_5 <- abs(dewma$dewma.all[w_subj_df$age_days == dd]) < 5

          comp_id <-
            if (sum(compare & less_5) == 1){
              w_subj_df$id[which(compare)]
            } else {
              ""
            }
          # the rest are to be excluded
          rem_ids <- w_subj_df$id[w_subj_df$age_days == dd &
                                    w_subj_df$id != comp_id]
        }
        # update criteria
        criteria <- w_subj_df$id %in% rem_ids
      }

      w_subj_keep[as.character(w_subj_df$id)][criteria] <-
        "Implausible"
      w_subj_reason[as.character(w_subj_df$id)][criteria] <-
        paste0("Implausible, Step ",step)

      w_subj_df <- w_subj_df[!criteria,]
    }

    # 11w, W distinct/moderate EWMA ----
    # 11w. Check remaining weight values depending on amount of values within
    # EWMA

    # check if there are two distinct, and if so, are they ordered?
    pair_distinct <- if (length(unique(w_subj_df$meas_m)) == 2){
      if (
        max(w_subj_df$age_years[w_subj_df$meas_m == unique(w_subj_df$meas_m)[1]]) <
        min(w_subj_df$age_years[w_subj_df$meas_m == unique(w_subj_df$meas_m)[2]])
      ){
        T
      } else {
        F
      }
    } else {
      F
    }

    # use all RV
    criteria <- rep(F, nrow(inc_df))

    # first, identify implausible for distinct pairs
    if (pair_distinct){
      # 11wa, W distinct ordered pairs ----
      # 11wa. Check pairs (2 distinct ordered values), where all first values are
      # of ages less than second values
      step <- "11wa, W distinct ordered pairs"

      # get first and last weight -- already ordered by age
      # all first are before all last, so these will be the two distinct values
      wt_first <- w_subj_df$meas_m[1]
      wt_last <- w_subj_df$meas_m[nrow(w_subj_df)]
      # weight difference at beginning and end of series
      wt_diff <- wt_last - wt_first

      # get age difference, in years
      ageyears_diff <-
        min(w_subj_df$age_years[w_subj_df$meas_m == unique(w_subj_df$meas_m)[2]]) -
        max(w_subj_df$age_years[w_subj_df$meas_m == unique(w_subj_df$meas_m)[1]])

      # compute "weight allow" how much change is allowed over time
      wta <- 4 + 18*log(1 + (ageyears_diff*12))
      # cap at 60
      wta[wta > 60] <- 60

      # is the difference outside the allowed weight change?
      exc_pairs <- abs(wt_diff) > wta

      # another criteria: consider the weight ratio
      wt_perc <-
        if (wt_last/wt_first < 1){
          wt_last/wt_first
        } else {
          wt_first/wt_last
        }

      # TODO: ASK WHAT THIS MEANS
      # set a limit for wts to be percent of other wts, focused on lower wts
      perc_limit <- .7
      if (any(w_subj_df$meas_m > 45)){
        perc_limit <- .4
      }

      # exclude if outside the limit (exc pairs will be a vector)
      exc_pairs <- exc_pairs | wt_perc < perc_limit

      criteria <- exc_pairs
    } else if (length(unique(w_subj_df$meas_m)) >= 2){ # (but not ordered)
      # 11wb, W moderate EWMA ----
      # 11wb. Check all other types, using a more moderate EWMA cutoff and other
      # criteria
      step <- "11wb, W moderate EWMA"

      inc_df <- w_subj_df

      # exclude the most extreme, then recalculate again and again
      rem_ids <- c()
      change <- T
      iter <- 1
      while (change){
        # TODO: ASK WHAT THIS MEANS -- VECTOR?
        # set a limit for wts to be percent of other wts, focused on lower wts
        perc_limit <- .7
        if (any(w_subj_df$meas_m > 45)){
          perc_limit <- .4
        }

        # figure out time difference between points
        ageyears_bef <- c(Inf, diff(inc_df$age_years))
        ageyears_aft <- c(diff(inc_df$age_years), Inf)
        minagediff <- ageyears_bef
        minagediff[ageyears_aft < ageyears_bef] <-
          ageyears_aft[ageyears_aft < ageyears_bef]
        # convert to days - rounded
        agedays_bef <- round(ageyears_bef*365.25)
        agedays_aft <- round(ageyears_aft*365.25)

        # figure out weight difference between points
        wt_bef <- c(NA, diff(inc_df$meas_m))
        wt_aft <- c(diff(inc_df$meas_m), NA)

        # 'polation (inter/extra-polation)
        # "interpolation" - between prior and next with error of 5 on either end
        binerr_interpol <- c(NA, sapply(2:(nrow(inc_df)-1), function(x){
          check_between(inc_df$meas_m[x],
                        inc_df$meas_m[x-1]-5, inc_df$meas_m[x+1]+5) |
            check_between(inc_df$meas_m[x],
                          inc_df$meas_m[x+1]-5, inc_df$meas_m[x-1]+5)
        }), NA)
        # extrapolation -- prior weights
        lepolate_p <- binerr_lepolate_p <- c(rep(NA,2))
        for (x in 3:nrow(inc_df)){
          slope <- (inc_df$meas_m[x-1] - inc_df$meas_m[x-2])/
            (inc_df$age_days[x-1] - inc_df$age_days[x-2])
          lepolate_p <- c(lepolate_p, round_pt(
            inc_df$meas_m[x-1] +
              slope*(inc_df$meas_m[x]-inc_df$meas_m[x-1]),
            .2
          ))

          # is current value between extrapolated and 2 previous
          binerr_lepolate_p <- c(
            binerr_lepolate_p,
            check_between(inc_df$meas_m[x],
                          inc_df$meas_m[x - 2] - 5, lepolate_p[x] + 5) |
              check_between(inc_df$meas_m[x],
                            lepolate_p[x] - 5, inc_df$meas_m[x - 2] + 5)
          )
        }
        # extrapolation -- next weights
        lepolate_n <- binerr_lepolate_n <- c()
        for (x in 1:(nrow(inc_df)-2)){
          slope <- (inc_df$meas_m[x+1] - inc_df$meas_m[x+2])/
            (inc_df$age_days[x+1] - inc_df$age_days[x+2])
          lepolate_n <- c(lepolate_n, round_pt(
            inc_df$meas_m[x+1] +
              slope*(inc_df$meas_m[x+1]-inc_df$meas_m[x]),
            .2
          ))

          # is current value between extrapolated and 2 next
          binerr_lepolate_n <- c(
            binerr_lepolate_n,
            check_between(inc_df$meas_m[x],
                          inc_df$meas_m[x + 2] - 5, lepolate_n[x] + 5) |
              check_between(inc_df$meas_m[x],
                            lepolate_n[x] - 5, inc_df$meas_m[x + 2] + 5)
          )
        }
        lepolate_n <- c(lepolate_n, rep(NA,2))
        binerr_lepolate_n <- c(binerr_lepolate_n, rep(NA,2))

        # compute "weight allow" how much change is allowed over time
        wta <- 4 + 18*log(1 + (minagediff*12))
        # cap at 60
        wta[wta > 60] <- 60

        # calculate ewma
        ewma_res <- ewma_dn(inc_df$age_days, inc_df$meas_m)

        dewma <- inc_df$meas_m - ewma_res
        colnames(dewma) <- paste0("d",colnames(ewma_res))

        # moderate EWMA exclusion criteria
        exc_mod_ewma <-
          (dewma$dewma.all > wta &
             dewma$dewma.before > (.75*wta) &
             dewma$dewma.after > (.75*wta)) |
          (dewma$dewma.all < -1*wta &
             dewma$dewma.before < (-1*.75*wta) &
             dewma$dewma.after < (-1*.75*wta))

        # alternate ewma criteria -- if difference with adjacent within wta and
        # difference in agedays <= 14
        # prior is close in age but far in weight
        alt_ewma_exc <-
          agedays_bef <= 14 &
          abs(wt_bef) > wta |
          (dewma$dewma.all > wta &
             dewma$dewma.after > (.75*wta)) |
          (dewma$dewma.all < -1*wta &
             dewma$dewma.after < (-1*.75*wta))
        # next is close in age but far in weight
        alt_ewma_exc <- alt_ewma_exc |
          (agedays_aft <= 14 &
             abs(wt_aft) > wta |
             (dewma$dewma.all > wta &
                dewma$dewma.before > (.75*wta)) |
             (dewma$dewma.all < -1*wta &
                dewma$dewma.before < (-1*.75*wta)))
        exc_mod_ewma <- exc_mod_ewma | alt_ewma_exc
        exc_mod_ewma[is.na(exc_mod_ewma)] <- F

        # identify binerr criteria -- edge ones are marked as true
        binerr_lepolate_n[is.na(binerr_lepolate_n)] <- T
        binerr_lepolate_p[is.na(binerr_lepolate_p)] <- T
        binerr_interpol[is.na(binerr_interpol)] <- F
        exc_binerr <- binerr_lepolate_p & binerr_lepolate_n & !binerr_interpol
        # alternate binerr crtieria -- if difference with adjacent within wta and
        # difference in agedays <= 14
        alt_exc_binerr <-
          agedays_bef <= 14 &
          abs(wt_bef) > wta &
          binerr_lepolate_n
        alt_exc_binerr <- alt_exc_binerr |
          agedays_aft <= 14 &
          abs(wt_aft) > wta &
          binerr_lepolate_p
        exc_binerr <- exc_binerr | alt_exc_binerr
        exc_binerr[is.na(exc_binerr)] <- F

        # ratio to ewma can also lead to exclusion
        percewma <- inc_df$meas_m/ewma_res
        exc_perc <- percewma$ewma.all < perc_limit &
          percewma$ewma.before < perc_limit &
          percewma$ewma.after < perc_limit

        criteria_new <- (exc_mod_ewma & exc_binerr) | exc_perc

        if (all(!criteria_new)){
          change <- F
        } else {
          # ewma ratio determines which to exclude -- highest ewmaratio
          ewmaratio <- abs(dewma$dewma.all)/wta
          # boost middle values -- be more strict
          ewmaratio[-c(1, length(ewmaratio))] <-
            ewmaratio[-c(1, length(ewmaratio))] + .2

          # figure out the most extreme value and remove it and rerun
          to_rem <- which.max(ewmaratio[criteria_new])

          # keep the ids that failed and remove
          rem_ids[length(rem_ids)+1] <- inc_df[criteria_new,][to_rem, "id"]
          inc_df <- inc_df[inc_df$id != rem_ids[length(rem_ids)],]

          # check if this is viable -- you need at least three points, otherwise
          # we're done
          if (nrow(inc_df) < 3){
            change <- F
          }
          # update iteration
          iter <- iter + 1
        }
      }

      # form results into a logical vector
      criteria <- rep(F, nrow(w_subj_df))
      criteria[w_subj_df$id %in% rem_ids] <- T
    }

    # implausible ids from the step
    impl_ids <- as.character(w_subj_df$id)[criteria]

    # update and remove
    w_subj_keep[c(impl_ids)] <- "Implausible"
    w_subj_reason[impl_ids] <- paste0("Implausible, Step ",step)

    w_subj_df <- w_subj_df[!w_subj_df$id %in% c(impl_ids, rv_impl_ids),]

    # 12, distinct 1 ----
    # 12.  determine if single values in height/weight fall within BMI criteria
    step <- "12, distinct 1"

    # only do this if there's 1 distinct in either height or weight
    if (length(unique(h_subj_df$meas_m)) == 1 |
        length(unique(w_subj_df$meas_m)) == 1){
      # match and create BMIs

      # possible removal of height/weights by bmi
      # h = height, w = weight
      comb_df <- comb_df_orig <-
        merge(h_subj_df, w_subj_df, by = "age_days", all = T,
              suffixes = c(".h", ".w"))
      # remove ones that don't match
      comb_df <- comb_df[!(is.na(comb_df$id.h) | (is.na(comb_df$id.w))),]

      # must be matches
      if (nrow(comb_df) > 0){
        # calculate bmi
        comb_df$bmi <- comb_df$meas_m.w/((comb_df$meas_m.h/100)^2)

        # check if it's in a more common range, then have wider limits
        h_exc_btw <-
          check_between(comb_df$bmi, 16, 60) &
          (comb_df$meas_m.h < 88 | comb_df$meas_m.h > 231) &
          length(unique(h_subj_df$meas_m)) == 1
        w_exc_btw <-
          check_between(comb_df$bmi, 16, 60) &
          (comb_df$meas_m.w < 30 | comb_df$meas_m.w > 271) &
          length(unique(w_subj_df$meas_m)) == 1

        # if there's a bmi for the subject, but outside the limits, AND none
        # within the limits AND it's 1D, we exclude
        h_bmi_out <-
          all(!check_between(comb_df$bmi, 16, 60)) &
          length(unique(h_subj_df$meas_m)) == 1
        w_bmi_out <-
          all(!check_between(comb_df$bmi, 16, 60)) &
          length(unique(w_subj_df$meas_m)) == 1

        # remove based on above criteria
        rem_ids_ht <- comb_df$id.h[h_exc_btw | rep(h_bmi_out, nrow(comb_df))]
        rem_ids_wt <- comb_df$id.w[w_exc_btw | rep(w_bmi_out, nrow(comb_df))]

        # update and remove
        h_subj_keep[rem_ids_ht] <- "Implausible"
        h_subj_reason[rem_ids_ht] <- paste0("Implausible, Step ",step)
        w_subj_keep[rem_ids_wt] <- "Implausible"
        w_subj_reason[rem_ids_wt] <- paste0("Implausible, Step ",step)

        h_subj_df <- h_subj_df[!h_subj_df$id %in% rem_ids_ht,]
        w_subj_df <- w_subj_df[!w_subj_df$id %in% rem_ids_wt,]
      }

      # combine again, then check if still 1D
      # h = height, w = weight
      comb_df <- comb_df_orig <-
        merge(h_subj_df, w_subj_df, by = "age_days", all = T,
              suffixes = c(".h", ".w"))
      # remove ones that don't match
      comb_df <- comb_df[!(is.na(comb_df$id.h) | (is.na(comb_df$id.w))),]

      # no bmis available -- no matches
      if (nrow(comb_df) == 0){
        exc_ht <-
          !check_between(h_subj_df$meas_m, 139, 206) &
          length(unique(h_subj_df$meas_m)) == 1
        exc_wt <-
          !check_between(w_subj_df$meas_m, 40, 225) &
          length(unique(w_subj_df$meas_m)) == 1

        # remove based on above criteria
        rem_ids_ht <- h_subj_df$id[exc_ht]
        rem_ids_wt <- w_subj_df$id[exc_wt]

        # update and remove
        h_subj_keep[rem_ids_ht] <- "Implausible"
        h_subj_reason[rem_ids_ht] <- paste0("Implausible, Step ",step)
        w_subj_keep[rem_ids_wt] <- "Implausible"
        w_subj_reason[rem_ids_wt] <- paste0("Implausible, Step ",step)

        h_subj_df <- h_subj_df[!h_subj_df$id %in% rem_ids_ht,]
        w_subj_df <- w_subj_df[!w_subj_df$id %in% rem_ids_wt,]
      }
    }

    # 13h, H error load ----
    # 13h. compute error load -- whether there are too many errors and all
    # should be excluded
    step <- "13h, H error load"

    # no need to do this if everything is already excluded
    if (nrow(h_subj_df) > 0){
      # compute errors -- without sde
      num_err <- sum(!grepl("same day", h_subj_reason) & h_subj_reason != "")
      # compute include -- without sde
      num_inc <- sum(!grepl("same day", h_subj_reason) & h_subj_reason == "")

      # if there are greater than 33.33% errors for a subject, fill everything
      # with errors
      if (num_err/num_inc > 1/3){
        criteria <- h_subj_reason == ""

        # remove based on above criteria
        rem_ids <- names(h_subj_reason)[criteria]

        # update and remove
        h_subj_keep[rem_ids] <- "Implausible"
        h_subj_reason[rem_ids] <- paste0("Implausible, Step ",step)

        # no need to update h, since we're done
      }
    }

    # 13w, W error load ----
    # 13w. compute error load -- whether there are too many errors and all
    # should be excluded
    step <- "13w, W error load"

    # no need to do this if everything is already excluded
    if (nrow(w_subj_df) > 0){
      # compute errors -- without sde
      num_err <- sum(!grepl("same day", w_subj_reason) & w_subj_reason != "")
      # compute include -- without sde
      num_inc <- sum(!grepl("same day", w_subj_reason) & w_subj_reason == "")

      # if there are greater than 33.33% errors for a subject, fill everything
      # with errors
      if (num_err/num_inc > 1/3){
        criteria <- w_subj_reason == ""

        # remove based on above criteria
        rem_ids <- names(w_subj_reason)[criteria]

        # update and remove
        w_subj_keep[rem_ids] <- "Implausible"
        w_subj_reason[rem_ids] <- paste0("Implausible, Step ",step)

        # no need to update h, since we're done
      }
    }

    # add to output ----

    # add results to full dataframe
    if (nrow(h_df) > 0){
      df[names(h_subj_keep), "result"] <- h_subj_keep
      df[names(h_subj_reason), "reason"] <- h_subj_reason
      df[names(h_subj_mean_sde), "mean_sde"] <- h_subj_mean_sde
    }
    if (nrow(w_df) > 0){
      df[names(w_subj_keep), "result"] <- w_subj_keep
      df[names(w_subj_reason), "reason"] <- w_subj_reason
      df[names(w_subj_mean_sde), "mean_sde"] <- w_subj_mean_sde
    }

    # complete ----

  }

  return(df)
}
