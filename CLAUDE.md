# CLAUDE.md — gc-github-latest (growthcleanr)

**Last updated:** 2026-04-19 (R-vs-R Apr/Jan Comparison Session 10 (Support functions — `ewma()`, `ewma_cache_*`, `as_matrix_delta()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `get_dop()`, `calc_and_recenter_z_scores()`, `.child_valid()`, `.child_exc()`, `sd_median()`, `read_anthro()`): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of support functions defined outside `cleanchild()` in reference (`pediatric_support.R` + `Infants_Main.R`) against current `child_clean.R`. **1 finding opened and closed (AJ16 — Bug fix).** **AJ16** Bug fix closed — `calc_and_recenter_z_scores()` CDC boundary `>= 4` (reference `Infants_Main.R:2485`) vs `> 5` (current `child_clean.R:2383`): reference helper used `ageyears >= 4` as the condition for pure CDC, so for 4–5 year olds the helper applied pure CDC while the main z-score pipeline (which uses `ageyears > 5`) applied the 2–5y smooth blend — an internal inconsistency in the reference. Current correctly uses `> 5` for both, matching the main pipeline. Current also adds NA fallbacks for the smooth zone (`smooth_val & is.na(cn.orig_cdc) → cn.orig_who` and vice versa) absent in reference. `Infants_Main.R:2485` vs `child_clean.R:2383`. Pitfall: **Boundary changes**. Items NOT flagged: `valid()` → `.child_valid()` (known rename); exclusion code string updates in `.child_valid()` (known — `Exclude-C-Temp-Same-Day`/`Exclude-C-Extraneous`/`Exclude-C-CF` vs reference's `Exclude-Temporary-Extraneous-Same-Day`/`Exclude-Extraneous`/`Exclude-Carried-Forward`); `na_as_false()` removal — only called from legacy `pediatric_clean.R`, never from `Infants_Main.R` (infant path); `temporary_extraneous()` removal — only for non-infant legacy path; `swap_parameters()` removal — only for non-infant legacy path; orphan `valid()` at `Infants_Main.R:4991` — exact duplicate of `pediatric_support.R:valid()`, R loaded whichever came last (both identical); `as_matrix_delta()` — identical logic in both; `ewma()` window 25 → 15 (AJ12), return data.frame → named list (Session 7b), O(n²) matrix → O(n) arithmetic (Session 7b) — all known intentional; `ewma_cache_init()`/`ewma_cache_update()` — new in current, no reference equivalent (known intentional performance feature); `sd_median()` — identical logic in both; `get_dop()` — identical logic in both; `identify_temp_sde()`: `id` → `internal_id` tiebreaker (known), rounding removal (per procedure), direct integer return (Session 8 D-b/D-c) — all known; `calc_otl_evil_twins()`: `oob` → `otl` column rename (known), rounding removal (per procedure) — known; `calc_and_recenter_z_scores()` optional `measurement.to.z`/`measurement.to.z_who` params (known intentional, `cf_detail` feature), dead `setkey()` removal (Session 9 D-a) — known; `read_anthro()`: `prelim_infants` removal — infant path always used `prelim_infants = TRUE`, which loaded `growthfile_who.csv.gz` — current always loads both WHO and CDC ext files directly; closure logic for infant path is identical to reference's `prelim_infants = TRUE` path. Tests unchanged 63/48/28/41/13; no tests re-run since no code change. Adult tests not re-run. **Session 10 status: 1 finding opened and closed, no code change. R-vs-R Apr/Jan Comparison series complete — all 10 sessions closed.** R-vs-R Apr/Jan Comparison Session 9 (Child Step 19 Pairs/Singles + Step 21 Error Load + Step 22 Output): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Steps 19+21+22 (~lines 4725–4960) against current `child_clean.R` Steps 19+21+22 (~lines 4735–4966). **1 finding opened and closed (AJ15 — Bug fix).** **AJ15** Bug fix closed — Step 21 `non_error_codes` omits `"Not cleaned"` in reference: reference lists `"Missing"` as a non-error code but NOT `"Not cleaned"`, so HC rows with agedays > 3 × 365.25 (assigned `"Not cleaned"` in reference preprocessing) are counted as `n_errors` per `(subjid, HEADCM)` group, inflating the error-load ratio and potentially triggering error-load escalation for subjects with many HC > 3y measurements. Current correctly adds `Exclude-Not-Cleaned` to `non_error_codes`, excluding these out-of-scope rows from both numerator and denominator — consistent with the existing intent of excluding `"Missing"` rows. `Infants_Main.R:4873–4884` vs `child_clean.R:4894–4898`. Pitfall: **Factor levels / exclusion codes**. Items NOT flagged — Step 19: sort `id`→`internal_id` (known); `valid()`→`.child_valid()` (known); keyed DOP lookup + `nomatch = NULL` (F109 walkthrough Session 13); `isTRUE()` pair/single defensive wrapping; `1:nrow` → `seq_len` (cosmetic); dead `abs_tbd.sd`/`abs_ctbd.sd`/`med_dop`/`med_cdop` removed (dead in reference too); rounding removed (per procedure); exclusion code rename (known); `.SDcols` trimmed (F113 walkthrough Session 13). Step 21: `valid_set` removal (F123 walkthrough Session 14); error.load.threshold parameter replaces `.4` (known intentional — documented in procedure); rounding removed (per procedure); old CF rescue codes removed from non_error_codes (rescued CFs are now Include rows). Step 22: dead `checkpoint_data` merge block removed (logically dead in reference); `cf_rescued` added to return_cols (known CF rescue scheme); `ewma1_it1.*` loop removed (known debug removal); `cf_detail` columns added (known CF rescue scheme). Tests unchanged 63/48/28/41/13; no tests re-run since no code change. Adult tests not re-run. **Session 9 status: 1 finding opened and closed, no code change.** R-vs-R Apr/Jan Comparison Session 8 (Child Step 17 Height/HC Velocity): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Step 17 (~lines 4210–4715) against current `child_clean.R` Step 17 (~lines 4115–4733). **1 finding opened and closed (AJ14 — Bug fix).** **AJ14** Bug fix closed — reference HC velocity section is effectively a no-op: for-loop condition `!is.na(df$whoagegrp.hc)` always FALSE because `whoagegrp.hc` is never set (reference sets `whoagegrp.ht` instead), making the loop body never execute; additionally the loop body references `sub_m_who_ht_vel` (HT variable, undefined in HEADCM context) and `.ht`-suffix column names; result: HC pairs only receive the fixed `mindiff = -1.5` fallback with no upper bound check. Current correctly implements HC velocity via pre-merged `.hc`-suffix columns and `fcase` lookup. `Infants_Main.R:4452–4523` vs `child_clean.R:4557–4609`. Pitfalls: **NA / empty-set handling** + **Parameter scope**. Items NOT flagged: vectorized pre-filter over all valid rows (safe efficiency optimization); WHO merge with `.hc`-suffix column copy (intentional HC bug fix); Tanner batch pre-computation (equivalent); dplyr→`shift()`; `id`→`internal_id`, `valid()`→`.child_valid()`, exclusion code rename, `cat()`→`message()`, rounding removal (all known per procedure); F94/F95/F97/F101/F102 (known from walkthrough Session 12); AJ12 `ewma_window`; dead `id_all`/`iter_count` removed (Session 12 F-bycatch); debug block removed (Session 12); HC 24-month endpoint boundary not flagged (reference HC entirely non-functional, no valid baseline exists). Tests unchanged 63/48/28/41/13; no tests re-run since no code change. Adult tests not re-run. **Session 8 status: 1 finding opened and closed, no code change.** Next session: Session 9 — Child Step 19 (Pairs/Singles) + Step 21 (Error Load) + Step 22 (Output). R-vs-R Apr/Jan Comparison Session 7 (Child Step 15 EWMA2 Moderate + Child Step 16 Birth HT/HC): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Steps 15+16 (~lines 3793–4201) against current `child_clean.R` Steps 15+16 (~lines 3673–4106). **0 findings opened.** Items NOT flagged: tbc_range pre-filter before the Step 15 while loop (range ≤ 1 groups cannot pass `|dewma| > 1`; safe efficiency optimization); DOP snapshot timing (`dop_snap` keyed at each iteration start vs. reference's live per-closure `data.df` scan — current is strictly more correct, eliminating implicit within-iteration ordering dependency); `include.temporary.extraneous = TRUE` removed from `include_counts` (Session 11 F74, confirmed intentional — temp SDEs resolved in Step 13 Phase B4 before Step 15 runs); sp_key-level Step 16 birth filter (Session 11 F85, confirmed intentional behavior-neutral tightening — avoids wasted iteration for subjects with birth HT but no birth HC); EWMA caching (`ewma_cache_init`/`ewma_cache_update`/`ewma2_caches`/`ewma2b_caches` — known intentional performance feature); `cols_to_drop_15_16` drops p_plus/p_minus/tbc.p_plus/tbc.p_minus/first_meas after Step 16 (Session 11 F81); all rounding-tolerance removals (per procedure, do NOT flag); exclusion code rename to `Exclude-C-Traj`; `id`→`internal_id`; `valid()`→`.child_valid()`; `cat()`→`message()`; additional reference-closure arguments to `calc_and_recenter_z_scores()`; `compare_df[!is.na(tbc.sd)]` NA guard after DOP snapshot lookup (defensive equivalent); `rm()` calls after loops; dplyr/case_when→data.table. Tests unchanged 63/48/28/41/13; no new test run since no code change. Adult tests not re-run. **Session 7 status: 0 findings, no code change.** Next session: Session 8 — Child Step 17 (Height/HC Velocity). R-vs-R Apr/Jan Comparison Session 6 (Child Step 13 Final SDE): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Step 13 (~lines 3434–3792) against current `child_clean.R` Step 13 (~lines 3375–3671). **1 finding opened and closed with no code change.** **AJ13** Bug fix closed — Phase B3 EWMA Inf guards: reference has no `is.infinite → NA` conversion for `ewma_fill` or `spa_ewma`; in the edge case where a `(subjid, param, agedays)` group contains only temp SDE rows (all stable Include rows excluded by earlier steps such as BIV/EWMA1), `max(ewma.all[!was_temp_sde], na.rm = TRUE)` returns `-Inf`, which propagates to `spa_ewma = -Inf`, `absdewma = Inf`, `min_absdewma > 1` — erroneously marking all rows SDE-All-Extreme (Extraneous). Current adds explicit `is.infinite → NA_real_` guards for both variables, preventing the erroneous marking and correctly keeping the best-tiebreaker row as Include. `Infants_Main.R:3722–3730` vs `child_clean.R:3593–3607`. Pitfall: **NA / empty-set handling**. Items NOT flagged: debug code (fwrite/stop) → warning() cleanup; dplyr → data.table refactor; all code/exclusion-code renames; double-rounding removal (per procedure); `ewma_window` (AJ12 — known intentional); `median_tbc_gt2` dead variable in reference removed; Phase B4 merge guard added defensively; `suppressWarnings` added to `max()` calls (cleanup, also prevents -Inf propagation that AJ13 addresses explicitly); `!is.na(min_absdiff_rel_to_median)` guard in Phase B2 SDE-All-Extreme (equivalent to reference's row-level NA check). Tests unchanged 63/48/28/41/13; no new test run since no code change. Adult tests not re-run. **Session 6 status: 1 finding opened and closed, no code change.** Next session: Session 7 — Child Step 15 (EWMA2 Moderate) + Step 16 (Birth HT/HC). R-vs-R Apr/Jan Comparison Session 5 (Child Step 11 EWMA1): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Step 11 (~lines 3257–3434) against current `child_clean.R` Step 11 (~lines 3206–3373). **1 finding opened and closed with no code change.** **AJ12** Intentional (other) closed — EWMA window default 25 → 15: reference `ewma()` has `window = 25` as default and Step 11 passes no explicit window argument; current `ewma()` has `window = 15` as default and all call sites (Steps 11, 13, 15/16, 17 and `ewma_cache_init` at Steps 15/16) pass `window = ewma_window` consistently. Window=15 confirmed intentional by Carrie; added to procedure's "Known intentional changes" list. `Infants_Main.R:2013–2015` vs `child_clean.R:1695`. Pitfall: **Boundary changes**. Items NOT flagged: `|tbc.sd| > 3.5` pre-filter (logically equivalent efficiency optimization — exclusion criteria already require `|tbc.sd| > 3.5`); vectorized `diff`/`pmax` exp_vals computation (equivalent to reference's `sapply`/data.frame approach); EWMA caching (`cache_env`) and ctbc.sd shortcut — known intentional; `ewma1_it1.*` debug capture block removed — confirmed dead code per walkthrough Session 6; rounding removal in pot_excl criteria and worst-row sort — per procedure do NOT flag; `.SDcols` changes (`index`/`sex`/`nnte` removed, `param` added) — none used in closure / known renames; `subj_with_sde` pre-computed once (same in both; final 11c refresh catches any residual drift); all function/code renames (`temporary_extraneous_infants` → `identify_temp_sde`, exclusion code renames, `valid()` → `.child_valid()`, `cat()` → `message()`) — known intentional. Tests unchanged 63/48/28/41/13; no new test run since no code change. Adult tests not re-run. **Session 5 status: 1 finding opened and closed, no code change.** Next session: Session 6 — Child Step 13 (Final SDE resolution, multi-phase B1/B2/B3). R-vs-R Apr/Jan Comparison Session 4 (Child Step 7 BIV + Child Step 9 Evil Twins): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference `Infants_Main.R` Steps 7 and 9 against current `child_clean.R`. **2 findings opened, both closed with no code change.** **AJ10** Intentional (other) closed — absolute BIV weight age boundary: reference uses `agedays == 0` / `agedays != 0` for the 0.2 kg / 1 kg floor; current uses `agedays <= 365` / `agedays > 365` (extends 0.2 kg floor to entire first year, fixing exclusion of legitimate preterm weights 0.7–1.0 kg). Intentional bug fix documented in CLAUDE.md pipeline history (2026-03-16); not listed in procedure's "Known intentional changes" but confirmed intentional. `Infants_Main.R:3073–3075` vs `child_clean.R:3007–3010`. Pitfall: **Boundary changes**. **AJ11** Bug fix closed — Step 9 early-exit check: current adds `na.rm = TRUE` to `any(start_df$otl)`, preventing `if (NA)` error when `otl` contains NAs (from rows with missing tbc.sd) and no TRUE values. Reference's `any(start_df$oob)` without `na.rm` would error in that edge case. `Infants_Main.R:3187` vs `child_clean.R:3123`. Pitfall: **NA / empty-set handling**. Items NOT flagged: BIV code unification (`Absolute-BIV` + `Standardized-BIV` → `Exclude-C-BIV`) — known code rename; `biv.z.*` parametrization — known intentional (in procedure list); `v == 0` explicit tagging removed at `Infants_Main.R:3154` — dead code (v==0 → NaN at wrapper → Exclude-Missing); `ageyears` column drop after Step 7 — cleanup; `oob` → `otl` rename — known; rounding-tolerance removal in `calc_otl_evil_twins` — per procedure, do NOT flag; `paste0/table` → `.N by-group` for `not_single_pairs` — equivalent, same paste-collision audit note as Sessions 1 and 3; single data.table closure → per-group for-loop in Step 9 — equivalent (Evil Twins groups are fully independent; OTL detection gates same_sp_next/same_sp_prev, so cross-group interactions are impossible). Tests unchanged 63/48/28/41/13; no new test run since no code change. Adult tests not re-run. **Session 4 status: 2 findings opened and closed, no code change.** Next session: Session 5 — Child Step 11 (EWMA1, complex closure and EWMA caching). R-vs-R Apr/Jan Comparison Session 3 (Child Step 6 CF detection + rescue): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference Jan 2026 Step 6 (`__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` ~lines 2677–3042) against current `R/child_clean.R` Step 6 (~lines 2673–2989). **0 findings opened.** Major intentional change (CF rescue scheme redesign — procedure's largest intentional surface) confirmed and logged briefly, not analyzed further per procedure. All non-intentional diffs are housekeeping / cleanup with logic-equivalent behavior, logged under "Items NOT flagged (audit trail)" in the walkthrough note: reference `cf_string_length` dead column (`Infants_Main.R:2960` created, never read downstream, dropped at `:2981`) removed from current; reference drops never-created `is_eligible_include` in cleanup list (silent no-op on non-existent column) vs current's clean drop list; reference `wholehalfimp` HT/HC uses redundant `%% 1 | %% 0.5` OR (`Infants_Main.R:2869–2870`, `:2874–2875`) vs current's `%% 0.5` alone (`child_clean.R:2782, :2784`) — whole integers satisfy `%% 0.5 == 0` too so OR is redundant; reference leaves Step 6 working columns (`v.orig`, `wholehalfimp`, `seq_win`, `cs`, `absdiff`, `ageday_has_include`) on `data.df` past Step 6 (grep-confirmed no downstream readers — `seq_win` at `:1116` and `absdiff_rel_to_median` at Step 13 are independent scopes) vs current's explicit drop at `child_clean.R:2984–2989`; reference's `not_single` via `paste0(subjid, "_", param) %in% names(tmp)[tmp > 1]` (`Infants_Main.R:2694–2695`) has theoretical underscore-concat collision risk impossible in practice given fixed `WEIGHTKG/HEIGHTCM/HEADCM` param vocabulary (isomorphic to Session 1's AJ4 — not re-flagged, kept as audit-trail note that same latent pattern was reviewed). **Known intentional changes** (logged briefly, not analyzed): CF rescue scheme redesign (lookup-table-based `cf_rescue = c("standard", "none", "all")` + new `cf_rescued` output column + optional `cf_detail` → `cf_status`/`cf_deltaZ`; rescued CFs become `Include` rather than getting different exclusion codes; reference's four rescue codes `Exclude-1-CF-deltaZ-<0.05` / `<0.1-wholehalfimp` / teen-2-plus variants gone; age / interval / param / rounding-specific thresholds 0.05 / 0.20 / 0.40 / NR replace fixed 0.05 / 0.1; reference's age / sex gating (teens only for multi-CF rescue) handled via lookup's age-bin axis instead; `cf_rescue = "all"` rescues every detected CF including shared-SPA CFs — Session 4a change); `include.carryforward` deprecated (no longer wraps Step 6; remains as `.child_valid()` argument only); exclusion code rename (`Exclude-Carried-Forward` → `Exclude-C-CF`, `Exclude-SDE-Identical` → `Exclude-C-Identical`, `Exclude-Temporary-Extraneous-Same-Day` → `Exclude-C-Temp-Same-Day`); `id` → `internal_id` in sort keys at `child_clean.R:2711` (CF detection sort) and `:2955` (SDE-Identical restore); `valid()` → `.child_valid()`; `temporary_extraneous_infants()` → `identify_temp_sde()` at temp-SDE re-run call (helper walked in Sessions 2 + 8, not re-walked); `cat()` → `message()` for progress messages; `not_single` rewrite from `table(paste0(...))` + `%in%` to data.table by-group `.N > 1` + right-join; `absdiff` rounding-tolerance (`janitor::round_half_up` double-rounding) removed per procedure; large commented-out dplyr/map_lgl and legacy data.table blocks at `Infants_Main.R:2762–2808` removed; `### CP DROP IN` / `### CP UP` Stata-era markers and commented-out alternate wholehalfimp formulations at `:2838–2850` removed; stale-changelog comments (`CP TOGGLED TO FALSE`, `Removed nnte filter`, `Fix:`, `Fix CF rescue logic`) removed; `#### CP Untested, replicate of above code in data.table()` orphan stub removed. Tests unchanged 63/48/28/41/13; no new test run since no code change. Adult tests not re-run — no adult files in scope. **Session 3 status: 0 findings, session closed with no code change.** Next session: Session 4 — Child Step 7 (BIV, known intentional change: eight `biv.z.*` per-cell parameters replacing single `sd.extreme` / `z.extreme`) + Child Step 9 (Evil Twins, low–medium complexity). R-vs-R Apr/Jan Comparison Session 2 (Child Step 5 Temp SDE + Early Child Step 13 SDE-Identicals + `identify_temp_sde()` helper): per `R-child-AprJanComparison-procedure.md`. Pure R-vs-R diff of reference Jan 2026 (`__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` + `pediatric_support.R`) against current `R/child_clean.R`. 2 findings batched (`R-child-AprJanComparison-findings.md`), both closed with no code change needed. **AJ8** Bug fix closed — Early Step 13 `keep_id` empty-group handling: current guards with explicit `if (length(incl_ids) == 0L) NA_integer_` (`child_clean.R:2638–2643`); reference's `ifelse(..., min(id[exclude=="Include"], na.rm=TRUE), max(...))` on all-excluded `(subjid, param, agedays, v)` groups returns `Inf`/`-Inf` with an R warning "no non-missing arguments to min/max" (`Infants_Main.R:2634–2636`). Behavior-equivalent downstream — in an empty-Include group, `has_dup = (n_same_value > 1)` is FALSE (because `n_same_value = sum(exclude == "Include") = 0`), so the filter `has_dup & exclude == "Include" & id != keep_id` cannot fire regardless of `keep_id`. Reference does not produce wrong results but leaks warnings to users on datasets with all-excluded same-day groups. Logic pitfall category: **NA / empty-set handling** (`min()`/`max()` on empty subset → `Inf`/`-Inf` with warning). **AJ9** Bug fix closed — reference `valid()` regex latent bug: reference assigns bare `"Missing"` and `"Not cleaned"` codes (no `Exclude-` prefix) at `Infants_Main.R:1207` (init from is.na(v) | agedays < 0), `:1212` (HC > 3y), `:1387` (post-recenter NA-tbc.sd safety), `:1391` (HC ≥ 5y). Reference's `!grepl("^Exclude", exclude)` at `pediatric_support.R:21` returns TRUE for both (neither starts with "Exclude"), silently including these rows in `valid.rows` at every step calling `valid()`. Most likely trigger in practice: a subject with two HC measurements in the 3–5y window (both `"Not cleaned"`, both carrying real `tbc.sd` from the recentering pass before the `:1387` safety mark overwrites NA tbc.sd) where Step 5's `extraneous.this.day := (.N > 1)` silently re-labels one as `"Exclude-Temporary-Extraneous-Same-Day"` (the Missing/Not-cleaned row with NA `absdmedian.spz` sorts last, so loser assignment depends on upstream row order). Current's 2026-04-14/16 exclusion-code rename (`Missing` → `Exclude-Missing`, `Not cleaned` → `Exclude-Not-Cleaned`) fixes this as a side effect — regex at `child_clean.R:5020` is identical, but code strings now start with `Exclude-` so `^Exclude` correctly returns FALSE. Logic pitfall categories: **Factor levels / exclusion codes** + **`.child_valid()` flag scope**. Rename was framed as cosmetic in the procedure's "Known intentional changes" section; AJ9 documents the bug-fix dimension explicitly on the findings log. **Session 2 status: both findings closed, no code change.** No new test run since no code changes; baseline unchanged 63/48/28/41/13. Adult tests not re-run — no adult files in scope. Next session: Session 3 — Child Step 6 (CF). Major intentional change (new CF rescue scheme with `cf_rescue` parameter + lookup thresholds); expect narrower comparison. R-vs-R Apr/Jan Comparison Session 1 (preprocessing — z-scores, GA correction Child Step 2b, recentering): per `R-child-AprJanComparison-procedure.md`; pure R-vs-R diff of `__Pipeline/current-infant-R-growthcleanr/growthcleanr/R/Infants_Main.R` (Jan 2026 reference) against current `R/child_clean.R`. 7 findings batched (`R-child-AprJanComparison-findings.md`): **AJ1** Bug fix open — `"HEIGHIN"` typo in reference's param-validation list silently rendered HEIGHTIN data Missing (current has correct `"HEIGHTIN"`); **AJ2** Bug fix open — `!is.na(sd.orig)` guard added to potcorr_wt computation in current; reference would propagate NA into potcorr_wt if first-WT row's measurement was Missing; **AJ3** Intentional (other) open — intwt low-weight floor bounds changed (250–560g→570g vs 100–499g→500g) tied to Fenton table swap (CSD method); current's lower bound `>= 100` and exclusive `< 500` upper documented in inline comment; **AJ4** Bug fix open — reference's `tmp[subjid %in% names(table(subjid) > 1),]` was effectively no-op (`names()` of a logical vector returns ALL names regardless of TRUE/FALSE); current uses `names(which(...))`. Practical impact nil (single-row subjects fail the downstream abssumdiff comparison naturally), but semantics now correct; **AJ5** Bug fix open — reference's `c(orig_colnames, id)` in column-keep filter used the function-parameter vector `id` rather than the string `"id"`; current uses `"id"`. Benign (orig_colnames already contained `"id"`) but cleaner; **AJ6** Bug fix FIXED — at exactly agedays=1461 (4 years exactly) for potcorr subjects, reference assigned `sd.corr := sd.orig` (the post-blend `>= 4` override step won), but current was assigning `sd.corr := sd.c` (smooth-blend branch's inclusive `<= 4` won). Resolved per Carrie: changed second `fcase()` branch upper bound from `ageyears_2b <= 4` to `ageyears_2b < 4` so the smooth-blend window is half-open `[2, 4)` and age=4 falls to default `sd.orig`. Affects only integer agedays=1461 (a single day per potcorr subject); ages 1460 and 1462 unchanged. Inline comment added naming the half-open convention; **AJ7** Intentional (other) open — HC ≥ 5y exclusion code recoded `'Missing'` → `'Exclude-Not-Cleaned'` (semantic change beyond the literal `Missing`→`Exclude-Missing` rename); current's comment cites cross-policy consistency with the HC > 3y `Exclude-Not-Cleaned` rule. Tests unchanged 63/48/28/41/13 pre- and post-fix. **Session 1 status: all 7 findings closed.** AJ6 fixed in current code; AJ1/AJ2/AJ4/AJ5 confirmed correct in current as-is (bug fixes already present); AJ3/AJ7 confirmed intentional, no change needed (per Carrie). Items NOT flagged (verified equivalent or already in procedure's known-changes list): WHO/CDC merge → closure refactor; sd.c assignment via `fcase()` collapse; sd.c_temp save/restore eliminated; pc-subset working pattern; LENGTHCM relabel; Fenton LMS→CSD formula; `prelim_infants` removal; `cat()`→`message()`; round_stata removal (per procedure: do NOT flag); cf_rescued init; recentering NHANES/derive paths removal; `id`→`internal_id` rename throughout; outer batching restructure (out of scope for preprocessing diff). Next session: Session 2 — Child Step 5 (Temp SDE) + Early Step 13 (SDE-Identicals). Session 14 Child Step 21 walkthrough: detailed code/comments/narrative/tests pass on Main Child Step 21 (Error Load escalation to all remaining Includes when per-(subjid, param) error ratio exceeds threshold) in `child_clean.R` (~lines 4860–4901 pre-edit). Adult has its own `eval_error_load()` in `adult_support.R` (Adult Step 14) with permissiveness-driven thresholds (0.41 / 0.29), param-specific error-code enumeration, and a `denom >= 3` sample-size floor instead of a `mincount` error-count floor — out of scope. All 21 checklist items applied (wrapper-only items 14–19 n/a). No behavior changes; tests unchanged 63/48/28/41/13. Fix-now F123–F132, all comment / narrative / small-cleanup: **F123** dead `valid_set <- rep(TRUE, nrow(data.df))` at `:4871` removed — unlike Child Steps 7/9/11/17/19 where `valid_set` is computed via `.child_valid()`, Step 21's was unconditionally TRUE since the initial v3.0.0 commit (`27bcc0f`, git-history-verified), making the two `valid_set &` / `valid_set,` uses at `:4881` / `:4897` no-ops; Step 21 is order-independent and correctly counts all rows (error-vs-non-error classification happens inside the counting expression via `non_error_codes`). **F124** F58/F75 stale changelog `# CF rescue codes removed — rescued CFs are now "Include" (stored in cf_rescued column)` at `:4874` rewritten as current-state naming each code's role: Identical/Extraneous (SDE housekeeping resolved in Child Step 13), CF (detected unrescued CF, data-entry artifact; rescued CFs are Include via cf_rescued and count as Includes in the denominator), Missing (no measurement), Not-Cleaned (HC > 3 × 365.25 days, out of scope). **F125** redundant + inaccurate short-form `# Denominator is errors + includes (excludes SDE/CF/Missing)` at `:4887` removed (outer comment at `:4868–4870` already explains the ratio structure, and the short form missed `Not-Cleaned`). **F126** roxygen `error.load.mincount` at `:211–212` rewritten — pre-edit "minimum count of exclusions on parameter before considering excluding all measurements" was too broad ("exclusions" lumped in SDE/CF/Missing/Not-Cleaned) and ambiguous ("on parameter" not "per (subject, param) group"); post-edit names the exact exclusion codes that do NOT count as real errors (Identical/Extraneous/CF/Missing/Not-Cleaned) and states per-(subject, param) scope. **F127** roxygen `error.load.threshold` at `:213–214` rewritten — pre-edit "percentage of excluded measurement count to included measurement count that must be exceeded before excluding all measurements of either parameter" was imprecise on formula (reads as errors : includes ratio but formula is `errors / (errors + includes)`), scope ("either parameter" suggests both params together but trigger is per-(subject, param)), and resulting code (not named); post-edit states the formula explicitly, names the per-(subject, param) scope, notes SDE/CF/Missing/Not-Cleaned exclusion from both numerator and denominator, cross-references `error.load.mincount`, and names `Exclude-C-Too-Many-Errors`. **F128** narrative Code-location cell (`child-algorithm-reference.md:1527`) expanded to name `.child_exc()` (F99/F115 pattern — exclusion-code constructor not previously mentioned). **F129** narrative Overview at `:1531` expanded — pre-edit "excluding SDEs, CFs, and Missing from the denominator" missed `Not-Cleaned` and understated the rule ("from the denominator" — the codes are excluded from BOTH numerator and denominator); post-edit names all four categories and both sides of the ratio, plus adds "for that subject-param" to the exclusion scope sentence. **F130** narrative Logic step 5 at `:1539` — "all Include rows" clarified to "all Include rows within that (subjid, param) group" to match the per-group trigger scope. **F131** narrative "Variables created and dropped" at `:1559` rewritten — pre-edit listed `n_includes` as a working column, but only `err_ratio` and `n_errors` are actually assigned as columns via `c("err_ratio", "n_errors") := {...}`; `n_includes` and `denom` are local scratch variables inside the `{}` block and never become columns. Post-edit names the two actual columns, notes their per-group constancy and cleanup-at-end, and describes `n_includes` / `denom` as local scratch inside the expression. **F132** narrative Checklist finding 3 at `:1565` rewritten — pre-edit was internally inconsistent (heading "Denominator excludes SDEs/CFs" vs body "from both numerator and denominator") and missed Missing/Not-Cleaned; post-edit reads "Numerator and denominator exclude SDEs/CFs/Missing/Not-Cleaned" with matching body. Parameter table row audit (checklist item 6) confirmed: `error.load.mincount` (default 2) and `error.load.threshold` (default 0.5) defaults match across signature `:343–344`, `cleanchild()` internal signature `:2526–2527`, CLAUDE.md Configurable Parameters table, and narrative Configurable parameters table at `:449–450`. Step 21 does NOT call `.child_valid()`, which is correct by design — counting needs to see ALL rows including BIV / Evil-Twins / EWMA / velocity / pair / single exclusions (must contribute to `n_errors`); filtering via `.child_valid()` first would collapse `n_errors` and defeat the step. Man pages regenerated via `roxygenise()`. Deferred items: **none** — all identified items were fix-now. Discussion points resolved during walkthrough: `valid_set` dead since initial v3.0.0 commit — never had a meaningful filter; `.child_valid()` skip rationale (real-error codes are "not currently Include" so `.child_valid()` would exclude them from the count); child vs adult design differences (single non-param-specific error list + `mincount` floor vs param-specific error enumeration + `denom >= 3` floor) documented but not normalized. Session 13 Child Step 19 walkthrough: detailed code/comments/narrative/tests pass on Main Child Step 19 (Pairs and Singles Evaluation) in `child_clean.R` (~lines 4722–4833 pre-edit). Adult has its own 1D-single / 2D-ordered-pair / 2D-non-ordered-pair paths (Adult Steps 10H, 11Wa, 11Wa2, 13) that operate on raw measurements with permissiveness-driven limits instead of DOP z-score comparison, so this pass is child-scoped. All 21 checklist items applied (wrapper-only items 14–19 n/a). No behavior changes; tests unchanged 63/48/28/41/13. Fix-now F103–F122, all comment / narrative / small-cleanup. F103 consolidated the two stale sort comments at `:4729–4730` (`# order just for ease later` + the F68/F73/F80/F90 `# Add id for consistent SDE order` pattern) into a single current-state block describing the (subjid, param, agedays, internal_id) sort's downstream use (`df[1]` / `df[2]` indexing in the 2-row pair branch; `internal_id` as a deterministic tiebreaker against a rare same-ageday case that does not normally reach Step 19). F104 `# Compute valid_set AFTER sort, not before` at `:4733` (F58/F75 stale-changelog pattern) rewritten as current-state rationale for what `valid_set` gates (Include rows whose sp has 1 or 2 Includes, counting only Includes so sp's whose excluded measurements leave them with 1 or 2 Includes still flow through). F105 the 5-line Problem/Solution DOP-snapshot comment at `:4743–4747` rewritten as current-state rationale — freezes Include rows (with both `tbc.sd` and `ctbc.sd`) before the per-(subjid, param) closure runs, so the DOP lookup for the second-processed param cannot miss rows that the first-processed param just excluded. F106 `# save initial exclusions to keep track` at `:4752` rewritten to name both vectors (`ind_all` for index-keyed back-writes; `exclude_all` as the closure's return vector). F107 six Stata-era letter sub-labels (`19A` at `:4756`, `19D` at `:4763`, `19B/C` at `:4779`, `19E` at `:4785`, `19F/G` at `:4795`, `19H` at `:4811`) rewritten as descriptive current-state comments — they appeared out of execution order (A → D → B/C → E → F/G → H), did not map onto the narrative's numbered list, and `voi` (in 19D) was undefined terminology. F108 bundled with F107 — the one `voi` occurrence (grep-confirmed no other `voi` / `VOI` anywhere in `R/*.R`) is now gone. F109 dead `else { df[, comp_diff := rep(NA, nrow(df))] }` branch at `:4775–4777` — the partial keyed lookup `dop_snapshot[.(subjid, dop_param)]` without `nomatch` returned 1 NA row when no match, so `nrow(dop) > 0` was always TRUE and the else branch was unreachable (behavior was still correct because `median(NA) = NA` and `abs(NA - x) = NA`, so `comp_diff` ended up NA either way). Added `nomatch = NULL` to the keyed lookup so the else branch is reachable and the no-DOP case skips the per-row for-loop entirely. F110 `# Use id tie-breaker for deterministic selection` at `:4786` (F80/F84 pattern) rolled into the F107 19E rewrite; new comment names `internal_id` explicitly and documents direction ("lowest internal_id wins — that row is the one excluded"). F111 vague `# save the results` at `:4808` rewritten as current-state description of why `exclude_all` is refreshed from `df$exclude` after the pair rule. F112 stale changelog `# Previous code at line 4529 unconditionally overwrote exclude_all, losing 2-meas results` at `:4825–4826` (F58/F75 pattern; stale line number) rewritten as current-state rationale for why the single-exclusion write is index-keyed (to preserve any pair exclusion already recorded above). F113 unused `'id'` removed from `.SDcols` at `:4833` (Session 7 F77 pattern; post-edit: `'index', 'internal_id', 'agedays', 'subjid', 'param', 'tbc.sd', 'ctbc.sd', 'exclude'`). F114 cosmetic `>=365.25` → `>= 365.25` spacing fix at `:4800`. F115 narrative Code-location cell (`child-algorithm-reference.md:1464`) expanded to name `.child_valid()` and `.child_exc()` alongside `get_dop()` (F99 pattern). F116 narrative "lowest id as tiebreaker" at line 1479 → "lowest internal_id" with direction clarified (F80/F84 pattern). F117 narrative "Variables created and dropped" at line 1505 listed `dop_tbc` and `dop_ctbc` as working columns, but grep confirmed those are notational shorthand in the comp_diff formula, not actual variables in the code; rewritten to distinguish local scalars (`diff_tbc.sd`, `diff_ctbc.sd`, `diff_agedays` in the pair branch) from the closure-local column (`comp_diff`) and note explicitly that no columns are persisted to `data.df`. F118 narrative Checklist finding 6 "All 3 codes exist" at line 1514 — miscount (Step 19 produces exactly two codes); rewritten to name `Exclude-C-Pair` and `Exclude-C-Single` explicitly and point at `exclude.levels.peds`. F119 narrative Overview "After all cleaning steps" at line 1468 — imprecise (Step 21 runs after Step 19); changed to "After the prior individual-value cleaning steps". F120 narrative Scope / Operates-on at line 1460 — "1–2 remaining measurements" tightened to "1–2 remaining Include measurements" to match the code's Include-only sp_count. F121 narrative Single evaluation at line 1488 and Checklist finding 4 at line 1512 — clarified that singles use only `tbc.sd` (no `ctbc.sd` check) but have an implicit potcorr guard via the `comp_diff` cross-check against DOP (both sides of the comparison are on the uncorrected scale). F122 narrative DOP snapshot subsection and Pair evaluation numbered list expanded — added explicit mention that `ctbc.sd` is also snapshotted, named `diff_agedays` in the pair numbered list, described `comp_diff` directly rather than via ambiguous `dop_tbc` notation, and reworded step 5 of the pair rule to match code. Callers-line in the `get_dop()` support-function subsection (`child-algorithm-reference.md:114`) updated to match current line numbers (3845 for Step 15/16 pre-loop, 4769 for Step 19). Deferred items: none; all identified items were fix-now. (The pair `> 4` / `> 2.5` and single `> 3` / `> 5` thresholds are deliberately strict `>` by design, not gaps to close.) Discussion points resolved during walkthrough: partial keyed lookup default returns 1 NA row, making the pre-edit else branch dead; non-potcorr `ctbc.sd` equals `tbc.sd` (because `sd.corr := sd.orig` for non-potcorr subjects), so `|diff_ctbc.sd|` check is equivalent to `|diff_tbc.sd|` and the `| is.na(diff_ctbc.sd)` disjunct is a defensive guard rather than the mechanism that exempts non-potcorr subjects; single rule has no explicit `ctbc.sd` check but the `comp_diff` cross-check against the DOP's `tbc.sd` is itself uncorrected-scale vs uncorrected-scale, so for a well-behaved potcorr subject the extremeness cancels in the disagreement check. Session 12 Child Step 17 walkthrough: detailed code/comments/narrative/tests pass on Main Child Step 17 (Height/HC Velocity) in `child_clean.R` (~lines 4094–4706 pre-edit). Adult has no velocity/raw-diff step (uses raw-measurement-based limits); out of scope. All 20 checklist items applied plus new item 21 "Look for age gaps" added to `algorithm-walkthrough-procedure.md` during the session (gap-search surfaced F101 and F102 part 2, so added as a formal item for future walkthroughs; Carrie plans a focused retrospective over already-reviewed steps with this lens). Six behavior-narrowing fixes (all test counts unchanged 63/48/28/41/13): F94 Tanner NA-out filter now keys on `tanner.months < 30` (midpoint-based, same key as the merge) instead of `agedays/30.4375 < 30` (this-row-only) — per Carrie "this has been wrong for a long time"; pairs whose midpoint age ≥ 30 months but starting age < 30 months now correctly use Tanner (pre-filter `:4197`, inner loop `:4441`). F95 paired `< 365.25` / `> 365.25` formulas for mindiff / maxdiff (both Tanner and Tanner-WHO-fallback variants) changed to `< 365.25` / `>= 365.25` to cover the `d_agedays == 365.25` boundary; formulas are continuous at gap = 1 year so behavior-equivalent for integer agedays (pre-filter `:4221`/`:4225`, inner loop `:4452`/`:4456`). F97 2-row branch strict `>` → `>=` in all eight `abs(tbc.sd) > shift(abs(tbc.sd), ...)` comparisons at `:4637–4652`; previously, tied `abs(tbc.sd)` pairs with a violating diff flagged NEITHER row, so the violation escaped Step 17 entirely — now both flag and `order(-absval, internal_id)` picks lowest internal_id. F101 pre-filter's dead `whoagegrp.ht > 24 | whoagegrp.ht.lead > 24` OR branch at `:4231–4234` — per Carrie the intent was NA-out WHO "if the second height is >24 months" (boundary-crossing pair), not the dead-as-written capped-column check; rewritten in pre-filter as direct `(agedays + d_agedays) / 30.4375 > 24` on `whoagegrp.ht`, and added per-iteration NA-out of `who_mindiff_ht` / `who_maxdiff_ht` in the inner loop after the fcase (whoagegrp.ht is pre-merged pre-loop on static data, but d_agedays changes per iteration). Behavior change: boundary-crossing HT pairs now fall through to Tanner or default mindiff=-3 rather than using a WHO reference extrapolated past the 24-month boundary. F102 (collapse) `d_agedays < 20 → 1` + `20 ≤ gap < 46 → 1` collapsed into `d_agedays < 46 → 1`; `153 ≤ gap < 199 → 6` + `d_agedays >= 200 → 6` collapsed into `d_agedays >= 153 → 6` — latter closes the former `d_agedays == 199` integer gap (pre-filter `:4229–4234`, inner loop `:4464–4469` pre-edit). whoinc.age.ht NA init added at top of inner-loop HT block (`df[, whoinc.age.ht := NA_integer_]` at `:4477`) as defense-in-depth so any future boundary change that reintroduces a gap cannot leave a stale value across iterations (HC block already had equivalent NA init). Comment/narrative fix-now F86–F91, F93, F98–F100, plus F-bycatch dead `iter_count` / `id_all` removed. F86 stale `# Bug fix: was redundantly re-reading from disk each batch.` at `:4104` rewritten as current-state. F87 two `# Chris updated this to >= from ==` stale-changelog lines at `:4465`/`:4470` removed (sanity grep confirms no remnants in `R/*.R`). F88 stale `# Fix Min-diff tie-breaking` at `:4663` rewritten as current-state rationale. F89 stale line-number reference `lines 4030-4031` at `:4499` replaced with section-name pointer "already transformed in 17G above". F90 `# Add id for consistent SDE order` comments at four sort sites (`:4152`, `:4421`, `:4519`, `:4570`) rewritten as current-state `internal_id` wording (F68/F73/F80 pattern; Session 11 deferred `:4148` site closed here, plus three more). F91 3-line tiebreaker comment cluster at `:4683–4685` collapsed to 2 lines (claim + rationale) matching Session 11 F78. F93 narrative HC interval table rewritten — pre-edit listed `gap < 20 → 1`, `20 ≤ gap < 46 → 1`, `gap ≥ 200 → 6` for HC and framed as "same gap-based mapping as HT," but HC has NO 1-month interval (smallest is 2-month, 46–75 days) and NO `≥ 200` fallback; corrected to four intervals only with explicit note that pairs outside the reference band use only the decrease-side default (`mindiff = -1.5`) with no upper bound — consistent with HC's narrower physiologic range. F98 narrative "Configurable parameters in scope: None" corrected — `ewma_window` is in scope at `:4593`; added single-row table. F99 narrative Code-location cell expanded to name specific support functions (`ewma()` direct call not cache API, `.child_valid()`, `.child_exc()`) and distinguish Tanner (pre-loaded in `cleangrowth()` and passed in) from WHO HT/HC (loaded inside `cleanchild()`). F100 vestigial tautology comment at `:4694` from Session 11 dead-branch collapse rewritten as current-state iteration rationale. F-bycatch dead `id_all <- copy(df$id)` at `:4390` removed (never read; `ind_all` is the parallel copy actually used). F-bycatch dead `iter_count` variable at `:4415`/`:4418` removed. Two items deferred (both previously decided or low-value): F96 Step 17 EWMA caching — Step 17 uses `ewma()` directly rather than `ewma_cache_*`, and 2-row branch invokes `ewma()` even though tiebreaker keys on `abs(tbc.sd)` directly; Carrie explicitly deferred. Known Issue (already in CLAUDE.md): WHO velocity tables loaded inside `cleanchild()` rather than via `gc_preload_refs()`. Structural/stylistic defers: `dewma.all` computed but not consumed in Step 17 closure (kept for EWMA-trio parallelism); redundant pre-filter `pf <- pf[order(...)]`; `not_single` via `table(paste0(...))` rather than data.table `by` aggregation. Tests: 63/48/28/41/13, identical to baseline even though six of the fixes were behavior changes in principle (test data does not exercise the affected boundary conditions — consistent with being rare edge cases). Adult tests not re-run. Session 11 Child Steps 15/16 walkthrough: detailed code/comments/narrative/tests pass on Main Child Steps 15 (EWMA2 moderate) and 16 (Birth HT/HC) in `child_clean.R` (~lines 3660–4088 pre-edit). Adult has no z-score pipeline and a separate moderate-EWMA path (Adult Step 11Wb) — out of scope. All 20 checklist items applied; wrapper-only items 14–19 n/a. Fix-now F73–F85 plus two behavior-neutral tightenings requested by Carrie. F73 `# Add id for consistent SDE order` at `:3671` (moderate EWMA, not SDE — and `id` vs `internal_id` loose wording; Session 10 deferred to this walkthrough) rewritten to name the downstream prior/next-neighbor + first_meas sort dependency. F74 `# Include temp SDEs in count so they get p_plus/p_minus calculated` at `:3675` rewritten — the comment was misleading on two counts (temp SDEs are all resolved in Main Step 13 Phase B4, and the count doesn't gate p_plus/p_minus which is filtered by `sp_key & param`, not exclude); at Carrie's request, `include.temporary.extraneous = TRUE` flag also removed (a straggler should manifest as "too few Include rows," not get counted toward the n ≥ 3 threshold). F75 `# Fix first_meas logic` stale-changelog line at `:3705` removed; surrounding per-param comments consolidated into a single current-state block describing the WT (agedays > 0) vs HT/HC (non-birth only) definitions. F76 stale `valid()` reference at `:3758` → `.child_valid()`; in-loop first_meas recalc comments consolidated to match F75 structure. F77 `temporarily excluded before iteration starts` at `:3911` rewritten — birth HT/HC rows are not temporarily excluded; they are filtered out of `step15_filter` at `:3725`/`:3762` and never enter `df`. F78 three lines of tiebreaker comments at `:3915–3917` collapsed to two (claim + rationale). F79 stale `valid()` at `:3981` → `.child_valid()`. F80 `# Use id tie-breaker` at `:4061` (Step 16) → `internal_id` wording (cross-refs Step 15). F81 drop-columns comment at `:4083–4085` rewritten — was "after Step 15" (drop runs after Step 16) and "+/-5% rule (Step 15)" (p_plus/p_minus/tbc.p_plus/tbc.p_minus feed the addcrit perturbation check in BOTH Steps 15 and 16; first_meas is Step 15 only); local variable renamed `cols_to_drop_15` → `cols_to_drop_15_16`. F82 narrative Code-location cell (`child-algorithm-reference.md:1190`) — `ewma()` removed (not directly called; Steps 15/16 use only the incremental `ewma_cache_init()` / `ewma_cache_update()` pair), `get_dop()` added (Step 15 last-ext DOP lookup). F83 narrative "±5% rule" (line 1205) → "addcrit perturbation check" (±5% is WT-only; HT/HC uses ±1 cm; and it's not a specific rule but the check wrapping every rule). F84 narrative "lowest id as tiebreaker" (line 1285) → "lowest internal_id — the row with the lowest internal_id among tied candidates becomes the exclusion" (direction clarified). F85 narrative "DOP snapshot refreshed each iteration" (line 1289) scoped to Step 15 only (Step 16 has no DOP lookup). Second behavior-neutral tightening: Step 16 birth filter at `:3948–3987` switched from subject-level (`subjid %in% subj_with_birth`) to per-sp_key (`sp_key %in% sp_with_birth`). A subject with birth HT but no birth HC (or vice versa) was previously running `ewma_cache_init()` on the non-birth-carrying param for one wasted iteration (all four Step 16 rules require `agedays == 0`, so the non-birth group produced no candidates and dropped from `sp_to_process` on the next iteration). Exclusion set is identical; narrative at `child-algorithm-reference.md:1281` already described the tightened intent correctly, with a one-line clarification added for the asymmetry case. Deferred items (behavior-neutral, low-value): `df[, exp_vals := cache$exp_vals]` at `:3801` / `:4009` (cache-hit branch writes `exp_vals` to local `df` but doesn't read it; kept for symmetry with the cache-miss branch); `.SDcols` contains unused `sex` and `v` in both Step 15 (`:3926–3927`) and Step 16 (`:4070–4071`) closures (analogous to Session 7 Step 11 trim); redundant `order()` at `:3672` (data.df already keyed at end of Phase B4); structural `15a`/`15b`/`15c` + `16a`/`16b` subsection headers matching Session 7 Step 11 convention; `child_clean.R:4148` `# Add id for consistent SDE order` comment at the Step 17 pre-loop sort (same F68/F73 pattern; belongs to the Step 17 walkthrough); file-header `ewma()` / `ewma_cache_*` bullet at `:94` (F66-style expansion possible but current wording general enough). Tests unchanged: 63/48/28/41/13. Adult tests not re-run. Session 10 Child Step 13 walkthrough: detailed code/comments/narrative/tests pass on Main Child Step 13 (Final SDE Resolution) in `child_clean.R` (~lines 3361–3653 pre-edit). Adult algorithm has its own separate SDE resolution (Adult Steps 3, 9H, 10W) and is out of scope; this pass is child-scoped. Early Child Step 13 was already walked in Session 2 and is not re-walked. Inventory: Phase A uses `identify_temp_sde()` with `exclude_from_dop_ids = temp_sde_ids_step13` (the only non-NULL caller, confirmed Session 8); Phase B3 uses `ewma()`. All 20 checklist items applied — the wrapper-only items 14–19 are n/a for Step 13. Fix-now F66–F72, no behavioral changes: F66 file-header bullet for `identify_temp_sde()` at `child_clean.R:89` expanded from `(Step 5)` to `Child Steps 5, 6, 7, 9, 11, 13` to match the function roxygen (F50); F67 stale `# Filter ... before dplyr chain` comment at `:3412` rewritten — the block below uses only data.table; F68 `# Include id for deterministic SDE order` comment at `:3450` clarified to name `internal_id` explicitly (the `setkey` below uses `internal_id`, not user `id`); F69 Stata-era `"Ba" value` parenthetical at `:3554` removed (single-hit stale token, no current-reader decodable meaning) and replaced with a current-state description of how `maxdiff` drives the EWMA exponent; F70 narrative Phase B1 description in `child-algorithm-reference.md:535–542` rewritten — the pre-edit "catch identicals that emerged after Child Steps 5–11 removed intervening values" framing was false (Early Child Step 13 already resolves same-day same-value Include pairs, `v` is never modified downstream, every temp-SDE rerun preserves one-Include-per-SPA, and `cf_rescue = "all"` restored values sit on different days than their originators; Phase B1 is a no-op safety net), and the "Only Include rows are candidates" claim was only half accurate (whole-day check at lines 3433–3435 has no `exclude == "Include"` guard; only the partial check at `:3446` does); narrative now describes Phase B1 as defensive with explicit invariant chain. Phase B1 code left in place — removing it would be behavior-adjacent (edge cases breaking the invariant would silently go through Phase B2/B3 as Extraneous rather than Identical, losing a diagnostic distinction). F71 narrative "Variables created and dropped" section rewritten — pre-edit list incorrectly attributed `median.spz`/`median.dopz`/`absdmedian.spz`/`absdmedian.dopz` to `data.sde` (they live on `identify_temp_sde()`'s local copy and are discarded on return); rewrite enumerates the actual per-phase Main Step 13 working variables (B1 defensive; B2 one-day; B3 EWMA) with code-location pointers and documents the `keep_cols_sde` / `setdiff` drop mechanism at `:3651–3652`. F72 Checklist finding 1 wording "HC → HT in all three DOP median calculations" clarified to match the Phase B2 narrative's `WT↔HT, HC→HT` convention. Out-of-scope items not opened as new deferreds: Phase B1 code removal (behavior-adjacent — safer as defensive); Step 15/16 line 3667 `# Add id for consistent SDE order` comment (belongs to the Step 15/16 walkthrough); `suppressWarnings(min/max)` + `is.infinite()` pattern at `:3583–3589` / `:3596–3601` / `:3473–3477` (already in `CLAUDE.md → Known Issues → Robustness audit → Low priority cleanup`). Tests unchanged: 63/48/28/41/13. Adult tests not re-run (no adult files touched). Session 9 `calc_and_recenter_z_scores()` support-function pass + narrative fixes + new "Support Functions" section. Part 1: dedicated comment-only + small-cleanup walkthrough of the shared z-score + recentering helper in `child_clean.R` (~lines 2243–2311 pre-edit). Inventory confirmed only two callers, both in the Child Step 15/16 pre-loop in `cleanchild()` (`:3617` with `cn = "p_plus"` and `:3619` with `cn = "p_minus"`). Adult algorithm has no z-score pipeline and does not call this helper, so the pass is child-scoped and adult tests were not re-run. Fix-now items F56–F62: unified roxygen header replacing the plain-`#` pre-function block with formal `@param` / `@return` and a description naming both callers and the per-batch closure-reuse optimization (F56); stale "# for infants, use z and who" line removed (F57); three "Fix ..." / "Bug ..." stale-changelog comment blocks rewritten as current-state prose with accurate invariants ("who_val / smooth_val mutually exclusive" etc.) derived during the walkthrough (F58–F60); stale "line 773" line-number pointer to the main z-score pipeline dropped in favor of the more general "matches cleangrowth() exactly" wording already added in F58 (F61); awkward "now recenter -- already has the sd.median ..." phrasing rewritten as a single current-state sentence (F62). Additional cleanup: dead `setkey(df, subjid, param, agedays)` at line 2303 removed — the only downstream ops (`tbc.cn := ...` and `setnames`) do not need a key and the caller's merge is by `on = "index"` (D-a). Walkthrough-triggered narrative correction F-bycatch1: wrapper-narrative `Age blending — calc_and_recenter_z_scores() (CF rescue)` subsection at lines 399–408 was attributed to Child Step 6 / CF rescue, but the helper is called only from Child Step 15/16; heading and opening sentence rewritten (formula table unchanged). F-bycatch2: `R/child_clean.R:91` file-header bullet attributed the helper to "Steps 11 and 15" — rewritten to "Child Step 15/16 pre-loop". Mid-session user request: check `.child_valid()`, `.child_exc()`, `get_dop()` carefully before commit — F63: `.child_valid()` roxygen title said "for cleanbatch" (stale — current function is `cleanchild()`) and description omitted two of three include flags; replaced with full roxygen naming all three flags (`include.temporary.extraneous`, `include.extraneous`, `include.carryforward`) and the Step 5 / 13 / 6 origins of their corresponding exclusion codes, plus formal `@param` / `@return`. F64: `.child_exc()` had a plain-`#` comment block with no `@keywords internal` / `@noRd`; converted to formal roxygen; call-site count updated ~40 → ~50 to match current Grep result. F65: `get_dop()` had an F52-pattern orphan `#` block above floating `#' @keywords internal`; unified into a full roxygen header noting the scalar-only contract (the `if / else if / else` chain requires a scalar test; both callers pass `df$param[1]`). Part 2: stale `child-gc-narrative-2026-04-13.md` → `child-algorithm-reference.md` rename applied across current-state docs (wrapper narrative 6 sites, this CLAUDE.md 1 site, `__Pipeline/CLAUDE.md` 1 site); historical `walkthrough-todo-*.md` logs intentionally left alone. Dedicated `sd_median()` subsection added to wrapper narrative `Shared Helpers`, between `gc_preload_refs()` and `get_dop()`, documenting the midyear-interpolation procedure (7 steps: year-of-age → 19-cap → pooled median → midyear-anchor at `floor((ageyears + 0.5) * 365.25)` → `approx()` linear interpolation → `rule = 2` endpoint clamp → sex-duplicate) and the not-sex-stratified caveat; the line-419 parenthetical was shortened to a pointer. `R/child_clean.R:91` file-header bullet for `calc_and_recenter_z_scores()` corrected from "(Steps 11 and 15)" to "(Child Step 15/16 pre-loop)". Part 3: new "Support Functions" top-level section inserted in `child-algorithm-reference.md` after Key Concepts and before Architecture, with nine subsections (`.child_valid()`, `.child_exc()`, `get_dop()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `ewma()`, `ewma_cache_*` / `as_matrix_delta()`, `calc_and_recenter_z_scores()`) following a fixed Purpose / Callers / Inputs+return / Key invariants / Code location format; plus a "Cross-references to wrapper helpers" block pointing at `read_anthro()` / `gc_preload_refs()` / `sd_median()` in the wrapper narrative. `.child_valid()` content migrated out of Key Concepts (shortened to a pointer); the TOC index at lines ~101–110 rewritten as a quick jump-off to the new section. Wrapper-narrative Shared Helpers intro updated to list `.child_exc()` and `get_dop()` and to point at the new "Support Functions" subsection explicitly. Man pages regenerated via `roxygenise()`. Tests unchanged: 63/48/28/41/13. Adult tests not re-run (no adult files touched). Session 8 `identify_temp_sde()` support-function pass: dedicated comment-only + small-cleanup walkthrough of the shared temp-SDE resolution helper in `child_clean.R` (~lines 2045–2200 pre-edit). Inventory confirmed seven call sites in `cleanchild()` (Steps 5, 6, 7, 9, 11 mid-loop, 11 end-of-step, 13); only the Step 13 caller passes non-NULL `exclude_from_dop_ids`. Adult algorithm has its own separate `temp_sde()` in `adult_support.R:449` and is out of scope. Fix-now items F50–F55: title + description roxygen no longer frames as "Step 5 and Step 13" (now names all 7 callers explicitly and describes the NULL-vs-non-NULL contract); two inline `Step 5 vs Step 13` comments rewritten; orphan "identicals precondition" `#` block between `@noRd` and the function def moved into roxygen `@details` as the opening paragraph; five Stata-era `# ----- STEP 1/2/3/4/5 -----` section headers rewritten as current-state prose (STEP 2's description corrected from "Distribute SP medians across all rows for each subject" to "Compact (subjid, param) -> median.spz lookup table" — it was actually creating a lookup table, not distributing); stale "Use by (not keyby) to preserve group order" design-note-to-self comment trimmed; formal `@param df` / `@param exclude_from_dop_ids` / `@return` tags added (precedent from Session 7b F44). Deferred items D-a–D-d also resolved: dead `is.null(df$subjid)` / `is.null(df$param)` defensive blocks removed (all 7 callers pass both columns, branches never fired); redundant `& valid.rows` mask in return dropped (`df$extraneous` already gated by guarded assignment); string-keyed named-vector lookup (`setNames(...); as.character(...)`) replaced with direct integer indexing (`result <- logical(nrow(df)); result[df$orig_row] <- df$extraneous`); two Step 5 / Step 6 call sites switched from base-R `data.df$exclude[identify_temp_sde(...)] <- ` to data.table `data.df[identify_temp_sde(...), exclude := ...]` idiom to match the other four sites (Step 11 mid-loop `sde_result <- identify_temp_sde(...)` separate-variable pattern intentionally left alone — it is not a base-R-vs-data.table style mismatch). Man pages regenerated via `roxygenise()`. Tests unchanged: 63/48/28/41/13. Adult tests not re-run (no adult code touched). Session 7b EWMA support-function pass: dedicated comment-only walkthrough of `ewma()`, `ewma_cache_init()`, `ewma_cache_update()`, and `as_matrix_delta()` in `child_clean.R`. Closes D33 — Stata-style `# 6./a./b./c./i./ii./iii.` block inside `ewma()` (~lines 1662–1682) rewritten as current-state R-idiom prose using current variable names (`tbc.sd`/`ctbc.sd`/`ewma.all`/`before`/`after`/`delta`); `@return` roxygen corrected from "Data frame with 3 variables" to named list (matches actual return + the data.table `:=` rationale on the existing inline comment); `@param ewma.exp` clarified to scalar-or-vector with note that internal callers always pass per-observation vectors; "self-weight is zero (replaces ifelse(delta == 0, …) approach)" parenthetical dropped; paired "O(n) arithmetic instead of matrix copy + multiply" comments rephrased as current-state Before/After block; "Bug fix: was only checking 1 position on each side; extended to 2 positions on each side" line removed in `ewma_cache_update()` (preceding rationale paragraph kept). Bycatch: Child Step 15/16 header "Restructured to use global iterations for efficiency / Key changes:" block rewritten as a single current-state paragraph. `ewma_cache_init()` and `as_matrix_delta()` reviewed and left unchanged. Scope correction: adult algorithm uses its own `adult_ewma_cache_init()` / `adult_ewma_cache_update()` in `adult_support.R` and does not call the child EWMA functions, so this pass is child-scoped. Man pages regenerated via `roxygenise()`. Tests: child unchanged 63/48/28/41/13; adult baseline (198 unit + 1508/1508 × 4 levels) confirmed pre-edit and not affected by edits. Session 7 Child Step 11 EWMA1 deeper pass: removed stale `# Return z-scores and EWMA1 iteration 1 values for comparison` comment in cleanchild output assembly (leftover from Session 6 F34 debug-parameter removal); trimmed Step 11 per-group closure `.SDcols` from 9 cols to 6 (dropped unused `index`/`id`/`sex`, kept `internal_id`/`param`/`agedays`/`tbc.sd`/`ctbc.sd`/`exclude`); clarified narrative "the closure resets those subjects' temp SDEs" → "the post-pass recalc block…"; added `identify_temp_sde()` to Step 11 narrative Code-location cell alongside `ewma()`. Tests unchanged: 63/48/28/41/13. Step naming convention (2026-04-17): all step references prefixed with "Child" or "Adult", headers use name-first format like "EWMA1: Extreme EWMA (Child Step 11)" — applied across both algorithm narratives and both CLAUDE.mds. Session 6 Child Step 11 EWMA1 walkthrough: removed dead `debug` parameter and `ewma1_it1.*` capture block — the block never populated any columns because EWMA fields lived inside a per-group closure that only returned `exclude`, confirmed empirically on stress data. Rewrote Child Step 11 code header block from changelog-style to current-state rationale; replaced Stata-notation `exc_*==0`/`exc_*==2` comments; removed "Fixed cdewma sign" historical comment; `lowest id` → `lowest internal_id` in worst-row comment. Added explicit 11a/11b/11c sub-sections (pre-filter / iteration loop / end-of-step temp SDE refresh) in code and narrative. Child Step 11 narrative rewritten for current-state only; Code-location cell now names file and functions; stale Checklist items 1–3 dropped and renumbered; added "Configurable parameters in scope for Child Step 11" subsection. Session 4b Child Step 7 BIV walkthrough: replaced dead `sd.extreme`/`z.extreme` parameters with 8 per-cell `biv.z.*` parameters; rewrote Child Step 7 narrative for "current state only"; fixed Stata-style 7d comment and inaccurate "overwrites non-temp codes" note; documented that full child permissiveness framework is deferred until after v3.0.0 and that `Child-growthcleanr-permissiveness-specs.md` is to be ignored during walkthroughs/reconciliation. Session 4a: removed run_cf_detection optimization; cf_rescue="all" now rescues every detected CF including shared-SPA ones, with Child Step 13 resolving multi-Include SPAs; HEADCM threshold placeholder moved to Known Issues; pre-session cleanup migrated D5/D6/D9+D10 to Known Issues → Open (wrapper) and fixed D4/D16 in child_clean.R)

## Overview

growthcleanr (gc) is an R package that identifies and flags implausible anthropometric measurements (height, weight, head circumference) in electronic health record data. It does not remove any values — each row gets an `exclude` code: `"Include"` or an exclusion code naming the step and reason.

Repository: https://github.com/carriedaymont/growthcleanr Branch: `efficiency-updates` (v3.0.0, not yet merged to main) CRAN: v2.2.0 (outdated — do not use for development)

This CLAUDE.md is for **Claude Code agents working on gc code**. For pipeline/Qual-AD context, see `__Pipeline/CLAUDE.md`.

**CRITICAL:** The only valid local working copy of growthcleanr is `__Pipeline/gc-github-latest/`. Stale copies exist in `/Users/Shared/` — do NOT use them. Never edit files outside `__Pipeline/` without Carrie's explicit written permission in chat.

**Standing permission** (2026-04-17): While working in `__Pipeline/gc-github-latest/`, it is OK to edit `__Pipeline/CLAUDE.md` without asking — typically to sync its "Last updated" / milestone summary with work that happened in `gc-github-latest/`. Edits to other files outside `gc-github-latest/` still require explicit per-change permission.

**Repo scope note:** Only the `gc-github-latest/` directory is a git repository (origin: `carriedaymont/growthcleanr`, branch `efficiency-updates`). The parent `__Pipeline/` directory is **not** under version control and `__Pipeline/CLAUDE.md` is not in any repo, so it does not need to be staged, committed, or pushed — edits to it are saved on disk only. When committing/pushing growthcleanr work, only stage files under `gc-github-latest/`.

---

## Algorithms

| Algorithm | Status | Default path | Notes |
|-----------|--------|--------------|-------|
| Child | Active, primary | `child_clean.R` | Default pediatric path in v3.0.0 |
| Adult | **Closed pending validation** | `adult_clean.R` + `adult_support.R` | Permissiveness framework, 4 exclusion levels. Do not modify without checking with Carrie first. |

**Legacy removed (2026-04-16):** The legacy pediatric algorithm (`pediatric_clean_legacy.R`, `pediatric_support_legacy.R`), `R/deprec/`, `R/modified_source_code/`, and associated parameters (`use_legacy_algorithm`, `prelim_infants`, `lt3.exclude.mode`, `ewma.exp`, `recover.unit.error`) have all been removed. Prior CRAN versions (v2.2.0) retain the legacy code for backward compatibility.

This CLAUDE.md covers both the child and adult algorithms.

---

## Code Structure

### Key files

| File | What it contains |
|------|------------------|
| `R/child_clean.R` | `cleangrowth()` entry point + exports (top of file); `gc_preload_refs()` (pre-loads reference closures for repeated calls); `cleanchild()` main algorithm + support functions (`.child_valid()`, `.child_exc()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `as_matrix_delta()` (internal), `ewma_cache_init()`/`ewma_cache_update()`, `get_dop()`, `read_anthro()`) |
| `R/utils.R` | Shared utilities (does NOT contain `.child_valid()`) |
| `R/adult_clean.R` | `cleanadult()` main algorithm — permissiveness framework, 14 steps |
| `R/adult_support.R` | Adult support functions — permissiveness presets, EWMA, BIV, height groups, evil twins, error load, etc. |
| `inst/extdata/` | Reference tables (growth charts, recentering file, velocity tables) |
| `tests/testthat/` | Test suite (see Testing section) |

### Preprocessing → main algorithm split

`cleangrowth()` (in `child_clean.R`) handles:
- Input validation, data.table construction
- Imperial → metric conversion (HEIGHTIN, WEIGHTLBS)
- LENGTHCM → HEIGHTCM reclassification (param only, no measurement adjustment — known limitation)
- Age cutpoint split (pediatric vs. adult)
- Z-score calculation (CSD method, WHO/CDC blending)
- GA correction (Child Step 2b) for potcorr subjects
- Recentering → tbc.sd and ctbc.sd
- Missing value identification
- Batching and dispatch to `cleanchild()`

`cleanchild()` processes one batch through all cleaning steps (see Complete List of Steps below).

### Batching

Two layered batching systems:

1. **Outer wrapper** (memory management): Divides subjects into batches (default 2,000, configurable via `batch_size`).

2. **Inner batching** (parallelism): Subdivides for parallel processing. With `parallel = FALSE` (current default), this is a no-op.

---

## Input/Output Format

### Required input (individual vectors, not a dataframe)

| Parameter | Type | Description |
|-----------|------|-------------|
| `subjid` | any | Subject identifier |
| `param` | character | `"WEIGHTKG"`, `"WEIGHTLBS"`, `"HEIGHTCM"`, `"HEIGHTIN"`, `"LENGTHCM"`, or `"HEADCM"` |
| `agedays` | numeric | Age in days at measurement |
| `sex` | any | `0`/`"m"`/`"M"` = male; `1`/`"f"`/`"F"` = female |
| `measurement` | numeric | Recorded value in units specified by `param` |
| `id` | any (required) | Unique row identifier; preserved in output. Can be numeric, character, UUID, etc. |

### Input handling

- `internal_id` (sequential `1:N`) is created at entry for all internal sorting and tiebreaking. The user's `id` is never used for algorithm logic — only for output.

- `measurement == 0` → replaced with `NaN` (treated as missing)
- Imperial converted to metric before processing
- LENGTHCM relabeled to HEIGHTCM (no measurement adjustment)
- Data split by `adult_cutpoint` (default 20 years)

### Sorting (critical for deterministic results)

```
setkey(data.df, subjid, param, agedays, id)
```

The `id` field breaks ties for same-day duplicates. Without it in the sort key, SDE resolution order is undefined.

### Output

`cleangrowth()` returns a data.table with:

| Column | Description |
|--------|-------------|
| `id` | User-provided row identifier (preserved as-is, any type) |
| `exclude` | Exclusion code or `"Include"` |
| `param` | Growth parameter |
| `cf_rescued` | CF rescue reason (empty if not rescued) |
| `sd.orig_who` | WHO CSD z-score |
| `sd.orig_cdc` | CDC CSD z-score |
| `sd.orig` | Blended CSD z-score |
| `tbc.sd` | Recentered blended z-score |
| `ctbc.sd` | Recentered corrected z-score |
| `final_tbc` | ctbc.sd for potcorr, tbc.sd for others |

**Adult-specific output columns** (NA for child rows):

| Column | Description |
|--------|-------------|
| `mean_ht` | Subject mean included height (used in adult 2D WT step) |
| `bin_result` | Binary `"Include"`/`"Exclude"` (adult rows only) |

**Child-specific output columns** (NA for adult rows): `cf_rescued`, `sd.orig_who`, `sd.orig_cdc`, `sd.orig`, `tbc.sd`, `ctbc.sd`, `final_tbc`.

**Breaking change from v2.2.0:** v3.0.0 returns a data.table with multiple columns (not a character vector of exclusion codes like CRAN v2.2.0). Any code consuming `cleangrowth()` output must expect a data.table, not a vector.

---

## Complete List of Exclusion Codes (child algorithm)

Child exclusion codes follow the `Exclude-C-{Reason}` format. Codes are NOT param-specific — the param for each row is in the data (the `param` column). The `.child_exc()` helper function (line ~142 of `child_clean.R`) generates these codes: `.child_exc(reason)` → `Exclude-C-{reason}`.

| Code | Child Step | Param | Description |
|------|------|-------|-------------|
| `Include` | — | All | Value passes all checks |
| `Exclude-Missing` | Init | All | NA, NaN, agedays < 0 |
| `Exclude-Not-Cleaned` | Init | HEADCM | HC with agedays > 5 × 365.25 (no WHO reference above 5 years) |
| `Exclude-C-CF` | 6 | All | Carried forward: identical to prior-day value (not rescued) |
| `Exclude-C-BIV` | 7 | All | Biologically implausible value (absolute or standardized) |
| `Exclude-C-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-C-Traj-Extreme` | 11 | All | Extreme EWMA outlier (EWMA1) |
| `Exclude-C-Identical` | Early 13, 13 | All | Same-day duplicate, identical value |
| `Exclude-C-Extraneous` | 13 | All | Same-day extraneous value (SDE loser by EWMA or other criteria) |
| `Exclude-C-Traj` | 15, 16 | All | Moderate EWMA outlier (EWMA2): covers middle, first, last, birth WT, birth HT/HC, and all sub-variants |
| `Exclude-C-Abs-Diff` | 17 | HT, HC | Height/HC velocity violation (increase > max or decrease > min allowed) |
| `Exclude-C-Pair` | 19 | All | Subject with 2 measurements, one excluded |
| `Exclude-C-Single` | 19 | All | Subject with 1 measurement, excluded |
| `Exclude-C-Too-Many-Errors` | 21 | All | Error ratio exceeds threshold |
| `Exclude-C-Temp-Same-Day` | 5 | All | **Internal only** — temporary SDE flag, must be resolved before output |

Codes are shared across all params — use the `param` column in the data to determine which parameter a code applies to.

**`Exclude-Missing` and `Exclude-Not-Cleaned`** are shared codes assigned in `cleangrowth()` before dispatch.

**`exclude.levels`** has been updated to include all child and adult codes (was previously missing adult codes, causing silent NAs in factor output).

**Simplified vs. detailed codes:** Unlike the adult algorithm which has distinct codes for each sub-step (e.g., separate moderate EWMA codes for firstRV vs allRV), the child algorithm uses simplified codes. For example, all EWMA2 exclusions use `Traj` regardless of whether the value was a middle, first, last, or birth measurement. All SDE resolutions use either `Identical` or `Extraneous`. The specific sub-step logic is in the code comments but not reflected in the exclusion code.

### CF rescue

**Parameter:** `cf_rescue = c("standard", "none", "all")`
- `"standard"` (default): age/interval/param-specific lookup thresholds
- `"none"`: no rescue (all CFs excluded)
- `"all"`: every detected CF rescued, including CFs on an SPA that also has another Include. Child Step 13 final-SDE resolution handles any resulting multi-Include SPAs. Use when the caller wants to treat CFs as plausible and let downstream SDE logic pick among candidates.

**cf_rescued column values:**

| Code | Meaning |
|------|---------|
| `""` | Not CF, or CF not rescued |
| `"Rescued"` | CF rescued by lookup threshold (standard mode) |
| `"Rescued-All"` | CF rescued because cf_rescue="all" |

**Optional cf_detail columns** (enabled by `cf_detail = TRUE`):
- `cf_status`: NA (not CF candidate), "CF-NR" (excluded), or "CF-Resc" (rescued)
- `cf_deltaZ`: absolute z-score difference from originator

CF rescue thresholds are defined in `cf-rescue-thresholds.md` — lookup tables by age bin × interval bin × param × rounding type. Three levels (0.05 / 0.20 / 0.40) plus NR (no rescue). Derived from 100K synthetic tracker population; full methodology in `__Pipeline/CF-exploration/cf-threshold-schemes.md`.

---

## Complete List of Steps (child algorithm)

Step numbers are not consecutive — they align with the Stata implementation. Child Steps 3, 4, 8, 10, 12, 14, 18, 20 either do not exist or are handled within other steps.

| Child Step | Name | Brief description |
|------|------|-------------------|
| (preprocessing) | Z-score calculation | WHO/CDC CSD z-scores, blending, rounding |
| (preprocessing) | GA correction (Child Step 2b) | Fenton-based correction for potcorr subjects |
| (preprocessing) | Recentering | Subtract population medians → tbc.sd/ctbc.sd |
| Early 13 | SDE-Identicals | Remove same-day identical values before CF detection |
| 5 | Temporary SDE | Temporarily flag same-day duplicates |
| 6 | Carried Forwards | Identify and optionally rescue carried-forward values |
| 7 | BIV | Exclude biologically implausible values |
| 9 | Evil Twins | Exclude adjacent extreme values |
| 11 | EWMA1 (Extreme) | Exclude extreme EWMA outliers |
| 13 | Final SDE | Resolve remaining same-day duplicates |
| 15 | EWMA2 (Moderate) | Exclude moderate EWMA outliers |
| 16 | Birth HT/HC | EWMA2 variant for birth values |
| 17 | Height/HC Velocity | Exclude values exceeding velocity limits |
| 19 | Pairs and Singles | Evaluate subjects with 1–2 remaining measurements |
| 21 | Error Load | Exclude all if error ratio too high |
| 22 | Output | Assemble return columns |

---

## Adult Algorithm

### Overview

The adult algorithm (`cleanadult()`) uses a permissiveness framework with 4 preset levels (`loosest`, `looser`, `tighter`, `tightest`) that control all thresholds simultaneously. Default is `looser`. Individual parameters can override presets. Unlike the child algorithm, the adult algorithm does not compute z-scores — it works directly with raw measurements. It uses BMI (computed internally) for some threshold decisions.

### cleanadult() Interface

Called internally by `cleangrowth()`. Input: data.table with `id`, `internal_id`, `subjid`, `sex`, `agedays`, `param`, `measurement`. `internal_id` is a sequential integer created by `cleangrowth()` for deterministic sorting and tiebreaking. Accepts HEIGHTCM/HEIGHTIN and WEIGHTKG/WEIGHTLBS (converts internally). Sex is not used by the algorithm but required by the package interface.

**Output:** Returns input df plus `result` (exclusion code), `mean_ht`, and optionally `bin_result` (default ON), `extraneous`, `loss_groups`, `gain_groups`. Internal columns (`meas_m`, `ageyears`, `age_days`) are dropped. The user's original `id` is preserved untouched. `internal_id` (integer) is used for all internal sorting and named vector indexing (via `as.character()` for names, `as.integer()` for join keys).

### Adult Algorithm Steps

1H/1W BIV → 2W RV markers → 3H/3W Temp SDE → 4W Weight Cap → 9Wa Evil Twins → 9H HT SDE → 9Wb Extreme EWMA → 10H HT Distinct → 10W WT SDE → 11H Mean HT → 11Wa 2D Ord WT → 11Wa2 2D Non-Ord WT → 11Wb Moderate EWMA → 13 Distinct 1D → 14 Error Load

**There is no Adult Step 12W.**

| Adult Step | Name | Brief description |
|------|------|-------------------|
| 1H/1W | BIV | Biologically implausible value exclusion (HT, WT, BMI) |
| 2W | RV markers | Mark repeated values for linked mode |
| 3H/3W | Temp SDE | Temporarily flag same-day duplicates |
| 4W | Weight Cap | Exclude weights at physical scale maximum |
| 9Wa | Evil Twins | Adjacent pair with implausible weight difference |
| 9H | HT SDE | Same-day height resolution (identical + extraneous) |
| 9Wb | Extreme EWMA | Extreme EWMA weight outliers |
| 10H | HT Distinct | 2D height pairs and 3+D height windows |
| 10W | WT SDE | Same-day weight resolution |
| 11H | Mean HT | Compute subject mean height (used by 2D WT) |
| 11Wa | 2D Ord WT | 2D ordered weight pairs (wtallow/perclimit) |
| 11Wa2 | 2D Non-Ord WT | 2D non-ordered weight pairs |
| 11Wb | Moderate EWMA | Moderate EWMA weight outliers + error load escalation |
| 13 | Distinct 1D | Single-measurement exclusion (HT, WT, BMI limits) |
| 14 | Error Load | Exclude all if error ratio exceeds threshold |

### Adult Exclusion Codes

#### Non-SDE Codes

| Code | Adult Step | Description |
|------|------|-------------|
| `Include` | — | Not excluded |
| `Exclude-A-BIV` | 1H/1W | Biologically implausible value (height or weight) |
| `Exclude-A-Scale-Max` | 4W | Weight at scale maximum |
| `Exclude-A-Scale-Max-Identical` | 4W | All weights identical at scale max |
| `Exclude-A-Scale-Max-RV-Propagated` | 4W | RV copy of scale-max exclusion (linked mode) |
| `Exclude-A-Evil-Twins` | 9Wa | Adjacent pair with implausible weight difference |
| `Exclude-A-Traj-Extreme` | 9Wb | Extreme EWMA outlier (independent mode) |
| `Exclude-A-Traj-Extreme-firstRV` | 9Wb | Extreme EWMA outlier (linked firstRV pass) |
| `Exclude-A-Traj-Extreme-firstRV-RV-Propagated` | 9Wb | RV copy of firstRV extreme exclusion (linked mode) |
| `Exclude-A-Traj-Extreme-allRV` | 9Wb | Extreme EWMA outlier (linked allRV pass) |
| `Exclude-A-Ord-Pair` | 10Ha | 2D height pair outside band (one excluded) |
| `Exclude-A-Ord-Pair-All` | 10Ha | 2D height pair outside band (all excluded) |
| `Exclude-A-Window` | 10Hb | 3+D height outside window (one excluded) |
| `Exclude-A-Window-All` | 10Hb | 3+D height outside window (all excluded) |
| `Exclude-A-2D-Ordered` | 11Wa | 2D ordered weight pair outside wtallow/perclimit |
| `Exclude-A-2D-Non-Ordered` | 11Wa2 | 2D non-ordered weight pair |
| `Exclude-A-Traj-Moderate` | 11Wb | Moderate EWMA outlier (independent or firstRV) |
| `Exclude-A-Traj-Moderate-RV-Propagated` | 11Wb | RV copy of firstRV moderate exclusion (linked mode) |
| `Exclude-A-Traj-Moderate-allRV` | 11Wb | Moderate EWMA outlier (linked allRV pass) |
| `Exclude-A-Traj-Moderate-Error-Load` | 11Wb | 4+ consecutive moderate EWMA candidates |
| `Exclude-A-Traj-Moderate-Error-Load-RV` | 11Wb | Error load escalation to entire patient (linked) |
| `Exclude-A-Single` | 13 | 1D measurement outside limits (height or weight) |
| `Exclude-A-Too-Many-Errors` | 14 | Error ratio > threshold |

Adult codes are NOT param-specific — the param for each row is in the data (the `param` column).

**Note:** The `-N` round suffix (iteration number) was removed from adult trajectory codes in the 2026-04-14 exclusion code rename. Trajectory codes no longer include the loop iteration.

#### SDE Codes

| Code | Adult Steps | Description |
|------|-------|-------------|
| `Exclude-A-Identical` | 9H/10W | Same-day identical values (keep one) |
| `Exclude-A-Extraneous` | 9H/10W | Same-day non-identical value (SDE loser) |

### Adult Permissiveness Presets

Default: `"looser"`

| Parameter | loosest | looser | tighter | tightest |
|-----------|---------|--------|---------|----------|
| BIV HT (cm) | 50–244 | 120–230 | 142–213 | 147–208 |
| BIV WT (kg) | 20–500 | 30–270 | 36–159 | 39–136 |
| BIV BMI | 5–300 | 12–65 | 16–45 | 18–40 |
| 1D limits | split (BMI/no-BMI) | same as BIV | same as BIV | same as BIV |
| wtallow formula | PW-H (piecewise) | PW-H (piecewise) | PW-L (piecewise-lower) | allofus15 |
| UW scaling | UW-based (see wtallow-formulas.md) | UW-based | UW-based | cap limited by PW-L |
| ET caps | wtallow cap + 20 | wtallow cap + 20 | wtallow cap + 20 | allofus15-cap-12m |
| perclimit ≤45 kg | 0.5 | 0.5 | 0.7 | 0.7 |
| perclimit 45–80 kg | 0.4 | 0.4 | 0.4 | 0.4 |
| perclimit >80 kg | 0 (disabled) | 0 (disabled) | 0.4 | 0.4 |
| error_load_threshold | 0.41 | 0.41 | 0.29 | 0.29 |
| mod_ewma_f | 0.75 | 0.75 | 0.60 | 0.60 |
| max_rounds | 100 | 100 | 100 | 100 |
| ht_band | 3" | 3" | 2" | 2" |
| allow_ht_loss | TRUE | FALSE | FALSE | FALSE |
| allow_ht_gain | TRUE | TRUE | TRUE | FALSE |
| repval_handling | independent | independent | linked | linked |

For the authoritative cross-algorithm parameter and threshold index (including the role and code location of each preset column — `mod_ewma_f`, `perclimit_*`, `error_load_threshold`, etc.), see [`parameters-reference.md`](parameters-reference.md).

### Key Differences: Adult vs. Child

| Aspect | Child | Adult |
|--------|-------|-------|
| Z-scores | CSD z-scores (WHO/CDC blend) | No z-scores — raw measurements |
| Data structure | Single `data.df`, never removes rows | Copies rows into shrinking dataframes |
| SDE tiebreaking at birth | Keep lowest `id` (pre-postnatal shift) | N/A (no births) — always keeps highest `id` |
| Permissiveness levels | Not yet implemented (planned) | 4 levels: loosest/looser/tighter/tightest |
| Parameters | param = HT, WT, HC | param = HT, WT only (no HC) |
| Rounding tolerance | None (removed) | 0.12 cm/kg on all threshold comparisons |
| Head circumference | Supported (WHO only, ≤3y cleaned) | Not applicable |
| `perclimit` scope | N/A | 11Wa: subject-level max wt; 11Wb: observation-level |
| Sort determinism | `setkey(data.df, subjid, param, agedays, internal_id)` | All sorts include `internal_id` (integer) as final tiebreaker; `as.character()` used for named vector keys |
| Missing-as-infinity | N/A | `ifelse(is.na(...), Inf, ...)` for edge EWMA values |

---

## Configurable Parameters (child)

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `biv.z.wt.low.young` | -25 | Child Step 7 | Lower unrecentered CSD z cutoff for WEIGHTKG at `ageyears < 1` |
| `biv.z.wt.low.old` | -15 | Child Step 7 | Lower unrecentered CSD z cutoff for WEIGHTKG at `ageyears >= 1` |
| `biv.z.wt.high` | 22 | Child Step 7 | Upper unrecentered CSD z cutoff for WEIGHTKG (all ages) |
| `biv.z.ht.low.young` | -25 | Child Step 7 | Lower unrecentered CSD z cutoff for HEIGHTCM at `ageyears < 1` |
| `biv.z.ht.low.old` | -15 | Child Step 7 | Lower unrecentered CSD z cutoff for HEIGHTCM at `ageyears >= 1` |
| `biv.z.ht.high` | 8 | Child Step 7 | Upper unrecentered CSD z cutoff for HEIGHTCM (all ages) |
| `biv.z.hc.low` | -15 | Child Step 7 | Lower unrecentered CSD z cutoff for HEADCM (all ages) |
| `biv.z.hc.high` | 15 | Child Step 7 | Upper unrecentered CSD z cutoff for HEADCM (all ages) |
| `error.load.mincount` | 2 | Child Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Child Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Default uses built-in rcfile; pass a data.table with columns `param`, `sex`, `agedays`, `sd.median` for custom recentering |
| `cf_rescue` | `"standard"` | Child Step 6 | CF rescue mode: `"standard"` (age/interval/param-specific lookup), `"none"` (all CFs excluded), `"all"` (every detected CF rescued; Child Step 13 resolves any multi-Include SPAs) |
| `cf_detail` | FALSE | Child Step 6 | If TRUE, add `cf_status` and `cf_deltaZ` columns to output |
| `ewma_window` | 15 | Child Steps 11, 13, 15/16, 17 | Max Include observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `quietly` | TRUE | All | Suppress progress messages |
| `ref_tables` | NULL | All reads | Pre-loaded closures from `gc_preload_refs()`; skips disk reads |
| `cached_results` | NULL | Partial run | data.table from prior `cleangrowth()` call. Auto-detects changed subjects if `changed_subjids` is NULL; uses explicit list if both provided |
| `changed_subjids` | NULL | Partial run | Optional vector of subject IDs to re-run. If NULL and `cached_results` provided, changed subjects are auto-detected |
| `batch_size` | 2000 | Batching | Number of subjects per processing batch |
| `length.adjust` | FALSE | Preprocessing | If TRUE, subtracts 0.7 cm from LENGTHCM measurements with `agedays > 2 × 365.25` before z-score calculation, converting post-infancy recumbent-labeled lengths to the standing equivalent assumed by the reference standards |

### gc_preload_refs() and partial runs (cached_results)

Canonical documentation for both features lives in `wrapper-narrative-2026-04-17.md` under "Partial runs and preloaded references" (in the Batching and Dispatch section). Quick reference:

- `gc_preload_refs()` returns a list of CDC/WHO reference closures (`list(mtz_cdc_prelim, mtz_who_prelim)`) that can be passed to `cleangrowth()` via `ref_tables = refs` to skip per-call disk reads (~0.9 sec/call saved).
- `cached_results` + optional `changed_subjids` let you re-process only subjects whose input rows changed vs. a prior run. Auto-detect mode (omit `changed_subjids`) compares subjid/param/agedays/sex/measurement against the cache.

Observed speedups (200 subjects, 50 000 calls):

| % subjects changed | Time/rep | vs standard | 50K saving |
|---|---|---|---|
| 100% (full run, refs only) | 3.40 s | 1.3× | 15 hrs |
| 50% changed | 1.93 s | 2.3× | 35 hrs |
| 25% changed | 1.16 s | 3.8× | 46 hrs |
| 10% changed | 0.58 s | 7.7× | 54 hrs |
| 5% changed | 0.40 s | 11.1× | 56 hrs |

Designed for two workflows: (1) error-injection pipelines where most subjects are unchanged per rep, and (2) secure-environment runs where only a subset of subjects has new data.

### Configurable Parameters (adult, via cleangrowth)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adult_permissiveness` | `"looser"` | Sets defaults for all adult sub-parameters |
| `adult_scale_max_lbs` | `Inf` | Physical scale upper limit in lbs |
| `ewma_window` | `15` | Max observations on each side for EWMA (Adult Steps 9Wb, 11Wb); shared with child |

All adult sub-parameters (BIV limits, 1D limits, wtallow formula, etc.) can be passed individually to `cleanadult()` to override the preset. See `adult_clean.R` roxygen for the full list. `cleangrowth()` exposes `adult_permissiveness`, `adult_scale_max_lbs`, and `ewma_window` for adult.

---

## Key Concepts for Code Work

### The `.child_valid()` function

Critical gatekeeper — returns a logical vector of which rows are eligible for each step. The `include.*` flags control whether temporarily excluded categories participate:

- `.child_valid(df)` — only non-excluded, non-Exclude-Missing, non-Exclude-Not-Cleaned
- `.child_valid(df, include.temporary.extraneous = TRUE)` — also temp SDEs
- `.child_valid(df, include.extraneous = TRUE)` — also permanent SDEs
- `.child_valid(df, include.carryforward = TRUE)` — also CFs

Works by checking the text of the `exclude` column, not factor ordering. Uses the new `Exclude-Missing` and `Exclude-Not-Cleaned` code names.

### CSD z-scores (not LMS)

All z-scores use the Conditional Standard Deviation method:
```
If measurement < M:  sd = (measurement - M) / csd_neg
If measurement >= M: sd = (measurement - M) / csd_pos
```
This is intentionally more sensitive to extreme high values than standard LMS z-scores. Intermediate z-score rounding has been removed (previously used Stata-style rounding for cross-platform validation).

### WHO/CDC age blending

| Age range | HT/WT formula | HC |
|-----------|---------------|-----|
| < 2 years | WHO only | WHO only |
| 2–5 years | `(cdc × (age-2) + who × (5-age)) / 3` | WHO only |
| > 5 years | CDC only | WHO only (through 5y) |

HC is WHO-only at all ages. HC > 5y is `Exclude-Not-Cleaned` (no WHO reference data above 5 years). HC measurements 3–5y are cleaned by all standard steps; Child Step 17 applies only `mindiff = -1.5` with no upper bound for these ages since the WHO HC velocity reference ends at 24 months.

### Age-dependent id tiebreaking

At birth (agedays == 0): keep lowest id (earliest measurement, before postnatal fluid shifts). All other ages: keep highest id (later measurement, likely more careful). This differs from the adult algorithm, which always keeps the highest id.

### Designated Other Parameter (DOP)

Weight's DOP is height; height's DOP is weight; HC's DOP is height. Used in Child Steps 5, 6, 11, 13, 15, 19 for cross-parameter plausibility checks.

### data.table reference semantics

The child algorithm keeps all rows in a single `data.df` data.table — rows are never physically removed. Exclusions work by setting the `exclude` column and using `valid()` to filter. This differs from the adult algorithm, which copies rows into separate shrinking dataframes.

`:=` modifies `data.df` in place. Be careful not to accidentally modify copies or joined tables.

---

## Testing

### testthat suite

Full details in `testing-reference.md` (test-by-test listing, coverage gaps, run instructions).

6 actively maintained test files in `tests/testthat/`:

| File | Tests | Assertions | Scope |
|------|-------|------------|-------|
| `test-cleangrowth.R` | 13 | ~80 | `cleangrowth()` API: adult integration (7 tests), child-adult spanning (3 tests) |
| `test-child-regression.R` | 8 | ~48 | Frozen counts, spot checks, cross-sample stability, Missing, HC, preload_refs, changed_subjids |
| `test-child-algorithms.R` | 24 | ~40 | CF rescue (modes + threshold cells), Evil Twins/OTL (1/2/3 unit errors + collateral), Error Load (threshold/mincount/constructed), HC boundary, cf_detail, parallel, GA correction (potcorr + near-potcorr), Birth EWMA2 (extreme/normal WT, extreme HT), LENGTHCM identity |
| `test-child-parameters.R` | 7 | ~15 | include.carryforward, biv.z.wt.high, ewma_window, error.load params, imperial units, LENGTHCM |
| `test-child-edge-cases.R` | 12 | 24 | Single subject, sparse data, all-NA, mixed NA, SDE-Identical, negative agedays, HEADCM >3yr, extreme values, density mix, CF, deterministic |
| `test-adult-clean.R` | 198 | 198 | All 14 adult steps, 4 permissiveness levels, edge cases (via `cleanadult()` directly) |
| **Total** | **262** | **~405** | |

Additional test files not actively maintained: `test-cdc.R` (not modified in v3.0.0); `test-utils.R` (10 failures — old API without `id` parameter, tests utility functions `splitinput`, `recode_sex`, `longwide`, `simple_bmi`).

Run with:
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'
```

### Adult regression harness (separate from testthat)

Regression test: `tests/test_harness.R` runs `cleanadult()` against `inst/testdata/adult-gc-test-ALL-PHASES.csv` — 1508 rows at all 4 permissiveness levels.

```bash
# Unit tests
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-adult-clean.R")'

# Regression test (any permissiveness level)
Rscript tests/test_harness.R [loosest|looser|tighter|tightest]
```

**Test CSV columns:** `id`, `subjid`, `param`, `agedays`, `sex`, `measurement`, plus `expected_loosest`, `expected_looser`, `expected_tighter`, `expected_tightest`, `test_description`, `test_category`, `precision_sensitive`.

**Deferred test gaps:** Error load with -5 exponent, weight scaling at permissiveness levels.

**IMPORTANT: Reinstall before testing.** Tests load `library(growthcleanr)`, which uses the *installed* package, not source files. After any code changes, always reinstall before running tests:
```r
devtools::install_local("gc-github-latest",
  force = TRUE, upgrade = "never")
```
Without this, tests run against stale installed code and results are misleading — tests may pass or fail for the wrong reasons. (`devtools::load_all()` loads source code in the current R session but does not affect fresh `Rscript` processes used by `testthat::test_file()`.)

### Runtime benchmarks

| Dataset | Size | Sequential | Parallel 2 | Parallel 4 |
|---------|------|-----------|------------|------------|
| syngrowth (full) | 77,721 rows, 2,697 subj | ~114 sec | ~65 sec (1.7×) | ~30 sec (3.9×) |
| syn.csv (1,036 subj subset) | 28,434 rows | ~18 sec | — | — |
| 500-subject subsample | ~14K rows | ~11 sec | ~8 sec | ~6 sec |

Parallel produces identical results to sequential (verified 77,721 rows, 0 mismatches). Requires installed package — see Known Issues. Run in background from Claude Code (`run_in_background: true`).

---

## Known Issues

### Open (child)

- **HEADCM CF rescue thresholds reuse HEIGHTCM matrices (confirmed intentional, 2026-04-19).** `.cf_rescue_lookup()` uses the same `other` and `imperial` matrices for HEADCM as for HEIGHTCM. Both parameters are WHO-only with similar reference-standard structure; the WHO HC velocity files used in Child Step 17 are in raw units (cm/interval) and not directly usable for z-score threshold derivation. HT thresholds are confirmed appropriate as-is.
- [ ] **F96: Step 17 EWMA caching** (performance, deferred 2026-04-19): Child Step 17 calls `ewma()` directly on every iteration rather than the `ewma_cache_*` API used by Steps 11/13/15/16; the 2-row branch also invokes `ewma()` even though its tiebreaker keys on `abs(tbc.sd)` rather than the EWMA result. No correctness impact. `child_clean.R:~4593`.
- **`parallel = TRUE` requires installed package.** Will fail with `load_all()` only — workers need `system.file()` access to extdata. Install with `devtools::install_local(".")` first.
- [ ] **Ensure recommendation to use all available data is prominent in user-facing documentation.** Using all available data (including outside the study age range) is intentional and recommended — surrounding context improves cleaning accuracy. As a consequence, the same measurement may receive a different exclusion outcome depending on what other data are available for that subject; this is expected behavior, not a bug. Users who run growthcleanr on a pre-filtered dataset (e.g., ages 2–10 only) may get less accurate results than if they ran on the full dataset and then restricted to the age range of interest for analysis. This recommendation should appear prominently in the vignette and/or user documentation.

### Robustness audit (2026-04-14)

Comprehensive audit of code patterns that could warn or degrade on large/unusual datasets. No blocking issues found.

**Investigated and confirmed safe:**
- **Division by zero in adult slope calculations** (`adult_support.R` ~lines 823, 836): `(ap1 - ap2)` and `(an1 - an2)` denominators cannot be zero because all same-day duplicates are resolved before these steps run (only one Include per subject-param-ageday). If zero somehow occurs, it indicates an upstream bug and should break loudly.
- **Division by zero in adult BMI** (`adult_support.R` ~line 1476): `ht_val` comes from included heights that already passed BIV (min 50cm at loosest). Cannot be zero.
- **HT+WT same-day pairing assumption:** Both algorithms handle unpaired measurements correctly. Adult uses `intersect()` for BMI days — skips if no shared days. Child DOP is a secondary tiebreaker with explicit NA→Inf fallback when the other param doesn't exist for a subject.

**Low priority cleanup (functional but could be cleaner):**
- [ ] **`suppressWarnings(min/max)` + `is.infinite()` pattern** (child_clean.R ~lines 4006, 4109, 4117, 4126): Suppresses warnings from `min()`/`max()` on empty subsets, then converts `-Inf`/`Inf` to `NA`. Works correctly but masks any unexpected warnings. Could be rewritten with explicit `if (length(...) == 0)` guards to match the style of the 5 locations fixed in commit `64c186b`.
- [ ] **`which.min`/`which.max` on empty vectors** (child, ~lines 4442, 5201): Returns `integer(0)` on empty input. Surrounding logic has row-count checks, so safe in practice.
- [ ] **Growing vectors with `c()` in loops** (adult, ~lines 291-292, 559): O(n²) pattern — `exc_ids <- c(exc_ids, new_id)`. Negligible for typical subject sizes but could slow down on subjects with thousands of measurements.
- [ ] **Step 15/16 low-value walkthrough deferrals** (2026-04-19): (a) dead `df[, exp_vals := cache$exp_vals]` in cache-hit branches (`:3801`/`:4009`) — written to local `df` but not read in that branch, kept for symmetry with cache-miss branch; (b) unused `sex`/`v` in `.SDcols` for Step 15 (`:3926–3927`) and Step 16 (`:4070–4071`) closures; (c) redundant `order()` at `:3672` — `data.df` already keyed at end of Step 13 Phase B4; (d) structural `15a/15b/15c` + `16a/16b` subsection headers (cosmetic); (e) file-header EWMA bullet at `child_clean.R:94` — general wording sufficient, Support Functions section in `child-algorithm-reference.md` has detail. All behavior-neutral.

### Open (wrapper)

Deferred from walkthrough sessions 1–3 (2026-04-16):

- [ ] **`cleangrowth()` `@return` roxygen incomplete.** Roxygen lists only part of the returned data.table's columns. The actual return also includes `internal_id, subjid, agedays, sex, line, bin_exclude, tri_exclude`, conditional columns (`cf_status`/`cf_deltaZ` under `cf_detail=TRUE`), and internal working columns (`v`, `v_adult`, `measurement`, `age_years`, `newbatch`, `sd.corr`, `sd.orig_uncorr`, `z.orig`, `uncorr`, `potcorr`). Before updating the roxygen, decide which internal columns should be dropped from output vs. kept as part of the public contract, then align `@return` to the final schema.
- [ ] **Velocity reference tables not in `gc_preload_refs()`.** `tanner_ht_vel` is loaded once per `cleangrowth()` call; WHO velocity tables (`who_ht_vel_for_age`, `who_hc_vel_for_age`) are loaded once per batch inside `cleanchild()`. Both are tiny (< 1 kB each gzipped) and the `fread()` calls are fast, so correctness is unaffected. Extension: add the three tables to `gc_preload_refs()` and plumb optional parameters through `cleangrowth()` / `cleanchild()`; keep existing `fread()` fallbacks.

### Open (adult)

- [ ] **Deferred test gaps:** Error load with -5 exponent, UW scaling edge cases (very low/high UW).
- [ ] **Performance:** `setkey(df, subjid)` optimization deferred.
- [ ] **`adult_ewma_cache_update()` not wired into iteration loops (2026-04-20).** `remove_ewma_wt()` (`adult_support.R:942`) and `remove_mod_ewma_wt()` (`:1071`) both call `adult_ewma_cache_init()` at the top of every round, rebuilding the full O(n²) delta matrix each iteration. `adult_ewma_cache_update()` (`:347`) is an O(n) incremental update that subtracts column j and trims the matrix; it was implemented but never wired in. **Integration plan for `remove_ewma_wt()`:** (1) Move `adult_ewma_cache_init()` call to before the while loop; (2) Remove it from inside the loop; (3) After `subj_df <- subj_df[-to_rem, ]`, add `cache <- adult_ewma_cache_update(cache, to_rem)`. Results should be regression-identical: the adult EWMA uses a fixed exponent (−5) so no neighbor-row recomputation is needed (unlike the child, whose variable exponent −1.5 to −3.5 triggers a 4-neighbor rebuild loop). The position-based windowing approximation (zero entries near the window boundary may not shift correctly after removal) is shared with the child cache and negligible in practice — at `ewma_window = 15`, weights at the boundary are `(gap_days + 5)^(−5)` ≈ 10⁻¹¹ to 10⁻¹⁴. `remove_mod_ewma_wt()` has the same structure and can be updated the same way; it is more complex (7-step loop) so should be done separately after `remove_ewma_wt()` is validated. Both callers are listed in `var_for_par` for parallel worker export. **Adult code is closed pending validation — do not implement without checking with Carrie.**

### Fixed (recent — condensed 2026-04-16)

See git history for full details on each fix.

- **Pre-walkthrough cleanup (2026-04-16):** Resolved deferred walkthrough items: Child Step 17 dead inner `if (count_exclude >= 1)` branch collapsed; adult dead `keeper_id` removed (2 locations); adult dead `"temp extraneous"` filter removed in `eval_2d_nonord()`; Child Step 15 Birth WT block now uses `.child_exc(param, "Traj")` for consistency with surrounding code (was 4 hardcoded literals). (The `debug` parameter added in this round was subsequently removed during the Child Step 11 walkthrough — the capture block was dead code that never populated the promised columns.)
- **Doc/dead-code cleanup (2026-04-16):** Deleted 6 orphaned man pages for removed extdata files. Removed dead `nnte` column assignment and unused NULL declarations (`sum_val`, `no_dup_val`, `no_outliers`, `no_bigdiff`, `nottoofar`, `nnte_full`). Fixed 3 stale comments (legacy file refs). Cleaned `_pkgdown.yml` (removed deleted functions/datasets). Updated vignettes (removed `prelim-infants-algorithm.Rmd`, fixed stale params in `configuration.Rmd`, `usage.Rmd`, `large-data-sets.Rmd`). Fixed `ewma()` roxygen (added missing `@param` tags, separated from `as_matrix_delta`). `as_matrix_delta` now internal (removed from exports). Regenerated all man pages via `roxygenise()`.
- **Legacy removal (2026-04-16):** Deleted `pediatric_clean_legacy.R`, `pediatric_support_legacy.R`, `R/deprec/`, `R/modified_source_code/`, 5 legacy extdata files. Removed parameters: `use_legacy_algorithm`, `prelim_infants`, `lt3.exclude.mode`, `ewma.exp`, `recover.unit.error`. Simplified `read_anthro()` (removed individual WHO txt path), `gc_preload_refs()` (2 closures, not 3), recentering (removed NHANES/derive string paths), dispatch (direct `cleanchild()` only), `var_for_par` (removed legacy functions).
- **Exclusion code work (2026-04-14–16):** Codes renamed to `Exclude-C-{Reason}` / `Exclude-A-{Reason}` format; param indicators removed; `Missing` → `Exclude-Missing`; adult `-N` round suffix removed; `.child_exc()` helper; `exclude.levels` updated with all adult codes.
- **CF rescue thresholds (2026-04-14):** Age/interval/param- specific lookup tables replace fixed thresholds.
- **internal_id (2026-04-12–16):** Sequential integer for all internal sorting; user `id` preserved untouched.
- **Walkthrough cleanups (2026-04-15–16):** 30+ items: `cat()` → `message()`, dead code removal, `batch_size` parameter, correctness fixes.
- **Earlier fixes:** BIV preterm threshold, error load bug (hardcoded .4 → parameter), missing data bug (duplicate `valid()` definitions), CF rescue re-inclusion, z-score blending bug, parallel processing fixes, outer batching wrapper fix.

---

## Next Priorities

1. **Child algorithm work** — current focus. Adult is closed pending clinician validation.
2. **Design extreme/clinical test patients from literature:** Create synthetic patients based on real clinical scenarios with literature-sourced growth values (e.g., hydrocephalus with shunt placement, craniosynostosis, failure to thrive, severe obesity, Turner syndrome, growth hormone deficiency). Goal: verify gc handles physiologically real but extreme trajectories correctly — not just random perturbations.

### Child permissiveness framework — deferred until after v3.0.0

A full child permissiveness framework (analogous to the adult 4-level system) is **deferred until after v3.0.0 is released**. For v3.0.0, child parameters are exposed as individual user-settable scalars (e.g., the eight `biv.z.*` Child Step 7 cutoffs added 2026-04-17), not as permissiveness presets.

`Child-growthcleanr-permissiveness-specs.md` contains draft scaffolding for the future framework and is kept in the repo for reference. **During walkthroughs, code reviews, and documentation reconciliation it should be ignored** — it is not part of the current implementation contract and is not expected to reflect current code or parameter names.

### Adult status (closed 2026-04-09)

**Do not modify adult code or tests without checking with Carrie first.** Adult algorithm is closed pending clinician validation. Walk-through complete (2026-04-09): 198/198 unit tests, 1508/1508 regression tests at all 4 permissiveness levels. Deferred items in `walkthrough-todo-2026-04-09.md`. Completed items:
- Adult integration into package (2026-04-03)
- wtallow redesign (2026-04-09)
- Algorithm walk-through with wtallow reconciliation (2026-04-09)

---

## CRAN Preparation Checklist

Identified 2026-04-03. Items marked [x] are done.

**Timeline note (updated 2026-04-04):** CRAN is not the immediate deadline. Priority order: (1) adult validation, (2) child validation, (3) CRAN cleanup (can overlap with validation as long as changes don't affect algorithm performance).

### Critical (R CMD check ERRORs/WARNINGs)

- [ ] **`parallel` not in DESCRIPTION Imports:** Imported in NAMESPACE but not declared in DESCRIPTION. Add to Imports. (~2 min)
- [ ] **Vignettes excluded via `.Rbuildignore`:** `^vignettes$` prevents vignette building. Stale parameter references cleaned up (2026-04-16): removed `prelim-infants-algorithm.Rmd`, updated `configuration.Rmd`, `usage.Rmd`, `large-data-sets.Rmd`. Still need to verify vignettes build cleanly or remove them.
- [ ] **Package tarball > 5 MB:** `inst/extdata/` alone is ~5.3 MB. `test_syngrowth_sas_output_compare.csv.gz` (1.8 MB) is the largest — move to a companion data package or remove if only needed for development. `stress_test_data.csv` (4.4 MB) in `tests/testthat/` also ships. (~1–2 hours to audit and reorganize)

### High (R CMD check NOTEs)

- [x] **Missing NULL declarations in `cleanadult()`:** Added `result <- mean_ht <- ... <- NULL` for data.table `:=` columns. (Fixed 2026-04-03)
- [x] **`cat()`/`print()` in `cleanadult()`:** Changed to `message()` for CRAN-preferred output handling. (Fixed 2026-04-03)
- [x] **`R/deprec/` and `R/modified_source_code/` removed (2026-04-16):** Entire directories deleted.
- [x] **`cat()`/`print()` in child algorithm (2026-04-15):** All `cat()`/`print()` calls converted to `message()`.

### Medium (pre-submission polish)

- [ ] **Verify `LICENSE` file exists:** DESCRIPTION says `MIT + file LICENSE`. Confirm plain `LICENSE` file (not just `LICENSE.md`) is at package root.
- [x] **Orphaned man pages cleaned up (2026-04-16):** Deleted `bmianthro.Rd`, `lenanthro.Rd`, `weianthro.Rd`, `nhanes-reference-medians.Rd`, `growth_cdc_ext.Rd`, `fentlms_forz.Rd`. Regenerated `cleangrowth.Rd`, `ewma.Rd`, `fenton2025_ms_lookup_smoothed.Rd` via `roxygenise()`. `as_matrix_delta` removed from exports (internal helper).
- [x] **Roxygen return value for `cleangrowth()`:** Updated to describe data.table with all columns. (Fixed prior to 2026-04-16; confirmed in generated `cleangrowth.Rd`.)
- [ ] **`plyr` dependency:** Used for `ddply` in both child and adult parallel dispatch (Phase 12, Phase 13). Long-term: replace with `data.table` / `foreach` equivalent. Replacing plyr also eliminates the two harmless `codetools` warnings emitted per `parallel = TRUE` call (see wrapper reference Phase 5). (~2 hours)
- [ ] **Test runtime:** Ensure total test suite completes in < 10 min. Long-running tests should use `skip_on_cran()`. (~30 min to audit and add skips)

---

## Coding Conventions

- **Sex:** 0 = male, 1 = female
- **Param values:** `"WEIGHTKG"`, `"HEIGHTCM"`, `"HEADCM"` (after input conversion)
- **Use `data.table`** throughout — consistent with gc internals
- **Rounding:** Stata-style rounding (`round_stata()`) has been removed from the child algorithm z-score pipeline. Z-scores are no longer rounded to 0.001 at intermediate steps. The function was only needed for Stata comparison.
- **Factor levels for `exclude`:** Assigning an unlisted string to a factor silently produces NA — always verify new codes exist in `exclude.levels`
- **Sort order:** `setkey(data.df, subjid, param, agedays, id)` — assumed by many steps. Re-sort after any modification that could change order.
- **`by`-group correctness:** `by = subjid` vs. `by = .(subjid, param)` — wrong grouping silently produces wrong results. Check carefully.
- **C++/Rcpp ruled out** — CRAN target, Rtools not always available on enterprise systems
- **Report weight before height** in summaries and output

---

## Pitfalls (for Claude agents)

Things that have caused bugs before or fail silently:

- **Do not create additional `.child_valid()` definitions.** `.child_valid()` is defined in `child_clean.R`. Do not duplicate it in other files — R loads all `.R` files alphabetically and the later definition wins silently.
- **`parallel = TRUE` requires installed package.** Will fail with `load_all()` only — workers need `system.file()` access to extdata. Install with `devtools::install_local(".")` first.
- **Do not modify sort order without re-sorting.** Many steps assume `setkey(data.df, subjid, param, agedays, id)`. If you add/modify rows, call `setkey()` again.
- **Do not assign unlisted strings to the `exclude` factor.** Unlisted values silently become NA. Always verify new exclusion codes exist in `exclude.levels` first.
- **Do not use `by = subjid` when you need `by = .(subjid, param)`.** Wrong grouping produces wrong results with no error or warning.
- **No intermediate z-score rounding.** `round_stata()` has been removed from the child algorithm. Do not reintroduce rounding at intermediate z-score steps.
- **Do not have duplicate support files in `R/`.** R loads all `.R` files in `R/` alphabetically. If two files define the same function, the later one wins silently. This caused breakage when both `02a_support.R` (new) and `adult_support.R` (old) coexisted — the old file's `check_between`/`round_pt` overwrote the new versions. **Rule:** Only one adult support file should exist in `R/`.
- **Adult `result` vs child `exclude`:** The adult algorithm uses `result` as its exclusion column name. `cleangrowth()` maps `res$result` → `exclude` in the combined output. Do not rename `result` to `exclude` inside `cleanadult()`.
- **Adult rounding tolerance:** 0.12 cm/kg on all threshold comparisons. This is intentional and should not be removed.

---

## GitHub and Distribution

- **GitHub:** `carriedaymont/growthcleanr`, branch `efficiency-updates`
- **CRAN:** v2.2.0 (does not include child algorithm or v3.0.0 changes)
- **Install from local:**
  ```r
  devtools::install_local("gc-github-latest",
    force = TRUE, upgrade = "never")
  ```
- **Install from GitHub:**
  ```r
  devtools::install_github("carriedaymont/growthcleanr",
    ref = "efficiency-updates",
    force = TRUE, upgrade = "never")
  ```
- v3.0.0 has not been released on GitHub or CRAN yet

---

## Algorithm Narrative Documents

Three companion narratives document the package. All live in `gc-github-latest/`.

| Narrative | Scope |
|---|---|
| `wrapper-narrative-2026-04-17.md` | `cleangrowth()` public interface, preprocessing (incl. Child Step 2b GA correction), z-score calculation pipeline (CSD, WHO/CDC blending, recentering), batching, dispatch, output reassembly, wrapper-level configurable parameters, canonical Code Review Checklist |
| `child-algorithm-reference.md` | Child algorithm cleaning steps (Child Steps 5, 6, 7, 9, 11, 13, 15/16, 17, 19, 21, 22), child-specific concepts (DOP, `.child_valid()`, EWMA, CF rescue), Support Functions section (algorithm-internal helpers), child-specific variables, child exclusion codes, brief Z-Score Summary pointing to wrapper |
| `adult-algorithm-narrative.md` | Adult algorithm cleaning steps (Adult Steps 1, 2W, 3, 4W, 9Wa/9H/9Wb, 10H/10W, 11H/11Wa/11Wa2/11Wb, 13, 14), permissiveness framework, rounding tolerance |

The two algorithm narratives reference the wrapper narrative for shared material (input format, z-score mechanics, batching, the code review checklist) so the same content does not need to be maintained in three places.

### Status

**Wrapper narrative: extraction pass complete (Pass 1b, 2026-04-17).** All wrapper-level material has been moved out of the child narrative and into a dedicated file. The `cleangrowth()` line-by-line walkthrough and output reassembly mechanics are deferred to a fill-in pass (Pass 2).

**Child narrative: COMPLETE for all walked steps.** Child Steps 5, 6, 7, 9, 11, 13, 15/16, 17, 19, 21, 22 all documented with summary tables, overview, logic details, rationale, and code review checklist findings.

**Adult narrative: COMPLETE.** Adult algorithm closed pending clinician validation; do not modify without checking with Carrie first.

**Open questions and issues:** Tracked inline in each narrative's "Open Questions" or "Checklist findings" sections. Items marked `[x]` are resolved; `[ ]` items are still open.

### Structure

Each step or phase section includes:
- Summary table (scope, prior/next step, exclusion codes)
- Overview
- Key terms and variable names
- Logic and implementation
- Rationale for selected decisions
- Code review checklist findings (12-point checklist; the canonical checklist itself lives in the wrapper narrative)

### Step naming convention

All step references are prefixed with "Child" or "Adult". Headers use name-first format: `## EWMA1: Extreme EWMA (Child Step 11)`, `## Moderate EWMA (Adult Step 11Wb)`. Inline references: `Child Step N`, `Adult Step NX`, `Child Steps 15/16`. See the 2026-04-17 entry in the "Last updated" header for the rationale.

---

## Terms Reference

This is the **canonical** terms reference (primary source). `__Pipeline/CLAUDE.md` has a convenience copy; update both if adding terms.

| Term | Full name |
|------|-----------|
| ewma | Exponentially weighted moving average |
| dewma | Deviation from EWMA |
| absdewma | Absolute value of dewma |
| spa | Subject-parameter-ageday combination |
| sde | Same-day extraneous value |
| dop | Designated other parameter |
| tbc | To-be-cleaned (recentered z-score) |
| ctbc | Corrected to-be-cleaned (recentered corrected z-score) |
| potcorr | Potentially correctable (prematurity) |
| biv | Biologically implausible value |
| csd | Conditional standard deviation |
| nnte | No need to evaluate (legacy; always FALSE in current code) |
| otl | Out of line (used in Evil Twins step; formerly "oob"/"out of bounds") |

---

## Running growthcleanr

### Standalone (from gc-github-latest directory)

```r
library(growthcleanr)
library(data.table)

dt <- fread("path/to/data.csv")
result <- cleangrowth(
  subjid = dt$subjid,
  param = dt$param,
  agedays = dt$agedays,
  sex = dt$sex,
  measurement = dt$measurement,
  parallel = FALSE
)
```

### From Pipeline (error-impact)

See `__Pipeline/CLAUDE.md` for pipeline-specific wrapper instructions. Working directory must be `error-impact/` so `here::here()` resolves correctly.
