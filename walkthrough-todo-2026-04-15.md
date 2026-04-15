# GC Full Walk-Through — Deferred To-Dos

Items identified during the 2026-04-15 walk-through that are deferred for later resolution.

## Summary

**Scope:** Full walk-through of `cleangrowth()` wrapper (chunks 1-8), support functions, and `cleanchild()` algorithm (all steps). Adult algorithm was walked through on 2026-04-09 (see `walkthrough-todo-2026-04-09.md`) and is closed pending clinician validation.

**Bugs fixed during walk-through:** 2
1. `cutpoint_update` dead code (chunk 3) — cutpoint clamping was computed but never applied
2. CLAUDE.md child exclusion code table was completely wrong (chunk 4) — listed ~30 nonexistent codes

**Documentation fixes during walk-through:** 2
1. `__Pipeline/CLAUDE.md`: "id defaults to NULL" → "id is required"
2. `gc-github-latest/CLAUDE.md`: Complete rewrite of exclusion code table

**Total deferred items:** 37
- CRAN compliance (`cat()`→`message()`): #4, #13, #15, #19, #20, #37 — bulk fix recommended
- Dead code / commented-out code: #10, #17, #25, #36
- Code clarity / style: #7, #8, #12, #14, #18, #23, #28, #31
- Minor correctness: #5, #16
- Potential improvements: #1, #2, #3, #6, #9, #11, #21, #22, #24, #26, #27, #29, #30, #32, #33, #34, #35

**No bugs found in the child algorithm.** The code is well-structured with clear step separations, appropriate safety checks (duplicate Include warnings), incremental EWMA caching for performance, and deterministic tiebreaking throughout.

---

## Pre-walk-through setup

No issues — all tests pass. Baseline counts:
- test-cleangrowth.R: 92 PASS
- test-child-regression.R: 54 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS, 2 WARN (codetools)
- test-child-parameters.R: 21 PASS, 1 WARN (expected deprecation)
- test-adult-clean.R: 198 PASS
- Adult regression: 1508/1508 at all 4 levels

---

## Wrapper walkthrough

### Chunk 1: Parameter deprecation & validation (lines ~333-425)

**Fixed:**
- __Pipeline/CLAUDE.md line 246-247: said "wrapper merges by row position when id = NULL" but code actually stop()s. Updated to "id is required".

**Deferred:**

1. **`ewma.exp` silently ignored for child algorithm** — `cleangrowth()` line 350, `ewma.exp = -1.5`. This parameter is passed to `cleanlegacy()` (lines 1485, 1538) but NOT to `cleanchild()` (lines 1497-1517). The child algorithm computes its own age-dependent exponents internally (-1.5 for ≤1y, -3.5 for ≥3y, linear between). No behavioral impact — child always uses correct values. **Fix:** Document `ewma.exp` as legacy-only in roxygen and CLAUDE.md. Consider adding a warning when `ewma.exp` is explicitly set with `use_legacy_algorithm = FALSE`.

2. **No conflict detection for `include.carryforward` + `cf_rescue`** — Lines 395-403. If user passes both `include.carryforward = TRUE` and `cf_rescue = "none"`, the deprecated parameter silently wins (overwrites to "all"). Unlikely edge case since both deprecated and new param would be used simultaneously. **Fix:** Add a warning when both are explicitly set to conflicting values (low priority).

3. **`sdmedian.filename` and `sdrecentered.filename` are legacy-only** — Lines 345-346. Used only by legacy algorithm (lines 1346, 1415). Not dead code. **Fix:** Document as legacy-only in roxygen (low priority).

### Chunk 2: Partial-run mode (lines ~428-507)

No issues found. Input vector alignment correct (lines 501-506). Auto-detect comparison uses identical transformations to main path. Hash ordering deterministic via setkeyv.

### Chunk 3: Data construction + input validation (lines ~509-558)

**Fixed:**
- **`cutpoint_update` dead code — cutpoint clamping broken.** Lines 551-558 computed a clamped cutpoint into `cutpoint_update` but never assigned it back to `adult_cutpoint`. Lines 772/778 used the unclamped value. Roxygen says values outside 18-20 are limited, but they weren't. Fixed: now directly assigns to `adult_cutpoint`.

**Deferred:**

4. **`cat()` not `message()` on unknown param** — Line 539. Uses `cat()` instead of `message()`, not gated behind `quietly`. CRAN compliance issue. **Fix:** Change to `message()` (CRAN polish).

5. **`|` vs `||` in `adult_scale_max_lbs` validation** — Line 534: `!is.numeric(adult_scale_max_lbs) | adult_scale_max_lbs < 0`. Should use `||` to short-circuit and avoid evaluating `< 0` on non-numeric input. **Fix:** Change to `||` (minor).

6. **Partial-run `line` and `internal_id` renumbering** — When partial-run filters input vectors before `data.all.ages` construction (lines 497-507), `line = seq_along(measurement)` and `internal_id = seq_len(.N)` are based on the filtered subset, not the original positions. After the partial-run merge (lines 1807-1813), the final `order(as.numeric(internal_id))` sort produces different ordering than a full run. All data is correct but output row order differs from a full run. Low impact since callers can sort by `id` or `line`. **Fix:** Preserve original `line` values for partial-run subjects (could pass original line/internal_id through, or renumber after merge).

### Chunk 4: Outer batching + batch-invariant setup (lines ~558-763)

**Fixed:**
- **CLAUDE.md child exclusion code table was completely wrong.** The table listed ~30 detailed codes (e.g., `Exclude-C-{p}-SDE-EWMA`, `Exclude-C-{p}-EWMA2-middle`, `Exclude-C-{p}-Min-diff`) that don't exist in the code. The actual child algorithm uses simplified suffixes via `.child_exc()`: CF, BIV, Traj, Traj-Extreme, Identical, Extraneous, Abs-Diff, Pair, Single, Too-Many-Errors, Evil-Twins. Updated the table to match the code. The `exclude.levels.peds` factor levels in the code are correct — they include the actual child codes plus legacy codes for backward compatibility.

**Deferred:**

7. **dplyr usage for batching** — Lines 560-565. Uses `%>%`, `select()`, `distinct()`, `mutate()`, `row_number()` instead of data.table equivalents. Inconsistent with the data.table convention used everywhere else. **Fix:** Convert to `unique(data.all.ages[, .(subjid)])` and `[, batch := (.I - 1L) %/% 2000L + 1L]`.

8. **`print("Batching Begin AAA")`** — Line ~557 (after cutpoint fix). Debug print statement, not gated behind `quietly`. **Fix:** Remove or gate behind `if (!quietly) message(...)`.

9. **Hardcoded batch size 2000** — Line ~565. Should be a `batch_size` parameter on `cleangrowth()`. Already noted in CLAUDE.md known issues.

10. **`tanner.fields` and `who.fields` are dead code** — Lines 724-725, 758-759. Computed from velocity table column names but never referenced anywhere (not passed to `cleanchild()` or `cleanlegacy()`). **Fix:** Remove the 4 lines.

11. **Velocity tables loaded for all-adult datasets** — Lines 714-759 load Tanner and WHO velocity tables before checking whether pediatric data exists (line 785). If data is all adults, these are loaded unnecessarily (~0.01 sec, negligible). **Fix:** Move inside the `if (nrow(data.all) > 0)` block (low priority).

### Chunk 5: Batch loop, ped/adult split, imperial/LENGTHCM (lines ~761-885)

No bugs. Batch boundary safety confirmed — subjects assigned to batches at subject level (whole subjects, never split). Imperial conversion factors correct (2.54 cm/in, 2.2046226 lbs/kg). LENGTHCM→HEIGHTCM relabel is a known limitation (no measurement adjustment).

**Deferred:**

12. **Stale `"valid"` in child `var_for_par`** — Line 802. The child algorithm uses `.child_valid()`, not `valid()`. The `"valid"` export sends the legacy function to workers, where it sits unused. **Fix:** Replace with `".child_valid"` (or remove if package loading handles it).

13. **Debug `print()` in batch loop** — Line 764: `print(paste("Loop " , i, " start.", sep = ""))`. Not gated behind `quietly`. **Fix:** Remove or gate.

14. **Stale Stata design comments** — Lines 839-878. Long block of Stata algorithm design notes referencing `exc_*`, `subjid_*` variables. Still informative as algorithm principles but the Stata-specific language is misleading. **Fix:** Condense to algorithm principles only, remove Stata references (low priority).

15. **`cat()` not `message()` on line 834** — "Begin processing pediatric data..." CRAN compliance. **Fix:** Change to `message()`.

### Chunk 6: Z-score calculation + Step 2b GA correction (lines ~885-1469)

No bugs found. Z-score blending boundaries correct and consistent between main calculation (lines 906-927) and `sd.c` calculation (lines 1098-1108). Fenton split-normal calculation correct. Potcorr optimization (skip Fenton I/O when no potcorr subjects) is sound.

**Deferred:**

16. **Harmless bug on line 1173** — `names(table(subjid) > 1)` returns ALL subject names instead of only those with count > 1 (should be `names(which(table(subjid) > 1))`). **Harmless:** single-measurement subjects have `abssumdiff=0` and never trigger the `sd.corr_abssumdiff > sd.orig_abssumdiff` replacement. **Fix:** Change to `names(which(table(subjid) > 1))` for correctness (trivial).

17. **Large commented-out debug code blocks** — Lines ~1383-1410. ~30 lines of commented-out `saveRDS` debug output and `debug_step` logic. Not harmful but clutters the code. **Fix:** Remove (low priority).

18. **Confusing `c(orig_colnames, id)` on line 1229** — `data.all[, colnames(data.all) %in% c(orig_colnames, id), with = FALSE]`. Here `id` is a variable holding the column name `"id"`, but it reads as if `id` is a bare column name string. Confusing but functionally correct since `id` is defined earlier as a parameter. **Fix:** Replace with `c(orig_colnames, "id")` for clarity (trivial).

19. **Multiple `cat()` calls in z-score section** — Lines 887, 897, 1196-1199, 1288, 1306, 1326, 1342, 1346, 1358, 1414, 1426-1434. Same CRAN compliance issue as #4, #13, #15. **Fix:** Change all to `message()`.

20. **Debug `cat()` for uncorr** — Lines 1196-1199. `cat(sprintf("DEBUG uncorr: %d subjects flagged..."))` — debug output gated behind `!quietly`, but labeled "DEBUG" so probably should be removed. **Fix:** Remove or change to informational `message()`.

### Chunk 7: Algorithm dispatch (lines ~1470-1700)

**Deferred:**

21. **`data.adult[, id := line]` overwrites user's `id`** — Line 1659. Sets `id` to `line` (sequential integer) before passing to `cleanadult()`. The user's original `id` is NOT lost — it's restored during output assembly via the `merge(data.batch, full_out, by = "line")` at line 1759, since `data.batch` retains the original `id`. However, the overwrite is confusing and unnecessary now that `cleanadult()` uses `internal_id` for all internal operations. **Fix:** Remove this line; pass `internal_id` through `data.adult` if needed (medium priority — risk of unintended side effects, test carefully).

22. **Non-deterministic adult batching** — Line 1651: `sample(num.batches, ...)` assigns subjects to batches randomly with no seed. With `num.batches = 1` (the default for `parallel = FALSE`), all subjects get batch 1 so this is a no-op. With parallel, batch assignment varies between runs. Child batching is deterministic (sequential). The existing code comment ("is the randomness necessary here?") already flags this concern. **Fix:** Use deterministic assignment like child path: `(seq_along(subjid.unique) - 1L) %/% batch_size + 1L` (low priority since parallel adult is rare).

23. **`data.adult[, measurement := v_adult]` rename** — Line 1658. Creates a `measurement` column from `v_adult` for `cleanadult()`. Not a bug but inconsistent: child path uses `v` (measurement after metric conversion); adult path has its own `v_adult` (the pre-conversion value kept separately at line ~810). The adult imperial conversion happens inside `cleanadult()` itself. This is fine but worth noting for the wrapper narrative.

24. **Missing `Exclude-Missing` for adults before dispatch** — Line 1683 applies `Exclude-Missing` AFTER `cleanadult()` returns (`res[is.na(measurement) | agedays < 0, result := "Exclude-Missing"]`). In the child path, `Exclude-Missing` is set BEFORE dispatch (line 1257). The adult algorithm's Step 1 (BIV) handles NAs by excluding them, so missing values don't cause errors, but the inconsistent placement could confuse readers. **Fix:** Document this asymmetry in the wrapper narrative (low priority).

25. **Commented-out old return block** — Lines 1699-1714. Old return logic that returned only `exclude` as a vector. Fully replaced by the data.table return below. **Fix:** Remove (trivial).

26. **`print()` in batch loop completion** — Line 1785: `print(paste("Loop ", i, " complete.", sep = ""))`. Debug output not gated behind `quietly`. Same class of issue as #8 and #13. **Fix:** Remove or gate behind `if (!quietly) message(...)`.

### Chunk 8: Output assembly + partial-run merge (lines ~1718-1854)

No bugs found. Output assembly logic is correct: child and adult results combined via `line`, merged with `data.batch` to restore original columns, checkpoint diagnostics merged by `id`, then sorted by `line`. Partial-run merge correctly excludes unchanged cached subjects by `subjid`, appends new results, and sorts by `internal_id`. Safety check for `Exclude-C-Temp-Same-Day` surviving to output is good defensive practice. `bin_exclude` and `tri_exclude` derivation is clean and correct.

**Deferred:**

27. **`checkpoint_data` merge by `id` may be fragile** — Lines 1768-1772. Merges checkpoint columns by `id`. For child rows, `id` is the user's original; for adult rows, `id` was overwritten to `line` (see #21). Since `checkpoint_data` only exists for the child algorithm path, adult rows just get NAs, so this works. But if checkpoint_data were ever extended to adults, it would silently fail. **Fix:** Document this dependency in code comment (trivial).

28. **`match()` reordering** — Line 1774: `all_results <- all_results[match(data.batch$id, all_results$id)]`. If duplicate IDs exist in `all_results` after the merge (shouldn't happen, but defensively), `match()` returns only the first match. No issue in practice since `id` is required to be unique. **Fix:** None needed, but a comment would clarify the uniqueness assumption.

### Support functions: read_anthro, gc_preload_refs, ewma, calc_and_recenter_z_scores

**read_anthro** (lines 1873-2143): Closure-based z-score calculator. Loads reference tables once, returns a function that computes z-scores by joining against the loaded data. Two modes: LMS z-scores (`csd = FALSE`) and CSD z-scores (`csd = TRUE`). WHO/CDC source selection inside closure uses `agedays < 5*365.25` boundary — consistent with the external blending window. No bugs.

**gc_preload_refs** (lines 2177-2183): Thin wrapper — calls `read_anthro` three times with different arguments. No issues.

**ewma** (lines 2223-2311): Core EWMA calculation. Uses matrix multiply for O(n²) initial computation. Age-dependent exponents applied row-wise via `sweep()`. Windowing (default 15 positions each side) zeroes out weights beyond window. Delta matrix caching via `cache_env` avoids recomputation when tbc and ctbc share same agedays. No bugs.

**ewma_cache_init/update** (lines 2320-2478): Incremental EWMA cache. `_init` does O(n²) first build; `_update` does O(n) subtraction when one observation is excluded. Neighbor exponent recalculation checks 2 positions each side (bug fix noted in code). Handles `skip_ctbc` optimization when tbc == ctbc. No bugs.

**calc_and_recenter_z_scores** (lines 2827-2891): Recalculates z-scores for replacement values (used in Steps 3, 4). Same WHO/CDC blending logic as main calculation.

**sd_median** (lines 2509-2546): Derives recentering medians from input data. Groups ages > 19 together. Uses `approx()` with `rule = 2` for extrapolation.

**Deferred:**

29. **`calc_and_recenter_z_scores` missing NA fallbacks** — Lines 2869-2875. The main z-score calculation (lines 920-927) includes fallback logic for smooth_val rows where one z-score source is NA: `(smooth_val & is.na(sd.orig_cdc))` falls back to WHO, and vice versa. `calc_and_recenter_z_scores` omits these fallbacks. For standard reference data this never triggers (WHO and CDC both cover 2-5y), but it's less robust. **Fix:** Add matching NA fallback logic for consistency (low priority).

30. **`as_matrix_delta` defined in legacy support file** — `pediatric_support_legacy.R` line 167. Both child EWMA (`ewma()`) and legacy algorithm use this. Works because R loads all `R/` files. Not a bug, but if the legacy file were ever removed, child EWMA would break. **Fix:** Move to `child_clean.R` or a shared utilities file (low priority).

31. **`sd.orig_who` can contain CDC data for ages ≥ 5y** — The WHO closure (`cdc.only = FALSE`) returns CDC z-scores for ages ≥ 5*365.25. The variable name `sd.orig_who` is therefore misleading for those rows. The external blending gives zero weight to WHO at age 5+ and `who_val` excludes age ≥ 2y (for non-HC), so the CDC values in `sd.orig_who` are never incorrectly used as WHO. Cosmetic only. **Fix:** Add comment explaining this (trivial).

---

## Child algorithm steps

### cleanchild() entry (lines 3019-3108)

No bugs. Entry re-creates `index` (correct — outer batching causes non-contiguous indices). Logging setup correct.

### Early Step 13: SDE-Identicals (lines 3109-3143)

No bugs. Age-dependent id selection (lowest at birth, highest otherwise). Correct guard against empty Include subsets in `keep_id` calculation. Cleanup of temp columns.

### Step 5: Temporary SDEs (lines 3145-3166)

No bugs. `temporary_extraneous_infants()` called with correct column subset including `internal_id`. Safety check for duplicate Includes after resolution uses `warning()` (correct).

### Step 6: Carried Forwards (lines 3168-3505)

No bugs. Complex step with multiple sub-phases:
- CF detection (pre-filtered by duplicate values for efficiency)
- Positional string detection (matches Stata: originators are Includes where next is CF)
- Forward propagation of string numbers and originator z-scores via iterative `shift()`
- CF rescue with 3 modes: "none", "all", "standard" (lookup tables)
- SDE-Identicals correctly removed before CF string calculation and restored after
- `wholehalfimp` flag for imperial rounding detection
- `cf_detail` output populates `cf_status` and `cf_deltaZ`

### Step 7: BIV (lines 3507-3607)

No bugs. Absolute limits (age-dependent for weight), then standardized z-score limits using `sd.orig_uncorr` (unrecentered). Correctly uses `!grepl(biv_pattern, exclude)` to avoid overwriting absolute BIV codes with standardized BIV. Temp SDEs redone after.

### Step 9: Evil Twins (lines 3609-3716)

No bugs. Restructured from single-closure to per-group processing. OTL calculation via `calc_otl_evil_twins()`. While-loop removes worst OTL one at a time with tiebreaker hierarchy (med_diff, abs(tbc.sd), internal_id). Correctly handles `otl` column removal before re-calling `calc_otl_evil_twins()` to avoid self-assignment.

### Step 11: EWMA1 Extreme (lines 3718-3891)

No bugs. Global iteration structure — only re-processes subject-params with new exclusions each iteration. Pre-filter: skip groups where no Include has |tbc.sd| > 3.5. Age-dependent exponents correct (linear interpolation -1.5 to -3.5 over 1-3 years). Both tbc.sd and ctbc.sd EWMA computed (with ctbc skip optimization). First/last values protected. EWMA1 iteration 1 values saved for debugging.

**Deferred:**

32. **EWMA1 debug columns saved unconditionally** — Lines 3844-3855. `ewma1_it1.*` columns are always saved for debugging, even in production. These add 6 columns per row. **Fix:** Gate behind a `debug` parameter or remove (low priority).

### Step 13: Final SDE (lines 3892-4180)

No bugs. Complex SDE resolution:
- Phase B2: One-day SDEs (all measurements on a single day)
- Phase B3: SDE-EWMA resolution using EWMA from non-temp-SDE values
- Phase B4: Merge results back to main data (correctly preserves permanent exclusions)

**Deferred:**

33. **`suppressWarnings(min/max)` + `is.infinite()` pattern** — Lines 4003-4007, 4106-4109, 4114-4116, 4123-4128. Same pattern as noted in CLAUDE.md robustness audit. Works correctly but masks unexpected warnings. Matches the 4 locations already documented; no new instances.

34. **Velocity tables re-read from disk per batch — by design** — Lines 4627-4668. `who_max_ht_vel`, `who_ht_vel_3sd`, `who_max_hc_vel`, `who_hc_vel_3sd` are read from disk inside `cleanchild()`. The raw WHO reference tables are static, but Step 17 merges them against the current data based on age intervals between *included* measurements. As measurements are excluded during cleaning (Steps 7, 9, 11, 13, 15, 16), the intervals between remaining included measurements change, so the velocity table merge must be recomputed with each batch's current inclusion state. The Tanner table is passed as a parameter from `cleangrowth()` but the WHO tables are read inside `cleanchild()` — both are re-merged against the data fresh within Step 17's while loop. **The disk reads themselves could potentially be cached** (read once in `cleangrowth()`, pass as parameters like Tanner), but the merge/application step cannot be cached because it depends on the current exclusion state. **Status:** Deferred — the disk read cost is negligible for typical datasets (1-2 batches). A future optimization could pre-load the WHO velocity tables once and pass them as parameters, analogous to `ref_tables`/`gc_preload_refs()`, while still re-merging per iteration.

### Step 15: EWMA2 Moderate (lines 4182-4461)

No bugs. Global iteration structure (same pattern as Step 11). Pre-calculates p_plus/p_minus and their z-scores once upfront. Incremental EWMA cache (`ewma_cache_init`/`ewma_cache_update`) persists across iterations. Complex exclusion rules for middle/first/last/birth-WT with addcrit thresholds. DOP lookup uses keyed snapshot table for O(log n) lookups.

### Step 16: Birth HT/HC (lines 4464-4604)

No bugs. Separate step for birth HT/HC (excluded from Step 15's filter). Same EWMA cache pattern as Step 15. Only processes subjects with birth measurements. Slightly simplified rules (birth only, no first/last/middle).

### Step 17: Raw Differences (lines 4613-5223)

No bugs. Complex velocity check with Tanner (HT ≥2.5y) and WHO (HT <2y, HC) reference data. Pre-filter identifies groups with violations before the expensive per-group while loop. Correct exponent formulas for Tanner: `^0.33` for <1y, `^1.5` for >1y. HC has separate tolerance (±1.5 cm vs ±3 cm for HT). Birth adjustments applied. While-loop removes one violation per iteration with DEWMA-based tiebreaking. Two-value pairs handled separately with simpler tbc.sd comparison.

### Step 19: 1 or 2 Measurements (lines 5225-5341)

No bugs. Counts only Include rows for singles/pairs determination (Stata match). DOP snapshot taken BEFORE by-group processing to avoid ordering dependency. Pair criteria: |diff| > 4 with ≥365 days or |diff| > 2.5 with <365 days. Single criteria: |tbc.sd| > 3 with comp_diff > 5, or |tbc.sd| > 5 with no DOP.

### Step 21: Error Load (lines 5343-5385)

No bugs. Error ratio = errors / (errors + includes), excluding SDE/CF/Missing from both numerator and denominator. Uses configurable `error.load.threshold` (default 0.5) and `error.load.mincount` (default 2). The hardcoded `.4` bug was previously fixed.

### Step 22: Output (lines 5387-5441)

No bugs. Creates `final_tbc` (ctbc.sd for potcorr, tbc.sd otherwise). Returns specified columns including z-scores, EWMA1 debug values, and cf_detail if requested.

**Deferred:**

35. **Unused variables computed in Step 19** — Already noted in CLAUDE.md known issues (`abs_tbd.sd`, `abs_ctbd.sd`, `med_dop`, `med_cdop`). **Confirmed resolved:** grep for these variables finds no matches in current `child_clean.R`. Remove this item from CLAUDE.md known issues.

36. **Commented-out old `valid()` function** — Lines 5444-5470. Old version of `.child_valid()` still present as a comment block. **Fix:** Remove commented-out code (trivial).

37. **Multiple `cat()` calls throughout `cleanchild()`** — Lines 3101, 3117, 3151, 3178, 3199, 3257, 3414, 3423, 3456, 3514, 3649, 3737, 3754, 3764, 3894, 3933, 4191, 4264, 4274, 4467, 4492, 4502, 4615, 4899, 5227, 5346, 5390. All gated behind `!quietly`. Same CRAN compliance issue as wrapper — should be `message()`. **Fix:** Bulk change (low priority, CRAN polish).

---

## Adult algorithm steps

Adult algorithm is **closed pending clinician validation** (per CLAUDE.md). Walk-through was completed on 2026-04-09 with findings in `walkthrough-todo-2026-04-09.md`. Not repeating here. See that file for deferred items.

---

## Post-walkthrough: test suite verification

All tests pass at baseline counts — no regressions from walkthrough.

| Test file | Result | Baseline |
|---|---|---|
| test-cleangrowth.R | 92 PASS | 92 PASS |
| test-child-regression.R | 54 PASS | 54 PASS |
| test-child-edge-cases.R | 28 PASS | 28 PASS |
| test-child-algorithms.R | 40 PASS, 2 WARN | 40 PASS, 2 WARN |
| test-child-parameters.R | 21 PASS, 1 WARN | 21 PASS, 1 WARN |
| test-adult-clean.R | 198 PASS | 198 PASS |
| Adult regression (looser) | 1508/1508 (100%) | 1508/1508 |

**Walk-through complete.**
