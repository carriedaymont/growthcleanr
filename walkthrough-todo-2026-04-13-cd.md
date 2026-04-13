# Child GC Walk-Through — Deferred To-Dos

Items identified during the 2026-04-13 walk-through that are deferred for later resolution.

---

## Pre-walk-through setup

Baseline tests all pass:
- Regression: 54 pass
- Edge cases: 23 pass, 8 warnings
- Parameters: 24 pass
- Algorithms: 6 known failures (unimplemented `cf_rescue_threshold`), 15 pass
- Stress: 38 pass, 3 warnings

---

## Discrepancies requiring Carrie's review

Items where code, comments, and narrative disagree and the correct behavior is unclear. These are NOT assumed to be bugs — the code, comments, or narrative could each be wrong.

### D1. `id` vs `internal_id` in SDE-Identical tiebreaking (4 locations)

**Issue:** Four SDE-Identical tiebreaking locations use `min(id)`/`max(id)` instead of `min(internal_id)`/`max(internal_id)`. All other tiebreakers (Step 5, Step 13 one-day SDE, Step 13 SDE-EWMA, Step 9, Step 11, Steps 15/16) correctly use `internal_id`. 

If `id` is numeric and sequential (as in test data), behavior is identical. If `id` is character, UUID, or non-sequential, `min(id)`/`max(id)` could select the wrong row.

**Locations:**
- Line 1101-1103: GA correction preprocessing SDE-Identical filter
- Line 2974-2976: Early Step 13 SDE-Identicals
- Line 3763-3765: Main Step 13 SDE-Identical (whole-day)
- Line 3773-3775: Main Step 13 SDE-Identical (partial)

**Question for Carrie:** Should these use `min(as.numeric(internal_id))`/`max(as.numeric(internal_id))` instead? The sort at each location already uses `internal_id`, suggesting the keeper selection should too. (The 2026-04-12 `internal_id` rename may have missed these.) 

**Answer**: Yes, these should all be internal_id and as.numeric. We must have missed those previously.

### D2. `Exclude-SDE-All-Exclude` — never assigned

**Issue:** This exclusion code is in `exclude.levels` (line 573), in the `child_sde_codes` list (line 1786), and in the `non_error_codes` error load list (line 5171), but no line of code ever assigns it to any row. It's listed in both the narrative exclusion code table and CLAUDE.md.

**Question for Carrie:** Is this code supposed to be assigned somewhere (e.g., when all same-day values were already excluded by prior steps), or is it dead and should be removed from the tables? The narrative lists it as a Step 13 code.

**Answer**I need to review this step in more detail. 

### D3. `Exclude-SDE-EWMA-All-Extreme` — never assigned

**Issue:** This exclusion code is in `exclude.levels` (line 620) and in `child_sde_codes` (line 1788), but never assigned. The EWMA-based SDE-All-Extreme check at line 3940 uses `Exclude-SDE-All-Extreme` (without "EWMA" in the name).

**Question for Carrie:** Should this code replace `Exclude-SDE-All-Extreme` at line 3940 (to distinguish EWMA-based from one-day), or is it dead? Also not in `non_error_codes` — if it were ever assigned, it would incorrectly count as an error in error load.

**Answer**I need to review this step in more detail. 

### D4. `Not cleaned` not in error load `non_error_codes`

**Issue:** HC rows with agedays > 3 × 365.25 are marked `Not cleaned` (line 1218). In Step 21 error load, `Not cleaned` is not in the `non_error_codes` list (lines 5170-5176), so these rows would count as errors in the error load ratio.

**Question for Carrie:** Should `Not cleaned` be in `non_error_codes`? It's an age-based scope limitation, not a cleaning error. (The `valid()` function already excludes `Not cleaned` rows, so they wouldn't be counted as Include either — they'd only affect the error count in the numerator.)

**Answer**: Yes

### D5. `height.tolerance.cm` — declared but never passed or used

**Issue:** `height.tolerance.cm` (default 2.5) is declared as a parameter of `cleangrowth()` (line 321) and documented in the narrative (line 848), but it is never passed to `cleanbatch_child()` and never referenced in Step 17's velocity logic. Step 17 uses hardcoded velocity reference tables instead.

**Question for Carrie:** Was this parameter intended to be used in Step 17 as a mindiff override? Or is it a legacy parameter from the Stata version that should be removed from the interface?

**Answer**: I think this was a new value that we were going to add, but this is going to be saved for a later version. Please remove. We'll leave this hardcoded for now.

### D6. `lt3.exclude.mode` — passed but never referenced

**Issue:** `lt3.exclude.mode` (default "default") is passed to `cleanbatch_child()` (line 1444, 1462, 1496, 1519) and accepted as a parameter (line 2876), but never referenced in Step 19 code. The "flag.both" mode described in the roxygen documentation is not implemented.

**Question for Carrie:** Should "flag.both" mode be implemented in Step 19, or should this parameter be removed (or documented as not yet implemented)?

**Answer**: I need to review this more. Let's discuss

---

## Fixes applied during walkthrough

### Narrative fixes

- Updated `calc_and_recenter_z_scores()` section to reflect that the bug is fixed
- Updated Open Questions: marked batching wrapper as fixed, marked `stop()` → `warning()` as fixed, marked parallel processing as fixed
- Updated references from `id` to `internal_id` throughout step descriptions
- Updated Step 6 checklist: removed stale findings (#1 old comments, #2 cf_string_length, #3 uncleaned variables, #4 redundant wholehalfimp)
- Updated Step 6: corrected wholehalfimp description (only mod 0.5, not mod 1 OR mod 0.5)
- Updated Step 9 checklist: removed stale findings (#1 rounded aliases, #2 NNTE comment)
- Updated Step 15 table: added c.dewma.all threshold notes for last-high and last-ext-high
- Updated Step 17 checklist: removed stale finding #10 about unused variables
- Noted that Step 5 input includes `internal_id` column

### Code fixes

- Removed stale CP markers in Step 17 (lines 4819, 4825, 4830)
- Removed dead code in Step 22 final_tbc (unreachable inner if block)

---

## Deferred items

### Code cleanup (no behavioral impact)

1. **Dead exclusion codes in `exclude.levels`**: `Exclude-SDE-All-Exclude` and `Exclude-SDE-EWMA-All-Extreme` are defined but never assigned. See D2/D3 above — depends on Carrie's decision.

**From Carrie**: We'll discuss


2. **`v == 0` handler in Step 7 (line 3405)**: Applies to ALL rows (no valid_set filter) and could theoretically overwrite a `Missing` code. In practice, `v == 0` is converted to NaN during preprocessing so this should never fire. Consider whether this is dead code.

**From Carrie**: Please discuss with me what v does.

3. **Stale Stata reference in Step 13**: Line 3873 comment `# h. no need to calculate for R.` — cryptic reference to a Stata step. Harmless.

**From Carrie**: Please fix

4. **Line 3785**: `setkey(data.sde, subjid, param, agedays, id)` — should this be `internal_id` instead of `id`? This is the setkey after SDE-Identical resolution in Step 13, before temp SDE tracking. Using `id` for sorting could produce non-deterministic order with character ids.

**From Carrie**: Yes, always internal_id


5. **Batch-invariant operations inside loop**: `exclude.levels` definition (lines 559-621), velocity reference reads (lines 743-793), `read_anthro()` calls (lines 847-852) are repeated per batch. Moving before the loop would improve performance for many-batch runs. Not bugs.

**From Carrie**: Please do this

### Documentation completeness (no behavioral impact)

6. **All narrative line numbers are stale** by ~230-300 lines due to code changes since 2026-03-18. A systematic line number update would improve usability but is low priority since code search works fine.

**From Carrie**: Please remove line numbers. Better absent than wrong.


7. **Narrative doesn't mention `internal_id` in input column lists**: Step 5 and other steps now take `internal_id` as input (added 2026-04-12). Multiple narrative descriptions of input columns are incomplete.

**From Carrie**: Please add to narrative


8. **Step 15 exclusion rule table missing c.dewma.all column**: The table shows dewma threshold and extra condition but doesn't have a column for `c.dewma.all` thresholds. For most rules c.dewma.all matches dewma, but for `last-high` and `last-ext-high` it's hardcoded at 3 (not matching the dewma threshold).

**From Carrie**Let's discuss this more
