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

### D2. `Exclude-SDE-All-Exclude` — never assigned

**Issue:** This exclusion code is in `exclude.levels` (line 573), in the `child_sde_codes` list (line 1786), and in the `non_error_codes` error load list (line 5171), but no line of code ever assigns it to any row. It's listed in both the narrative exclusion code table and CLAUDE.md.

**Question for Carrie:** Is this code supposed to be assigned somewhere (e.g., when all same-day values were already excluded by prior steps), or is it dead and should be removed from the tables? The narrative lists it as a Step 13 code.

### D3. `Exclude-SDE-EWMA-All-Extreme` — never assigned

**Issue:** This exclusion code is in `exclude.levels` (line 620) and in `child_sde_codes` (line 1788), but never assigned. The EWMA-based SDE-All-Extreme check at line 3940 uses `Exclude-SDE-All-Extreme` (without "EWMA" in the name).

**Question for Carrie:** Should this code replace `Exclude-SDE-All-Extreme` at line 3940 (to distinguish EWMA-based from one-day), or is it dead? Also not in `non_error_codes` — if it were ever assigned, it would incorrectly count as an error in error load.

### D4. `Not cleaned` not in error load `non_error_codes`

**Issue:** HC rows with agedays > 3 × 365.25 are marked `Not cleaned` (line 1218). In Step 21 error load, `Not cleaned` is not in the `non_error_codes` list (lines 5170-5176), so these rows would count as errors in the error load ratio.

**Question for Carrie:** Should `Not cleaned` be in `non_error_codes`? It's an age-based scope limitation, not a cleaning error. (The `valid()` function already excludes `Not cleaned` rows, so they wouldn't be counted as Include either — they'd only affect the error count in the numerator.)

### D5. `height.tolerance.cm` — declared but never passed or used

**Issue:** `height.tolerance.cm` (default 2.5) is declared as a parameter of `cleangrowth()` (line 321) and documented in the narrative (line 848), but it is never passed to `cleanbatch_child()` and never referenced in Step 17's velocity logic. Step 17 uses hardcoded velocity reference tables instead.

**Question for Carrie:** Was this parameter intended to be used in Step 17 as a mindiff override? Or is it a legacy parameter from the Stata version that should be removed from the interface?

### D6. `lt3.exclude.mode` — passed but never referenced

**Issue:** `lt3.exclude.mode` (default "default") is passed to `cleanbatch_child()` (line 1444, 1462, 1496, 1519) and accepted as a parameter (line 2876), but never referenced in Step 19 code. The "flag.both" mode described in the roxygen documentation is not implemented.

**Question for Carrie:** Should "flag.both" mode be implemented in Step 19, or should this parameter be removed (or documented as not yet implemented)?

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

### Code fixes (2026-04-13 session 1 — narrative only saved)

- Updated `calc_and_recenter_z_scores()` section to reflect bug is fixed
- Updated Open Questions: marked batching wrapper, `stop()`→`warning()`, parallel as fixed
- Updated `id` → `internal_id` references throughout step descriptions
- Updated Step 6 checklist, wholehalfimp description, Step 9 checklist
- Updated Step 15 table: added c.dewma.all threshold notes
- Updated Step 17 checklist: removed stale finding #10
- Noted Step 5 input includes `internal_id`

### Code fixes (2026-04-13 session 2 — code + narrative saved)

- **D1 FIXED:** `id` → `internal_id` (with `as.numeric()`) at all 4 SDE-Identical tiebreaking locations (GA correction preprocessing, Early Step 13, Main Step 13 whole-day, Main Step 13 partial)
- **D4 FIXED:** Added `"Not cleaned"` to `non_error_codes` in Step 21 error load
- **D5 FIXED:** Removed `height.tolerance.cm` parameter (roxygen + declaration)
- **#3 FIXED:** Removed stale Stata comment `# h. no need to calculate for R.` in Step 13
- **#4 FIXED:** `setkey(data.sde, subjid, param, agedays, id)` → `internal_id` in Step 13
- **Stale CP markers FIXED:** Replaced 3 stale CP comments in Step 17 with descriptive comments
- **Dead code FIXED:** Removed unreachable inner if block in Step 22 final_tbc

---

## Remaining items

### Still to do (from Carrie's review)

1. **#5: Move batch-invariant operations before the loop** — `exclude.levels` definition, velocity reference reads, `read_anthro()` calls are repeated per batch. Moving before the outer loop would avoid redundant re-computation. IN PROGRESS — was interrupted by API errors.

2. **#6: Remove line numbers from narrative** — All narrative line numbers are stale. Better absent than wrong.

3. **#7: Add `internal_id` to narrative input column lists** — Multiple step descriptions list input columns without `internal_id`.

### Needs discussion with Carrie

4. **D2: `Exclude-SDE-All-Exclude`** — never assigned. Dead code or missing implementation? Carrie needs to review Step 13 in detail.

5. **D3: `Exclude-SDE-EWMA-All-Extreme`** — never assigned. Should it replace `Exclude-SDE-All-Extreme` at the EWMA check, or is it dead? Also not in `non_error_codes`. Carrie needs to review Step 13.

6. **D6: `lt3.exclude.mode`** — passed to `cleanbatch_child()` but never referenced. "flag.both" mode not implemented. Needs discussion.

7. **#2: `v == 0` handler in Step 7** — applies to ALL rows (no valid_set filter). In practice, `v == 0` is converted to NaN during preprocessing so this should never fire. Carrie wants to discuss what `v` does.

8. **#8: Step 15 exclusion rule table** — missing `c.dewma.all` column. For `last-high` and `last-ext-high`, `c.dewma.all` is hardcoded at 3. Needs discussion.

### Code cleanup (low priority, no behavioral impact)

9. **Dead exclusion codes in `exclude.levels`**: `Exclude-SDE-All-Exclude` and `Exclude-SDE-EWMA-All-Extreme` — depends on D2/D3 decisions.
