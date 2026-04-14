# Adult GC Walk-Through — Deferred To-Dos

Items identified during the 2026-04-03 walk-through that are deferred for later resolution.

---

## Pre-walk-through setup

### `ewma_cache_init` name collision (pre-existing, not introduced by walk-through)

**File:** `R/adult_support.R` and `R/child_clean.R`

**Issue:** Both files define `ewma_cache_init` with different signatures:
- Adult (`adult_support.R`): `ewma_cache_init(agedays, meas, ewma_exp = -5, ewma_window = 15)`
- Child (`child_clean.R`): `ewma_cache_init(agedays, z_tbc, z_ctbc, exp_vals, ids, window = 15)`

R loads files alphabetically; `child_clean.R` loads after `adult_support.R` and overwrites the adult version. When adult code calls `ewma_cache_init(..., ewma_window = ewma_window)`, it hits the child version, which does not accept `ewma_window`. This causes 8 test failures in `test-adult-clean.R`.

**Tests affected:** 8 failures involving `remove_mod_ewma_wt` (moderate EWMA step), specifically tests that exercise tightest/loosest comparison, `mod_ewma_f`, and `allofus15` formula.

**Fix needed:** Rename the adult version (e.g., `adult_ewma_cache_init`) in `adult_support.R` and update all adult call sites. Update `child_clean.R`'s parallel export string list accordingly.

**Fix:** Renamed adult versions to `adult_ewma_cache_init` / `adult_ewma_cache_update` in `adult_support.R`, updated 3 call sites, added new names to `child_clean.R` clusterExport list. All 208 adult tests now passing. **RESOLVED.**

---

## Steps walk-through (populated as walk-through proceeds)

### Steps 1H/1W (BIV) and 2W (Repeated Values)

**4. Bug: BMI BIV check not implemented**

`overall_bmi_min` and `overall_bmi_max` are resolved from presets and unpacked into local variables but never used. `biv_df` only contains height and weight rows — no BMI row is constructed and no BMI comparison is executed in Step 1. For looser/tighter/tightest, `overall_bmi` equals `single_bmi` so there is no practical impact. For loosest, `overall_bmi_min = 5` vs `single_bmi_min = 10`, so BMI values 5–10 are not caught until Step 13.

Fix requires: construct same-day ht+wt pairs, compute BMI, compare against `overall_bmi_min`/`overall_bmi_max` with no tolerance (per rounding tolerance spec), exclude the pair if BMI is outside limits.

Also update: narrative configurable parameters table (add clarifying note), permissiveness spec (note that BMI BIV check is applied at Step 1, add implementation detail).

**Priority:** Medium — only affects loosest, and Step 13 catches it later. But the code and docs claim it happens in Step 1.

**Fix:** Added Step 1B BMI BIV check block in `adult_clean.R` after 1W, before 2W. Guard condition: `overall_bmi_min < single_bmi_min || overall_bmi_max > single_bmi_max` — only runs at loosest (overall BMI [5,300] vs single [10,250]); skipped at other levels where the two limits are equal (applying the check earlier than intended would disrupt later steps). Pairs first ht and first wt (by id) on each shared ageday; excludes both with `Exclude-A-HT-BIV`/`Exclude-A-WT-BIV` if BMI is outside `[overall_bmi_min, overall_bmi_max]` (no tolerance on BMI). Four tests added (low BMI excluded, high BMI excluded, valid pair not excluded, different-day pair not checked). Narrative updated: scope row, parameter table notes, logic steps including guard explanation. 220/220 unit tests passing, 1483/1483 regression tests at all 4 permissiveness levels. **RESOLVED.**

---

### Pre-step sections (Key Concepts, Rounding Tolerance, Permissiveness Framework)

**1. Key concepts stubs — narrative**

The Key Concepts section has three blank stub entries: EWMA, Iterations, and Exclusion rescue. Content needed; deferred for user to fill in.

**Fix:** Filled all three stubs: EWMA describes the weighted average, exponent e=−5, three variants (all/before/after), window parameter, and edge-value fallback. Iterations describes the single-worst-out-per-round dynamic loop pattern used throughout. Exclusion rescue describes the moderate EWMA (11Wb) error-load escalation as the one rescue mechanism in the adult algorithm. **RESOLVED.**

**2. Permissiveness spec: wtallow tables incomplete**

The wtallow tables section shows the loosest/looser (piecewise) table only. Tables for `piecewise-lower` (tighter) and `allofus15` (tightest) are missing. Defer until those formulas are finalized and stable.

**Fix:** Updated all three formula descriptions with corrected detail (piecewise 0–1m segment is log-curved, not linear; piecewise-lower now includes starting values and 12-month formula) and added breakpoint tables for all three formulas. **RESOLVED.**

**3. Rounding tolerance — Step 4W dewma_bef/aft detail — RESOLVED**

Fixed during Step 4W walk-through. Rounding tolerance row now reads: "`dewma_all` vs 50 kg; `dewma_before`/`dewma_after` vs 45 kg (0.9 factor applied to 50 kg base); adjacency diffs vs 50 kg."

---

### Step 4W (Scale Max / Weight Cap)

**4. Step 4W: Defensive NA→Inf for dewma_bef/aft are dead code**

In `adult_clean.R`, lines ~365–368, four `ifelse(is.na(...), ±Inf, ...)` conversions guard against NA in `dewma_bef` and `dewma_aft`. However, `adult_ewma_cache_init()` never produces NA for before/after EWMA — edge cases (first/last observation) fall back to `ewma_all`. The safe conversions are harmless but unused.

**Priority:** Low — no behavioral impact. Address if cleaning up dead code before CRAN.

**Fix:** Removed the four `ifelse(is.na(...), ±Inf, ...)` conversions and their associated comment. `dewma_bef`/`dewma_aft` used directly. Added a brief comment explaining why NA guards are not needed. **RESOLVED.**

---

### Step 9Wa (Evil Twins)

**5. Step 9Wa: Redundant pairs guard inside `evil_twins()`**

In `adult_support.R`, `evil_twins()` checks `n < 3` at the top of each round (line ~748, breaks early), then checks `n_inc < 3` again at line ~790 after computing deviations. No rows are removed between those two checks, so the second check is always redundant.

**Priority:** Low — no behavioral impact. Address if cleaning up dead code before CRAN.

**Fix:** Removed redundant `n_inc < 3` guard as part of item 6 fix. **RESOLVED.**

---

**6. Bug: Evil twins uses 3 fixed rounds instead of dynamic iteration**

`evil_twins()` uses `for (round in 1:3)` — exactly 3 rounds regardless of whether OOB pairs remain after round 3. The EWMA steps (`remove_ewma_wt`, `remove_mod_ewma_wt`) correctly use `while (change)` and iterate until no candidates remain.

With only 3 rounds, a patient with 4+ evil twin pairs will have the remaining pairs missed by this step; they may or may not be caught by later steps (EWMA, 2D).

**Fix needed:** Replace `for (round in 1:3)` with a `while`-based loop that breaks when no OOB pairs are found or `n < 3`. Pattern to follow: `remove_ewma_wt()` in `adult_support.R`.

Also update: narrative (remove "3 fixed rounds" language; describe dynamic iteration with early exit when no OOB pairs remain).

**Priority:** Medium — affects patients with many large weight swings. Straightforward fix.

**Fix:** Changed `for (round in 1:3)` to `while (TRUE)` in `evil_twins()`. Removed redundant `n_inc < 3` guard (items 5 and 6 addressed together). Added test with 9 alternating weights (4 OOB pairs) to verify all 4 excluded. Narrative updated (3 locations). All tests passing. **RESOLVED.**

---

### Step 10H (Height Distinct Values)

**7. Step 10H: 3D re-sort missing id tiebreaker**

In `adult_clean.R`, the 10Hb branch re-sorts `h_subj_df` at line ~716 with `order(h_subj_df$ageyears)` — no id tiebreaker. In practice safe because Step 9H resolves same-day duplicates before reaching 10Hb. Could become non-deterministic if two height measurements on different days map to identical floating-point `ageyears` values.

**Fix needed:** `h_subj_df <- h_subj_df[order(h_subj_df$ageyears, as.numeric(h_subj_df$id)), ]`

**Priority:** Low — no known behavioral impact. Address as part of general sort determinism cleanup.

**Fix:** Added `as.numeric(h_subj_df$id)` as tiebreaker. Consistent with all other sorts in the file. **RESOLVED.**

---

### Step 14 (Error Load)

**8. Step 14: Dead codes in `eval_error_load()` `sde_codes` vector**

In `adult_support.R`, `eval_error_load()` defines:
```r
sde_codes <- c("Same-day-identical", "Same-day-extraneous",
               "Exclude-A-HT-Identical", "Exclude-A-HT-Extraneous",
               "Exclude-A-WT-Identical", "Exclude-A-WT-Extraneous")
```
`"Same-day-identical"` and `"Same-day-extraneous"` are not produced by the adult algorithm (they are child algorithm artifacts). These two entries never match any result value and can be removed.

**Priority:** Low — no behavioral impact.

**Fix:** Removed `"Same-day-identical"` and `"Same-day-extraneous"` from `sde_codes`. **RESOLVED.**

---

**9. Step 14: `Exclude-A-WT-Scale-Max-Identical` in `error_codes_ht`**

In `adult_support.R`, `eval_error_load()` includes `"Exclude-A-WT-Scale-Max-Identical"` in `error_codes_ht`. This is a weight-only exclusion code that can never appear on a height row, so it never matches in the height pass. Dead code. Narrative corrected to remove it from the "Height errors" list.

**Priority:** Low — no behavioral impact.

**Fix:** Removed `"Exclude-A-WT-Scale-Max-Identical"` from `error_codes_ht`. **RESOLVED.**
