# GC Full Walk-Through — Deferred To-Dos

Items identified during the 2026-04-16 (second) walk-through that are deferred for later resolution.

---

## Pre-walk-through setup

All tests pass. Baseline counts:
- test-cleangrowth.R: 92 PASS
- test-child-regression.R: 50 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS (2 warnings expected)
- test-child-parameters.R: 21 PASS (1 warning expected)
- test-adult-clean.R: 198 PASS
- Adult regression: 1508/1508 at all 4 levels

---

## Fixed during walkthrough

1. **Adult `evil_twins()` sort determinism** (`adult_support.R` lines 865, 874, 913): Three `order()` calls used character `internal_id` without `as.numeric()`, causing lexicographic sort ("10" < "2"). Fixed by adding `as.numeric()` to all three.

2. **Header DOP comment** (`child_clean.R` line 62): Said "height/HC→weight" but HC's DOP is height, not weight. Fixed to "height→weight, HC→height".

3. **Stale `### CP REPLACE D/U` markers** (`child_clean.R` lines 2753, 2767): Dead commented-out code and stale marker tags in `calc_and_recenter_z_scores()`. Removed.

4. **Stale `# Oriignal Valid Toggled off` comment** (`child_clean.R` line 5334): Typo and meaningless. Removed.

---

## Wrapper walkthrough

No deferred items. Wrapper code is clean after prior walkthrough fixes.

---

## Child algorithm steps

### RESOLVED: `setkey(data.df, ..., id)` vs `internal_id` consistency
- **File/lines:** `child_clean.R` lines 2975, 3035, 4077 (and others)
- **Issue:** `setkey` used user's `id` for tiebreaking, while `order()` calls used `as.numeric(internal_id)`.
- **Resolution (2026-04-16):** Changed `internal_id` to integer type, assigned after `setkey(..., id)` so it reflects id-sorted order. All `setkey()` calls now use `internal_id` instead of `id`. Adult output joins use `as.integer(names(...))` for type matching.

### Defer: Step 17 velocity disk reads per batch
- **File/lines:** `child_clean.R` lines 4524-4565
- **Issue:** WHO HT and HC velocity CSV files are read from disk every time `cleanchild()` is called (once per batch). The wrapper only loads HC velocity files (for child algorithm), so Step 17 re-reads HT velocity files AND re-reads HC velocity files that are different from those passed by the wrapper.
- **Why deferred:** Performance only (~0.01 sec per batch), not correctness. For simulation loops, savings would be negligible.
- **Exact fix:** Load all 4 velocity tables (HT maxvel, HT vel, HC maxvel, HC vel) once in the wrapper and pass as parameters to `cleanchild()`.

### Defer: Step 17 dead code branch
- **File/lines:** `child_clean.R` lines 5080/5104
- **Issue:** Inside `if (count_exclude > 0)` block, the inner `if (count_exclude >= 1)` is always TRUE. The `else` branch (line 5108: `testing <- FALSE`) is dead code.
- **Why deferred:** No behavioral impact — dead branch never executes.
- **Exact fix:** Remove the inner `if/else` and keep only the body of the `>= 1` branch.

### Defer: Step 15 hardcoded exclusion code string
- **File/lines:** `child_clean.R` lines 4287-4294
- **Issue:** Uses hardcoded `"Exclude-C-Traj"` instead of `.child_exc(param, "Traj")`. Functionally identical but inconsistent with rest of the code.
- **Why deferred:** No behavioral impact.
- **Exact fix:** Replace `"Exclude-C-Traj"` with `.child_exc(param, "Traj")`.

### Defer: `ewma()` sort order (line 2155)
- **File/lines:** `child_clean.R` line 2155
- **Issue:** `order(agedays)` without tiebreaker. Same-day values would get arbitrary order. In practice, EWMA is called after SDE resolution, so same-day values shouldn't exist. Also, same-day values would get identical EWMA values regardless of order.
- **Why deferred:** No behavioral impact in practice.
- **Exact fix:** Add `index` or positional tiebreaker, or document why it's safe.

### Defer: Incomplete comments
- **File/lines:** `child_clean.R` lines 3436, 3447
- **Issue:** Sentences trail off ("Min/max HT based on analysis in do file from", "Also, 13 is z=-6 for 22 0/7 in Fenton and").
- **Why deferred:** Comments only, no behavioral impact.
- **Exact fix:** Either complete the sentences or remove them.

### Defer: Step 11 debug columns
- **File/lines:** `child_clean.R` lines 3738-3748
- **Issue:** Six `ewma1_it1.*` columns are unconditionally saved. Should be gated behind a debug parameter.
- **Why deferred:** No behavioral impact — extra columns in output, but callers don't use them.

---

## Adult algorithm steps

### Defer: Dead `keeper_id` variable
- **File/lines:** `adult_clean.R` lines 558, 955
- **Issue:** `keeper_id <- s_df$internal_id[1]` is computed in both 9H and 10W SDE loops but never used. Exclusion logic uses `s_df$internal_id[-1]` for `rem_ids` instead.
- **Why deferred:** Dead code, no behavioral impact.
- **Exact fix:** Remove both `keeper_id` assignments.

### Defer: Stale `"temp extraneous"` string in `eval_2d_nonord()`
- **File/lines:** `adult_support.R` line 1407
- **Issue:** `w_subj_keep[all_wt_ids] != "temp extraneous"` — `w_subj_keep` never contains this string. The condition is always TRUE (dead logic).
- **Why deferred:** No behavioral impact — condition passes through harmlessly.
- **Exact fix:** Remove the `!= "temp extraneous"` condition from the filter.

### Defer: Growing vectors with `c()` in loops
- **File/lines:** Various locations in `adult_clean.R` and `adult_support.R`
- **Issue:** O(n^2) pattern — `exc_ids <- c(exc_ids, new_id)`. Already noted in CLAUDE.md Known Issues.
- **Why deferred:** Negligible for typical per-subject sizes.

### Defer: `sd_median()` exported but internal-only
- **File/lines:** `child_clean.R` line 2406
- **Issue:** Function is `@export`ed but only used internally within `cleangrowth()`.
- **Why deferred:** CRAN cleanup, no behavioral impact.
- **Exact fix:** Change `@export` to `@keywords internal` and `@noRd`, remove from NAMESPACE.
