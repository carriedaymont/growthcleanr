# GC Full Walk-Through — Deferred To-Dos

Items identified during the 2026-04-16 (third) walk-through. Focus: batching and parallel processing.

---

## Pre-walk-through setup

All baseline tests pass:
- test-cleangrowth.R: 92 PASS
- test-child-regression.R: 50 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS (2 expected warnings)
- test-child-parameters.R: 21 PASS (1 expected warning)
- test-adult-clean.R: 198 PASS
- Adult regression harness: 1508/1508 at all 4 permissiveness levels

---

## Issues found and fixed

### F1. Unguarded `stopCluster(cl)` in child parallel path (BUG)
- **File/line:** `child_clean.R` line 1473
- **Issue:** `stopCluster(cl)` was inside the `num.batches > 1` branch but NOT guarded by `if (parallel)`. If a user passes `num.batches > 1` with `parallel = FALSE`, `cl` was never created, causing an error. The adult path at line 1588 correctly uses `if (parallel) { stopCluster(cl) }`.
- **Fixed:** Wrapped with `if (parallel) stopCluster(cl)` to match adult pattern.

### F2. `rbindlist` missing `fill = TRUE` (BUG)
- **File/line:** `child_clean.R` line 1686
- **Issue:** `rbindlist(results_list, use.names = TRUE)` without `fill = TRUE`. When some outer batches have zero pediatric data (all-adult subjects) while others have pediatric data, the result data.tables have different column sets. Specifically:
  - **Checkpoint columns** (`sd.orig`, `sd.corr`, `tbc.sd`, `ctbc.sd`, etc.) are only added when the batch has pediatric data and `checkpoint_data` exists (line 1657-1661).
  - **cf_detail columns** (`cf_status`, `cf_deltaZ`) are only added when `cf_detail = TRUE` and the batch has pediatric data (line 1637-1640).
  - `rbindlist` with `use.names = TRUE` but without `fill = TRUE` would error on mismatched columns.
- **Trigger condition:** Dataset with > `batch_size` subjects where some batches (by subjid sort order) contain only adult-age measurements.
- **Fixed:** Added `fill = TRUE` to `rbindlist` call.

### F3. `var_for_par` missing `.child_exc` and `.cf_rescue_lookup`
- **File/lines:** `child_clean.R` lines 806-814
- **Issue:** `var_for_par` listed 12 internal functions for parallel worker export, but omitted `.child_exc()` (used in ~40 data.table j expressions throughout `cleanchild()`) and `.cf_rescue_lookup()` (used in Step 6 CF rescue). Not a functional bug because `.paropts = list(.packages = c("data.table", "growthcleanr"))` loads the full package namespace on workers, making all internal functions accessible. But inconsistent with the belt-and-suspenders approach of the other listed functions.
- **Fixed:** Added both to `var_for_par`.

---

## Wrapper walkthrough (batching + parallel focus)

### Architecture summary (verified correct)

Two layered batching systems:
1. **Outer batching** (lines 568-684): Divides all subjects into groups of `batch_size` (default 2000). Each group goes through the full pipeline (z-scores, recentering, algorithm dispatch, output assembly). Deterministic assignment: `(seq_len(.N) - 1L) %/% batch_size + 1L`.
2. **Inner batching** (lines 828-834, 1546-1551): For parallel processing, subdivides subjects within an outer batch into `num.batches` sub-batches for dispatch to `cleanchild()`/`cleanadult()` via `ddply`. Random assignment: `sample(num.batches, n, replace = TRUE)`.

### Verified correct
- **Batch boundary safety:** Both outer and inner batching assign whole subjects to batches (lines 572-573, 830-834). No subject is ever split.
- **Recentering consistency:** Child algorithm always uses fixed reference medians file (`rcfile-2023-08-15_format.csv.gz`), loaded once per outer batch. The `sd.recenter` data.table persists across loop iterations, so the file is only read on the first outer batch (line 1263 checks `!is.data.table(sd.recenter)`).
- **`internal_id` across batches:** Assigned on the full dataset (line 540), so values are non-contiguous within batches. `cleanchild()` recreates `index` (line 2972) for contiguous within-batch indexing. No code assumes `internal_id` is contiguous.
- **Z-score blending boundaries consistent:** Wrapper (lines 884-905) and `calc_and_recenter_z_scores()` (lines 2738-2764) use identical boundaries: WHO < 2y, blend 2-5y (weight = (5-age)/3), CDC > 5y, HC always WHO.
- **Partial-run filtering:** Happens before `data.all.ages` is built (lines 505-515). Input vectors are filtered in lockstep. No alignment issues.
- **Output assembly:** `full_out` combines child+adult by `line` (unique per row). Merge with `data.batch` by `line` with `all.x = TRUE` correctly handles all cases. `<` and `>=` split at `adult_cutpoint * 365.25` is complementary — every row goes to exactly one path.

---

## Child algorithm steps (batching-sensitive code)

### Verified correct
- **`cleanchild()` re-sorts on entry:** Line 2968 `setkey(data.df, subjid, param, agedays, internal_id)` ensures correct sort order regardless of how `ddply` passes data.
- **`index` recreation:** Line 2972 `data.df[, index := 1:.N]` creates contiguous indices within the batch. Comment at lines 2969-2971 correctly documents why.
- **`ewma.fields` parameter shadowed:** Line 3983 redefines `ewma.fields` locally with the same value as the parameter passed from the wrapper. Harmless but redundant.
- **Velocity tables in Step 17:** `cleanchild()` loads its own WHO HT and HC velocity files (lines 4515-4558), shadowing the wrapper's `who.ht.vel` parameter. The local files are correct (HT velocity for HT step, HC velocity for HC step). The wrapper's `who.ht.vel` actually contains HC data for the child algorithm path (misleading naming). Already noted as deferred in previous walkthrough.

---

## Adult algorithm steps

No batching-specific issues found. Adult parallel dispatch (`ddply` at line 1567) and `stopCluster` (line 1588, guarded by `if (parallel)`) are correct.

---

## Deferred items

### FIXED: `as.numeric(internal_id)` calls removed
- **Files:** `child_clean.R` (~30 occurrences), `adult_clean.R` (~9 occurrences), `adult_support.R` (~6 occurrences)
- **Fix:** Removed all `as.numeric()` wrappers around `internal_id` since it is now integer.

### FIXED: Parallel cluster moved outside outer batch loop
- **Fix:** Unified parallel setup block (cluster creation, `var_for_par` export, `registerDoParallel`) moved before the outer loop. Single `stopCluster(cl)` after the loop. Both child and adult parallel setup blocks removed from inside the loop. Eliminates ~2-5 sec overhead per batch per cluster cycle.

### FIXED: Inner batch assignment now deterministic
- **Fix:** Replaced `sample(num.batches, n, replace = TRUE)` with `(seq_along(subjid.unique) - 1L) %% num.batches + 1L` in both child and adult inner batching.

### Defer: Legacy `sd.recenter = "derive"` per-batch inconsistency
- **File/lines:** `child_clean.R` lines 1293-1297
- **Issue:** In the legacy algorithm path, when `sd.recenter` is not specified and N >= 5000, recentering medians are derived from the batch data: `data.all[exclude < 'Exclude', sd_median(...)]`. Since `data.all` is one outer batch, different batches get different medians. This could give slightly inconsistent recentering across a large dataset.
- **Why deferred:** Legacy algorithm only (deprecated). Child algorithm always uses fixed reference file.
- **Exact fix:** Derive medians from `data.all.ages` (full dataset) before the outer loop, or always use NHANES reference for legacy.

### FIXED: Walkthrough procedure stale `as.numeric(internal_id)` guidance
- **Fix:** Updated `algorithm-walkthrough-procedure.md` line 184 to describe `internal_id` as integer with no conversion needed.

### FIXED: Wrapper WHO velocity loading removed for child algorithm
- **Fix:** Removed the `else` branch that loaded HC velocity files as `who.ht.vel` for the child algorithm path. Wrapper WHO velocity loading now only happens for the legacy algorithm (`!use_child_algorithm`). Removed `who.ht.vel` parameter from both `cleanchild()` dispatch calls and from the `cleanchild()` function signature. `cleanchild()` loads its own HT and HC velocity files internally in Step 17.

---

## Post-walk-through test results (after all fixes including deferred items)

All tests pass after all fixes (bugs + deferred cleanup):
- test-cleangrowth.R: 92 PASS
- test-child-regression.R: 50 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS (2 expected warnings)
- test-child-parameters.R: 21 PASS (1 expected warning)
- test-adult-clean.R: 198 PASS
- Adult regression harness: 1508/1508 at all 4 permissiveness levels
