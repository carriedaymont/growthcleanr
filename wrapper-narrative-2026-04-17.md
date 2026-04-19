# Wrapper, Preprocessing, and Output Technical Reference (growthcleanr)

# Reference structure

growthcleanr is an R package for automated data cleaning of anthropometric data. It consists of 3 main components: a child algorithm, an adult algorithm, and a function cleangrowth() that performs pre- and post-processing. Each main component has its own narrative, which provides technical detail about the algorithm and code.

This document provides a technical reference for `cleangrowth()` (the package's public entry point), including the preprocessing, reference-table loading, batching, dispatch, and output-assembly infrastructure that surrounds the two algorithm cores. The child algorithm was developed first, and the adult algorithm was implemented in R with a different structure and integrated later. Pre-processing is not entirely consistent between the two algorithms, which is reflected in cleangrowth(). For example, some child-specific processing happens in cleangrowth(). The algorithms have different optional parameters and optional output because of different features of the algorithm and differences in implementation.

Reference Last Updated: 2026-04-18 Initial draft by: Claude (Opus), with Carrie Daymont.

Code Last Updated: 2026-04-18

This document describes what the current R code does in the package wrapper. The two algorithm cores are documented separately:

- `child-algorithm-reference.md` — child algorithm (Child Steps 5, 6, 7, 9, 11, 13, 15/16, 17, 19, 21, 22)
- `adult-algorithm-narrative.md` — adult algorithm (Adult Steps 1, 2W, 3, 4W, 9Wa/9H/9Wb, 10H/10W, 11H/11Wa/11Wa2/11Wb, 13, 14)

Child Step 2b (Gestational Age Correction) is documented here because it runs as preprocessing inside `cleangrowth()`, before dispatch to either algorithm.

---

## How to Read This Document

This narrative is structured to mirror the algorithm narratives where useful. Child and adult algorithm steps are all named and numbered within the code and the reference. This wrapper reference includes **named phases** rather than numbered steps to describe the function of `cleangrowth()`. 

---

## Scope and Contents

This narrative covers everything that runs in the package infrastructure:

- The `cleangrowth()` public interface (input vectors, output data.table)
- Input validation, unit conversion, LENGTHCM relabeling
- `internal_id` assignment and the user-`id` preservation contract
- Z-score calculation pipeline (CSD method, WHO/CDC blending, recentering, Child Step 2b GA correction)
- Adult/child age split and dispatch to `cleanchild()` / `cleanadult()`
- Reference-table loading, `gc_preload_refs()`, partial-run caching (`cached_results` + `changed_subjids`)
- Outer-batch memory management and inner parallel batching
- Output reassembly (combining child + adult results, NA semantics for cross-algorithm columns)
- Wrapper-level configurable parameters
- Cross-cutting conventions and pitfalls

The two algorithm narratives focus on what happens once a batch of pediatric or adult data is dispatched to `cleanchild()` or `cleanadult()`.

---

## Key Concepts (wrapper-level)

***Data retention:*** Although we refer to values as being excluded, growthcleanr does not remove any values from the dataset returned to users. It adds a column to the user's input recommending inclusion or exclusion. It does exclude values from further consideration in cleaning. For example, a value of 690 kg will be excluded early in the algorithm and is ignored thereafter. As shorthand, we refer to this process as excluding a value.

***Three growth parameters:*** growthcleanr cleans weight (WEIGHTKG), height/length (HEIGHTCM/LENGTHCM), and head circumference (HEADCM). LENGTHCM is relabeled to HEIGHTCM — only the `param` value changes; the measurement value is **not** altered, i.e. no adjustment is made for the ~0.7 cm difference between supine length and standing height because of inconsistency in measurement and recording practices. Z-scores use WHO's length reference <2y and CDC's height reference ≥2y. HEADCM is only cleaned through age 3 years (marked `Exclude-Not-Cleaned` for `agedays > 3 × 365.25`); WHO HC reference data ends at 5 years, and HC rows at `agedays ≥ 5 × 365.25` also get `Exclude-Not-Cleaned` since they cannot be z-scored. The adult algorithm cleans HT and WT only (no HC).

***Child-adult split***: Rows are assigned to the child or adult algorithm based on row-level age. By default, rows with `agedays < 20*365.25` are evaluated by the child algorithm, and rows with `agedays >= 20*365.25` are evaluated by the adult algorithm. If a subject has rows that span both ages, their child and adult rows will be evaluated separately and then combined in the output.

***CSD z-scores, not LMS:*** Z-scores for the child algorithm are calculated in the wrapper. growthcleanr z-score calculations use the Conditional Standard Deviation (CSD) method, not the standard LMS method. The CSD method uses the absolute value of half the distance between the median and +2SDs (SDpos) and the median and -2SDs (SDneg) as stand-ins for a standard deviation. The CSD z-score then equals (value - median)/SD, using SDpos for values above the median and SDneg for values below the median. The LMS method accounts for skewness in the distribution of weights, which is often useful, but in cleaning it reduces the ability of z-scores to distinguish large absolute differences in weight for heavier children, reducing cleaning accuracy. The child algorithm uses several variants of CSD z-scores (`sd.orig`, `tbc.sd`, `sd.corr`, `ctbc.sd`, `final_tbc`) — see "Z-Score Infrastructure" below. The adult algorithm does not use z-scores at all; it works with raw measurements and BMI only.

***Dual z-scores for prematurity correction:*** For child subjects whose data are consistent with prematurity or SGA, a corrected z-score (`sd.corr`) is calculated using values derived from the Fenton 2025 percentile reference (CSD method). Corrected z-scores are used as an extra check to prevent exclusion of values that are plausible for a baby with very low birth weight. See the Child Step 2b section below for the full GA-correction logic.

***Recentering:*** After calculating CSD z-scores, the algorithm subtracts a population-level median z-score for each parameter, sex, and age to create "to-be-cleaned" SD scores (`tbc.sd` from `sd.orig`, `ctbc.sd` from `sd.corr`). Because the median z-score changes with age, using an unrecentered z-score can lead to bias in exclusions when comparing z-scores across long time periods. The recentering medians come from a fixed reference file (`rcfile-2023-08-15_format.csv.gz`) derived from outpatient records from the Children's Hospital of Philadelphia.

***User `id`*** is preserved and returned to the user. ***`internal_id` is the algorithm's internal identifier used for sorting:*** `cleangrowth()` assigns `internal_id` as a sequential integer `1:N` in id-sorted order at session start. All internal sorting and tiebreaking uses `internal_id`. The user's original `id` column is preserved as-is into output and is never used for algorithm logic. Both algorithms share this contract.

***Age-dependent `internal_id` tiebreaking (child):*** When resolving same-day duplicates, the child algorithm uses an age-dependent rule: at birth (`agedays == 0`), keep the lowest `internal_id` (earliest measurement, before postnatal fluid shifts and interventions); at all other ages, keep the highest `internal_id` (later measurement, which may represent a more careful re-measurement). The adult algorithm always keeps the highest `internal_id`.

---

## Architecture: `cleangrowth()`

`cleangrowth()` (in `child_clean.R`, despite the file name) is the package's public entry point. It accepts five parallel input vectors plus `id` and a long parameter list, builds a single working data.table, splits it by age into pediatric and adult strata, computes z-scores and recentering for the pediatric branch, dispatches each stratum to its specialized cleaning routine (`cleanchild()` / `cleanadult()`), and reassembles outputs in the original input row order.

The two algorithm cores (`cleanchild()`, `cleanadult()`) live in `child_clean.R` and `adult_clean.R` respectively. The wrapper guarantees them clean, validated, batched input and combines their results into a single output data.table.

### Phase overview

`cleangrowth()` runs as 15 sequential phases. Phases 6–14 execute once per outer batch; all other phases run once per call.

| Phase | Purpose |
|---|---|
| 1 | Parameter validation, NULL-binding block |
| 2 | Partial-run filter to changed subjects (if `cached_results` set) |
| 3 | Build `data.all.ages`; assign `internal_id` |
| 4 | Outer-batch construction (memory); load shared constants |
| 5 | Parallel cluster setup (no-op if `parallel = FALSE`) |
| 6 | Per-batch pediatric/adult row split |
| 7 | Imperial → metric conversion (pediatric) |
| 8 | `LENGTHCM` → `HEIGHTCM` relabel |
| 9 | Inner-batch assignment + CSD z-scores (WHO, CDC, blend) |
| 10 | GA correction for `potcorr` subjects (Child Step 2b) |
| 11 | Recentering; `Exclude-Missing` / `Exclude-Not-Cleaned` assignment |
| 12 | Checkpoint capture + dispatch to `cleanchild()` |
| 13 | Dispatch to `cleanadult()`; post-dispatch `Exclude-Missing` |
| 14 | Per-batch output reassembly (child + adult harmonization) |
| 15 | Batch teardown, partial-run merge, `bin_exclude`/`tri_exclude` |

Most phases are pre- and post-processing bookkeeping. Phases 9 and 10 (z-scores) are the only wrapper-level phases that directly affects child algorithm behavior.

### Phase-by-phase walkthrough of `cleangrowth()`

#### Phase 1: Parameter validation

`cf_rescue` is constrained to `"standard" | "none" | "all"` via `match.arg()`. A long block of `<- NULL` assignments declares the data.table column names that will be created later; these exist purely to silence `R CMD check`'s "no visible binding for global variable" notes and have no runtime effect.

Quality checks run shortly afterward (after the data.table is built): `adult_cutpoint` must be numeric and is clamped to `[18, 20]`; `adult_scale_max_lbs` must be numeric and non-negative; rows whose `param` is not one of the six accepted strings have their measurement set to `NA` (with a console message), so they survive into the output but flow through as `Exclude-Missing`.

#### Phase 2: Partial-run filtering (`cached_results` / `changed_subjids`)

If `cached_results` is non-NULL, `cleangrowth()` enters partial-run mode. With `changed_subjids = NULL`, it auto-detects which subjects to re-run by hashing each subject's `(param, agedays, sex, v)` tuples (with the `measurement == 0 → NaN` transform applied so the comparison matches the cached form) and comparing against the cache. Added, removed, and modified subjects are computed via `setdiff` / `intersect` plus a per-subject `paste(..., collapse = "\n")` hash join. With `changed_subjids` supplied explicitly, auto-detection is skipped.

If the changed set is empty, the cache is returned untouched. Otherwise the input vectors (including `id`) are filtered to changed subjects only, and the rest of the function processes that subset. Full-vs-changed reassembly happens in Phase 15. See **Partial runs and preloaded references** in the Batching and Dispatch section for use-case context.

#### Phase 3: data.table construction and `internal_id` assignment

The five input vectors plus a `line` index (the original input row position, which is what the output will eventually be re-sorted to) are bound into `data.all.ages`. Two measurement columns are stored: `v` is the metric-target value with `0 → NaN` applied (used by the child algorithm), and `v_adult` is the original measurement preserved in its caller-supplied units (used by `cleanadult()`, which performs its own imperial conversion). `sex` is normalized: `0`/`"m"`/`"M"` → `0L`, `1`/`"f"`/`"F"` → `1L`, anything else → `NA`.

`id` is required, and a length mismatch between `id` and `measurement` raises a hard `stop()`. `id` is stored as-is on `data.all.ages`. The table is then keyed by `(subjid, param, agedays, id)` and `internal_id := seq_len(.N)` is assigned. `internal_id` is the sequential integer used for all internal sorting and tiebreaking; the user's `id` is never touched by algorithm logic and is only used for output joining.

#### Phase 4: Outer-batch construction and shared constants

Subjects (not rows) are partitioned into outer batches of `batch_size` (default 2000); the `batches` data.table maps `subjid → batch`. Each outer-batch iteration operates on a fresh `data.all` / `data.adult` / `data.batch` copy that gets garbage-collected before the next iteration, which caps peak memory regardless of input size. The inner parallel batching (Phase 9) is layered on top.

Two batch-invariant constants are also built here: the `exclude.levels` factor (concatenation of `exclude.levels.peds` and `exclude.levels.adult`, used everywhere the `exclude` factor is constructed) and the Tanner height-velocity reference (loaded once from `inst/extdata/tanner_ht_vel.csv.gz`, used by Child Step 17 of `cleanchild()`). The Tanner table is loaded even for all-adult datasets — the cost is negligible (~0.01 s).

#### Phase 5: Parallel cluster setup

When `parallel = TRUE`, a single `makeCluster()` is created here and torn down once after the outer loop (Phase 15). All algorithm-internal functions are exported via `clusterExport()` — both child (`.child_valid`, `ewma`, `read_anthro`, `as_matrix_delta`, `identify_temp_sde`, `get_dop`, `calc_otl_evil_twins`, `calc_and_recenter_z_scores`, `.child_exc`, `.cf_rescue_lookup`, `ewma_cache_init` / `update`) and adult (`cleanadult`, `permissiveness_presets`, `resolve_permissiveness`, `compute_et_limit`, `compute_perc_limit`, `compute_wtallow`, `evil_twins`, `propagate_to_rv`, `remove_ewma_wt`, `remove_mod_ewma_wt`, `eval_2d_nonord`, `eval_1d`, `eval_error_load`, etc.). `num.batches` defaults to `getDoParWorkers()` if not set; in serial mode it defaults to `1`. The cost of cluster creation is amortized across all outer batches.

**Known harmless warnings.** With `parallel = TRUE`, `plyr::ddply` emits two `codetools` warnings per `cleangrowth()` call of the form `<anonymous>: ... may be used in an incorrect context: '.fun(piece, ...)'` — one from the child dispatch (Phase 12), one from the adult dispatch (Phase 13). These are false positives: `codetools` flags plyr's internal use of `...` as potentially misused, but the usage is correct. They do not affect results. They will disappear when `plyr::ddply` is replaced with a `data.table` / `foreach` equivalent (tracked as a CRAN-prep item in `gc-github-latest/CLAUDE.md`).

#### Phase 6: Outer-batch loop — pediatric/adult split

For each outer batch `id_batch`, three working copies are made from `data.all.ages`:

- `data.all` — pediatric rows (`agedays < adult_cutpoint * 365.25`)
- `data.adult` — adult rows (`agedays >= adult_cutpoint * 365.25`)
- `data.batch` — all rows for this batch's subjects (used for output assembly)

The split is row-level, not subject-level: a subject who has measurements before and after the cutpoint will appear in both `data.all` and `data.adult`. There is no cross-stratum logic — each side is processed independently and their results are concatenated in Phase 14. If either side is empty, the corresponding processing block is skipped and an empty data.table is used as the placeholder.

#### Phase 7: Imperial → metric conversion (pediatric only)

Before z-score calculation, `data.all` is converted in place: `HEIGHTIN → HEIGHTCM` (`v := v * 2.54`) and `WEIGHTLBS → WEIGHTKG` (`v := v / 2.2046226`). The `param` strings are then overwritten with the metric labels. The adult algorithm performs its own conversion internally on `v_adult`, which is why `v_adult` is preserved unchanged on `data.all.ages`.

#### Phase 8: LENGTHCM → HEIGHTCM relabel

`LENGTHCM` rows are relabeled `HEIGHTCM`. The relabel happens here so that all downstream pediatric logic can treat HEIGHTCM as a single param. The measurement is **not** adjusted for the supine-versus-standing offset (~0.7cm) because of inconsistent measurement and recording practices across health systems. All linear growth measurements <2y are evaluated as length and all >=2y are evaluated as height.

#### Phase 9: Inner-batch assignment and z-score calculation

`data.all` is keyed by `subjid`, then a per-subject batch column is added (`(seq_along(unique(subjid)) - 1L) %% num.batches + 1L`). When `num.batches > 1` this drives the `ddply(..., .parallel = parallel)` call later; when `num.batches == 1` it is a no-op.

Z-score calculation uses WHO and CDC reference closures from `read_anthro()`. If `ref_tables` was supplied via `gc_preload_refs()`, the preloaded closures are used and the disk read is skipped (~0.9 s saved per call); otherwise `read_anthro()` is invoked twice (once with `cdc.only = TRUE`, once with `cdc.only = FALSE`). Each closure is applied to `(param, agedays, sex, v, TRUE)` to populate `sd.orig_cdc` and `sd.orig_who`. The age-blending logic writes the final `sd.orig`; `sd.orig_uncorr` is a snapshot captured before GA correction. Full CSD and blending formulas are in the **Z-Score Infrastructure** section.

#### Phase 10: GA correction (Child Step 2b)

Subjects whose first weight is `< -2 SD` at `< 10 months` are flagged `potcorr` and re-scored against Fenton 2025 curves using postmenstrual age; the corrected z-score (`sd.corr`) acts as an extra check against inappropriate exclusion for plausibly-preterm infants. If no subjects qualify, the block short-circuits and skips all Fenton I/O, assigning `sd.corr := sd.orig` and `uncorr := 0L`. The result columns (`sd.corr`, `uncorr`, `potcorr`) are merged back onto `data.all` by `id`; `data.all` is then trimmed to its original column set plus those three. Full mechanics (`potcorr_wt` flagging, two Fenton merges, CSD z-score formula, corrected/blended `sd.c`, abssumdiff revert rule) are in the **Gestational Age Correction (Child Step 2b, preprocessing)** section.

#### Phase 11: Recentering and Missing/Not-Cleaned assignment

`data.all` is rekeyed `(subjid, param, agedays, internal_id)`, and a fresh row index `index := 1:.N` is added. The `exclude` factor column is created with two pre-step assignments:

- `Exclude-Missing` for `is.na(v) | agedays < 0`
- `Exclude-Not-Cleaned` for `param == "HEADCM" & agedays > 3 * 365.25`

`cf_rescued` is initialized to `""` (Child Step 6 populates it). Recentering loads the built-in `rcfile-2023-08-15_format.csv.gz` (or a user-supplied `sd.recenter` data.table if provided), joins on `(param, sex, agedays)`, and computes `tbc.sd := sd.orig - sd.median` and `ctbc.sd := sd.corr - sd.median`; the table is then rekeyed to `(subjid, param, agedays, internal_id)`. Two safety checks follow: rows where `tbc.sd` is `NA` (recentering miss) get `Exclude-Missing`, and HC rows with `agedays >= 5 * 365.25` get `Exclude-Not-Cleaned` (WHO HC reference only goes to 5 years). Combined with the pre-recentering `agedays > 3 * 365.25 → Exclude-Not-Cleaned` assignment, this keeps a single consistent HC exclusion code across both "we don't clean HC >3y" and ">=5y has no reference data" regions. Recentering rationale and the `sd_median()` derivation are in the **Z-Score Infrastructure** section.

#### Phase 12: Checkpoint capture and dispatch to `cleanchild()`

Before dispatch, `checkpoint_data` is captured: a copy of `(id, subjid, param, agedays, sd.orig, sd.corr, sd.orig_uncorr, z.orig, tbc.sd, ctbc.sd)` columns that exist on `data.all`. These are merged back onto the output by `id` after dispatch, so the user can inspect intermediate z-scores even though `cleanchild()` only returns `(id, exclude, cf_rescued, ...)`.

Dispatch follows the `num.batches` switch: serial calls `cleanchild()` once on `data.all`; parallel calls `ddply(data.all, .(batch), cleanchild, .parallel = parallel, ...)`. Both paths pass the same parameters: BIV cutoffs (eight `biv.z.*`), `cf_rescue`, `cf_detail`, `error.load.*`, `ewma_window`, the preloaded Tanner table, the WHO/CDC `measurement.to.z` closure, and `ref_tables` (for any internal z-score recomputation inside `cleanchild()`). The result is `ret.df`. See **Batching and Dispatch** for inner-batch and parallel mechanics.

#### Phase 13: Adult dispatch to `cleanadult()`

`data.adult` gets its own inner-batch column (`newbatch`) and an `age_years` column. `measurement := v_adult` is set so `cleanadult()` reads from the original-units column (it does its own imperial conversion). Dispatch mirrors the child path: serial → single `cleanadult()` call, parallel → `ddply(..., .(newbatch), cleanadult, ...)`. Only `adult_permissiveness` and `adult_scale_max_lbs` are exposed at the wrapper level; all other adult parameters fall back to the permissiveness preset.

After `cleanadult()` returns, `Exclude-Missing` is applied **post-dispatch** (`res[is.na(measurement) | agedays < 0, result := "Exclude-Missing"]`). This differs from the child path, where `Exclude-Missing` is set *before* dispatch (because `.child_valid()` needs it pre-set to filter rows). The adult algorithm handles NAs internally via its BIV step, so post-dispatch assignment is correct on that side.

#### Phase 14: Output reassembly

`ret.df` (child) and `res` (adult) are concatenated into `full_out` keyed by `line`, with column harmonization: `cf_rescued` is empty string for adult rows; `mean_ht` and `bin_result` are `NA` for child rows. If `cf_detail = TRUE` and `cf_status` is present in the child output, the detail columns are added with adult rows getting `NA`.

`full_out` is then merged onto `data.batch` by `line` (so all original input columns are preserved on the output), then back-merged with `checkpoint_data` by `id` (adult rows get `NA` for the checkpoint columns). The result is reordered to match `data.batch$id` (input row order for this batch) and finally `setorder(all_results, line)` restores the original input row order. The batch result is appended to `results_list`.

A safety guard: any `Exclude-C-Temp-Same-Day` surviving to the final output triggers a `stop()` — that code is internal to Child Steps 5/13 and must be resolved before assembly. The adult algorithm handles temp SDEs differently: it uses an internal `extraneous` boolean column on its working dataframe (set by `temp_sde()` in `adult_support.R`, consumed and resolved inside Adult Steps 9H/10W) rather than an `Exclude-A-Temp-Same-Day` factor code. Adult output therefore never carries a temp-SDE exclusion code, and the Phase 14 safety guard is specific to the child path.

#### Phase 15: Outer-batch teardown, partial-run merge, and bin/tri exclude

After the outer loop, the parallel cluster (if any) is stopped, and `results_list` is `rbindlist`-ed with `use.names = TRUE, fill = TRUE` (handles missing `cf_status` / `cf_deltaZ` columns when `cf_detail = FALSE`). One more `setorder(all_results, line)` restores the input row order.

In partial-run mode, `all_results` is `rbind`-ed with unchanged-subject rows from `cached_results` (filtered by `!subjid %in% changed_subjids`) and ordered by user `id`. User `id` is preserved untouched through both the cached and re-processed paths and is contract-guaranteed unique, so it produces a stable, meaningful whole-dataset order. (The full-run path is unaffected: it still restores by `line`.)

Two derived columns are added last so they reflect both full-run and partial-run paths consistently:

- `bin_exclude`: `"Include"` vs `"Exclude"` (always present)
- `tri_exclude` (opt-in via `tri_exclude = TRUE`): `"Include"` vs `"Same-Day"` (for any `Exclude-{C,A}-Identical` or `Exclude-{C,A}-Extraneous`) vs `"Exclude"` (all other Exclude values)

The final `all_results` data.table is returned. See **Output Format** for the column inventory and contract.

---

## Batching and Dispatch

### Purpose

The outer batching wrapper (added by Chris Penney) divides subjects into groups for memory management. Without batching, preprocessing (z-score calculation, recentering, reference table merges) for the full dataset must fit in memory at once. For large datasets (100s of thousands or millions of subjects), this can cause memory pressure. The wrapper processes one batch at a time, then combines results.

### Subject-level batching

Batching occurs at the **subject level**, not the row level. All of a given subject's rows — across all parameters (WEIGHTKG, HEIGHTCM, HEADCM) and all ages — are assigned to the same batch. This is a logistical choice: batching by subject before the child/adult age split keeps both halves of a subject's data within the same outer-batch iteration, so output reassembly can combine child and adult results for that subject without cross-batch bookkeeping. The child and adult algorithms themselves are independent and do not share trajectory data across the age split.

### Interaction with the child/adult age split

After batching by subject, `cleangrowth()` splits each batch's data by `adult_cutpoint` (default 20 years, range 18–20) into pediatric (`data.all`) and adult (`data.adult`) subsets. A subject whose rows span the cutpoint will have some rows sent to `cleanchild()` and others to `cleanadult()`. This is intentional: the child and adult algorithms are independent pipelines with no interaction between them. A subject's child measurements are cleaned using only their other child measurements, and their adult measurements are cleaned using only their other adult measurements.

Both sets of rows for such a subject are guaranteed to be in the same batch because batching is done on `subjid` before the age split. The result assembly at the end of each batch iteration recombines the child and adult results for the batch's subjects, and `rbindlist()` combines across batches.

### Two batching systems in the code

The code contains two layered batching systems with different purposes:

**System 1 — Outer wrapper (for memory management):**
- In `child_clean.R`, inside `cleangrowth()`.
- Creates batches of up to `batch_size` subjects (default 2000; configurable via the `batch_size` parameter): `patients[, batch := (seq_len(.N) - 1L) %/% batch_size + 1L]`.
- Iterates the batches with a `for` loop. Each iteration filters the input to the current batch's `subjid` values (both the pediatric `data.all` and the adult `data.adult` are limited to those subjects), runs preprocessing and the child/adult algorithms on just that batch, and pushes the per-batch result onto `results_list`.
- After the loop, `rbindlist(results_list)` combines per-batch results and the merge against `data.batch` restores user input columns.

**System 2 — Inner batching (for parallelism):**
- Also in `child_clean.R`, inside each outer-batch iteration.
- Assigns round-robin batch numbers to subjects in the already-filtered `data.all`: `(seq_along(subjid.unique) - 1L) %% num.batches + 1L`.
- When `parallel = FALSE`, `num.batches = 1`, making this a no-op: all subjects in the outer batch go to a single `cleanchild()` call.
- When `parallel = TRUE`, dispatches via `ddply(..., .parallel = TRUE)` to process inner batches in parallel across CPU cores.

The inner batching subdivides whatever data the outer wrapper provides into parallel work units. With `parallel = FALSE`, it has no effect. The adult algorithm has its own inner parallel batching system, structured identically.

### Parallel processing

`parallel = TRUE` runs the inner batches concurrently. The worker processes are separate R sessions that resolve libraries from `.libPaths()` and need `system.file()` access to the `inst/extdata/` reference tables. This means growthcleanr must be installed (via `devtools::install_github()`, `devtools::install_local()`, or CRAN), not merely source-loaded via `devtools::load_all()`.

**Practical note for developers:** There is a limitation that applies only to package developers. Most uses of `devtools` — including `devtools::install_github("carriedaymont/growthcleanr", ...)` and `devtools::install_local("gc-github-latest")` — actually **install** the package into your R library, and `parallel = TRUE` works normally after those calls. The issue only appears when a developer uses `devtools::load_all()` to load the source tree into the current R session without installing (a workflow used for rapid iteration on package code). The workaround in that case is to run `devtools::install_local(".")` once, which takes only a few seconds.

Results are bit-identical to sequential processing (verified on the 77,721-row `syngrowth` dataset). Typical observed speedups on 2,697 subjects with many errors (which results in longer run-time): ~114 s sequential, ~65 s at `num.batches = 2` (1.7×), ~30 s at `num.batches = 4` (3.9×).

Batch-invariant setup (the `exclude.levels` factor and the Tanner / WHO velocity reference files) is performed once outside the outer loop so the cost is not repeated per batch.

### Partial runs and preloaded references

Two `cleangrowth()` parameters support workflows that call the function many times on overlapping data (simulation loops, error-injection pipelines, secure-environment workflows with updated data extracts). Both are optional and independent; either or both can be used.

**`ref_tables` + `gc_preload_refs()`.** `gc_preload_refs()` reads the CDC and WHO growth-chart reference files once and returns a named list of closures (`mtz_cdc_prelim` and `mtz_who_prelim`). Passing that list as the `ref_tables` argument to subsequent `cleangrowth()` calls skips the per-call `read_anthro()` disk reads (≈0.9 seconds per call). The return value is a session-lifetime object; rebuild only if the reference files themselves change.

```r
refs <- gc_preload_refs()
res  <- cleangrowth(..., ref_tables = refs)
```

**`cached_results` (+ optional `changed_subjids`).** When rerunning `cleangrowth()` on input that differs from a prior run for only a subset of subjects, pass the prior result as `cached_results`. Only subjects whose input rows differ are re-processed; unchanged subjects are copied from the cache. Subjects are independent in all cleaning operations (by-subject grouping, fixed-reference recentering), so partial runs produce the same output a full run would for the subjects that changed, and bit-identical output for the ones that didn't.

There are two modes for cached-results:

- ***Auto-detect*** (preferred — `cached_results` provided, `changed_subjids = NULL`): `cleangrowth()` compares the incoming `subjid`, `param`, `agedays`, `sex`, and `measurement` against the cache (with the same `0 → NaN` transformation applied internally) and re-processes any subjects with added, removed, or modified rows.
- ***Explicit*** (`cached_results` + `changed_subjids` provided): only the listed subjects are re-processed; the comparison step is skipped.

With `quietly = FALSE`, the wrapper prints a one-line summary: `Auto-detected N changed subjects (X added, Y modified, Z removed, W unchanged)`. Output row order matches input order (sorted by the session-assigned `internal_id`).

Observed runtimes (200 subjects, 50 000 calls, 5 % of subjects changed per call): ≈0.40 s/call vs ≈3.4 s/call for a full run with preloaded refs — roughly an 11× speedup and ~56 hours saved over the simulation. See `gc-github-latest/CLAUDE.md` for the full benchmark table.

---

## Data Requirements and Input Format

### Required input columns

`cleangrowth()` accepts individual vectors (not a dataframe):

| Parameter | Type | Description |
|-----------|------|-------------|
| `subjid` | any | Subject identifier |
| `param` | character | `"WEIGHTKG"`, `"WEIGHTLBS"`, `"HEIGHTCM"`, `"HEIGHTIN"`, `"LENGTHCM"`, or `"HEADCM"` |
| `agedays` | numeric | Age in days at measurement |
| `sex` | any | `0`, `"m"`, or `"M"` for male; `1`, `"f"`, or `"F"` for female |
| `measurement` | numeric | Recorded value in units specified by `param` |
| `id` | numeric or character (optional) | Unique row identifier; if NULL, generated as `1:N` |

### Input handling

- `measurement == 0` is replaced with `NaN` (treated as missing) in the internal variable `v`
- Imperial measurements are converted to metric before processing: HEIGHTIN multiplied by 2.54 and relabeled HEIGHTCM; WEIGHTLBS divided by 2.2046226 and relabeled WEIGHTKG
- LENGTHCM is relabeled to HEIGHTCM. (No length-height correction is applied because of variability in measurement and recording practices across health systems.) 
- `sex` is recoded to integer: 0 = male, 1 = female
- `agedays` is cast to integer
- Data are split by `adult_cutpoint` (default 20 years): rows with `agedays < cutpoint * 365.25` go to the child algorithm; rows at or above go to the adult algorithm

### Sorting

The primary sort order is:

```
setkey(data.df, subjid, param, agedays, id)
```

This sort is critical for deterministic results. The `id` field breaks ties when multiple measurements share the same `subjid`, `param`, and `agedays` (same-day duplicates). If `id` is not included in the sort key, SDE resolution order is undefined and results may differ between sequential and parallel runs.

### Internal identifiers

- `line` — row number in original input order; used to restore output order
- `index` — sequential row number after sorting; used as internal key within `cleanchild()`
- `id` — user-provided or auto-generated; used as tiebreaker in sorting and SDE resolution
- `v` — working measurement value (metric, with 0 → NaN)
- `v.orig` — copy of `v` before any transformations (unit error recovery, swapped parameter correction)

---

## Z-Score Infrastructure

The child algorithm computes and uses multiple z-score variants. Understanding which is used where is essential for debugging. The adult algorithm does not use z-scores; it works with raw measurements (and BMI computed internally).

### Calculation method: CSD (Conditional Standard Deviation)

All z-scores in the cleaning algorithm use the CSD method:

```
If measurement < M:  sd = (measurement - M) / csd_neg
If measurement >= M: sd = (measurement - M) / csd_pos
```

Where:
- `M` = age/sex-specific median from reference table
- `csd_pos` = half the distance from M to the value at z = +2, precomputed in reference tables
- `csd_neg` = half the distance from M to the value at z = -2, precomputed in reference tables

This is NOT the standard LMS z-score used in clinical practice. CSD z-scores are intentionally more sensitive to extreme high values for weight, which helps the algorithm detect implausible measurements.

### Reference tables

| Table file | Source | Used for |
|------------|--------|----------|
| `growthfile_who.csv.gz` | WHO 2006 | HT/WT age < 5y; all HEADCM |
| `growthfile_cdc_ext_infants.csv.gz` | CDC 2000 extended | HT/WT age ≥ 2y |
| `fent_foraga.csv.gz` | Fenton | Weight → est. gestational age |
| `fenton2025_ms_lookup_smoothed.csv` | Fenton 2025 | Approximate values derived from published percentile charts: M, S_upper (half-distance between median and +2SDs), S_lower (half-distance between median and -2 SDs) for CSD z-scores |
| `rcfile-2023-08-15_format.csv.gz` | Derived | Recentering medians |

### Z-score variants

| Variable | Source | Description |
|----------|--------|-------------|
| `sd.orig_who` | WHO CSD | WHO-only CSD z-score |
| `sd.orig_cdc` | CDC CSD | CDC-only CSD z-score |
| `sd.orig` | Blended | Age-blended WHO/CDC CSD z-score (see below) |
| `sd.orig_uncorr` | Copy | Copy of `sd.orig` before GA correction; used in CF rescue (Child Step 6) and BIV (Child Step 7) |
| `sd.corr` | GA-corrected | Fenton-corrected z-score for potcorr subjects; equals `sd.orig` for others |
| `tbc.sd` | Recentered | `sd.orig - sd.median`; primary score used by most algorithm steps |
| `ctbc.sd` | Recentered corrected | `sd.corr - sd.median`; used alongside `tbc.sd` in Evil Twins and EWMA steps |
| `final_tbc` | Output | `ctbc.sd` for potcorr subjects, `tbc.sd` for others |


### Age blending — main z-score (`sd.orig`)


| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | `sd.orig = sd.orig_who` |
| 2–5 years | HT, WT | `(sd.orig_cdc × (age-2) + sd.orig_who × (5-age)) / 3` |
| > 5 years | HT, WT | `sd.orig = sd.orig_cdc` |
| any age | HEADCM | `sd.orig = sd.orig_who` |

Note: The blending window is 2–5 years with divisor 3. At exactly 2 years, the formula yields 100% WHO. At exactly 5 years, it yields 100% CDC.

***Code-verified detail:*** The code assigns WHO for `ageyears < 2` and CDC for `ageyears > 5` (strict inequality). For ages exactly 2.0 or 5.0, the smooth_val condition (`>= 2 & <= 5`) applies, so they go through the blending formula. This means ages exactly at the boundaries use the blending formula, which gives 100% WHO at 2.0 and 100% CDC at 5.0.

### Age blending — corrected z-score (`sd.corr`)

For potcorr subjects, `sd.corr` is calculated in two stages:

**Stage 1 — Compute `sd.c` (the corrected z-score before smoothing into original):**

First, corrected WHO and CDC z-scores (`sd.c_who`, `sd.c_cdc`) are calculated using corrected age (`cagedays`) with an adjustment for standing vs. supine position when chronological age > 730 days but corrected age ≤ 730 days (+0.8 cm for WHO, +0.7 cm for CDC).

Then `sd.c` is assigned using WHO/CDC blending with the same 2–5 year window as `sd.orig`:
- HEADCM: WHO at all ages
- HT/WT age ≤ 2: WHO
- HT/WT age 2–5: weighted blend `(who × (5-age) + cdc × (age-2)) / 3`
- HT/WT age ≥ 5: CDC

For potcorr subjects age ≤ 2 with corrected WHO available and post-menstrual age ≥ 350 days, `sd.corr` is set to `sd.c_who`.

Then Fenton is preferred over corrected WHO for ages ≤ 2: if Fenton (`unmod_zscore`) is available, it overwrites `sd.c`; corrected WHO is used only as fallback.

**Stage 2 — Smoothing corrected into original:**

| Age range | Formula |
|-----------|---------|
| ≤ 2 years (potcorr) | `sd.corr = sd.c` (fully corrected) |
| 2–4 years | `(sd.orig × (4 - age) + sd.c × (age - 2)) / 2` |
| > 4 years | `sd.corr = sd.orig` (no correction) |

Note: This smoothing window is 2–4 years with divisor 2, which is narrower than the main WHO/CDC blending window (2–5 years, divisor 3).

### Age blending — `calc_and_recenter_z_scores()` (Child Step 15/16 p_plus / p_minus)

This helper function, used in the Child Step 15/16 pre-loop (EWMA2 / moderate trajectory outliers), recalculates recentered z-scores for each row's perturbed `p_plus` and `p_minus` measurement values. It uses the same WHO/CDC blending formula as the main z-score calculation (2–5 year window) so that the perturbed z-scores are on the same footing as `tbc.sd`.

| Age range | Parameter | Formula |
|-----------|-----------|---------|
| < 2 years | any | WHO only |
| 2–5 years | HT, WT | `(cdc × (age-2) + who × (5-age)) / 3` |
| > 5 years | HT, WT | CDC only |
| any age | HEADCM | WHO only |

### Recentering

After z-score calculation, the algorithm subtracts age-, sex-, and parameter-specific median z-scores to produce tbc.sd:

```
tbc.sd = sd.orig - sd.median
ctbc.sd = sd.corr - sd.median
```

Recentering uses a precomputed reference file (`inst/extdata/rcfile-2023-08-15_format.csv.gz`) whose medians were built once offline — see **Shared Helpers → `sd_median()`** below for the midyear-interpolation procedure. Callers can override by passing a `sd.recenter` data.table to `cleangrowth()`.

The recentering median file is indexed by param, sex, and agedays. It is merged into the data by a rolling join on these keys.

**Required schema for a custom `sd.recenter` data.table.** A user-supplied override must be a data.table with exactly these columns:

| Column | Type | Meaning |
|---|---|---|
| `param` | character | `"WEIGHTKG"`, `"HEIGHTCM"`, or `"HEADCM"` |
| `sex` | integer | 0 (male) or 1 (female) |
| `agedays` | integer | Age in days; rolling join matches the nearest lower agedays in the table to each row |
| `sd.median` | numeric | Population median (in CSD z-score units) to subtract from `sd.orig` / `sd.corr` |

`setkey(sd.recenter, param, sex, agedays)` is called on the supplied table before the rolling join; the caller does not need to set the key but the columns must be named exactly as above. Values outside the (param, sex, agedays) range covered by the table receive the clamped endpoint median.

---

## Cross-Algorithm Differences

The two algorithms differ in many cross-cutting respects that are worth understanding when reading either codebase or when combining results.

| Aspect | Child algorithm | Adult algorithm |
|--------|-----------------|-----------------|
| Data structure | Single `data.df` data.table; rows are never physically removed | Copies rows into separate working dataframes (`h_subj_df`, `w_subj_df`) that shrink as values are excluded |
| Exclusion mechanism | Set `exclude` column; filter with `.child_valid()` | Remove from working dataframe |
| Output column | `exclude` | `result` (mapped to `exclude` by `cleangrowth()` in Phase 14) |
| `Exclude-Missing` assignment | Pre-dispatch in `cleangrowth()` | Post-dispatch in `cleangrowth()` |
| Z-scores | CSD z-scores (WHO/CDC blend) | None — raw measurements (BMI computed internally for some thresholds) |
| Param scope | HT, WT, HC | HT, WT only (no HC) |
| Rounding tolerance | None | 0.12 cm/kg on all threshold comparisons |
| SDE tiebreaking at birth | Lowest `internal_id` (pre-postnatal shift) | N/A (no births) — always keeps highest `internal_id` |
| Permissiveness levels | Not yet implemented (planned post-v3.0.0) | 4 levels: loosest/looser/tighter/tightest |
| Iteration exponent (EWMA) | -1.5 | -5 |
| `perclimit` scope | N/A | Adult Step 11Wa: subject-level max wt; Adult Step 11Wb: observation-level |
| Reference data | WHO 2006, CDC 2000 extended, Fenton 2025 | None (no z-scores); uses `inst/extdata/adult_columns.csv` for column-level scaling only |
| Adult-specific output cols | NA on child rows | `mean_ht`, `bin_result` |
| Child-specific output cols | `cf_rescued`, `sd.orig*`, `sd.corr`, `tbc.sd`, `ctbc.sd`, `final_tbc` | NA on adult rows |
| Sort key (within algorithm) | `setkey(data.df, subjid, param, agedays, internal_id)` | All sorts include `internal_id` (integer) as final tiebreaker; `as.character(internal_id)` used for named-vector keys |
| Missing-as-infinity edge handling | N/A | `ifelse(is.na(...), Inf, ...)` for edge EWMA values |

For the full child working-dataframe mechanics (how `.child_valid()` works, the temp-SDE soft-flag mechanism, SDE-Identical removal during CF detection), see `child-algorithm-reference.md → Working Dataframe and Output`.

---

## Variable Glossary (wrapper-level)

For child-specific variables (CF working columns, EWMA working columns, exclusion/status columns), see `child-algorithm-reference.md → Variable Glossary`.

### Input and identification

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `line` | integer | Preprocessing | Original row number; used to restore output order |
| `id` | any (required) | Input | User-provided row identifier. Required (errors if `NULL`). Preserved untouched into output; not used as a sort key or tiebreaker by the algorithm. |
| `internal_id` | integer | `cleangrowth()` | Sequential integer `1:N` assigned in id-sorted order at session start. Used throughout as the final sort-tiebreaker and for named-vector indexing where applicable (age-dependent tiebreaking keeps the lowest `internal_id` at `agedays == 0`, the highest at `agedays > 0`). |
| `index` | integer | `cleanchild` | Sequential `1:.N` within each outer-batch dataframe, reassigned at the top of `cleanchild()` so merges and subset operations within the batch have contiguous keys. |
| `subjid` | any | Input | Subject identifier (any type; preserved as-is into output). |
| `param` | character | Input (post-conversion) | `"WEIGHTKG"`, `"HEIGHTCM"`, or `"HEADCM"` after imperial→metric relabeling. |
| `agedays` | integer | Input | Age in days at the measurement. |
| `sex` | integer | Input | 0 = male, 1 = female. |
| `v` | numeric | Preprocessing | Working measurement value (metric; 0 → NaN). |
| `v.orig` | numeric | Child Step 5 | Copy of `v` captured at the start of Child Step 5, used in Child Step 6 CF detection. |

### Z-score variables

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `sd.orig_who` | numeric | Preprocessing | WHO-only CSD z-score |
| `sd.orig_cdc` | numeric | Preprocessing | CDC-only CSD z-score |
| `sd.orig` | numeric | Preprocessing | Blended WHO/CDC CSD z-score |
| `sd.orig_uncorr` | numeric | Preprocessing | Copy of sd.orig before GA correction |
| `sd.corr` | numeric | Child Step 2b | GA-corrected CSD z-score |
| `sd.median` | numeric | Recentering | Population median z-score for recentering |
| `tbc.sd` | numeric | Recentering | Recentered blended z-score (`sd.orig - sd.median`) |
| `ctbc.sd` | numeric | Recentering | Recentered corrected z-score (`sd.corr - sd.median`) |
| `final_tbc` | numeric | Output | ctbc.sd for potcorr subjects, tbc.sd for others |

### Prematurity correction variables

| Variable | Type | Created | Description |
|----------|------|---------|-------------|
| `potcorr` | logical | Child Step 2b | TRUE if subject may need GA correction (subject-level) |
| `potcorr_wt` | logical | Child Step 2b | TRUE for the qualifying first weight row (row-level) |
| `uncorr` | integer | Child Step 2b | 1 if correction was reverted (made z-scores more extreme) |
| `fengadays` | numeric | Child Step 2b | Estimated gestational age in days (from Fenton) |
| `pmagedays` | numeric | Child Step 2b | Post-menstrual age in days |
| `cagedays` | numeric | Child Step 2b | Corrected age in days |

---

## Output Format

`cleangrowth()` returns a data.table with all original input columns plus additional columns. The columns returned from `cleanchild()` are:

| Column | Type | Description |
|--------|------|-------------|
| `id` | original | User-provided row identifier |
| `line` | integer | Original row order |
| `exclude` | factor→character | Exclusion code or "Include" |
| `param` | character | Growth parameter |
| `cf_rescued` | character | CF rescue reason (empty if not rescued) |
| `sd.orig_who` | numeric | WHO CSD z-score (if available) |
| `sd.orig_cdc` | numeric | CDC CSD z-score (if available) |
| `sd.orig` | numeric | Blended CSD z-score |
| `tbc.sd` | numeric | Recentered blended z-score |
| `ctbc.sd` | numeric | Recentered corrected z-score |
| `final_tbc` | numeric | ctbc.sd for potcorr, tbc.sd for others |

After merging with adult results, `cleangrowth()` also includes checkpoint diagnostic columns (`sd.corr`, `potcorr`, `uncorr`, `sd.orig_uncorr`) when the child algorithm was used, and adult-specific columns (`mean_ht`, `bin_result`) when the adult algorithm was used. Cross-algorithm columns are NA on rows that were processed by the other algorithm (e.g., `tbc.sd` is NA on adult-processed rows).

For the **child-specific exclusion code list** and CF rescue reason codes, see `child-algorithm-reference.md → Output Format`.

For the **adult-specific exclusion code list**, see `adult-algorithm-narrative.md` (or `gc-github-latest/CLAUDE.md → Adult Exclusion Codes`).

---

## Configurable Parameters (wrapper-level)

The parameters in the table below refer to parameters used in the wrapper code. Some of these are child- or adult-specific but are handled in the wrapper.

For parameters used only in the child algorithm code (BIV cutoffs, CF rescue modes, EWMA window, error load, etc.), see `child-algorithm-reference.md → Configurable Parameters`.

For parameters used only in the adult algorithm code (permissiveness presets, wtallow/ET caps, `mod_ewma_f`, `ht_band`, etc.), see `adult-algorithm-narrative.md` (and the `adult_permissiveness` / `adult_scale_max_lbs` rows below for the wrapper-level adult parameters).

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `id` | NULL | Preprocessing | Required row identifier preserved into output. Can be any type (numeric, character, UUID). If `NULL`, `cleangrowth()` errors. |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing the pediatric and adult paths. Subjects with measurements spanning the cutpoint have their pediatric and adult rows cleaned by separate algorithms. |
| `sd.recenter` | NA | Preprocessing (recentering) | If left `NA`, the built-in rcfile is used (`inst/extdata/rcfile-2023-08-15_format.csv.gz`). Pass a data.table to override with custom recentering medians. |
| `ref.data.path` | `""` | Preprocessing | Directory override for reference files (CDC/WHO z-score tables, rcfile, velocity tables, Fenton data). Default empty string means use the installed package's `inst/extdata/`. |
| `batch_size` | 2000 | Batching (outer wrapper) | Number of subjects per outer-wrapper batch. Subjects are never split across batches. |
| `parallel` | FALSE | Batching (inner) | If `TRUE`, inner batches dispatch via `ddply(..., .parallel = TRUE)` for concurrent processing. Requires the package to be installed (not just `load_all()`-ed). |
| `num.batches` | NA | Batching (inner) | Number of inner (parallel) batches. Auto-set when `parallel = TRUE`; ignored when `parallel = FALSE` (acts as 1). |
| `log.path` | NA | Batching (inner) | If non-NA and `parallel = TRUE`, each worker sinks messages to `{log.path}/cleangrowth-{date}-batch-{n}.log`. |
| `ref_tables` | NULL | Performance | Pre-loaded reference closures from `gc_preload_refs()`. Skips repeated disk reads across many `cleangrowth()` calls; see "Partial runs and preloaded references" above. |
| `cached_results` | NULL | Performance | Prior `cleangrowth()` output. Unchanged subjects are copied from cache; changed subjects are re-processed. Either auto-detected or listed explicitly via `changed_subjids`. See "Partial runs and preloaded references" above. |
| `changed_subjids` | NULL | Performance | Optional explicit list of subject IDs to re-process when using `cached_results`. If `NULL`, changed subjects are auto-detected by comparing input to cache. |
| `tri_exclude` | FALSE | Output | If `TRUE`, adds a `tri_exclude` column with a three-level summary (Include / Exclude / N/A). Convenience for downstream filtering. |
| `quietly` | TRUE | All | Suppress progress messages. When `FALSE`, the wrapper prints batch and step progress via `message()`. |
| `adult_permissiveness` | `"looser"` | Adult dispatch | Permissiveness preset passed to `cleanadult()`. See `adult-algorithm-narrative.md`. |
| `adult_scale_max_lbs` | `Inf` | Adult dispatch | Physical scale upper limit in lbs; passed to `cleanadult()`. |

---

## Gestational Age Correction (Child Step 2b, preprocessing)

| | |
|---|---|
| **Scope** | Subjects whose first weight suggests prematurity (potcorr subjects) |
| **Operates on** | All rows of potcorr subjects |
| **Prior step** | Z-score calculation (WHO/CDC) |
| **Next step** | Recentering |
| **Output columns** | `sd.corr`, `uncorr`, `potcorr` (merged back onto all rows) |
| **Exclusion code** | None — computes corrected z-scores only |
| **Code location** | `cleangrowth()` in `child_clean.R` (preprocessing block, before dispatch to `cleanchild()`) |
| **Controlled by** | No user parameters; uses built-in Fenton 2025 reference |

### Overview

Very preterm subjects whose anthropometric references should be preterm growth standards are identified as "potentially correctable" (`potcorr`) and have their z-scores recomputed against the Fenton 2025 reference using postmenstrual age. If the corrected values show a larger change in z-scores than the uncorrected (sum of corrected first-to-second, first-to-third, and first-to-fourth weight deviations sum larger than uncorrected), the correction is reverted for that subject (`uncorr := 1`, `sd.corr := sd.orig`). All three columns (`sd.corr`, `uncorr`, `potcorr`) are written back to `data.all` for every row; downstream steps read them through `final_tbc` and related recentered fields.

### Phase 1: potcorr identification

Subjects are flagged by a two-level rule.

- **Sort.** Sort by `(subjid, param, agedays, id_sort)` where `id_sort = internal_id` at `agedays == 0` and `-internal_id` otherwise. At birth this keeps the earliest-recorded row as the "first" weight; at later ages it keeps the most recent internal_id — the age-dependent tiebreaker used throughout the algorithm.
- **`potcorr_wt`** (row-level): a weight row is `potcorr_wt = TRUE` iff (a) it is the first weight row for that subject under the sort above (`seq_len(.N) == 1L` with `by = subjid`), (b) `sd.orig < -2`, and (c) age is less than 10 months. NA-safe: `sd.orig` is `NA` for `Exclude-Missing` rows so they never flag.
- **`potcorr`** (subject-level): `any(potcorr_wt, na.rm = TRUE)` per subject. Propagated to all rows of that subject.

### Phase 2: fast path when no potcorr

`has_potcorr <- any(data.all$potcorr, na.rm = TRUE)`. If FALSE, Child Step 2b short-circuits: `sd.corr := sd.orig`, `uncorr := 0`. Fenton reference files are not read. This typical-case fast path matters because potcorr subjects are often not present or comprise < 1 % of a dataset.

### Phase 3: Fenton correction (potcorr subjects only)

All correction work happens on a copy `pc <- data.all[subjid %in% pc_ids]`.

1. **Integer weight.** `intwt := trunc(v * 100) * 10` (grams rounded to 10 g); values in `[100, 500)` are floored to 500 to hit the Fenton table's minimum.
2. **Fenton merge #1 — weight → GA.** Merge `pc` with `fent_foraga` (from `inst/extdata/fent_foraga.csv.gz`) on `(sex, intwt)` to get `fengadays` (Fenton gestational age in days for which that weight (rounded to 10g) is the median).
3. **Propagate.** Use the qualifying first-weight row's `fengadays` for every row of the subject (`fengadays_subj := min(fengadays[potcorr_wt], ...)`, by subjid).
4. **Ages.** `pmagedays := agedays + fengadays` (postmenstrual age); `cagedays := pmagedays - 280` (corrected age). Both computed only when `fengadays < 259` (term cutoff).
5. **Fenton merge #2 — GA → M/S.** Map gc param names (`WEIGHTKG`/`HEIGHTCM`/`HEADCM`) to Fenton param names (`weight`/`length`/`headcirc`), then merge with `fenton2025_ms_lookup_smoothed.csv` on `(sex, fengadays, param)` to pull `M`, `S_upper`, `S_lower`.
6. **Fenton CSD z-score.** `v_fenton` is `v * 1000` for weight and `v` otherwise. `unmod_zscore := (v_fenton − M) / (S_upper * M)` if `v_fenton >= M`, else `/ (S_lower * M)`. This uses the same CSD method as WHO/CDC elsewhere in the algorithm.
7. **Correction failure.** If Fenton merge failed (no M), reset `potcorr_wt` to FALSE for that row and recompute the subject-level `potcorr`. A subject whose weight is outside the algorithms' Fenton range silently falls out of the correction. 
8. **Initial `sd.corr`.** Rows with `unmod_zscore` non-NA get `sd.corr := unmod_zscore`. Rows still `NA` fall back to `sd.orig`.
9. **Corrected WHO/CDC z-scores.** Compute `sd.c_who` and `sd.c_cdc` via the in-memory `measurement.to.z_who` / `measurement.to.z` closures at `cagedays`. For HEIGHTCM between ages 2 and 2+ at Fenton-corrected-age-≤2, apply the supine/standing offsets (+0.8 cm for WHO, +0.7 cm for CDC) before the z-score call. This is done because the formal recommendation is for linear growth to be measured as length <2y and height >=2y age. Since we are converting chronologic to corrected age, if recommendations are followed we need to account for the length-height discrepancy (length is higher than height). We need a conistent way of handling this issue, but we recognize that practices vary widely regarding length-height measurement and documentation. Most health systems do not provide metadata re length vs. height in extracted EHR data, and metadata that are provided are not necessarily reliable. 
10. **`sd.c` (corrected blend).** HEADCM uses WHO at all ages. Other params: WHO for `ageyears_2b ≤ 2`, smoothed WHO→CDC blend for `2 < ageyears_2b < 5`, CDC for `ageyears_2b ≥ 5`.
11. **Fenton takes precedence under 2.** For `ageyears_2b ≤ 2` with Fenton available, `sd.c := unmod_zscore`.
12. **Final `sd.corr`.** Age-dependent:
    - `ageyears_2b ≤ 2` with potcorr: `sd.c`
    - `2 < ageyears_2b ≤ 4`: smoothed blend `(sd.orig * (4 − ageyears_2b) + sd.c * (ageyears_2b − 2)) / 2`
    - otherwise: `sd.orig`
    - Note that smoothing age intervals are intentionally different: CDC/WHO: 2-5 years; corrected/uncorrected: 2-4 years

### Phase 4: uncorr check (revert when correction worsens trajectory)

Only weight rows of potcorr subjects participate. A temporary `tmp` copy is built:

- First, all but one copy of identical values for the same subject/parameter/ageday (aka SDE-Identicals) are filtered via the age-dependent `keep_id` rule (match Early Child Step 13 behavior).
- `seq_win := sequence(.N)` per subject; keep rows with `seq_win ≤ 4 AND ageyears_2b < 2`.
- Drop subjects that have only one value or no non-NA z-scores.
- Compute `sd.corr_abssumdiff := |sum(sd.corr[1] − sd.corr)|` and `sd.orig_abssumdiff := |sum(sd.orig[1] − sd.orig)|` per subject.
- A subject's correction is reverted if `sd.corr_abssumdiff > sd.orig_abssumdiff` (and the first weight exists): in `pc`, `sd.corr := sd.orig`, `uncorr := 1`.

### Phase 5: merge results back

`data.all[, sd.corr := sd.orig]` and `uncorr := 0L` baseline, then merge the pc results back by `id`:

```r
data.all[pc_result, on = "id", `:=`(
  sd.corr = i.sd.corr,
  uncorr = i.uncorr,
  potcorr = i.potcorr
)]
```

Temporary columns (`agemonths`, `ageyears_2b`, `potcorr_wt`, etc.) are dropped before the data flows into recentering.

Walkthrough checklist findings for the GA-correction section (Phase 10 / Child Step 2b) are recorded in `walkthrough-todo-2026-04-17.md → GA correction walkthrough checklist findings`.

---

## Shared Helpers

These helpers are defined in `child_clean.R`, in the same file as `cleangrowth()`, and are used by the wrapper or exposed to external callers. Child-algorithm-internal helpers (`.child_valid()`, `.child_exc()`, `get_dop()`, `identify_temp_sde()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `ewma_cache_*`, `as_matrix_delta()`) are documented in `child-algorithm-reference.md → Support Functions`.

### `read_anthro(path = "", cdc.only = FALSE)`

Reads CDC and WHO growth-chart reference tables from `inst/extdata/` (`growthfile_cdc_ext_infants.csv.gz` and `growthfile_who.csv.gz`) and returns a closure that, given `(param, agedays, sex, measurement, return.sd)`, computes the appropriate CSD z-score for that observation. Two flavors are loaded based on `cdc.only`: TRUE returns a CDC-only `mtz_cdc_prelim` closure; FALSE returns a WHO-augmented `mtz_who_prelim` closure that uses WHO under 2 years and CDC above. `path` is an optional override for the reference directory (default `""` → installed package's `inst/extdata/`).

The closure approach lets the reference tables be resolved at load time (one disk read) and reused across many calls (one function call per row). It's the foundation of `gc_preload_refs()`.

### `gc_preload_refs(path = "")`

Preloads both reference closures via two `read_anthro()` calls and returns them as `list(mtz_cdc_prelim, mtz_who_prelim)`. Pass the returned list as `cleangrowth(..., ref_tables = refs)` to skip the per-call disk reads (~0.9 sec saved per call). The closures are session-lifetime — rebuild only if the reference files themselves change. See **Batching and Dispatch → Partial runs and preloaded references** above for benchmark numbers.

Note that velocity reference tables (`tanner_ht_vel`, `who_ht_vel_for_age`, `who_hc_vel_for_age`) are NOT included in `gc_preload_refs()` and are loaded separately (the Tanner table once per `cleangrowth()` call, the WHO velocity tables once per batch inside `cleanchild()`). All three are tiny so the loads are fast; a possible future extension to include them in `gc_preload_refs()` is tracked under `gc-github-latest/CLAUDE.md → Known Issues → Open (wrapper)`.

### `sd_median(param, sex, agedays, sd.orig)`

Exported helper that derives population median CSD z-scores by day of age. **Not called at runtime.** It was used once offline to build the packaged recentering file `inst/extdata/rcfile-2023-08-15_format.csv.gz`; `cleangrowth()` reads that precomputed file directly and does not invoke `sd_median()`. The function is retained as an exported API so that users who want to rebuild the recentering file against their own reference distribution can reproduce the procedure.

Returns a data.table keyed on `(param, sex, agedays)` with columns `param`, `sex`, `agedays`, and `sd.median`. The same `sd.median` value is replicated across `sex = 0` and `sex = 1` for each (param, agedays) — the procedure pools sexes before taking the median and does not sex-stratify.

Procedure:

1. Assign each input row its integer year-of-age (`ageyears = floor(agedays / 365.25)`).
2. Cap `ageyears` at 19 so all observations at age 19 or older are pooled into a single terminal year.
3. For each `(param, ageyears)` cell, compute `median(sd.orig)` pooled across sexes (rows with `is.na(sd.orig)` are dropped first).
4. Treat each per-year median as anchored at midyear-age in days — `floor((ageyears + 0.5) * 365.25)`, i.e. year 0 → day 182, year 1 → day 547, …, year 19 → day 7122.
5. Linearly interpolate `sd.median` by day of age between midyear anchors (`approx()`) onto the full day range from `floor(min(ageyears) * 365.25)` through `floor(max(ageyears + 1) * 365.25)`.
6. For days outside the interpolation range, clamp to the nearest midyear median (`approx(..., rule = 2)` — constant extrapolation to both endpoints).
7. Replicate the resulting (param, agedays) → sd.median map for both `sex = 0` and `sex = 1`.

Sex is carried as a column on the output for downstream join compatibility with the main algorithm's `setkey(sd.recenter, param, sex, agedays)` rolling join, but the underlying medians are not sex-stratified. A caller who wants sex-specific recentering can run the same procedure manually per sex.

### `get_dop(param_name)`

Returns the designated other parameter (DOP) for a given param. WEIGHTKG → HEIGHTCM, HEIGHTCM → WEIGHTKG, HEADCM → HEIGHTCM. Used by both algorithms in same-day duplicate resolution and cross-parameter plausibility checks. Cheap; called inline wherever a DOP is needed.

---

## Common Pitfalls

Cross-cutting wrapper-level pitfalls that have caused real bugs. The two algorithm narratives have their own algorithm-specific pitfalls.

### Sorting and identifiers

- **Sort-key contract.** `setkey(data.df, subjid, param, agedays, internal_id)` is assumed throughout the pediatric branch (and the analogous adult sort key inside `cleanadult()`). If you add or modify rows in preprocessing, call `setkey()` again or downstream tiebreakers will be undefined.
- **`internal_id` is the algorithm identifier; user `id` is not.** Algorithm logic must use `internal_id` for sorting and tiebreaking. The user's `id` is preserved as-is into output and must not be renumbered or repurposed. Length mismatch between `id` and `measurement` raises a hard `stop()`.
- **Output row order in partial runs is deterministic but semantically meaningless.** `cleangrowth(cached_results = ..., ...)` returns rows sorted by `internal_id`, which is reassigned `1..K` over the changed-subject input only and collides with the cached `1..N`. Per-row data is correct; sort by user `id` if you need a stable order.

### Parallel execution

- **`parallel = TRUE` requires installed package.** Will fail with `load_all()` only — workers need `system.file()` access to `inst/extdata/`. Install with `devtools::install_local(".")` first.
- **Cluster export list must include all algorithm-internal helpers.** When adding a new helper that's called from inside `cleanchild()` or `cleanadult()`, add it to the `clusterExport()` list in Phase 5 — otherwise parallel runs will fail with "object not found" inside workers but serial runs will pass.
- **Batch-invariant setup belongs outside the outer loop.** Loading `exclude.levels`, the Tanner table, and any other per-run-but-not-per-batch constants once before the loop saves real time on multi-batch runs.

### Factor and code conventions

- **Factor pitfall.** Assigning an unlisted string to the `exclude` factor silently produces `NA`. Always verify new exclusion codes exist in `exclude.levels` first. The `exclude.levels` factor is built from `exclude.levels.peds` + `exclude.levels.adult` in Phase 4.
- **Duplicate function definitions.** R loads all `.R` files in `R/` alphabetically. If two files define the same function, the later one wins silently. Do not duplicate `.child_valid()` (defined in `child_clean.R`) elsewhere. Only one adult support file should exist in `R/` for the same reason — duplicate files have caused silent overwrites of `check_between` / `round_pt` in the past.
- **Adult `result` vs child `exclude` column names.** `cleanadult()` uses `result` as its exclusion column; `cleangrowth()` maps it to `exclude` in the combined output (Phase 14). Do not rename `result` inside `cleanadult()`.
- **Adult rounding tolerance** (0.12 cm/kg) is intentional and should not be removed. It allows for rounding to nearest half-cm and floating point precision.

### Cross-algorithm assembly

- **`Exclude-Missing` is set at different times for child vs adult.** Child: pre-dispatch in `cleangrowth()` (because `.child_valid()` needs it pre-set). Adult: post-dispatch. Don't move either side without understanding why.
- **`Exclude-C-Temp-Same-Day` must be resolved before output.** Phase 14 has a hard `stop()` guard. If a new step introduces or preserves temp-SDE codes, ensure Child Step 13 (or an explicit cleanup) resolves them.
- **Cross-algorithm columns are NA-filled.** Adult rows have `NA` for `tbc.sd`, `ctbc.sd`, `cf_rescued`, etc.; child rows have `NA` for `mean_ht`, `bin_result`. Downstream consumers must handle this.

### Reference data

- **Reference-table loading is gated by `ref_tables`.** When set, disk reads are skipped — but the closures must still be valid for the current input. Rebuild `gc_preload_refs()` if the reference files change.
- **HC reference data only goes to 5 years.** Phase 11 sets `Exclude-Not-Cleaned` for HC at `agedays > 3 * 365.25` (pre-recentering) and also at `agedays >= 5 * 365.25` (post-recentering), so all non-cleaned HC rows share a single consistent exclusion code. The 3y–5y band has valid z-score data for reporting but is not cleaned by the algorithm.

---

## Cross-References

- Child algorithm steps (5, 6, 7, 9, 11, 13, 15/16, 17, 19, 21, 22): `child-algorithm-reference.md`
- Adult algorithm steps: `adult-algorithm-narrative.md`
- CF rescue lookup table methodology: `cf-rescue-thresholds.md`
- Wtallow formulas (adult): `wtallow-formulas.md`
- Z-score reference: `z-score-calculation-for-growthcleanr-reference.md`
- Testing reference: `testing-reference.md`
- Walk-through procedure and the cross-cutting code-review checklist applied to each step or phase: `algorithm-walkthrough-procedure.md`
