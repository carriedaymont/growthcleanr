# GC Full Walk-Through — 2026-04-16 (fourth walk)

Extremely detailed walk-through after the recent cleanup wave
(legacy removal, Fenton switch, `debug` param, orphan cleanup,
pre-walkthrough resolved items). Goal: find any regressions or
missed issues introduced by those changes.

Recent changes in scope:
- `bc86770` Pre-walkthrough cleanup: add `debug` parameter, resolve deferred items
- `8a7f92c` Remove `prelim_infants` parameter from `read_anthro()`, update stress test frozen counts
- `4a2c232` Post-legacy cleanup: delete orphaned docs, remove dead code, regenerate roxygen
- `2dd07a2` Rename `fentlms_foraga` → `fent_foraga`, delete old `fentlms_forz` files
- `c9209bf` Switch Fenton z-scores to `fenton2025_ms_lookup_smoothed.csv` (CSD method)
- `c93bb34` Remove legacy pediatric algorithm and archived code
- `f02923e` Batching/parallel walkthrough: 3 bug fixes, 5 deferred cleanups

This walk has been split into sessions per
`algorithm-walkthrough-procedure.md`. This file tracks **Session 1
(wrapper, partial)** only — remaining wrapper chunks, all child
steps, and all adult steps will be covered in subsequent sessions.

---

## Pre-walkthrough setup

Baseline tests pass:
- test-cleangrowth.R: 65 PASS, 1097 warnings (data.table recycle in cleanadult)
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS (2 codetools lint warnings, plyr/foreach)
- test-child-parameters.R: 13 PASS (1 expected deprecation warning)
- test-adult-clean.R: 198 PASS, 130 warnings (data.table recycle)
- Adult regression harness: 1508/1508 at loosest, looser, tighter, tightest

Observations:
- Test count drop vs. 2026-04-16c walk (92→65 cleangrowth, 21→13 child-parameters) is explained by commit c93bb34 removing legacy algorithm tests.
- data.table recycle warnings (1097 in cleangrowth tests, 130 in adult tests) are new or newly-surfaced. Flagged as W1 below; the user has asked for a dedicated follow-up session focused on test warnings.

---

## Session 1 — wrapper parameter handling and file-level docs

### F1. Stale top-of-file algorithm structure block — FIXED
- **File/lines:** `child_clean.R` lines 26–52 (pre-edit)
- **Issue:** Block listed `STEP 3: Unit error recovery (optional)` and
  `STEP 4: Swapped parameters detection` — both legacy Stata-era top-level
  steps that have never been part of the R child algorithm (Step 3 was
  controlled by the now-removed `recover.unit.error` parameter; Step 4 was
  never implemented). Block also omitted the real algorithm-step entries
  for Early 13 (SDE-Identicals) and STEP 16 (Birth HT/HC EWMA2).
- **Fix:** Rewrote the block so it matches CLAUDE.md's canonical step list.
  Split preprocessing (handled in `cleangrowth()`) from algorithm steps
  (inside `cleanchild()`). Explicitly calls out which step numbers do not
  exist.

### F2. Stale version stamp in file header — FIXED
- **File/line:** `child_clean.R` line 10 (pre-edit)
- **Issue:** `VERSION: 2026-03-20` was set to the date of an earlier
  commit and went stale. File version is tracked by git; embedding a
  date in the source creates a maintenance burden.
- **Fix:** Removed the `VERSION:` line entirely.

### F3. Phantom `clean_infants_with_sde()` reference in supporting-functions
       block — FIXED
- **File/lines:** `child_clean.R` lines 94–102 (pre-edit)
- **Issue:** Supporting-functions comment block listed
  `clean_infants_with_sde()` as handling Step 13 final SDE resolution.
  No such function exists in the package. Block also listed
  `temporary_extraneous_infants()` twice (once in the primary list and
  again in "Additional support functions defined in this file"), and
  omitted `calc_and_recenter_z_scores()`, the EWMA helpers, `.child_exc()`,
  `.cf_rescue_lookup()`, `get_dop()`, `read_anthro()`, and `gc_preload_refs()`.
- **Fix:** Replaced the block with an accurate list of the supporting
  functions that actually live in this file, with one-line purpose
  annotations.

### F4. `include.carryforward` deprecation branch is NA-unsafe — FIXED
- **File/lines:** `child_clean.R` line 350 (pre-edit)
- **Issue:** `if (include.carryforward == TRUE)` evaluates to `NA` if the
  user passes `include.carryforward = NA`, and `if (NA) { ... }` throws
  "missing value where TRUE/FALSE needed." A user-supplied NA is never a
  valid value for a deprecated logical flag, but the right thing to do is
  skip the deprecation branch and proceed with `cf_rescue` — not error.
- **Fix:** Changed to `if (isTRUE(include.carryforward))`. Added a comment
  explaining the NA-safety rationale.

---

## Findings reviewed and retracted

### R1. Line 88 "assign higher ids to later values" — NOT a bug
- Initial reading flagged this as contradicting the age-dependent
  tiebreaking rule (at agedays == 0 the algorithm keeps the lowest
  internal_id; at agedays > 0 it keeps the highest).
- On review: the line-88 guidance is a **user-side ordering convention**
  (higher id = later measurement), while the tiebreaking rule is an
  **algorithm behavior** (which id wins at birth vs. non-birth). They are
  complementary, not in conflict. The birth rule depends on the user
  following this convention to get semantically correct results
  ("earliest measurement at birth, before postnatal fluid shifts").
- Action: no change. Leave line 88 comment as-is.

---

## Session 1 — deferred items

### D1. `data.table()` column-length recycle warnings in `cleanadult()`
- **Surfaced by:** `test-cleangrowth.R` (1097 warnings) and
  `test-adult-clean.R` (130 warnings). Message:
  `Item N has K rows but longest item has M; recycled with remainder.`
  Backtrace: `cleanadult` → `data.table(...)` → `as.data.table.list`.
- **Why deferred:** User has requested a dedicated follow-up session on
  test warnings. Also: adult algorithm is closed pending clinician
  validation, so fixes require explicit approval.
- **Next-session plan:** Locate the specific `data.table()` call(s) in
  `cleanadult` that produce the recycle; determine whether columns are
  meant to be aligned (silent misalignment bug) or scalars being
  broadcast (harmless but should be explicit via `list()`-construction
  or `rep(..., length.out = ...)`).

### D2. `cleanchild()` accepts but never uses `include.carryforward`
- **File/lines:** `child_clean.R` line 2463 declares parameter; body
  never references it; wrapper passes it at lines ~1152 and ~1183.
- **Why deferred:** Clean removal touches the `cleanchild()` signature
  and two dispatch sites in the wrapper. Safe but out of scope for
  Session 1. Flag it during Session 5 (child finishing / support
  functions) when we're already in `cleanchild()` structure.
- **Exact fix:** Remove `include.carryforward` from `cleanchild()`
  signature; remove the two matching arguments in the wrapper's
  `cleanchild(...)` and `ddply(..., cleanchild, ...)` dispatch calls.

### D3. `fent_foraga.csv` uncompressed file in `inst/extdata`
- **File:** `inst/extdata/fent_foraga.csv` (~11 KB)
- **Observation:** Only `inst/extdata/fent_foraga.csv.gz` is loaded
  by `child_clean.R` line ~806. The uncompressed `.csv` sits alongside
  but is not referenced anywhere in R code or docs.
- **Why deferred:** Need to check with Carrie whether the
  uncompressed CSV is intentionally shipped for human inspection.
  If not, remove it (contributes to CRAN tarball size which is
  already flagged in the CRAN prep checklist).
- **Exact fix:** Either delete `fent_foraga.csv` or add a comment in
  `inst/extdata/README.md` (if one exists) explaining why both formats
  are shipped.

### D4. Legacy "SD-score recentering" comment block describes old logic
- **File/lines:** `child_clean.R` lines 1066–1077 (post-edit may have
  shifted). The comment block walks through the old per-dataset median
  interpolation procedure (midyear-age interpolation, etc.) that was
  used when recentering medians were derived from the data. Current
  code (lines ~1083–1108) just reads a precomputed reference file
  (`rcfile-2023-08-15_format.csv.gz`).
- **Why deferred:** The comment is still useful as context for how the
  rcfile was constructed, but it reads as if the code were doing that
  interpolation live. A clarifying rewrite is worth doing, but it's
  larger than a fix-now and should go into a wrapper-narrative pass.
- **Exact fix:** Rewrite the block as "this code reads a precomputed
  recentering file, which was itself built by the procedure described
  below: ..." so the distinction between construction and use is clear.

### D5. `cleangrowth()` roxygen `@return` block is incomplete
- **File/lines:** `child_clean.R` lines 255–267.
- **Issue:** Lists only `id, exclude, param, cf_rescued, sd.orig_who,
  sd.orig_cdc, sd.orig, tbc.sd, ctbc.sd, final_tbc, mean_ht, bin_result`
  but the actual returned data.table also carries `internal_id, subjid,
  agedays, sex, line, bin_exclude, tri_exclude (when requested),
  cf_status/cf_deltaZ (when cf_detail = TRUE), debug columns (when
  debug = TRUE), sd.corr, sd.orig_uncorr, z.orig, uncorr, potcorr, v,
  v_adult, measurement, age_years, newbatch`.
- **Why deferred:** Needs careful enumeration and audit to decide which
  columns are part of the public contract vs. internal leakage that
  should be dropped before return.
- **Exact fix:** Before updating the roxygen, decide the intended
  public output schema; drop internal working columns from the final
  return (e.g., `v`, `v_adult`, `measurement`, `newbatch`,
  `age_years`), then update `@return` to match.

### D6. Partial-mode output ordering uses meaningless `internal_id` key
- **File/lines:** `child_clean.R` lines 1365–1373.
- **Observation:** In partial-run mode, `internal_id` is re-assigned as
  `seq_len(.N)` over only the changed-subject input (line 490), so
  values 1..K collide with values 1..N in `cached_results`. The rbind
  + `order(internal_id)` at lines 1366–1372 produces a deterministic
  but semantically meaningless row order. Per-row data is correct.
  Existing comment at 1362–1364 acknowledges it.
- **User note (2026-04-16):** internal_id is only meaningful within a
  subject, not across subjects, so the current behavior is acceptable.
  Defer.
- **If revisited:** sort by user `id` (contract-guaranteed unique) in
  the partial-merge branch, or renumber `internal_id` after rbind, or
  drop the sort.

### D7. Top-of-file comment at line ~50 still says "The algorithm
       proceeds through 22 steps" implicitly
- After F1's rewrite the block no longer counts to 22, but the
  notion of "22 steps" is a Stata-lineage convention. Nothing
  broken — just flagging that future docs should not reintroduce
  the "22 steps" framing.

---

## Session 1 — post-fix test results

All tests pass after F1–F4:

- test-cleangrowth.R: 65 PASS (1097 warnings — D1, unchanged, expected)
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 40 PASS (2 codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS (1 deprecation warning, unchanged)
- test-adult-clean.R: 198 PASS (130 warnings — D1, unchanged)
- Adult regression harness:
  - loosest: 1508/1508
  - looser: 1508/1508
  - tighter: 1508/1508
  - tightest: 1508/1508

No changes to the warning counts are expected at this stage; the
warnings live in `cleanadult()` (D1) and are the subject of the
dedicated follow-up session.
