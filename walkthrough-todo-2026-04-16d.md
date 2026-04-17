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

### D1. `data.table()` column-length recycle warnings in `cleanadult()` — RESOLVED 2026-04-16 (commit deb26a9)
- **Surfaced by:** `test-cleangrowth.R` (1097 warnings) and
  `test-adult-clean.R` (130 warnings). Message:
  `Item N has K rows but longest item has M; recycled with remainder.`
  Backtrace: `cleanadult` → `data.table(...)` → `as.data.table.list`.
- **Root cause (found in follow-up session):** Silent misalignment bug.
  Lines 245 and 322 indexed `h_extraneous`/`w_extraneous` with integer
  `internal_id`, which R treats as positional subscripting rather than
  name-based. Because `internal_id` is globally sequential (not 1..N
  within the subject's vector), this extended the vector with NAs and
  set TRUE at the wrong positions while the intended named entries
  stayed FALSE. The `extraneous` column assembled in `h_out`/`w_out`
  was therefore semantically wrong at every row but silently dropped
  from output by default (only surfaces with `include_extraneous=TRUE`,
  which `cleangrowth()` never exposes).
- **Fix:** Wrapped `internal_id` in `as.character()` at both sites so R
  uses named indexing. Two-character-level change.
- **Verification:** test-cleangrowth.R 65 PASS / 0 warnings (was 1097);
  test-adult-clean.R 198 PASS / 0 warnings (was 130); adult regression
  1508/1508 at all 4 permissiveness levels.
- **Audit scope:** Confirmed only two `data.table()` call sites exist
  in adult code (`adult_clean.R:1217` h_out, `adult_clean.R:1232` w_out);
  zero in `adult_support.R`. Also confirmed no stray `cat()`, `print()`,
  or `browser()` — only the intentional `message(capture.output(print(
  table(df$result))))` at line 1279 inside `!quietly`.

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

### D3. `fent_foraga.csv` uncompressed file in `inst/extdata` — RESOLVED 2026-04-16 (commit deb26a9)
- **File:** `inst/extdata/fent_foraga.csv` (~11 KB)
- **Observation:** Only `inst/extdata/fent_foraga.csv.gz` is loaded
  by `child_clean.R` line ~806. The uncompressed `.csv` sits alongside
  but is not referenced anywhere in R code or docs.
- **Resolution:** Carrie deleted the uncompressed CSV; only
  `fent_foraga.csv.gz` remains in `inst/extdata/`. Deletion included
  in commit deb26a9.

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

---

## Session 2 — cleangrowth() preprocessing + start of cleanchild()

Scope: input validation, batching setup, imperial→metric, LENGTHCM,
z-score setup, Step 2b GA correction, parallel setup, output
assembly, plus the opening of `cleanchild()` (Early Step 13 SDE-
Identicals + Step 5 temporary SDE flagging).

### F5. Stale `cleanbatch()` reference — FIXED
- **File/lines:** `child_clean.R` lines 1132, 1134 (pre-edit).
- **Issue:** Banner comment `# pediatric: cleanbatch (most of steps)`
  and in-line note `# NOTE: the rest of cleangrowth's steps are done
  through cleanbatch().` both referred to a function named
  `cleanbatch()`. No such function exists — the dispatcher is
  `cleanchild()`.
- **Fix:** Replaced both references with `cleanchild`.

### F6. Stale "NNTE removed" comment block — FIXED
- **File/lines:** `child_clean.R` lines 1130–1131 (pre-edit).
- **Issue:** Two-line comment `# NNTE removed: deterministic pre-
  filters in CF Step 6 and SDE Step 13 / are used instead.` sat
  immediately before the pediatric dispatch. `nnte` is a legacy
  Stata concept long removed from the R algorithm (per CLAUDE.md
  Terms Reference "always FALSE in current code"); the "removed"
  message has no present-day meaning and just adds noise.
- **Fix:** Deleted the two-line block.

### F7. Stale TODO-like comment `# NOTE: MAY WANT TO SUBSET HERE` — FIXED
- **File/lines:** `child_clean.R` line 756 (pre-edit).
- **Issue:** Dangling "may want to" note that had never been
  actioned and no longer carries useful context after the Step 2b
  fast-path subset (`pc <- copy(data.all[subjid %in% pc_ids])`)
  was introduced for potcorr subjects.
- **Fix:** Deleted.

### F8. Stale reference to "INFANT DOCS" — FIXED
- **File/line:** `child_clean.R` line 751 (pre-edit).
- **Issue:** Comment `# NOTE: SD SCORES IN CODE ARE Z IN INFANT
  DOCS -- USE sd.orig ONLY` referred to the legacy `prelim_infants`
  documentation set, which was removed with the legacy pediatric
  algorithm (commit `c93bb34`, 2026-04-16). Comment is stale.
- **Fix:** Deleted.

### F9. Stale "Changed smoothing from 2-4y to 2-5y to match Stata" — FIXED
- **File/line:** `child_clean.R` line 722 (pre-edit).
- **Issue:** Historical marker for a change made years ago to bring
  the R blending window into agreement with Stata. Stata cross-
  comparison is no longer an active concern; the remaining two
  comments (`smooth z-scores/SD scores between ages 2-5yo` +
  `older uses cdc, younger uses who`) already describe the current
  behavior.
- **Fix:** Deleted the one-line historical note. Kept the two
  descriptive lines below it.

### F10. Stale Stata `.do` line reference in Early Step 13 comment — FIXED
- **File/lines:** `child_clean.R` lines 2552–2556 (pre-edit).
- **Issue:** Header block for Early Step 13 SDE-Identicals included
  `# Stata removes identicals early (lines 179-215 in do file)
  before Steps 5/6` and `# Updated to use age-dependent id
  preference`. The specific Stata line-number pointer is a
  cross-language maintenance burden that will only drift further
  out of sync; "Updated to use age-dependent id preference" is a
  commit-message-style note that doesn't describe current behavior.
- **Fix:** Rewrote the header to describe the rationale in
  present tense ("Same-day identical values are excluded before CF
  detection so that CF logic only sees distinct repeated values"),
  retained the age-dependent-id rule description, dropped both
  stale lines.

### F11. Stale `(line 1022)` reference in cleanchild index comment — FIXED
- **File/lines:** `child_clean.R` line 2527 (pre-edit).
- **Issue:** Comment explaining why `index` is reassigned inside
  `cleanchild()` said `# index was created on full dataset before
  batching (line 1022), so batches have non-contiguous indices.`
  The actual `data.all[, index := 1:.N]` assignment in cleangrowth
  is at line 1048 (post-legacy-removal) and will drift again on
  future edits. Embedded line numbers are fragile.
- **Fix:** Rewrote the comment to reference the function
  (`cleangrowth()`) and reason rather than a line number.

---

## Session 2 — findings reviewed and retracted

### R2. `adult_scale_max_lbs` NaN edge case — NOT a bug worth fixing
- `!is.numeric(adult_scale_max_lbs) || adult_scale_max_lbs < 0` in
  the validation block could produce `if (NA)` for a NaN input
  (NaN is numeric, so first clause is FALSE; `NaN < 0` is NA).
- Review: a user supplying NaN as a scale cap is not a real
  workflow. The error thrown (`missing value where TRUE/FALSE
  needed`) surfaces the bad input loudly enough. Not worth
  wrapping with `isTRUE()`.
- Action: no change.

### R3. `data.all$ageyears <- data.all$agedays/365.25` uses $<- — NOT a bug
- At line 725 the column is assigned via `$<-` rather than the
  idiomatic `data.all[, ageyears := agedays/365.25]`. Works
  correctly for a data.table, just inconsistent with surrounding
  `:=` style.
- Review: this is purely stylistic; no correctness issue, no
  performance issue at this scale. Fix-as-you-go if the block
  gets rewritten for other reasons.
- Action: no change.

### R4. `v / 2.2046226` vs `v * 0.45359237` imperial→metric — NOT a bug
- The weight conversion at line ~668 uses a divisor instead of a
  multiplier; worth confirming the constants agree. 1/0.45359237
  = 2.20462262… so `v / 2.2046226` and `v * 0.45359237` produce
  effectively identical values to within floating-point
  precision.
- Action: no change.

---

## Session 2 — deferred items

### D8. Dead-code fallback in GA-correction ordering preference
- **File/lines:** `child_clean.R` lines 931–935 (post-F9 line
  numbers may have shifted by −1).
- **Observation:** Block:
  ```r
  pc[, sd.c_temp := sd.c]
  pc[potcorr == TRUE & ageyears_2b <= 2 & !is.na(unmod_zscore),
     sd.c := unmod_zscore]
  pc[potcorr == TRUE & ageyears_2b <= 2 & is.na(sd.c) & !is.na(sd.c_temp),
     sd.c := sd.c_temp]
  pc[, sd.c_temp := NULL]
  ```
  The third statement is unreachable. Trace:
  - If the second statement's predicate fails (unmod NA),
    `sd.c` is unchanged, so `sd.c == sd.c_temp` in all rows.
    Then `is.na(sd.c) & !is.na(sd.c_temp)` is FALSE & any →
    FALSE, or any & TRUE=FALSE → FALSE. Never fires.
  - If the second statement's predicate succeeds, it assigns
    `sd.c := unmod_zscore` which is non-NA (predicate requires
    `!is.na(unmod_zscore)`), so `is.na(sd.c)` is FALSE.
    Never fires.
- **Why deferred:** Dead code with no behavioral impact. Worth
  removing as a cleanup during the dedicated Step 2b narrative
  pass, alongside D4 (legacy recentering comment) and the rest
  of the 2b housekeeping.
- **Exact fix:** Delete the `sd.c_temp` snapshot/restore block
  entirely. The Fenton-preferred assignment on line ~932 already
  leaves pre-existing corrected WHO/CDC values in place for rows
  where Fenton is unavailable (because the predicate doesn't
  match, so `:=` is a no-op).

### D9. Velocity reference tables (Tanner, WHO) not in `gc_preload_refs()`
- **File/lines:** `child_clean.R` `gc_preload_refs()` at line
  ~1624 returns only `list(mtz_cdc_prelim, mtz_who_prelim)`.
  Tanner HT velocity is loaded once per `cleangrowth()` call at
  lines ~590–597; WHO velocity tables (HT and HC) are loaded
  inside `cleanchild()` once per batch.
- **Why deferred:** Correctness is fine. Efficiency gain is
  real but small — the velocity tables are tiny (< 1 kB each)
  and `fread` on gzipped inst/extdata files is fast. Bundling
  them into `gc_preload_refs()` would only save a few ms per
  call. Low priority relative to the 50K-rep simulation use
  case, which is already dominated by recentering / z-score
  closure work.
- **Exact fix:** (a) Extend `gc_preload_refs()` to also read and
  return `tanner_ht_vel`, `who_ht_vel_for_age`, and
  `who_hc_vel_for_age`. (b) Add matching optional parameters
  to `cleangrowth()`/`cleanchild()` so these can be passed in.
  (c) Keep the existing fallback `fread()` calls for
  backwards compatibility when the refs aren't pre-loaded.

### D10. WHO velocity tables reloaded per batch
- **File/lines:** `cleanchild()` — see Step 17 area for
  `fread()` calls on `who_ht_vel_for_age` / `who_hc_vel_for_age`.
- **Observation:** These files are read inside `cleanchild()`,
  so with `num.batches > 1` each batch re-reads the same
  files. Hoisting them into `cleangrowth()` (like Tanner) and
  passing as arguments would eliminate the redundant reads.
- **Why deferred:** Minor efficiency only. In practice
  `parallel = FALSE` is the current default, which means
  `num.batches = 1`. And the read is fast.
- **Exact fix:** Hoist the two `fread()` calls into
  `cleangrowth()` before the batch loop; add arguments to
  `cleanchild()` to receive them; update the `var_for_par`
  export list and `ddply(.paropts = ...)` if needed. Natural
  to combine with D9 (preload-ref extension).

### D11. Redundant `with(data.all, …)` in exclude initialization
- **File/lines:** `child_clean.R` line ~1054:
  ```r
  data.all[, exclude := factor(with(data.all, ifelse(
    is.na(v) | agedays < 0, 'Exclude-Missing', 'Include'
  )), levels = exclude.levels, ordered = TRUE)]
  ```
- **Observation:** Inside `data.all[, …]`, column names are
  already resolved against the data.table — the `with(data.all,
  …)` wrapper is redundant. Style inconsistency with the rest
  of the file, not a correctness issue.
- **Why deferred:** Purely cosmetic. Simplify during the
  narrative pass for the pre-dispatch preprocessing section.
- **Exact fix:** `data.all[, exclude := factor(ifelse(is.na(v)
  | agedays < 0, 'Exclude-Missing', 'Include'), levels =
  exclude.levels, ordered = TRUE)]`.

### D12. D4 re-confirmed from Session 1
- The legacy "SD-score recentering" prose block
  (`child_clean.R` ~lines 1066–1077 post-F-fixes) still reads as
  if the code were doing on-the-fly midyear-age interpolation,
  which it isn't. Already tracked as D4; no change in
  Session 2.

---

## Session 2 — post-fix test results

All tests pass after F5–F11:

- test-cleangrowth.R: 65 PASS, 0 warnings (D1 was resolved in
  commit deb26a9; the 1097-warning count from the pre-Session-1
  baseline is now gone once the package is reinstalled)
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 40 PASS (2 codetools warnings, baseline,
  unchanged)
- test-child-parameters.R: 13 PASS (1 deprecation warning, baseline,
  unchanged)

Adult regression harness not re-run for this session (no adult
code touched).
