# GC Full Walk-Through — 2026-04-17 (continuation)

Continues the walk started in `walkthrough-todo-2026-04-16d.md` (Sessions
1–3 covered wrapper preprocessing, Early 13, and Step 5). This file
covers **Sessions 4a and 4b**:

- Session 4a — Step 6 (Carried Forwards): CF detection, `cf_rescue`
  modes (`"standard"` / `"none"` / `"all"`), `.cf_rescue_lookup()`
  helper, optional `cf_status` / `cf_deltaZ` output columns, chained
  CF propagation.
- Session 4b — Step 7 (BIV): absolute and standardized BIV thresholds,
  preterm weight fix, `sd.extreme` / `z.extreme` parameter coverage.

Completed fixes and the Stata-reference bulk pass from Sessions 1–3
live in `walkthrough-todo-2026-04-16d.md`.

---

## Pre-session cleanup (2026-04-17, before Session 4a)

Closed the open deferreds carried forward from the prior walk. Fixed
in place where small and localized; moved to `CLAUDE.md → Known
Issues` where they required product decisions or were low-priority
efficiency work; dropped one item that no longer represented a
to-do.

**Fixed in code (this session):**

### D4 — Legacy "SD-score recentering" comment block — FIXED
- **File/lines:** `child_clean.R` ~1038–1049 (inside `cleangrowth()`,
  just above the `if (!is.data.table(sd.recenter))` branch).
- **Issue:** The prose block walked through the old per-dataset
  midyear-age interpolation procedure (3.a–3.e), reading as if the
  code immediately below were performing that procedure live. Current
  code just reads `rcfile-2023-08-15_format.csv.gz` and merges the
  precomputed medians. Item (f) on measprev/measnext was orphaned
  context about later-step conventions, unrelated to recentering.
- **Fix:** Rewrote the block so it names the file being read, states
  what the read+merge does, and then describes the construction
  procedure used to build that file (pointing at `sd_median()` for
  the current-day implementation). Dropped the orphaned (f) item.

### D16 — Stata-style Roxygen prose in `identify_temp_sde()` — FIXED
- **File/lines:** `child_clean.R` ~1995–2008 (Roxygen above the
  function) and the redundant `#` comment block at ~2012–2019 after
  `@noRd`.
- **Issue:** Roxygen used Stata notation (`exc_*==0`, `median_tbc*sd`,
  `median_tbcOsd`). The `#` block below `@noRd` duplicated what the
  Roxygen already described, minus one useful precondition note.
- **Fix:** Rewrote the Roxygen in R idiom (`tbc.sd`, `median.spz`,
  `median.dopz`, `Exclude-C-Temp-Same-Day`) describing the selection
  rule in present tense. Trimmed the `#` block below `@noRd` to just
  the identicals-precondition note; dropped the redundant summary.

**Moved to `gc-github-latest/CLAUDE.md → Known Issues → Open (wrapper)`
(new subsection):**

- **D5** (`cleangrowth()` `@return` roxygen incomplete) — moved
  verbatim; needs a product decision about which internal columns
  belong in the public output contract before the roxygen can be
  finalized.
- **D6** (partial-run output row ordering is deterministic-but-
  meaningless) — moved verbatim; already accepted as-is, noted
  mitigation options.
- **D9 + D10** (velocity reference tables not in `gc_preload_refs()`;
  WHO velocity tables reloaded per batch) — consolidated into a single
  entry since they'd be fixed together.

**Dropped:**

- **D7** ("Don't reintroduce the 22-steps framing"). F1 in the prior
  walk already rewrote the file-header step block; this was a
  note-to-self, not a to-do. No code or doc action needed.

**Also noted (not fixed):**

- `CLAUDE.md → Known Issues → Open (adult)` has two duplicated bullet
  pairs (deferred test gaps + performance `setkey`, each listed
  twice). Not touching — adult code is closed pending clinician
  validation, and this is a pure doc cosmetic issue. Worth a two-
  line cleanup during the next adult pass.
- `sd_median()` internal comment block (`child_clean.R` ~1932–1943)
  contains item (f) on measprev/measnext that is equally orphaned
  there. Left alone; can be handled alongside a narrative pass for
  `sd_median()` if/when that function's docs are rewritten.

---

## Pre-session test baseline

Reinstalled the package (`devtools::install_local(".", force=TRUE,
upgrade="never")`). All child test suites pass after the D4 and D16
edits:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 40 PASS (2 codetools warnings, baseline,
  unchanged)
- test-child-parameters.R: 13 PASS (1 expected deprecation warning,
  baseline, unchanged)

Adult test suites and regression harness not re-run (no adult code
touched).

---

## Session 4a — Step 6 (Carried Forwards)

Scope: CF detection, `cf_rescue` modes (`"standard"` / `"none"` /
`"all"`), `.cf_rescue_lookup()` / `.cf_get_thresholds()` helpers,
optional `cf_status` / `cf_deltaZ` output columns, chained CF
propagation, `cf_rescued` output column semantics.

### Fix-now items

### F17. Dead `cf_threshold` store on `data.df` in Step 6 — FIXED
- **File/lines:** `child_clean.R` lines ~2839–2840 (pre-edit),
  plus the `cf_threshold` entry in `cols_to_drop_6` at
  ~lines 2879–2884.
- **Issue:** After the rescue decision was made against a
  local `cf_dt` data.table, the standard-mode branch stored
  `cf_threshold` back onto `data.df` with a comment "Store
  thresholds on data.df for cf_detail output." The store was
  never read — `cf_detail` writes `cf_deltaZ` from `absdiff`,
  not from `cf_threshold`, and `cols_to_drop_6` dropped the
  column at the end of Step 6. Dead assignment with a
  misleading comment.
- **Fix:** Removed the assignment and comment; removed
  `"cf_threshold"` from the `cols_to_drop_6` string vector.
  Kept the NULL declaration at ~line 2456 (needed for
  `cf_dt[, cf_threshold := ...]` NSE in the local data.table).
- **Verification:** Reinstalled package; test-child-regression.R
  48 PASS, test-child-algorithms.R 40 PASS (2 baseline warnings).

### F18. Step 6 narrative rewrite — FIXED
- **File/sections in `child-gc-narrative-2026-04-13.md`:**
  - Intro "CF rescue" summary paragraph (lines ~117–125).
  - Carried-forward variables table (lines ~689–699).
  - CF rescue reason codes table (lines ~766–779).
  - Step 6 section (formerly lines ~1172–1366).
- **Issue:** The narrative described the pre-2026-04-14 CF
  implementation (fixed `< 0.05` / `< 0.1` thresholds, single-CF
  vs. adolescent multi-CF rescue, 4 rescue codes, control via
  deprecated `include.carryforward`, a `cf_string_length`
  variable that no longer exists, and stale Checklist findings
  that had already been resolved). `map_lgl` dplyr blocks the
  narrative warned about have been removed from code.
- **Fix:** Rewrote all four sections to describe the current
  lookup-based rescue logic:
  - Header table updated: rescue codes are now `Rescued` /
    `Rescued-All`, control parameter is `cf_rescue` +
    `cf_detail`, code location points at the actual line
    ranges and the `.cf_rescue_lookup()` / `.cf_get_thresholds()`
    helpers.
  - Overview enumerates the three modes and the optional
    `cf_detail` diagnostic output.
  - Added "Detection/rescue optimization" block describing
    the `run_cf_detection` skip path and the `cf_rescued`
    inconsistency it introduces (with pointer to the Known
    Issues entry).
  - Phase numbering updated: 1 Detection, 2 Temp SDE
    re-evaluation, 3 String construction, 4 Rescue (three
    modes), 5 SDE-Identical restoration + cf_detail + cleanup.
  - Variables and rescue-code tables updated to current
    columns (with `orig_ageday`, `cf_interval` added;
    `cf_string_length` dropped) and current code values
    (`Rescued`, `Rescued-All` only).
  - Checklist findings rewritten against 2026-04-17 code
    state; retracted earlier findings that had been
    resolved (uncleaned columns, `cf_string_length`,
    `mod 1`/`mod 0.5` redundancy, stale `map_lgl` blocks).

### Step 6 deferreds — all closed this session

### F19. Remove `run_cf_detection` optimization (closes D17) — FIXED
- **File/lines:** `child_clean.R` formerly lines 2570–2573
  (variable + outer `if`) and line ~2842 (outer `}`); ~270
  lines of detection/rescue block inside.
- **Issue:** The optimization skipped the entire CF block
  when `cf_rescue == "all" && !cf_detail`. Inclusion outcome
  was correct (nothing flagged → everything stays Include)
  but the `cf_rescued` column was `""` instead of
  `"Rescued-All"` in that path, creating an inconsistency
  dependent on `cf_detail`.
- **Fix:** Removed the `run_cf_detection` variable and the
  outer `if (...) {...}` wrapper via `sed` unindent of the
  ~270-line body by 2 spaces (one indent level). CF
  detection, string construction, and rescue now always run
  in all three modes.

### F21. `cf_rescue = "all"` now truly rescues every CF — FIXED
- **File/lines:** `child_clean.R` ~lines 2791–2800 (the
  `"all"` branch of the rescue dispatch).
- **Background:** After F19, `cf_rescue = "all"` still
  restricted rescue to CFs with `seq_win > 0` and
  `ageday_has_include == FALSE`, on the reasoning that CFs
  sharing an SPA with a non-CF Include should stay excluded
  to preserve the one-Include-per-SPA SDE invariant. Carrie
  pointed out that this preferentially excludes CFs even
  when they are more consistent with the trajectory than the
  non-CF same-day value — which defeats the "ignore CFs"
  intent of `"all"` mode.
- **Fix:** Changed the `"all"` branch to rescue every row
  with `exclude == "Exclude-C-CF"` regardless of
  `ageday_has_include` / `seq_win`. Multi-Include SPAs that
  result are resolved by Step 13 final-SDE logic
  downstream. `"standard"` mode is unchanged — it still
  gates rescue on the `ageday_has_include == FALSE` +
  `seq_win > 0` eligibility (only some CFs get rescued in
  that mode, and a conservative stance on multi-Include
  SPAs is appropriate there).
- **Test updates** reverted to the pre-F19 expectations:
  - `tests/testthat/test-child-algorithms.R` test 4 "CF
    rescue: all mode rescues all CFs" — restored to
    `expect_equal(n_cf, 0)` plus a new `Rescued-All > 0`
    check. Section 1 banner comment updated to describe
    the "every detected CF rescued" behavior with a note
    about Step 13 resolving multi-Include SPAs.
  - `tests/testthat/test-child-parameters.R` test 3
    "include.carryforward = TRUE keeps CFs, FALSE excludes
    them" — restored the `n_cf_kept == 0` assertion;
    comment updated to note Step 13's role.

- **Verification:** All 5 child test files pass after
  reinstall: test-cleangrowth.R 65 PASS, test-child-
  regression.R 48 PASS, test-child-edge-cases.R 28 PASS,
  test-child-algorithms.R 41 PASS (up from 40; +1
  Rescued-All assertion), test-child-parameters.R 13 PASS.
  Only the expected 2 codetools warnings (baseline) and 1
  deprecation warning (expected).

### D18. HEADCM CF rescue thresholds are HEIGHTCM-derived — CLOSED (moved to Known Issues)
- Added to `gc-github-latest/CLAUDE.md → Known Issues →
  Open (child)` as a documented design decision rather than
  a walkthrough deferred. HT thresholds are a reasonable
  starting point; HC-specific derivation would need a
  dedicated synthetic-HC population study. Entry points at
  `CF-exploration/cf-threshold-schemes.md` for future work.

### F20. CLAUDE.md wording for `cf_rescue = "all"` (closes D19) — FIXED
- Updated in two places in `gc-github-latest/CLAUDE.md`:
  - "CF rescue" overview block: `"all"` description now
    reads "every detected CF rescued, including CFs on an
    SPA that also has another Include. Step 13 final-SDE
    resolution handles any resulting multi-Include SPAs."
  - Configurable Parameters (child) table: `cf_rescue` row
    updated to match.
- Narrative (`child-gc-narrative-2026-04-13.md`) Phase 4
  "all" description and Step 6 overview also updated;
  the obsolete "Detection/rescue optimization" block and
  "Known issue: cf_rescued under 'all'" block were removed
  since F19 eliminated the underlying optimization.

### Deferreds surfaced while in Step 6 territory (not Step 6 scope)

These narrative issues were found while navigating around
Step 6 but belong to other steps / shared sections. Logged
here so we don't lose track, but will be closed in later
sessions (likely during a narrative-level cleanup pass).

### D20. Exclusion codes table in narrative has stale Step 3/4 rows
- **File/lines:** `child-gc-narrative-2026-04-13.md` lines
  ~744–746 and the surrounding note at ~762–764.
- **Issue:** Table includes `Unit-Error-High` (Step 3),
  `Unit-Error-Low` (Step 3), and `Swapped-Measurements`
  (Step 4). These steps never existed in the R child
  algorithm — confirmed by the file-header rewrite (F1,
  Session 1). The tail note says `exclude.levels` includes
  legacy codes for backward compatibility, which is no
  longer true after the 2026-04-16 legacy removal.
- **Exact fix:** Delete the three stale rows. Rewrite the
  trailing note to just explain that exclusion codes are
  not param-specific.

### D21. Narrative uses `valid()` where code now uses `.child_valid()`
- **File/lines:** `child-gc-narrative-2026-04-13.md`
  lines ~139–160 (the `valid()` function description section).
- **Issue:** Function was renamed `.child_valid()` during
  legacy cleanup to avoid colliding with the legacy `valid()`
  name. Narrative still says `valid(df)`, `valid(df,
  include.temporary.extraneous = TRUE)`, etc.
- **Exact fix:** Global rename `valid()` → `.child_valid()`
  within the narrative's function description and nearby
  prose. Check whether the flag-name conventions
  (`include.temporary.extraneous`, `include.extraneous`,
  `include.carryforward`) are still accurate in the current
  `.child_valid()` signature.

### D22. `valid()` description lists stale non-Exclude codes
- **File/lines:** `child-gc-narrative-2026-04-13.md` line
  ~159 — "Non-Exclude codes like 'Include',
  'Unit-Error-High', 'Unit-Error-Low', 'Swapped-Measurements'
  pass through."
- **Issue:** Those three non-`Include` codes never existed
  in the R code. Misleading.
- **Exact fix:** Reduce the sentence to just "Non-Exclude
  codes like `"Include"` pass through."

### D23. Narrative references "known bug where the outer wrapper is currently defeated"
- **File/lines:** `child-gc-narrative-2026-04-13.md` lines
  ~193–196 (Preprocessing step 9 callout), ~242, ~315, ~367
  (Batching and Dispatch section).
- **Issue:** The outer batching wrapper bug was fixed
  (per `gc-github-latest/CLAUDE.md → Known Issues → Fixed`
  history; batching now works with `parallel = TRUE`
  producing identical results to sequential on 77,721
  rows). The narrative still points readers at a "known
  bug" for something that's resolved.
- **Exact fix:** Verify status against current code, then
  either remove the "known bug" language or reduce it to
  a historical footnote.

---

## Session 4a — post-fix test results

Reinstalled package; ran the child test suites:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (up from 40 — +1
  Rescued-All assertion in Test 4; 2 codetools warnings,
  baseline, unchanged)
- test-child-parameters.R: 13 PASS (1 expected deprecation
  warning, baseline, unchanged)

Adult test suites and regression harness not re-run (no adult
code touched).

---

## Session 4b — Step 7 (BIV)

Scope:
- Absolute BIV limits (per-param, per-age)
- Standardized BIV limits
- Preterm weight fix (<0.2 kg first year, <1 kg after)
- Fate of the `sd.extreme` / `z.extreme` parameters

### Pre-session baseline tests

After Session 4a, re-ran the child test suites against the
installed package:

- test-cleangrowth.R: 65 PASS
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings)
- test-child-parameters.R: 13 PASS (1 expected deprecation warning)

Adult suites not re-run; no adult code touched in this session.

### Fix-now items

### F22. `sd.extreme` / `z.extreme` replaced by eight per-cell `biv.z.*` parameters — FIXED
- **Background.** `sd.extreme = 25` and `z.extreme = 25` were
  declared in `cleangrowth()` and `cleanchild()` signatures
  and documented in Roxygen as Step 7 BIV thresholds, passed
  through the wrapper/dispatch, but never referenced anywhere
  in `cleanchild()`. All standardized-BIV cutoffs were
  hardcoded (`-25`, `-15`, `22`, `8`, `15`). They were live
  in the legacy pediatric algorithm (removed 2026-04-16) and
  never rewired. Test 4 in `test-child-parameters.R` asserted
  `expect_gte(n_ewma_strict, n_ewma_lenient)` on two values
  of `sd.extreme`; because the parameter was dead, both runs
  produced identical counts and the test passed vacuously.
- **Design decision (this session).** Remove `sd.extreme` and
  `z.extreme` outright; replace with one scalar per
  param/direction/age cell (eight in total). Defaults preserve
  the prior hardcoded behavior. Standardized-z scope only —
  absolute-BIV raw-value limits remain hardcoded.
- **Parameters added** (all to `cleangrowth()` and
  `cleanchild()` signatures, all Roxygen-documented):

  | Parameter | Default | Applied at |
  |---|---|---|
  | `biv.z.wt.low.young` | -25 | WEIGHTKG, `ageyears < 1` |
  | `biv.z.wt.low.old` | -15 | WEIGHTKG, `ageyears >= 1` |
  | `biv.z.wt.high` | 22 | WEIGHTKG, all ages |
  | `biv.z.ht.low.young` | -25 | HEIGHTCM, `ageyears < 1` |
  | `biv.z.ht.low.old` | -15 | HEIGHTCM, `ageyears >= 1` |
  | `biv.z.ht.high` | 8 | HEIGHTCM, all ages |
  | `biv.z.hc.low` | -15 | HEADCM, all ages |
  | `biv.z.hc.high` | 15 | HEADCM, all ages |

- **Files changed:**
  - `R/child_clean.R` — removed `sd.extreme`/`z.extreme` from
    `cleangrowth()` signature + Roxygen (lines 163–166,
    299–300); added the 8 new params in the same slots with
    new Roxygen. Replaced two pass-through blocks
    (lines ~1160–1166, ~1190–1196) with the 8 new
    pass-throughs. Removed `sd.extreme, z.extreme` from
    `cleanchild()` signature (lines 2472–2473); added the 8
    new params in the same slot. Step 7 standardized-BIV
    block now uses the param variables instead of hardcoded
    numeric literals.
  - `tests/testthat/test-child-parameters.R` — Test 4
    rewritten to exercise `biv.z.wt.high = 2` vs default
    `22` with `expect_gte(n_biv_strict, n_biv_default)`.
    Test title/comment updated accordingly.
  - `testing-reference.md` line 190 — row updated from
    "sd.extreme parameter" to "biv.z.wt.high parameter".
  - `gc-github-latest/CLAUDE.md` — two rows in the
    Configurable Parameters (child) table replaced with
    eight new rows; `test-child-parameters.R` scope row
    updated (sd.extreme → biv.z.wt.high).
  - `child-gc-narrative-2026-04-13.md` — parameters-table
    rows 842–843 replaced with eight new rows; Step 7
    "Standardized BIV thresholds" tables updated to show
    parameter + default per cell; "Configurable parameters
    in scope for Step 7" subsection rewritten.
  - `man/cleangrowth.Rd` — regenerated via
    `devtools::document()`.
- **Verification.** After reinstall, all 5 child test suites
  pass (65 / 48 / 28 / 41 / 13). Regression suite
  (test-child-regression.R) confirms defaults preserve prior
  behavior — no count or id mapping shifted.

### F23. Step 7 narrative rewritten for "current state only" — FIXED
- **Background.** `child-gc-narrative-2026-04-13.md` Step 7
  section (lines ~1593–1724 pre-edit) contained history-
  flavored content: a "BUG FIXED (2026-03-20)" block about
  the prior `v < 0.2` gap, Checklist items describing
  "BUG FIXED", "Cleaned up", and "v == 0 handler removed
  (2026-04-13)"; Checklist numbering skipped 3 (2 → 4);
  text referred to a non-existent `abs_biv` variable.
- **Fix:** Rewrote the section top-to-bottom in present
  tense. Removed all "BUG FIXED" / "Cleaned up" / "Removed"
  prose. Replaced `abs_biv` mentions with accurate
  descriptions of `biv_pattern <- "^Exclude-C-BIV$"`.
  Rephrased the "Step 5 rerun" wording (Overview and
  temp-SDE section) so it accurately describes
  `identify_temp_sde()` being rerun, not all of Step 5.
  Updated the summary table's "Code location" from "See
  code" to a concrete line range. Added a
  "Configurable parameters in scope for Step 7" subsection
  listing the eight new `biv.z.*` parameters (paired with
  F22). Checklist renumbered 1–8 consecutively.

### F24. Standardized-BIV guard description corrected — FIXED
- **File/lines:** `child-gc-narrative-2026-04-13.md`
  pre-edit lines 1653–1655 and 1711 referenced an
  "`exclude != abs_biv` guard" that does not exist in
  code. Code uses `!grepl(biv_pattern, exclude)` where
  `biv_pattern <- "^Exclude-C-BIV$"`.
- **Fix:** Rewrote the Standardized BIV thresholds
  preamble and Checklist item 4 to describe the actual
  guard, its rationale (valid_set was computed before the
  absolute block), and which rows can or cannot be
  overwritten by Step 7 (only `Exclude-C-Temp-Same-Day`
  rows are both in `valid_set` and have a non-`Include`
  code at entry; no other non-Include code can be
  overwritten).

### F25. "Step 5 rerun" phrasing in narrative Overview — FIXED
- **File/lines:** `child-gc-narrative-2026-04-13.md`
  pre-edit line 1615 said "After both, temp SDEs are
  re-evaluated (Step 5 rerun)." Only
  `identify_temp_sde()` is rerun; the rest of Step 5 is
  not.
- **Fix:** Overview now reads "After both blocks, the
  temp-SDE identification is rerun: all rows currently
  flagged `Exclude-C-Temp-Same-Day` are reset to
  `Include`, then `identify_temp_sde()` is called
  against the post-BIV state."

### F26. Stata-style 7d comment rewritten in R idiom — FIXED
- **File/lines:** `R/child_clean.R` pre-edit line 2970 —
  `# 7d.  Replace exc_*=0 if exc_*==2 & redo step 5
  (temporary extraneous)` (Stata notation: `exc_*==0`
  for Include, `exc_*==2` for Temp-Same-Day).
- **Fix:** Replaced with "`# 7d. Re-evaluate temp SDEs
  after BIV exclusions: reset all Exclude-C-Temp-Same-Day
  rows to Include, then rerun identify_temp_sde().
  Rationale: absolute or standardized BIV may have
  excluded the prior temp-SDE "keeper" on an SPA;
  another value in that SPA should now be flagged
  instead.`"

### F27. Inaccurate "overwriting non-temporary codes" note replaced — FIXED
- **File/lines:** `R/child_clean.R` pre-edit lines
  2932–2935 — inline comment claimed "This is the only
  step where an exclusion code can overwrite a non-
  temporary exclusion code (e.g., this could overwrite
  CF codes)." In practice, `valid_set` is built from
  `.child_valid(data.df, include.temporary.extraneous =
  TRUE)`, which only admits `Include` rows and
  `Exclude-C-Temp-Same-Day` rows. CF-excluded rows (and
  any other non-temp Exclude) are not in `valid_set`
  and cannot be overwritten by Step 7.
- **Fix:** Replaced with a two-line note explaining why
  the `!grepl(biv_pattern, exclude)` guard exists (it
  skips rows just assigned `Exclude-C-BIV` by the
  absolute block, because `valid_set` was computed
  beforehand) without making the spurious "overwrite CF"
  claim.

### Session 4b deferreds (not addressed)

### D24. Step 7 test coverage gaps — DEFERRED
- **Scope.** Standardized-BIV cells have no targeted unit
  tests (the only direct BIV tests in `test-child-edge-
  cases.R` all exercise absolute-BIV limits). Uncovered:
  - Preterm minimum weight (`v < 0.2` at
    `agedays <= 365`)
  - Post-365d minimum (`v < 1` at `agedays > 365`)
  - Birth maxima (`v > 10.5` WT, `v > 65` HT,
    `v > 50` HC at `agedays == 0`)
  - `v > 35` WT at `ageyears < 2`
  - Standardized WT `sd.orig_uncorr < biv.z.wt.low.young`,
    `< biv.z.wt.low.old`, `> biv.z.wt.high`
  - Analogous HT cells, HT `> biv.z.ht.high` (default 8)
  - HC standardized `< biv.z.hc.low`, `> biv.z.hc.high`
  - 7d temp-SDE re-evaluation after BIV (constructed case
    where an SPA keeper gets BIV'd and the other value
    picks up the temp-SDE flag)
- **Why deferred.** Coverage gaps; medium effort; no
  current accuracy issue (regression tests confirm
  behavior unchanged under defaults, and the new Test 4
  exercises at least one `biv.z.*` parameter path).
- **Where fixed later.** Likely in a dedicated test-
  coverage session for Step 7, alongside similar passes
  for other walked steps.

---

## Session 4b — post-fix test results

After reinstalling the package and regenerating
`man/cleangrowth.Rd` via `devtools::document()`:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 codetools warnings,
  baseline, unchanged)
- test-child-parameters.R: 13 PASS (1 expected deprecation
  warning from the `include.carryforward` test, baseline,
  unchanged)

Adult suites and regression harness not re-run (no adult
code touched).

---

## Session 5 — Step 9 (Evil Twins)

Scope: `child_clean.R` ~lines 3035–3147 (main logic including
9d temp-SDE rerun) and ~lines 2242–2271 (`calc_otl_evil_twins`
support function); narrative `child-gc-narrative-2026-04-13.md`
~lines 1775–1899.

### Pre-session baseline tests

After Session 4b, re-ran the child test suites against the
installed package:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings)
- test-child-parameters.R: 13 PASS (1 expected deprecation warning)

Adult suites not re-run; no adult code touched in this
session.

### Fix-now items

### F28. Step 9 narrative rewritten for "current state only" — FIXED
- **Background.** `child-gc-narrative-2026-04-13.md` Step 9
  section (pre-edit lines ~1775–1899) had two problems: (1)
  two Checklist findings described already-removed code —
  item 1 flagged `_rounded` aliases in `calc_otl_evil_twins`
  that grep confirms do not exist in the current function
  body (lines 2260–2263 use bare `tbc_next_diff` /
  `tbc_prev_diff` / `ctbc_next_diff` / `ctbc_prev_diff`);
  item 2 flagged a "non NNTE values" comment that is not
  present anywhere in `child_clean.R`. (2) The narrative
  still referenced `id` as the sort tiebreaker in the
  Setup pseudocode (line 1830) and the Checklist Boundaries
  entry (line 1881), but the code has used `internal_id` for
  sort tiebreaking since the 2026-04-12 internal_id rename.
  In addition, the Code-location cell was the same
  `"See code"` placeholder pattern Session 4b F23 replaced
  with concrete line ranges, and the Overview ("This step is
  now added…") framed the step as a change from a prior
  implementation rather than a current-state description.
- **Fix.** Rewrote the entire Step 9 section:
  - Code-location cell now reads "`child_clean.R`
    ~lines 3035–3142 (main logic, including 9d temp-SDE
    rerun); support function `calc_otl_evil_twins()` at
    ~lines 2242–2271".
  - Overview rewritten in current-state prose: describes
    why Evil Twins is needed (adjacent extremes distort
    each other's EWMA so EWMA-based exclusions alone can
    miss them) without the "original pediatric algorithm"
    framing.
  - 9a–9d subsections made explicit: 9a (global early-exit
    check via one `calc_otl_evil_twins()` call), 9b
    (per-group while-loop), 9c (bulk apply exclusions), 9d
    (cleanup + temp-SDE rerun). Removed the old "9D" inner
    label from the narrative's while-loop description and
    rephrased as numbered steps a–f within 9b.
  - Setup pseudocode updated to
    `order(subjid, param, agedays, internal_id)`.
  - Tiebreaker documented as "Lowest `internal_id`"
    (replaces "Lowest `id`") in both 9b step c and the
    Checklist Boundaries entry.
  - Added a "Configurable parameters in scope for Step 9"
    subsection noting there are none (the `> 5` OTL
    threshold is hardcoded in `calc_otl_evil_twins()`).
  - Checklist renumbered 1–8 consecutively after removing
    stale items 1 and 2; no other Checklist content
    removed.

### F29. Stale "param-specific code" inline comment removed — FIXED
- **File/lines:** `R/child_clean.R` pre-edit line 3043 —
  `# Evil Twins exclusion — param-specific code applied at
  final assignment`.
- **Issue.** Since the 2026-04-14 exclusion-code rename,
  child codes are not param-specific. `.child_exc()` (line
  137) accepts a `param_val` argument but ignores it and
  returns `paste0("Exclude-C-", suffix)`. The `param`
  argument passed at line 3137's bulk assignment
  `.child_exc(param, "Evil-Twins")` is carried for call-site
  consistency with ~40 other call sites, but the value is
  discarded. The comment misleadingly implied otherwise.
- **Fix.** Removed the comment as part of the broader
  rewrite of the Step 9 header block in F30 below; no
  replacement needed since the header block now describes
  the current step contents rather than the individual
  assignment mechanics.

### F30. "Session 16" changelog note rewritten as current-state rationale — FIXED
- **File/lines:** `R/child_clean.R` pre-edit lines 3036–3041
  — the Step 9 header comment block. Block contained (a)
  a past-tense description framing Evil Twins as a fix for
  a weakness in the "original pediatric growthcleanr
  algorithm" that "is now added", and (b) a "Optimization
  (Session 16): Restructured from single-closure operating
  on ALL valid rows to by-(subjid, param) group processing"
  changelog entry.
- **Issue.** Per walkthrough procedure item 10 (Stale
  content — references to prior versions), session-numbered
  changelog entries and "this step is now added" prose
  should not live in the implementation file. The rationale
  each conveyed is still useful — Evil Twins' reason-to-
  exist and why per-group processing is used — but needs
  to be expressed in current-state prose.
- **Fix.** Replaced the block with two current-state
  paragraphs: one describing what Evil Twins does and why
  (adjacent extremes distort each other's EWMA so
  EWMA-based exclusions alone can miss them; Step 9 runs
  before EWMA processing to catch them) and one describing
  the per-group processing strategy and its rationale
  (groups are small and independent; local-copy processing
  avoids repeatedly copying the full valid dataset per
  while-loop iteration; same pattern in Steps 11, 15, 17).
  Also updated the preceding sort-key comment block
  (pre-edit lines 3045–3051) to reference `internal_id` in
  its rationale prose; the actual `order()` call already
  used `internal_id`.

### F31. 9A/B/C/D/F step-letter labels unified to lowercase 9a/9b/9c/9d — FIXED
- **File/lines:** `R/child_clean.R` pre-edit line 3062
  (`# 9A/B/C`), line 3081 (no label; added `9b.`), line
  3098 (`# 9D:`), line 3142 (`# 9F.  redo temp
  extraneous`).
- **Issue.** Uppercase labels inherited from the Stata
  implementation; the gap at `9E` made the sequence
  self-inconsistent. Narrative now uses contiguous
  `9a / 9b / 9c / 9d`.
- **Fix.** Re-labeled as: `9a.` global early-exit check
  (line 3062), `9b.` per-group loop (added on line 3081),
  `9c.` median and distance-from-median inside the while
  loop (line 3098, was `9D`), `9d.` temp-SDE re-evaluation
  (line 3142, was `9F`). `9d`'s comment also gained an
  explicit rationale paragraph mirroring the Step 7 F26
  pattern.

### F32. `Lowest id` comment in worst-OTL tiebreaker updated — FIXED
- **File/lines:** `R/child_clean.R` pre-edit line 3106 —
  the tiebreaker hierarchy comment listed "3. Lowest id
  (deterministic)".
- **Issue.** The actual `order()` expression on pre-edit
  line 3109 uses `otl_rows$internal_id` (has since the
  2026-04-12 internal_id rename); comment was stale.
- **Fix.** Updated to "3. Lowest internal_id
  (deterministic)". Comment and code now agree.

### F33. All code line-number references removed from narrative — FIXED
- **Background.** Carrie flagged that line-number
  references in `child-gc-narrative-2026-04-13.md` drift
  every time code is edited; most had already been purged
  in earlier passes, but the Code-location cells in step
  summary tables still used ranges like `~lines 2565–2880`.
  Session 4b F23 had just added one such range for Step 7,
  and Session 5 F28 added another for Step 9 — both
  premised on the F23-era convention that has since been
  reconsidered.
- **Issue.** A grep for line-range patterns found six
  occurrences, all in Code-location table cells:
  - Step 2b (GA correction)
  - Early Step 13 (SDE-Identicals)
  - Step 5 (Temporary SDE)
  - Step 6 (Carried Forwards)
  - Step 7 (BIV)
  - Step 9 (Evil Twins)
- **Fix.** Rewrote all six cells to reference the
  containing function and file only (e.g., "Inline in
  `cleanchild()` in `child_clean.R`; support function
  `calc_otl_evil_twins()` defined earlier in
  `child_clean.R`"). No line numbers. Historical line
  references in session logs (this file and earlier
  walkthrough-todo files) are left intact — those are
  archival records of what the code looked like at the
  time of each session, not living narrative.
- **Convention going forward.** Narrative Code-location
  cells name the file and function(s); no line numbers.

### Session 5 deferreds (not addressed)

### D25. Step 9 targeted test coverage gaps — DEFERRED
- **Scope.** `test-child-algorithms.R` Section 2 (Tests
  5–8) exercises Evil Twins detection through realistic
  unit-error injection (HT×2.54). No targeted unit tests
  for the specific code paths:
  - OTL threshold boundary (diff exactly 5.0 → NOT OTL)
  - Cross-`(subjid, param)` adjacency guard (a HT row
    immediately followed by a WT row in the sorted table
    must not be compared as adjacent)
  - Tiebreaker hierarchy (constructed case where two OTL
    rows tie on `med_diff`, and then on `|tbc.sd|`, to
    exercise the `internal_id` final tiebreaker)
  - `sp_count_9 <= 2` pre-filter (subject-param with 2
    total rows, both OTL — must not be flagged)
  - `calc_otl_evil_twins()` edge case for `nrow(df) < 2`
    returning all `otl = FALSE`
  - 9d temp-SDE re-evaluation after an Evil Twins
    exclusion removes the prior SPA "keeper"
- **Why deferred.** Coverage gaps; medium effort; no
  current accuracy issue. Regression tests already
  confirm the behavior under defaults, and the existing
  Evil Twins/OTL tests exercise realistic end-to-end
  scenarios. Same rationale as Session 4b D24 (Step 7
  coverage gaps).
- **Where fixed later.** Dedicated test-coverage session
  for walked steps.

---

## Session 5 — post-fix test results

After reinstalling the package:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 codetools warnings,
  baseline, unchanged)
- test-child-parameters.R: 13 PASS (1 expected deprecation
  warning, baseline, unchanged)

Adult suites and regression harness not re-run (no adult
code touched).

---

## Session 6 — Step 11 (EWMA1 Extreme)

Scope: `child_clean.R` Step 11 block (main logic, per-group
closure, global iteration loop, end-of-step temp SDE refresh)
and the surrounding `debug` parameter plumbing; narrative
`child-gc-narrative-2026-04-13.md` Step 11 section plus the
parameters table and output-columns list; `gc-github-latest/
CLAUDE.md` and `NEWS.md`; `__Pipeline/CLAUDE.md` header.

### Pre-session baseline tests

After Session 5, re-ran the child test suites against the
installed package:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings)
- test-child-parameters.R: 13 PASS (1 expected deprecation warning)

Adult suites not re-run; no adult code touched in this session.

### Fix-now items

### F34. `debug` parameter + `ewma1_it1.*` capture block removed — FIXED
- **Background.** Commit `bc86770` (Pre-walkthrough cleanup,
  2026-04-16) added `debug = FALSE` to `cleangrowth()` and
  `cleanchild()` to gate a pre-existing capture block in
  Step 11 (originally unconditional — flagged as
  walkthrough-todo-2026-04-15 item 32 / 2026-04-16b). The
  documented behaviour: when `debug = TRUE`, emit six
  `ewma1_it1.*` columns (`ewma_all`, `ewma_before`,
  `ewma_after`, `dewma_all`, `dewma_before`, `dewma_after`)
  from the first EWMA1 iteration to aid diagnosis.
- **Issue.** The capture block (pre-edit `child_clean.R`
  lines 3275–3287) and the final-output column filter
  (pre-edit lines 4844–4853) are both dead code. The EWMA
  columns live on a local `df = copy(.SD)` inside the Step 11
  per-group closure; only `df$exclude` is returned to
  `data.df`. The block's gate `col %in% names(data.df)` is
  therefore always FALSE, and `get(col)` would resolve to the
  NULL declaration at `child_clean.R` line 2494 even if the
  gate ever passed. Verified empirically: running
  `cleangrowth(..., debug = TRUE)` on a 200-subject stress-
  data sample (29 real `Exclude-C-Traj-Extreme` results)
  produced no `ewma1_it1.*` columns in the output.
- **Design decision (this session).** Remove the `debug`
  parameter outright rather than fix the capture block.
  Rationale: the feature was freshly added, never verified
  to work, has no current consumer, and the walkthrough is
  focused on current-state correctness. A systematic
  diagnostic mechanism can be reintroduced later if needed
  (e.g., return-list restructuring from the per-group
  closure, or a separate diagnostic entry point).
- **Files changed:**
  - `R/child_clean.R` — removed `debug = FALSE` from the
    `cleangrowth()` signature and the corresponding
    `@param debug` Roxygen block; removed from both Roxygen
    pass-throughs (sequential `cleanchild()` call and
    `ddply(..., cleanchild, ...)` call); removed `debug = FALSE`
    from the `cleanchild()` signature; removed the Step 11
    capture block (`if (debug && iteration == 1) { ... }`) and
    the final-output `ewma1_it1.*` column filter in
    `cleanchild()`.
  - `child-gc-narrative-2026-04-13.md` — removed the
    parameters-table row for `debug`; removed the
    `ewma1_it1.*` diagnostics bullet from the Step 22 output
    section; removed the Step 11 "Debug columns preserved"
    subsection (folded into the overall Step 11 rewrite
    below).
  - `gc-github-latest/CLAUDE.md` — removed `debug` row from
    the Configurable Parameters (child) table; removed
    `ewma1_it1.*` reference from the `cleangrowth()` `@return`
    Known-Issues entry; rewrote the 2026-04-16 Pre-walkthrough
    cleanup entry to note that the `debug` parameter added
    in that round was subsequently removed here.
  - `NEWS.md` — removed the 3.0.0 "Added" bullet for the
    `debug` parameter.
  - `__Pipeline/CLAUDE.md` — refreshed "Last updated" header
    to note the debug-parameter removal (edited with
    standing permission — see gc CLAUDE.md header block).
  - `man/cleangrowth.Rd` — regenerated via
    `devtools::document()` after the Roxygen edits.
- **Verification.** After reinstall, all 5 child test suites
  pass: 65 / 48 / 28 / 41 / 13.

### F35. Step 11 header block rewritten for current-state — FIXED
- **Background.** Pre-edit `child_clean.R` lines 3141–3148
  contained a "Restructured to use global iterations for
  efficiency / Key changes: 1. …" changelog block describing
  what changed relative to an earlier implementation. Pre-edit
  lines 3150–3157 contained Stata-derived enumeration
  (`# 11. / a. / b. / i.`) with `exc_*==0` / `exc_*==2`
  notation. Same patterns as Session 5 F30 (Step 9 "Session
  16" changelog) and Session 4b F26 (Stata notation).
- **Fix.** Replaced both blocks with a single current-state
  description at the top of Step 11: what the step does,
  which rows participate, why per-pass single-exclusion + a
  global while loop, and a pointer to the 11a/11b/11c
  sub-sections introduced in F37.

### F36. Historical "Fixed cdewma sign" comment and `lowest id` wording — FIXED
- **File/lines:** pre-edit `child_clean.R` lines 3234–3235
  ("# Fixed cdewma sign for negative outliers / # Was:
  c.dewma.all > 3.5 for negative case (wrong - should be
  < -3.5)") and pre-edit line 3256 ("# Select worst: highest
  abs(tbc.sd + dewma.all), lowest id as tiebreaker").
- **Fix.** The "Fixed / Was:" block was replaced with a
  current-state two-line note describing what the guard does
  (`c.dewma.all` is tested in the same direction as
  `dewma.all`). The "lowest id" wording was corrected to
  "lowest internal_id" to match the `-internal_id` expression
  actually used in the `order()` call — same F32 correction
  pattern as Step 9.

### F37. Explicit 11a / 11b / 11c sub-sections added — FIXED
- **Background.** Step 11 previously had no labelled
  sub-sections; its natural structure (pre-filter +
  while-loop + end-of-step temp SDE refresh) was visible in
  the code only via spacing and comment cues. Session 5 F31
  established an `Xa` / `Xb` / `Xc` / `Xd` convention for
  Step 9; Carrie asked for that to extend to Step 11.
- **Fix.** Added three labelled sections in `child_clean.R`:
  - `11a. Pre-filter / setup` — count filter, extreme
    filter, and `subj_with_sde` capture.
  - `11b. Iteration loop` — global while loop (the per-group
    closure + per-pass temp-SDE recalc live inside).
  - `11c. End-of-step temp SDE refresh` — global reset and
    rerun of `identify_temp_sde()` before Step 13. The
    existing comment already noted the `exclude_from_dop_ids`
    distinction vs. Step 13; kept that wording.
  Narrative Step 11 section rewritten to mirror these
  labels, with the EWMA computation / exclusion criteria /
  worst-value selection / iteration bookkeeping as
  sub-sections under 11b.

### F38. Step 11 narrative rewritten for current-state only — FIXED
- **Background.** Pre-edit `child-gc-narrative-2026-04-13.md`
  Step 11 section (lines ~1937–2095):
  - `Code location` cell said "See code" — same placeholder
    F28/F33 replaced with a function+file reference in
    earlier sessions.
  - Worst-value-selection paragraph said "Lowest `id` breaks
    ties" (code uses `-internal_id`). F32 pattern.
  - "Debug columns preserved" subsection described a feature
    that doesn't work (retired here with F34).
  - Checklist items 1–3 described already-completed
    cleanup actions ("Removed `nnte` from `.SDcols`",
    "Removed stale comments", "Removed debug exit block") —
    same F28 pattern. Verified against current code: `nnte`
    does not appear in the Step 11 block or any `.SDcols`
    list in `cleanchild()`; "Round to 3 decimals" is not
    present; `saveRDS` only appears in the `child_clean.R`
    file-header documentation comment, not in Step 11.
  - No "Configurable parameters in scope for Step 11"
    subsection existed — same gap F23 / F28 added for
    Steps 7 and 9.
- **Fix.** Full rewrite of the Step 11 narrative:
  - `Code location` → "Inline in `cleanchild()` in
    `child_clean.R`; support function `ewma()` defined
    earlier in `child_clean.R`".
  - Added an Overview paragraph that also enumerates the
    11a / 11b / 11c sub-sections.
  - Reorganised content under the 11a / 11b / 11c labels;
    EWMA computation / exclusion criteria / worst-value
    selection / iteration bookkeeping as sub-sections under
    11b.
  - Worst-value selection paragraph updated to "Lowest
    `internal_id` breaks ties."
  - Removed the "Debug columns preserved" subsection.
  - Dropped stale Checklist items 1–3 and renumbered the
    remaining items 1–8.
  - Added "Configurable parameters in scope for Step 11"
    subsection listing `ewma_window` (default 15) with a
    note that the threshold constants are not
    user-configurable.

### Session 6 deferreds (not addressed)

### D26. Step 11 targeted test coverage gaps — DEFERRED
- **Scope.** `test-child-stress.R` line 134 asserts
  `catcount("Exclude-C-Traj-Extreme") == 58` as an
  aggregate regression count on the stress dataset;
  `test-child-regression.R` has no Step 11-specific
  assertions. No targeted unit tests exist for:
  - Count-filter branch (subject-param with exactly 2
    Include rows must not be processed)
  - Extreme-filter branch (subject-param with max
    `|tbc.sd| <= 3.5` must not be processed)
  - Positive- vs. negative-outlier branches at the exclusion
    gate
  - `ctbc.sd == tbc.sd` fast path vs. ctbc recompute path
  - `ctbc.sd` NA-escape on the `c.dewma.all` check
  - First/last Include protection (worst candidate at an
    endpoint must not be flagged)
  - Worst-row tiebreaker (multi-candidate → sort by
    `abs(tbc.sd + dewma.all)` → `internal_id` final
    tiebreaker)
  - Global iteration loop (constructed case where
    iteration 1 exclusion opens a new candidate that fires
    in iteration 2)
  - Targeted per-iteration temp-SDE recalc (subject with
    both an existing temp SDE and a new EWMA1 exclusion)
- **Why deferred.** Coverage gaps; medium effort; no
  current accuracy issue (regression counts unchanged on
  stress data). Same D24 / D25 precedent.
- **Where fixed later.** Dedicated test-coverage session
  for walked steps.

---

## Session 6 — post-fix test results

After reinstalling the package and regenerating
`man/cleangrowth.Rd` via `devtools::document()`:

- test-cleangrowth.R: 65 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 codetools warnings,
  baseline, unchanged)
- test-child-parameters.R: 13 PASS (1 expected deprecation
  warning, baseline, unchanged)

Adult suites and regression harness not re-run (no adult
code touched).

---

## Wrapper-level deferreds surfaced during narrative review (post-Session 6)

### D27. Deprecated-parameter translation shims in `cleangrowth()` — FIXED (2026-04-18)
**Resolution.** Removed `weight_cap` and `include.carryforward` from `cleangrowth()` argument list; deleted both deprecation-translation blocks and the associated warnings. Updated `test-child-parameters.R` Test 3 to use `cf_rescue = "all"` / `cf_rescue = "standard"` directly (baseline deprecation warning dropped from 1 → 0). Removed the `weight_cap` deprecation test from `test-cleangrowth.R`. Regenerated `cleangrowth.Rd`. Phase 1 header in wrapper reference renamed from "Deprecation handling and parameter validation" to "Parameter validation"; `weight_cap` row removed from Configurable Parameters table. Tests post-fix: `test-cleangrowth.R` 63 PASS 0 warnings; `test-child-parameters.R` 13 PASS 0 warnings.

*Original deferred item:*
- **File/lines:** `child_clean.R`, top of `cleangrowth()` (Phase 1 "Deprecation handling and parameter validation" in the wrapper technical reference).
- **Issue.** `cleangrowth()` translates two deprecated parameters to their current equivalents: `weight_cap` → `adult_scale_max_lbs` (deprecation warning + value copy) and `include.carryforward = TRUE` → `cf_rescue = "all"` (deprecation warning, with a stronger warning if the caller also sets `cf_rescue` explicitly). `isTRUE()` is used so `NA` slips past the deprecation branch silently. The shims exist for backward compatibility with pre-v3.0.0 callers.
- **Why defer/address.** v3.0.0 is already a breaking release (output changed from character vector to data.table, multiple parameters removed, exclusion codes renamed). The deprecation shims preserve only a small slice of the pre-v3.0.0 surface while adding entry-block complexity and a test-baseline deprecation warning. Dropping them simplifies Phase 1 and removes a category of parameter-validation branching.
- **Proposed fix.** Remove both translation blocks from `cleangrowth()`; drop `weight_cap` from the argument list; drop any `include.carryforward = TRUE` handling (leaving `cf_rescue` as the only control). Update tests and roxygen to reference current parameter names only; expect the "1 expected deprecation warning" baseline in `test-child-parameters.R` to drop to 0.
- **Where fixed later.** Next wrapper-level code pass / v3.0.0 parameter-surface cleanup.

### D28. Phase 11 HC boundary: reference says `>= 3*365.25 → Exclude-Not-Cleaned`, code says `>= 5*365.25 → Exclude-Missing` — FIXED (2026-04-18, option 2)
**Resolution.** Per Carrie's pick, both the 3y and 5y HC regions now use `Exclude-Not-Cleaned` (consolidated to a single HC exclusion code). Code change: `child_clean.R` line 1132 now reads `data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Exclude-Not-Cleaned']` (was `'Exclude-Missing'`). The pre-recentering `agedays > 3*365.25 → Exclude-Not-Cleaned` assignment at line 1057 is unchanged; the post-recentering check now overrides the generic "NA tbc.sd → Exclude-Missing" rule for the HC-above-5y case, keeping a consistent code across the 3y+ region. Ripple edits: wrapper reference Phase 11 (5*365.25 + rationale paragraph), wrapper reference Key Concepts (HEADCM bullet), wrapper reference Common Pitfalls (HC reference-data bullet), `child-gc-narrative-2026-04-13.md` Key Concepts, `gc-github-latest/CLAUDE.md` Exclusion Codes table and HC summary line. Tests post-fix: `test-cleangrowth.R`, `test-child-regression.R`, `test-child-algorithms.R` all pass (child-algorithms retains 2 baseline codetools warnings, unchanged).

*Original deferred item:*
- **File/lines:** `child_clean.R` line 1130 (`data.all[param == "HEADCM" & agedays >= 5*365.25, exclude := 'Exclude-Missing']`) vs. the wrapper technical reference Phase 11 post-recentering safety check, which now describes `agedays >= 3 * 365.25` → `Exclude-Not-Cleaned`. Related pre-recentering assignment at line 1057 uses `agedays > 3*365.25` → `Exclude-Not-Cleaned` (listed two bullets earlier in the same Phase 11 section).
- **Issue.** The wrapper reference and the code disagree on both the age threshold (3y vs 5y) and the exclusion code (`Exclude-Not-Cleaned` vs `Exclude-Missing`) for the post-recentering HC safety check. Carrie's stated preference: "Missing-Not-Cleaned" — intent is that the post-recentering safety check should also use `Exclude-Not-Cleaned` (not `Exclude-Missing`), and possibly at the 3y boundary rather than 5y.
- **Why defer.** Needs a product decision on (a) the threshold (3y or 5y), (b) the exclusion code (`Not-Cleaned` vs `Missing`), and (c) whether the post-recentering check is redundant with the pre-recentering 3y assignment. Depending on the answers, either the code or the reference needs to change. Text is left as-is in the reference for now; a walkthrough session should reconcile.
- **Where fixed later.** Child algorithm walkthrough session covering Phase 11 / HC handling; coordinate code change in `child_clean.R` line 1130 with reference text.

### D29. Check whether `parallel = TRUE` emits warnings — INVESTIGATED, DOCUMENTED (2026-04-18, no code change)
**Resolution.** Confirmed: `parallel = TRUE` emits exactly two `codetools` warnings per `cleangrowth()` call of the form `<anonymous>: ... may be used in an incorrect context: '.fun(piece, ...)'` — one each from the child and adult `ddply()` dispatches. These are false positives: `codetools` flags `plyr::ddply`'s internal use of `...` as potentially misused, but plyr's usage is correct. Documented in wrapper reference Phase 5 ("Known harmless warnings") and cross-linked to the CRAN-prep item on replacing `plyr`. No suppression applied — suppressing would mask future real warnings and create cleanup debt; the warnings go away naturally when `ddply` is replaced with `data.table`/`foreach`. See CLAUDE.md → CRAN Preparation Checklist for the plyr-replacement item.

*Original deferred item:*
- **Scope.** Carrie noted (during review of wrapper reference Phase 5) that `parallel = TRUE` may emit ~2 harmless warnings; cause unknown. Unverified.
- **Investigation needed.** Run `cleangrowth(..., parallel = TRUE)` on a small sample, capture any warnings with `withCallingHandlers()`, identify the root cause. Possible sources: `clusterExport()` on unexported helpers, `ddply()` deprecation messages, reference-file loading race conditions in workers.
- **Resolution options.** (a) If warnings are cosmetic, document them in Phase 5 of the wrapper reference and suppress with `suppressWarnings()` locally. (b) If they point to a real issue (e.g., missing `clusterExport` of a function), fix the root cause.
- **Where fixed later.** Next parallel-path review, or alongside D27/wrapper-level Phase 1 cleanup.

### D30. Partial-run output row order should use user `id` — FIXED (2026-04-18)
**Resolution.** Changed `child_clean.R` Phase 15 partial-run rbind reorder from `all_results[order(internal_id)]` to `all_results[order(id)]`. User `id` is contract-guaranteed unique and preserved untouched through both the cached and re-processed paths, producing a stable and semantically-meaningful output row order. Full-run path is unaffected (it already restores by `line`). Wrapper reference Phase 15 updated accordingly; CLAUDE.md Open (wrapper) entry removed. Tests post-fix: `test-child-regression.R` 48 PASS (partial-run tests unaffected).

*Original deferred item:*
- **File/lines:** `child_clean.R` Phase 15 rbind + reorder step (the `rbind` with `cached_results[!subjid %in% changed_subjids]` followed by `setorder(..., internal_id)`).
- **Issue.** Already noted in `gc-github-latest/CLAUDE.md → Known Issues → Open (wrapper)`: output is currently ordered by `internal_id`, which is reassigned `1..K` over only the changed-subject input and collides with the cached `1..N`, producing deterministic-but-semantically-meaningless row order. Carrie has confirmed that sorting by user `id` (contract-guaranteed unique and preserved untouched through both paths) is the preferred fix.
- **Proposed fix.** Replace `setorder(..., internal_id)` with `setorder(..., id)` for the final output-assembly step in partial-run mode. Verify that the full-run path is unaffected (it restores by `line` already, which is still correct there).
- **Where fixed later.** Next wrapper-level code pass. Related CLAUDE.md entry should be updated to reflect the chosen fix direction.

### D31. Lower Fenton floor from 250 g → 100 g in GA-correction integer-weight step — FIXED (2026-04-18)
**Resolution.** Changed `child_clean.R` line 833 from `pc[intwt >= 250 & intwt < 500, intwt := 500]` to `pc[intwt >= 100 & intwt < 500, intwt := 500]`. Inline comment updated to explain the rationale. Wrapper reference already described the target state (`[100, 500)`). Tests post-fix: `test-cleangrowth.R` 63 PASS, `test-child-regression.R` 48 PASS, `test-child-algorithms.R` 41 PASS.

*Original deferred item:*
- **File/lines:** `child_clean.R` line 833: `pc[intwt >= 250 & intwt < 500, intwt := 500]` (inside `cleangrowth()`, Phase 10 / Child Step 2b Fenton integer-weight block).
- **Issue.** The current floor only rounds `intwt` values in `[250, 500)` up to the Fenton table's minimum of 500. Values below 250 g currently fall out of the Fenton merge silently. In practice most <250 g weights are BIVs and would be excluded anyway, but Carrie has decided to lower the floor to 100 g so that users who configure lower BIV thresholds can still get Fenton-corrected z-scores for those weights. The wrapper technical reference Phase 10 / Phase 3 already describes the target state (`[100, 500)`).
- **Proposed fix.** Change the condition to `pc[intwt >= 100 & intwt < 500, intwt := 500]`. Add targeted test cases covering a weight in `[100, 250)` to confirm Fenton merge succeeds and produces a sensible `fengadays`/`unmod_zscore`.
- **Where fixed later.** Next wrapper-level or Child Step 2b code pass; small, localized change.

### D32. Rename "narrative" → "reference" everywhere — DEFERRED

- **Scope.** Rename the three algorithm/wrapper documents and all references to them so that the word "reference" is used consistently in filenames, body text, and cross-references; "narrative" should be retired as a term for these documents. The child document has already been renamed (`child-gc-narrative-2026-04-13.md` → `child-algorithm-reference.md`) and its in-body "narrative" references updated on 2026-04-18; the other two files and all cross-references to them still use the old pattern.
- **Files to rename (proposed names TBD).**
  - `wrapper-narrative-2026-04-17.md` → e.g., `wrapper-reference.md`
  - `adult-algorithm-narrative.md` → e.g., `adult-algorithm-reference.md`
- **Files to update after rename.**
  - Body of each renamed file: replace in-body uses of "narrative" referring to the document itself with "reference"; update cross-references to the other two renamed files.
  - [`child-algorithm-reference.md`](child-algorithm-reference.md): 10 remaining `wrapper-narrative-2026-04-17.md` cross-references to update.
  - [`gc-github-latest/CLAUDE.md`](CLAUDE.md): multiple references to `wrapper-narrative-2026-04-17.md`, `child-gc-narrative-2026-04-13.md`, `adult-algorithm-narrative.md`, and their descriptions (e.g., "Algorithm narrative documents" section, "Narrative status" section, per-section cross-references).
  - `__Pipeline/CLAUDE.md`: references to `wrapper-narrative-2026-04-17.md`, `child-gc-narrative-2026-04-13.md`, `adult-algorithm-narrative.md` in the "Algorithm narratives" section and elsewhere.
  - `gc-github-latest/algorithm-walkthrough-procedure.md`: any "narrative" mentions that refer to these documents.
  - Existing walkthrough-todo files: leave as-is (they are historical logs; renaming would obscure the record of what was called what at the time).
- **Why defer.** The three files should be renamed together so cross-references can be updated in a single pass. Splitting the rename risks broken links and inconsistent terminology during the interim.
- **Proposed approach.** Use `git mv` for each rename (preserves history). After all renames, sweep each updated file for (a) cross-references to the renamed files and (b) in-body uses of "narrative" that refer to the document itself (leaving uses that refer to "narrative" as a general term — e.g., "narrative description" — alone).
- **Where fixed later.** Next walkthrough session focused on documentation reconciliation, or whenever the next file-level rename is being done (e.g., when the wrapper or adult file is next substantively updated).

---

## GA correction (Phase 10 / Child Step 2b) — walkthrough checklist findings (2026-04-17)

Moved from the wrapper technical reference on 2026-04-18. These are walkthrough findings applied to the GA-correction section of `cleangrowth()`.

1. **Sort determinism.** Age-dependent sort key (`id_sort`) is assigned explicitly before `potcorr_wt` is marked.
2. **Fenton reference files.** Loaded only under `has_potcorr`; `fent_foraga.csv.gz` and `fenton2025_ms_lookup_smoothed.csv` both in `inst/extdata/`.
3. **Fast-path correctness.** When no subjects flag, the fast path writes `sd.corr := sd.orig` and `uncorr := 0L` — every row still has a non-missing `sd.corr` for downstream use.
4. **NA safety.** The `potcorr_wt` predicate uses `!is.na(sd.orig) & sd.orig < -2`, so `Exclude-Missing` rows (where `sd.orig` is `NA`) cannot flag.
5. **Fenton merge failure.** Correctly resets `potcorr_wt` to FALSE and recomputes subject-level `potcorr` so the subject is excluded from downstream Fenton logic without silently using partial results.
6. **Correction revert rule.** The abssumdiff comparison is computed only over weight rows with `ageyears_2b < 2` and `seq_win ≤ 4`. The check is first-weight-only (`is_first`), which matches its purpose (the correction is driven by the first qualifying weight).
7. **HEIGHTCM supine/standing offsets** (+0.8 WHO, +0.7 CDC) are applied only when `agedays > 730 & cagedays ≤ 730` — the band where chronological age is post-2y but corrected age is still under 2y.

---

## Session 7 — Step 11 (EWMA1 Extreme) deeper pass (2026-04-18)

Second pass over Child Step 11 after Session 6 (2026-04-17). Session 6 covered the main cleanup (F34–F38: debug parameter removal, 11a/11b/11c sub-sectioning, current-state narrative rewrite, `lowest id` → `lowest internal_id`). This pass surfaced four narrow items, all fixed in place, plus one deferred.

### Pre-session baseline tests

After the 2026-04-18 wrapper-reference / D27–D31 work, the installed package was retested:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings from parallel test)
- test-child-parameters.R: 13 PASS

### Fix-now items

### F39. Stale `# Return z-scores and EWMA1 iteration 1 values for comparison` comment removed — FIXED
- **File/lines:** `child_clean.R` pre-edit line 4784 — header comment for the output-assembly column-build block in `cleanchild()`.
- **Issue.** The `debug` parameter and `ewma1_it1.*` capture block were removed in Session 6 F34; the current code returns only essentials (`id`, `line`, `exclude`, `param`, `cf_rescued`), optional z-score columns, `final_tbc`, and (when `cf_detail = TRUE`) `cf_status` / `cf_deltaZ`. The "EWMA1 iteration 1 values for comparison" clause described a capture mechanism that no longer exists.
- **Fix.** Replaced the comment with current-state prose: `# Assemble return columns: essentials + z-scores + final_tbc + optional cf_detail columns.`

### F40. Dead columns in Step 11 `.SDcols` removed — FIXED
- **File/lines:** `child_clean.R` pre-edit line 3235 — `.SDcols` for the per-group closure at the end of Step 11b.
- **Issue.** The `.SDcols` list included `index`, `id`, and `sex`, none of which are referenced inside the per-group closure body. The closure reads `exclude` (for `.child_valid(df)`), `agedays`/`tbc.sd`/`ctbc.sd` (for EWMA and the pot_excl criteria), `param` (for `.child_exc(param, "Traj-Extreme")`), and `internal_id` (worst-row tiebreaker). `id`, `index`, and `sex` are copied into `.SD` per group with no consumer.
- **Fix.** Reduced the list to `c('internal_id', 'param', 'agedays', 'tbc.sd', 'ctbc.sd', 'exclude')`. No expected behavioral change (the removed columns were never read).
- **Verification.** All 5 child test suites pass after reinstall: 63 / 48 / 28 / 41 / 13, identical to baseline.

### F41. Narrative "the closure resets those subjects' temp SDEs" phrasing clarified — FIXED
- **File/lines:** `child-algorithm-reference.md` Step 11b iteration-bookkeeping item 2 (pre-edit line 967).
- **Issue.** The reset + `identify_temp_sde()` rerun happens in the outer iteration-loop block (`child_clean.R` ~lines 3243–3254), not inside the per-group closure. The word "closure" misleadingly suggested the per-group `copy(.SD)` closure.
- **Fix.** Replaced "the closure resets those subjects' temp SDEs to Include and reruns `identify_temp_sde()` on that subset" with "the post-pass recalc block resets those subjects' temp SDEs to Include and reruns `identify_temp_sde()` on that subset."

### F42. Narrative Code-location cell now names `identify_temp_sde()` — FIXED
- **File/lines:** `child-algorithm-reference.md` Step 11 summary table (pre-edit line 881).
- **Issue.** Code-location cell listed only `ewma()` as a support function used by Step 11, but Step 11b (per-iteration targeted recalc at `child_clean.R` ~lines 3250–3254) and Step 11c (end-of-step refresh at ~lines 3273–3274) also call `identify_temp_sde()`.
- **Fix.** Cell now reads: "Inline in `cleanchild()` in `child_clean.R`; support functions `ewma()` and `identify_temp_sde()` also defined in `child_clean.R`."

### Session 7 deferreds

### D33. Legacy Stata-style `# 6. / a. / b. / c. / i. / ii. / iii.` comment block in `ewma()` — FIXED (Session 7b, 2026-04-18)
- **File/lines:** `child_clean.R` ~lines 1662–1682 — numbered block inside the `ewma()` function body describing the original EWMA formula in pseudo-math notation (`EWMASDi`, `SDj…SDn`, `ΔAgej`, `EWMAZ`, `EWMAZbef`, `EWMAZaft`).
- **Issue.** Per walkthrough convention (F26/F30/F35), Stata-style enumeration and algorithm-design prose in the implementation file should be rewritten as current-state R-idiom comments. The block predates the R implementation and does not use current variable names (`tbc.sd`, `ctbc.sd`, `ewma.all`/`before`/`after`, `delta`).
- **Why deferred.** `ewma()` is a shared support function used by Child Steps 11, 13, 15/16, and 17, and also by the adult algorithm indirectly. The block is not inside Step 11's main logic; it's inside the `ewma()` function definition. Rewriting it belongs with a broader `ewma()` / `ewma_cache_init()` / `ewma_cache_update()` support-function pass rather than inside a single-step walkthrough. The technical content is still accurate (the formula described matches what the code computes).
- **Related.** `child_clean.R` ~lines 3571–3576 has a similar "Restructured to use global iterations for efficiency / Key changes:" block at the top of Child Steps 15/16 — Step 15/16 scope, not Step 11. To be handled during the Steps 15/16 walkthrough.
- **Where fixed later.** Dedicated `ewma()` / support-function cleanup pass (prompt proposed separately).

### Session 7 — post-fix test results

After reinstalling the package:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All counts identical to baseline. No regressions.

---

## Session 7b — EWMA support-function pass (2026-04-18)

Dedicated walkthrough of the shared EWMA machinery in `R/child_clean.R`:
`ewma()` (~lines 1657–1744), `ewma_cache_init()` (~1753–1804),
`ewma_cache_update()` (~1809–1911), and `as_matrix_delta()` (~1607). Closes
D33. Comment-only changes; no behavioral edits.

### Pre-session baseline tests

- test-cleangrowth.R: 63 PASS
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 41 PASS (2 pre-existing codetools warnings)
- test-child-parameters.R: 13 PASS
- test-adult-clean.R: 198 PASS
- adult regression test_harness.R: 1508/1508 at all 4 permissiveness levels

### Scope correction

- The original prompt described `ewma()` / `ewma_cache_init()` / `ewma_cache_update()` as "called from Child Steps 11, 13, 15/16, and 17 (and any adult counterparts)." On inventory, **only the child algorithm calls these functions.** The adult algorithm has its own parallel implementations `adult_ewma_cache_init()` and `adult_ewma_cache_update()` in `R/adult_support.R` and uses them in `cleanadult()` (and adult support functions). `as_matrix_delta()` is only called by the child `ewma()` and `ewma_cache_init()`.
- Net effect: this pass is child-scoped by code. Adult tests were still run as a baseline confirmation but are not affected by the edits.

### Fix-now items

### F43. Stata-style `# 6. / a. / b. / c. / i. / ii. / iii.` block in `ewma()` rewritten — FIXED (closes D33)
- **File/lines:** `R/child_clean.R` ~lines 1662–1682 (pre-edit).
- **Issue.** Numbered block predating the R implementation, written in pseudo-math notation (`EWMASDi`, `SDj…SDn`, `ΔAgej`, `EWMAZ`, `EWMAZbef`, `EWMAZaft`) and Stata-style enumeration. Did not use current R variable names (`tbc.sd`, `ctbc.sd`, `ewma.all` / `before` / `after`, `delta`).
- **Fix.** Replaced with a current-state R-idiom prose block at the top of the `ewma()` body that describes:
  - the weight formula `w_ij = (5 + |agedays_i - agedays_j|) ^ ewma.exp_i` with `w_ii = 0`;
  - why "+5" and a negative exponent are used;
  - what each of the three return variants (`ewma.all` / `ewma.before` / `ewma.after`) excludes and the first/last/`n <= 2` conventions;
  - the per-observation vector form of `ewma.exp` and how `sweep()` applies it row-wise;
  - the `window` parameter and its disable form (`Inf`);
  - the role of `cache_env` for paired tbc/ctbc calls.

### F44. `@return` roxygen corrected from "Data frame with 3 variables" to named list — FIXED
- **File/lines:** `R/child_clean.R` `@return` block in `ewma()` roxygen (lines ~1631–1638 pre-edit).
- **Issue.** Roxygen described a data frame, but the function explicitly returns a named list (line ~1739 inline comment: "data.table := compatible, avoids data.frame overhead"). The list form is the version current callers depend on (`df[..., (ewma.fields) := ewma(...)]` patterns at child_clean.R lines 3186, 3191, 3482, 4500). The `ewma.adjacent = FALSE` case was also undocumented.
- **Fix.** Roxygen now describes a named list of numeric vectors, with each variant's contents and the `n <= 2` / first / last conventions, plus the `ewma.adjacent = FALSE` one-element-list case.

### F45. `@param ewma.exp` clarified to scalar-or-vector — FIXED
- **File/lines:** `R/child_clean.R` `@param ewma.exp` in `ewma()` roxygen.
- **Issue.** Said only "Exponent to use for weighting." All four internal callers in `cleanchild()` pass a per-observation vector (`exp_vals`), not a scalar; the roxygen example happens to pass a scalar (works because `sweep()` recycles).
- **Fix.** Roxygen now states scalar or per-observation vector usage and notes that the algorithm's main internal callers pass a per-observation vector that varies by widest neighbor age gap.

### F46. "self-weight is zero (replaces ifelse(delta == 0, …) approach)" parenthetical dropped — FIXED
- **File/lines:** `R/child_clean.R` line ~1704 (pre-edit).
- **Issue.** Trailing parenthetical referenced a prior `ifelse` implementation that no longer exists — stale-content style.
- **Fix.** Comment now reads `# self-weight is zero so obs i is excluded from its own EWMA`.

### F47. "O(n) arithmetic instead of matrix copy + multiply" comments rephrased — FIXED
- **File/lines:** `R/child_clean.R` lines ~1722–1724 and ~1729–1731 (pre-edit, two paired comment blocks for `pred_weights` / `succ_weights`).
- **Issue.** Phrased as a comparison to a hypothetical alternative implementation rather than describing what the current code does.
- **Fix.** Consolidated into one current-state Before/After comment block that says the predecessor/successor contribution is subtracted directly from the cached weighted sum and row sum, and that this avoids rebuilding a separate weight matrix per adjacent-neighbor exclusion. No "instead of …" framing.

### F48. "Bug fix: was only checking 1 position on each side; extended to 2 positions on each side." removed in `ewma_cache_update()` — FIXED
- **File/lines:** `R/child_clean.R` lines ~1848–1849 (pre-edit).
- **Issue.** Stale "Bug fix: was X" comment style. The preceding paragraph already explains *why* up to 2 positions on each side need to be checked after removing `pos_j` (the shift in max-gap for `pos_j-1` and `pos_j-2`, and symmetrically on the other side).
- **Fix.** Two-line "Bug fix" tail removed; rationale paragraph kept.

### F49. Step 15/16 header "Restructured to use global iterations for efficiency / Key changes:" block rewritten — FIXED (bycatch)
- **File/lines:** `R/child_clean.R` ~lines 3570–3576 (pre-edit).
- **Issue.** Same stale-changelog style as D33: bullets numbered "1./2./3." describing what the restructure changed, in the implementation file. Step 15/16 scope, but trivially close in style and 5 lines long.
- **Fix.** Replaced with a single current-state paragraph describing the global-iteration design (per-batch p_plus / p_minus computed once; one shared while loop; only subject-param groups with new exclusions re-processed each iteration; converged groups are not re-walked).
- **Note.** This is bycatch from the EWMA support-function pass and was approved in conversation. The Steps 15/16 walkthrough should still happen; this just resolves one cosmetic comment block in advance so the convention stays clean.

### `ewma_cache_init()` — no fixes needed
- Comments are current-state. `skip_ctbc` shortcut, weight matrix build, windowing, and EWMA computation blocks all match what the code does. Leave as-is.

### `as_matrix_delta()` — no fixes needed
- Single-purpose internal helper (lines ~1605–1611): `Convert vector of ages (days) into pairwise absolute-difference matrix. Used by ewma() for age-gap weighting.` Comment matches behavior; implementation is correct. The construction `abs(matrix(rep(agedays, n), n, byrow = TRUE) - agedays)` could equivalently be written with `outer()` but neither correctness nor performance is affected. Leave as-is. (Note: `as_matrix_delta` was previously exported and is now internal — already done in earlier 2026-04-16 cleanup; man page deletion already documented in CLAUDE.md.)

### Session 7b deferreds

None. The two related items called out in D33 are both resolved this session:
- D33 itself (F43) — closed in `ewma()`.
- Step 15/16 header block (F49) — closed as bycatch.

Items intentionally not in scope and not opened as new deferreds:
- The roxygen `@examples` block for `ewma()` still names its result `e_df`, which reads as "ewma data frame". Now that the return contract is documented as a list, the variable name is mildly misleading but does not break the example. Cosmetic; not worth a separate deferred entry.

### Session 7b — post-fix test results

After reinstalling the package and regenerating man pages with `roxygen2::roxygenise()`:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions. Adult tests not re-run by request — edits do not touch `adult_clean.R` or `adult_support.R`.

---

## Session 8 — `identify_temp_sde()` support-function pass (2026-04-18)

Dedicated walkthrough of the shared temp-SDE resolution support function in `R/child_clean.R`: `identify_temp_sde()` (~lines 2045–2200 pre-edit; ~lines 2055–2200 post-edit). Comment-only plus small dead-code cleanup; no behavioral changes.

### Pre-session baseline tests

- test-cleangrowth.R: 63 PASS
- test-child-regression.R: 48 PASS
- test-child-edge-cases.R: 28 PASS
- test-child-algorithms.R: 41 PASS (2 pre-existing codetools warnings)
- test-child-parameters.R: 13 PASS
- Adult tests not run (see Scope correction).

### Scope correction

- The original prompt suggested `identify_temp_sde()` might be called from both child and adult algorithms. Inventory confirmed: called only from `cleanchild()` (7 call sites). The adult algorithm has its own separate `temp_sde()` function in `R/adult_support.R:449` with a different signature and semantics; it is not a shared support function.
- Net effect: this pass is child-scoped by code. Adult tests were not re-run — `adult_clean.R` and `adult_support.R` were not touched.
- Call-site inventory:

| Line (pre-edit) | Step | `exclude_from_dop_ids` |
|---|---|---|
| 2571 | Step 5 (initial pass) | NULL |
| 2682 | Step 6 (recalc after CF identification) | NULL |
| 2996 | Step 7d (recalc after BIV) | NULL |
| 3114 | Step 9d (recalc after Evil Twins) | NULL |
| 3260 | Step 11 mid-loop (subset recalc for affected subjects) | NULL |
| 3283 | Step 11c end-of-step (full-dataset recalc) | NULL |
| 3306–3307 | Step 13 (final SDE) | `temp_sde_ids_step13` |

Only the Step 13 caller passes a non-NULL `exclude_from_dop_ids`.

### Fix-now items

### F50. Title + `Step 5 and Step 13` roxygen simplification rewritten — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2015–2034 (pre-edit).
- **Issue.** Roxygen title said the function was for "Step 5 and Step 13." In practice the function is called from 7 sites across Steps 5, 6, 7, 9, 11 (mid- and end-of-step), and 13. Only Step 13 passes non-NULL `exclude_from_dop_ids`. The description paragraph on `median.dopz` repeated the same "Step 5 / Step 13" simplification.
- **Fix.** Title rewritten as `Identify temporary SDE (same-day extraneous) "losers" on each (subjid, param, agedays) group.` The `median.dopz` description paragraph no longer carries the step-specific framing. Added an explicit "Called from Steps 5, 6, 7, 9, 11 (mid-loop and end-of-step), and 13 of \\code{cleanchild()}" paragraph that names all callers and notes the NULL-vs-non-NULL contract.

### F51. Inline `Step 5 vs Step 13` comments updated — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2046–2048 (parameter-contract comment inside function body) and ~lines 2113–2117 (block header above the DOP-median branch) pre-edit, plus the two `# Step 5:` / `# Step 13:` comments inside the `if` / `else` branches themselves.
- **Issue.** Same oversimplification — these comments described the parameter contract as "Step 5 vs Step 13" rather than the actual contract "NULL vs non-NULL", and the inline `# Step 5:` / `# Step 13:` labels at the top of each branch framed the paths as if only those two callers existed.
- **Fix.** Parameter-contract comment rewritten to list all non-NULL callers explicitly: "Only the Step 13 caller passes ids; all earlier callers (Step 5 initial pass and the recalc callers in Steps 6, 7, 9, and 11) leave it NULL." The block header above the DOP-median branch now describes "When `exclude_from_dop_ids` is NULL (all callers except Step 13) …" and "When non-NULL (Step 13) …". The two branch-label comments were changed to "Default path (all callers except Step 13)" and "Step 13 path".

### F52. Orphan invariant comment between `@noRd` and function def moved into roxygen — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2041–2043 (pre-edit).
- **Issue.** A 3-line `#` comment block sat between the closing roxygen `@noRd` tag and the `identify_temp_sde <- function(...)` definition, describing an important invariant ("Identical same-day values are removed earlier (Early Step 13 SDE-Identicals), so by the time identify_temp_sde() runs, all same-day values in a given (subjid, param, agedays) group are dissimilar."). This is an important precondition but lived in an awkward transitional spot.
- **Fix.** Moved the invariant into the roxygen description as the opening paragraph of `@details` (immediately after the title), where it reads as the natural precondition for the selection rule described right below it. The orphan `#` block was deleted.

### F53. `# ----- STEP 1/2/3/4/5 -----` section headers rewritten as current-state prose — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2099, 2104, 2113, 2153, 2173 (pre-edit).
- **Issue.** Five Stata-era numbered section headers (`# ----- STEP 1: Calculate SP medians ... -----` style) inside the function body. Same overall style concern as Session 7b F43. Additionally, STEP 2's description ("Distribute SP medians across all rows for each subject") was inaccurate — the distribution happens in STEP 1 via data.table `:=`; STEP 2 actually creates a compact `(subjid, param) -> median.spz` lookup table used later by STEP 3 for cross-parameter anchoring.
- **Fix.** Replaced each of the five numbered headers with a short current-state prose comment:
  - STEP 1 → "Assign each valid row its SP median: median of tbc.sd within (subjid, param), computed over all valid rows (not just SDE days)."
  - STEP 2 → "Compact (subjid, param) -> median.spz lookup table. Used below to anchor each row's DOP median to the median of its other parameter." (accuracy fix)
  - STEP 3 → folded with F51; now describes the NULL vs non-NULL branch structure as the decision point.
  - STEP 4 → "For rows on an SDE day, compute the absolute distance of tbc.sd from the two median anchors. These distances drive the selection below."
  - STEP 5 → "Within each SDE group, pick the winner by sorting on absdmedian.spz (asc), then absdmedian.dopz (asc), then age-dependent internal_id. The first row after sorting is the keeper; all others are marked extraneous."

### F54. Stale "Use by (not keyby)" comment replaced (trim-only; code unchanged) — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2070–2072 (pre-edit).
- **Issue.** Three-line comment block read `# make a small copy of df with fields we need / # Include internal_id for age-dependent tiebreaker / # Use by (not keyby) to preserve group order, then sort explicitly`. The last line in particular reads as a design-note-to-self contrasting with an earlier alternative rather than describing what the code does. The `by = .(subjid, param, agedays)` on the next line is effectively a no-op grouping (same row count, grouped arrangement) and the real sort is done by the `order(...)` on the line after that — so the "Use by" hint was misleading.
- **Fix.** Replaced with a single current-state block: "Narrow df to the columns we need, grouped by (subjid, param, agedays) and then sorted with internal_id as the tiebreaker so same-day rows have a deterministic order for the winner selection below." Code (the `by = ...` column-subset on line 2073 and explicit `order(...)` on 2074) left unchanged per user instruction.

### F55. Added formal `@param` and `@return` tags — FIXED
- **File/lines:** `R/child_clean.R` roxygen for `identify_temp_sde()`.
- **Issue.** The description mentioned the return value and parameter semantics in prose but there were no formal `@param` or `@return` tags. Same precedent as Session 7b F44 for `ewma()`: adding the tags is useful even though the function is `@noRd`, because the source-level roxygen is still the contract for readers.
- **Fix.** Added three tags:
  - `@param df` — "data.table column-subset containing id, internal_id, subjid, param, agedays, tbc.sd, and exclude. Not mutated."
  - `@param exclude_from_dop_ids` — "optional vector of \\code{id} values to remove from the DOP median calculation only (they still contribute to the SP median). NULL by default; non-NULL only at Step 13."
  - `@return` — "A logical vector aligned to the caller's input row order. TRUE marks rows the caller should flag \\code{'Exclude-C-Temp-Same-Day'}; FALSE means leave the exclude column unchanged."

### Deferred items also resolved (per user request)

### D-a. Dead defensive code removed — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2058–2062 (pre-edit).
- **Issue.** `if (is.null(df$subjid)) df[, subjid := NA]` / `if (is.null(df$param)) df[, param := NA]` defensive blocks. All 7 callers pass both columns, so these branches never fire.
- **Fix.** Removed both `if` blocks (5 lines) and the leading "may be missing depending on where this is called from" comment.

### D-b + D-c. Return block simplified — FIXED
- **File/lines:** `R/child_clean.R` ~lines 2189–2199 (pre-edit) plus the earlier unused `original_rows <- df$orig_row` at ~line 2068 (pre-edit).
- **Issue (D-b).** `extraneous_result <- df$extraneous & valid.rows` had a redundant `& valid.rows` — `df$extraneous` is only ever set TRUE by the guarded `extraneous := ...` assignment above (which restricts to `valid.rows & extraneous.this.day`), so masking again on `valid.rows` was unnecessary.
- **Issue (D-c).** The return used a string-keyed named-vector lookup (`setNames(..., df$orig_row)` + `as.character(original_rows)`) to permute the result back into the caller's row order, with `as.character()` on integer values. Correct but inefficient and indirect.
- **Fix.** Rewrote the final block as:
  ```r
  result <- logical(nrow(df))
  result[df$orig_row] <- df$extraneous
  return(result)
  ```
  Also removed the unused `original_rows <- df$orig_row` assignment higher up (no longer referenced). Comment rewritten to describe what the mapping does and to note that `df$extraneous` can be used without re-masking because the guarded assignment above already restricts to valid rows on SDE days.

### D-d. Call-site style unified to data.table idiom (limited scope) — FIXED
- **File/lines:** `R/child_clean.R` lines 2571 (Step 5 initial) and 2682 (Step 6 recalc), pre-edit.
- **Issue.** Two call sites used base-R vector indexing: `data.df$exclude[identify_temp_sde(...)] <- 'Exclude-C-Temp-Same-Day'`. Four other call sites (lines 2996, 3114, 3283, 3306–3307) used the data.table `:=` idiom: `data.df[identify_temp_sde(...), exclude := 'Exclude-C-Temp-Same-Day']`. Both produce the same result, but the mixed style was a readability issue.
- **Fix.** Converted the two base-R-style sites (2571 and 2682) to the `:=` idiom to match the others. The Step 11 mid-loop call at line 3260 uses a different separate-variable pattern (`sde_result <- identify_temp_sde(...)` followed by an index-based apply) and was intentionally left unchanged — that pattern serves a different purpose (subsetted input, index-join apply to `data.df`) and is not a base-R vs data.table style mismatch.

### Session 8 deferreds

None. All items flagged during the pass (F50–F55, D-a, D-b, D-c, D-d) were approved and resolved.

Items intentionally not in scope and not opened as new deferreds:
- `dop_map` is defined inside `identify_temp_sde()` (lines ~2098 post-edit) and also inside `get_dop()` (lines ~2002–2012 post-edit). Two different shapes — `get_dop()` returns a scalar for scalar input; the in-function `dop_map` is a named character vector used for `dop_map[p]` lookups. Not worth collapsing.
- The `by = .(subjid, param, agedays)` column-subset on line 2078 (post-edit) is redundant in effect (the real sort is done by the explicit `order(...)` on the next line) but the pattern is harmless and preserving it avoids a code change during a comment-only pass. Left as-is per F54 trim-only scope.

### Session 8 — post-fix test results

After reinstalling the package and regenerating man pages with `roxygen2::roxygenise()`:

- test-cleangrowth.R: 63 PASS, 0 warnings
- test-child-regression.R: 48 PASS, 0 warnings
- test-child-edge-cases.R: 28 PASS, 0 warnings
- test-child-algorithms.R: 41 PASS (2 baseline codetools warnings, unchanged)
- test-child-parameters.R: 13 PASS

All child counts identical to baseline. No regressions. Adult tests not re-run — adult algorithm has its own separate `temp_sde()` in `adult_support.R`; `identify_temp_sde()` is child-only.
