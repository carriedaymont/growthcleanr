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
- Standardized BIV limits (`sd.extreme`, `z.extreme` parameters)
- Preterm weight fix (<0.2 kg first year, <1 kg after)

[to be populated as the walk proceeds]
