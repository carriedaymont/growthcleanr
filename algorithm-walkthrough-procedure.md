# Algorithm Walk-Through Procedure

Procedure for conducting a systematic review of a growthcleanr algorithm
(adult or child) before major milestones (clinician validation, CRAN submission,
significant algorithm changes). Covers code correctness, inline comments,
parameter specs, and narrative documentation.

---

## Purpose

A walk-through serves four goals simultaneously:

1. **Find code issues** — bugs, sort-order non-determinism, threshold inconsistencies, dead code
2. **Verify inline comments** — comments match what the code actually does
3. **Reconcile all documentation** — code, comments, parameter specs, and narrative tell the same story
4. **Remove stale content** — prior version references, edit dates, author notes

---

## Pre-Walk-Through Setup

### 1. Confirm baseline tests pass

**Adult:**
```bash
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-adult-clean.R")'

Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest
```

**Child:**
```bash
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-regression.R")'

NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-stress.R")'

NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-edge-cases.R")'

NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-algorithms.R")'

NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-parameters.R")'
```

All tests must pass before starting. If any fail, fix first.
(Child `test-child-algorithms.R` has 6 known failures for unimplemented
`cf_rescue_threshold` — these are pre-existing and acceptable.)

### 2. Create a dated todo file

Name: `walkthrough-todo-YYYY-MM-DD.md`

Template:
```markdown
# [Adult/Child] GC Walk-Through — Deferred To-Dos

Items identified during the YYYY-MM-DD walk-through that are deferred for later resolution.

---

## Pre-walk-through setup

[any issues found before starting]

---

## Steps walk-through (populated as walk-through proceeds)
```

### 3. Identify documents to review

**Adult:**
- `R/adult_clean.R` — main algorithm
- `R/adult_support.R` — support functions
- `adult-algorithm-narrative.md` — step-by-step documentation
- `gc-adult-permissiveness-spec-YYYY-MM-DD.md` — permissiveness framework
- `wtallow-formulas.md` — canonical specification for wtallow base formulas, UW adjustment, ET caps, and ceilings. Must be reconciled with code and narrative for every step that uses wtallow, ET caps, or related limits (Steps 9Wa, 9Wb, 11Wa, 11Wa2, 11Wb).

**Child:**
- `R/child_clean.R` — main algorithm and all support functions (single file, unlike adult's two-file split)
- `child-gc-narrative-2026-03-18.md` — step-by-step documentation (check for current version)
- No permissiveness spec (not applicable for 3.0.0)
- No wtallow equivalent (child uses z-score thresholds, not raw measurement limits)

---

## Walk-Through Procedure

Review 1–3 steps per session. For each step:

### Step A: Read the narrative section

Find and read the step's section in the narrative document. Note what it claims
the code does.

### Step B: Read the code

Find the corresponding code in the main algorithm file and any support functions.

- **Adult:** The main loop in `adult_clean.R` dispatches to support functions in
  `adult_support.R`.
- **Child:** Everything is in `child_clean.R` — both the main algorithm
  (`cleanchild()`) and all support functions (`valid()`,
  `temporary_extraneous_infants()`, `calc_otl_evil_twins()`,
  `calc_and_recenter_z_scores()`, `ewma()`, etc.).

### Step C: Apply the checklist

For each step, systematically check all items that apply. Items are marked
**[adult]**, **[child]**, or **[both]**.

1. **[both] Sort order determinism** — Every `order()` or sort call that could
   have ties must include an `as.numeric(internal_id)` tiebreaker. `internal_id`
   is character (created in `cleangrowth()`), so numeric conversion is required
   for correct sort order. Look for any call like `order(ageyears)` or
   `order(agedays)` without an internal_id tiebreaker.

   **[child only] Birth tiebreaking** — At agedays == 0, the algorithm keeps
   the *lowest* internal_id (earliest measurement, before postnatal fluid
   shifts). At all other ages, it keeps the *highest* internal_id (later, more
   careful measurement). Verify any age-dependent tiebreaking logic gets this
   right.

2. **[adult] Rounding tolerance** — All threshold comparisons should include
   the standard tolerance (`+ 0.12` for adult ht/wt). Check for bare
   comparisons like `> wta` that should be `> wta + 0.12`.

   **[child] Z-score correctness** — Verify that EWMA calculations and z-score
   comparisons use the correct z-score (`tbc.sd`, `ctbc.sd`, `sd.orig`, etc.)
   for their purpose. Confirm that the correct values are included in EWMA
   calculations (e.g., that values meant to be excluded aren't inadvertently
   contributing to the EWMA, and that recentered vs. non-recentered z-scores
   are used in the right contexts).

3. **[both] Dead code** — Conditions that can never be TRUE, variables computed
   but never used, codes included in lists that can never match any real value.

4. **[both] Exclusion code names** — Every exclusion code string in the code
   must exactly match the canonical list in `CLAUDE.md`. Check for typos, old
   names, or codes that were renamed.

5. **[both] Column names** — Inline comments and narrative must use the correct
   column names. Common sources of stale documentation:
   - Adult: `meas_m` not `meas_orig`; `cap_params` not `caps`
   - Child: `tbc.sd`, `ctbc.sd`, `sd.orig`, `data.df`

6. **[adult] Permissiveness table accuracy** — Narrative parameter tables must
   show correct values for all permissiveness levels, not just one. Cross-check
   against the preset code in `adult_support.R` (the `adult_presets()` function).

   **[child] Configurable parameter defaults** — Verify that configurable
   parameter defaults match the code and narrative. The child algorithm has
   ~13 configurable parameters (e.g., `sd.extreme`, `error.load.threshold`,
   `ewma_window`, `lt3.exclude.mode`, `recover.unit.error`). No permissiveness
   levels for 3.0.0.

7. **[both] Step linkage** — "Prior Step" and "Next Step" in each narrative
   summary table must match actual algorithm flow. Renumbered or reordered
   steps leave stale step references. Child step numbers are non-consecutive
   (Early 13, 5, 6, 7, 9, 11, 13, 15, 16, 17, 19, 21, 22), so stale
   references are especially likely.

8. **[both] Grepl vs exact matching** — Any `grepl(pattern, codes)` used to
   classify exclusion codes risks matching unintended codes. Prefer exact
   `%in%` matching when the set of codes is known. Especially dangerous:
   substring patterns that match multiple codes with different semantics.

9. **[both] Inline comment accuracy** — Comments that describe what the code
   does should be read skeptically. The code changed; the comments may not have.

10. **[both] Stale content** — References to prior versions, edit dates
    ("as of v2.1"), author attributions, "TODO" comments, and "check whether"
    notes that have already been checked. For child, also check for references
    to the legacy pediatric algorithm or Stata that no longer apply.

11. **[adult] wtallow/limit reconciliation** — For any step that uses `wtallow`,
    `compute_wtallow()`, `compute_et_limit()`, ET caps, or UW-based scaling,
    verify that the code implementation matches the canonical specification in
    `wtallow-formulas.md`. Check: base formula values, UW adjustment logic
    (highUW, lowUW, ceiling), cap values at 6m/12m, ET cap derivation, and
    the allofus15 cap-limiting rule. Also verify the narrative description of
    these limits matches both code and spec. Relevant steps: 9Wa (Evil Twins),
    9Wb (Extreme EWMA), 11Wa (2D Ordered), 11Wa2 (2D Non-Ordered), 11Wb
    (Moderate EWMA).

    **[child] EWMA/BIV/velocity threshold reconciliation** — For steps that
    use EWMA thresholds, BIV limits, or velocity tables, verify the code
    matches documentation:
    - EWMA extreme/moderate thresholds and `ewma.exp` weighting
    - BIV absolute and standardized limits (including age-dependent preterm
      weight fix: <0.2 kg first year, <1 kg after)
    - `height.tolerance.cm` and the Tanner/WHO velocity tables in Step 17

12. **[child] `valid()` flag correctness** — Each step uses `valid()` with
    specific flag combinations (`include.temporary.extraneous`,
    `include.extraneous`, `include.carryforward`). Verify that the correct
    flags are used for each step — using the wrong flags means the step
    operates on the wrong set of rows.

13. **[child] DOP (Designated Other Parameter) logic** — Weight's DOP is
    height, height's DOP is weight, HC's DOP is height. Many steps use DOP
    for cross-parameter plausibility checks. Verify the DOP assignments are
    correct.

### Step D: Classify findings

**Fix now** if:
- Behavioral bug (wrong result for any input)
- Incorrect exclusion code or threshold value
- Documentation that contradicts the code (misleads future work)
- Stale references that need removal

**Defer** if:
- Dead code with no behavioral impact
- Sort determinism issue that can't affect results in practice
- Documentation completeness (e.g., missing parameter table rows for edge cases)
- Refactoring/cleanup with no correctness impact

### Step E: Apply fixes

Fix-now items immediately. Run tests after each fix to confirm no regression.

### Step F: Document deferred items

Add each deferred item to the todo file with:
- Step number and file
- Description of the issue
- Why it's deferred (priority + behavioral impact)
- Exact fix needed (enough detail to implement without re-reading the code)

---

## Post-Walk-Through

### 1. Run full test suite

**Adult:**
```bash
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-adult-clean.R")'
Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest
```

**Child:**
```bash
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-regression.R")'
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-stress.R")'
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-edge-cases.R")'
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-algorithms.R")'
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-parameters.R")'
```

All tests must pass.

### 2. Write a summary

Record in session notes or CHANGELOG:
- Date of walk-through
- Steps covered
- Bugs found and fixed
- Deferred items (with link to todo file)
- Final test counts

---

## Common Patterns Found in Adult Walk-Through (2026-04-03)

These are the most common issues found; check for these specifically:

| Pattern | What to look for |
|---------|-----------------|
| Character id sort | `order(..., internal_id)` without `as.numeric()` — lexicographic sort is wrong |
| Missing tolerance | `> threshold` without `+ 0.12` on ht/wt comparisons |
| Grepl matching wrong codes | `grepl("Identical", ...)` matched `Scale-Max-Identical` (non-SDE) |
| Stale column name | `meas_orig` should be `meas_m`; `caps` should be `cap_params` |
| Stale step link | "Next Step: 12W" when Step 12 doesn't exist |
| Single-level param table | Tables showing only one permissiveness level |
| Dead grepl | `grepl("Same-day", ...)` on adult codes (no adult code contains "Same-day") |
| BMI BIV not implemented | Parameters resolved from presets but never used in the BIV check |
| Fixed iteration loop | `for (round in 1:3)` should be `while (change)` |

---

## Child-Specific Patterns to Watch For

| Pattern | What to look for |
|---------|-----------------|
| Wrong `valid()` flags | Missing `include.temporary.extraneous = TRUE` when temp SDEs should participate |
| Birth vs non-birth tiebreaking | `agedays == 0` should keep lowest internal_id; all other ages keep highest |
| Wrong grouping | `by = subjid` when it should be `by = .(subjid, param)`, or vice versa |
| DOP assignment error | HC's DOP is height (not weight); verify cross-parameter logic |
| Factor level silent NA | Exclusion code not in `exclude.levels` produces NA instead of error |
| Z-score blend boundaries | WHO-only < 2y, blend 2-5y, CDC-only > 5y; HC is always WHO |
| Wrong z-score variable | Using `sd.orig` where `tbc.sd` or `ctbc.sd` is needed, or vice versa |
| Stale Stata references | Comments referencing Stata line numbers or variables that no longer apply |
