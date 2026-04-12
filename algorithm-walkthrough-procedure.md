# Algorithm Walk-Through Procedure

Procedure for conducting a systematic review of a growthcleanr algorithm
(adult or child) before major milestones (clinician validation, CRAN submission,
significant algorithm changes). Covers code correctness, inline comments,
permissiveness/parameter specs, and narrative documentation.

---

## Purpose

A walk-through serves four goals simultaneously:

1. **Find code issues** — bugs, sort-order non-determinism, threshold inconsistencies, dead code
2. **Verify inline comments** — comments match what the code actually does
3. **Reconcile all documentation** — code, comments, permissiveness spec, and narrative tell the same story
4. **Remove stale content** — prior version references, edit dates, author notes

---

## Pre-Walk-Through Setup

### 1. Confirm baseline tests pass

```bash
# Adult unit tests
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-adult-clean.R")'

# Adult regression tests (all 4 permissiveness levels)
Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest

# Child tests
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-regression.R")'
```

All tests must pass before starting. If any fail, fix first.

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
- `R/child_clean.R` — main algorithm + support functions
- `child-gc-narrative-YYYY-MM-DD.md` — step-by-step documentation
- `Child-growthcleanr-permissiveness-specs.md` — if permissiveness framework exists

---

## Walk-Through Procedure

Review 1–3 steps per session. For each step:

### Step A: Read the narrative section

Find and read the step's section in the narrative document. Note what it claims
the code does.

### Step B: Read the code

Find the corresponding code in the main algorithm file and any support functions.
For adult, the main loop in `adult_clean.R` dispatches to support functions in
`adult_support.R`.

### Step C: Apply the 10-point checklist

For each step, systematically check:

1. **Sort order determinism** — Every `order()` or sort call that could have ties
   must include an `internal_id` tiebreaker. For adult, `internal_id` is stored
   as character inside `cleanadult()`, so use `as.numeric(internal_id)` not
   `internal_id`. Look for any call like `order(ageyears)` or
   `order(age_days)` without an internal_id tiebreaker.

2. **Rounding tolerance** — All threshold comparisons should include the standard
   tolerance (`+ 0.12` for adult ht/wt; none for child z-scores). Check for
   bare comparisons like `> wta` that should be `> wta + 0.12`.

3. **Dead code** — Conditions that can never be TRUE, variables computed but never
   used, codes included in lists that can never match any real value.

4. **Exclusion code names** — Every exclusion code string in the code must exactly
   match the canonical list in `CLAUDE.md`. Check for typos, old names, or codes
   that were renamed.

5. **Column names** — Inline comments and narrative must use the correct column
   names (`meas_m` not `meas_orig`; `cap_params` not `caps`; etc.). These are
   common sources of stale documentation.

6. **Permissiveness table accuracy** — Narrative parameter tables must show correct
   values for all permissiveness levels, not just one. Cross-check against the
   preset code in `adult_support.R` (the `adult_presets()` function).

7. **Step linkage** — "Prior Step" and "Next Step" in each narrative summary table
   must match actual algorithm flow. Renumbered or reordered steps leave stale
   step references.

8. **Grepl vs exact matching** — Any `grepl(pattern, codes)` used to classify
   exclusion codes risks matching unintended codes. Prefer exact `%in%` matching
   when the set of codes is known. Especially dangerous: substring patterns that
   match multiple codes with different semantics.

9. **Inline comment accuracy** — Comments that describe what the code does should
   be read skeptically. The code changed; the comments may not have.

10. **Stale content** — References to prior versions, edit dates ("as of v2.1"),
    author attributions, "TODO" comments, and "check whether" notes that have
    already been checked.

11. **wtallow/limit reconciliation** — For any step that uses `wtallow`,
    `compute_wtallow()`, `compute_et_limit()`, ET caps, or UW-based scaling,
    verify that the code implementation matches the canonical specification in
    `wtallow-formulas.md`. Check: base formula values, UW adjustment logic
    (highUW, lowUW, ceiling), cap values at 6m/12m, ET cap derivation, and
    the allofus15 cap-limiting rule. Also verify the narrative description of
    these limits matches both code and spec. Relevant steps: 9Wa (Evil Twins),
    9Wb (Extreme EWMA), 11Wa (2D Ordered), 11Wa2 (2D Non-Ordered), 11Wb
    (Moderate EWMA).

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

```bash
# Adult
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-adult-clean.R")'
Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest

# Child
NOT_CRAN=true Rscript -e \
  'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-child-regression.R")'
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

## Common Patterns Found in 2026-04-03 Adult Walk-Through

These are the most common issues found; check for these specifically:

| Pattern | What to look for |
|---------|-----------------|
| Character id sort | `order(..., id)` where id is character — use `as.numeric(id)` |
| Missing tolerance | `> threshold` without `+ 0.12` on ht/wt comparisons |
| Grepl matching wrong codes | `grepl("Identical", ...)` matched `Scale-Max-Identical` (non-SDE) |
| Stale column name | `meas_orig` should be `meas_m`; `caps` should be `cap_params` |
| Stale step link | "Next Step: 12W" when Step 12 doesn't exist |
| Single-level param table | Tables showing only one permissiveness level |
| Dead grepl | `grepl("Same-day", ...)` on adult codes (no adult code contains "Same-day") |
| BMI BIV not implemented | Parameters resolved from presets but never used in the BIV check |
| Fixed iteration loop | `for (round in 1:3)` should be `while (change)` |

---

## Notes for Child Walk-Through

The child algorithm has different conventions from adult:

- **No rounding tolerance** — threshold comparisons are exact (no `+ 0.12`)
- **`valid()` gatekeeper** — understand `include.temporary.extraneous`, `include.extraneous`, `include.carryforward` flags before reading any step
- **Z-scores not raw measurements** — comparisons use `tbc.sd`/`ctbc.sd`, not `meas_m`
- **Single `data.df`** — rows never removed; exclusions set `exclude` column
- **Age-dependent internal_id tiebreaking** — at birth (agedays == 0): keep lowest internal_id; all other ages: keep highest internal_id
- **No permissiveness levels yet** — parameter table review is simpler
- **DOP (Designated Other Parameter)** — weight's DOP is height, and vice versa; HC's DOP is height. Many steps use DOP for cross-parameter plausibility

Narrative document: `child-gc-narrative-2026-03-18.md` (may be dated; check for current version).
