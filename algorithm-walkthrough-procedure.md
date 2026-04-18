# Algorithm Walk-Through Procedure

Procedure for conducting a systematic review of growthcleanr code before
major milestones (clinician validation, CRAN submission, significant algorithm
changes). Covers the full code path: `cleangrowth()` wrapper (preprocessing,
dispatch, output assembly), the child algorithm, and the adult algorithm.
Checks code correctness, inline comments, parameter specs, and narrative
documentation.

---

## Purpose

A walk-through serves four goals simultaneously:

1. **Find code issues** — bugs, sort-order non-determinism, threshold inconsistencies, dead code
2. **Verify inline comments** — comments match what the code actually does
3. **Reconcile all documentation** — code, comments, parameter specs, and narrative tell the same story
4. **Remove stale content** — prior version references, edit dates, author notes

---

## Cross-cutting code-review checklist

Apply this general checklist to each step or phase during a walkthrough. The detailed implementation-specific checks in **Phase 1 → Step C** below are additional step-level items (rounding tolerance, wtallow reconciliation, permissiveness tables, etc.).

1. **Comments and narrative match code.** No stale comments describing removed logic; narrative accurately reflects what the code does.
2. **Exact boundaries.** `<` vs `<=`, `>` vs `>=`, and exact numeric tolerances (e.g., 0.4 vs 0.41) verified against spec.
3. **Uninitialized or unnecessary variables.** No variables used before assignment; no leftover variables from prior refactors.
4. **Efficiency.** Caching opportunities, unnecessary merges, redundant computations, operations that could be vectorized.
5. **Edge cases.** 0, 1, or 2 valid rows; empty groups; all-NA values — does code handle gracefully or crash?
6. **`.child_valid()` call correctness** (child algorithm only). Which `include.*` flags are passed? Consistent with step's purpose?
7. **`by`-group correctness.** `by = subjid` vs `by = .(subjid, param)` etc. — wrong grouping silently produces wrong results.
8. **Sort order assumptions.** Does code assume a particular sort? Is that guaranteed by a prior `setkey` / `setorder`?
9. **data.table reference semantics.** `:=` on `data.df` (intended) vs. accidentally modifying a copy or joined table.
10. **Factor level issues.** Exclusion codes assigned must exist in `exclude.levels`; assigning an unlisted string to a factor silently produces NA.
11. **Parameter scope.** Does the step handle all 3 params (HT, WT, HC) or only some? Is HC exclusion intentional and documented?
12. **Interaction with later steps.** Variables stored for downstream use — correctly scoped in `data.df`?

---

## Pre-Walk-Through Setup

### 1. Confirm baseline tests pass

**Important:** Reinstall the package first so tests run against current code:
```bash
Rscript -e 'devtools::install_local(".", force=TRUE, upgrade="never")'
```

**All test files:**
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-cleangrowth.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-edge-cases.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-algorithms.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-parameters.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-adult-clean.R")'

Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest
```

All tests must pass before starting. If any fail, fix first.

### 2. Create a dated todo file

Name: `walkthrough-todo-YYYY-MM-DD.md`

Template:
```markdown
# GC Full Walk-Through — Deferred To-Dos

Items identified during the YYYY-MM-DD walk-through that are deferred for later resolution.

---

## Pre-walk-through setup

[any issues found before starting]

---

## Wrapper walkthrough (populated as walk-through proceeds)

[cleangrowth() preprocessing, dispatch, output assembly]

---

## Child algorithm steps (populated as walk-through proceeds)

---

## Adult algorithm steps (populated as walk-through proceeds)
```

### 3. Identify documents to review

**Wrapper (`cleangrowth()` + support functions):**
- `R/child_clean.R` lines ~333–1855 — `cleangrowth()` entry point, preprocessing, dispatch, output assembly
- `R/child_clean.R` support functions — `read_anthro()`, `gc_preload_refs()`, `ewma()`, `calc_and_recenter_z_scores()`, `.cf_rescue_lookup()`
- `R/utils.R` — shared utilities
- **No narrative document yet.** Findings go in the walkthrough todo file. A wrapper narrative should be created later (TODO).

**Adult:**
- `R/adult_clean.R` — main algorithm
- `R/adult_support.R` — support functions
- `adult-algorithm-narrative.md` — step-by-step documentation
- `gc-adult-permissiveness-spec-YYYY-MM-DD.md` — permissiveness framework
- `wtallow-formulas.md` — canonical specification for wtallow base formulas, UW adjustment, ET caps, and ceilings. Must be reconciled with code and narrative for every step that uses wtallow, ET caps, or related limits (Steps 9Wa, 9Wb, 11Wa, 11Wa2, 11Wb).

**Child:**
- `R/child_clean.R` — main algorithm and all support functions (single file, unlike adult's two-file split)
- `child-gc-narrative-2026-03-18.md` — step-by-step documentation (check for current version)
- `Child-growthcleanr-permissiveness-specs.md` — **ignore during walkthroughs.** Draft scaffolding for a future child permissiveness framework, deferred until after v3.0.0 is released. It is not part of the current implementation contract and is not expected to reflect current code or parameter names. Do not reconcile code, comments, or other narratives against it.
- No wtallow equivalent (child uses z-score thresholds, not raw measurement limits)

---

## Walk-Through Procedure

### Suggested session plan (full walkthrough)

| Session | Scope |
|---------|-------|
| 1 | Wrapper: parameter handling, partial-run mode, data construction, input validation, batching setup, exclude.levels, velocity table loading |
| 2 | Wrapper: imperial→metric, LENGTHCM, z-score calculation, Step 2b (GA correction), parallel setup |
| 3 | Wrapper: algorithm dispatch, output assembly, partial-run merge, bin/tri_exclude, support functions (`read_anthro`, `gc_preload_refs`) |
| 4+ | Algorithm steps: child (Early 13, 5, 6, 7, 9, 11, 13, 15, 16, 17, 19, 21, 22), then adult (1–14). Review 1–3 steps per session. |

---

### Phase 0: Wrapper walkthrough

The `cleangrowth()` wrapper (~lines 333–1855 of `child_clean.R`) handles
everything before and after the algorithm steps. Review it in chunks before
moving to algorithm steps.

**For each chunk:**
1. Read the code
2. Apply relevant checklist items (see below — wrapper-specific items are marked **[wrapper]**)
3. Classify findings as fix-now or defer (same criteria as algorithm steps)
4. Log deferred items in the walkthrough todo file

**Wrapper chunks to review:**

| Chunk | Approx lines | What it does |
|-------|-------------|--------------|
| Parameter deprecation & validation | 333–365 | `weight_cap`→`adult_scale_max_lbs`, `include.carryforward`→`cf_rescue` |
| Partial-run mode | 428–507 | `cached_results` / `changed_subjids` auto-detect and explicit filtering |
| Data construction + input validation | 509–558 | `data.all.ages` build, `internal_id` creation, sex recoding, param validation, cutpoint clamping |
| Outer batching + batch-invariant setup | 559–763 | `exclude.levels` (child + adult), velocity table loading (Tanner, WHO), batching by 2000 subjects |
| Inner batch loop: ped/adult split | 764–782 | Age cutpoint split into `data.all` vs `data.adult`, `data.batch` reference copy |
| Imperial→metric, LENGTHCM | 789–885 | Unit conversion, param relabeling |
| Z-score calculation (CSD) | 888–936 | WHO/CDC closures, blending (WHO<2y, blend 2–5y, CDC>5y), HC always WHO |
| Step 2b: GA correction (potcorr) | 940–1469 | Fenton merge, split-normal z, corrected WHO/CDC, recentering, `sd.median` file, EWMA fields |
| Parallel setup | 796–825 | `var_for_par`, cluster creation/export |
| Algorithm dispatch | 1470–1700 | Params passed to `cleanchild()`/`cleanadult()`, adult-specific setup |
| Output assembly | 1719–1854 | Child+adult merge, `cf_rescued`/`mean_ht`/`bin_result`, `data.batch` join, partial-run merge, `bin_exclude`/`tri_exclude`, safety check |

---

### Phase 1: Algorithm steps walkthrough

Review 1–3 steps per session. For each step:

### Step A: Read the narrative section

Find and read the step's section in the narrative document. Note what it claims
the code does. (For wrapper chunks reviewed in Phase 0, there is no narrative
yet — review code and comments only.)

### Step B: Read the code

Find the corresponding code in the main algorithm file and any support functions.

- **Adult:** The main loop in `adult_clean.R` dispatches to support functions in
  `adult_support.R`.
- **Child:** Everything is in `child_clean.R` — both the main algorithm
  (`cleanchild()`) and all support functions (`.child_valid()`,
  `identify_temp_sde()`, `calc_otl_evil_twins()`,
  `calc_and_recenter_z_scores()`, `ewma()`, etc.).

### Step C: Apply the checklist

For each step, systematically check all items that apply. Items are marked
**[adult]**, **[child]**, or **[both]**.

1. **[both] Sort order determinism** — Every `order()` or sort call that could
   have ties must include an `internal_id` tiebreaker. `internal_id` is an
   integer assigned in id-sorted order by `cleangrowth()`. No numeric conversion
   needed. Look for any call like `order(ageyears)` or `order(agedays)` without
   an internal_id tiebreaker.

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
   configurable parameters (e.g., the eight `biv.z.*` Step 7 cutoffs,
   `error.load.threshold`, `ewma_window`, `cf_rescue`, `batch_size`). No
   permissiveness levels for 3.0.0.

   **[both] Parameters table row audit per session** — Each session must
   audit the rows of the narrative's Configurable Parameters table that
   correspond to the step(s) in scope. Cross-check against the current
   `cleangrowth()` / `cleanadult()` signature and defaults; remove rows for
   parameters that no longer exist (e.g., `recover.unit.error`,
   `use_legacy_algorithm`, `ewma.exp` were removed during 2026-04 cleanup).
   Table-wide cleanup should accumulate across sessions rather than be
   deferred to a separate pass.

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

14. **[wrapper] Input vector alignment** — All input vectors (`subjid`,
    `param`, `agedays`, `sex`, `measurement`, `id`) must stay aligned
    through filtering (e.g., partial-run `keep` subsetting). Check that no
    operation reorders one vector without the others.

15. **[wrapper] Batch boundary safety** — Subjects split across batches
    would break by-subject operations. Verify batching assigns whole
    subjects to the same batch, never splitting a subject's rows.

16. **[wrapper] Column survival across merge** — The `data.batch` →
    `full_out` merge by `line` must not drop or duplicate columns. Check
    `all.x = TRUE` behavior when adult or child results are empty
    (`data.table()` with 0 rows).

17. **[wrapper] Adult id overwrite** — Line ~1662: `data.adult[, id := line]`
    overwrites the user's original `id` for adult data with the `line`
    index. Investigate whether this loses the user's id or whether
    `internal_id` handles it. Compare with how child data preserves `id`.

18. **[wrapper] Z-score blending boundary consistency** — The WHO/CDC
    blending window in `cleangrowth()` preprocessing (lines ~908–930) must
    match `calc_and_recenter_z_scores()` used later in Steps 11/15. Check
    for boundary disagreements (e.g., `>=` vs `>` at 2 or 5 years).

19. **[wrapper] Partial-run correctness** — Auto-detect comparison
    (lines ~445–485) must apply the same transformations as the main path
    (0→NaN, sex recoding). Verify the hash-based comparison is
    order-independent within each subject (i.e., sorted before hashing).

20. **[wrapper] CRAN output compliance** — `cat()`/`print()` should be
    `message()` for CRAN. Flag any `cat()` or `print()` calls that are
    not gated behind `if (!quietly)`.

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

**Important:** Reinstall the package before testing so tests run against
current code, not a stale installed version:
```bash
Rscript -e 'devtools::install_local(".", force=TRUE, upgrade="never")'
```

**Integration + child + adult:**
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-cleangrowth.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-edge-cases.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-algorithms.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-parameters.R")'

NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-adult-clean.R")'

Rscript tests/test_harness.R loosest
Rscript tests/test_harness.R looser
Rscript tests/test_harness.R tighter
Rscript tests/test_harness.R tightest
```

All tests must pass.

### 2. Write a summary

Record in session notes or CHANGELOG:
- Date of walk-through
- Sections covered (wrapper chunks, algorithm steps, or both)
- Bugs found and fixed
- Deferred items (with link to todo file)
- Final test counts

### 3. Future: Create wrapper narrative

A narrative document for the `cleangrowth()` wrapper (preprocessing,
dispatch, output assembly) does not yet exist. The child and adult
narratives only cover what happens inside `cleanchild()` and
`cleanadult()`. After completing the first full walkthrough, create
a wrapper narrative document following the same format as the algorithm
narratives. This will serve as a reference for future walkthroughs and
for anyone reading the code.

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

## Wrapper-Specific Patterns to Watch For

| Pattern | What to look for |
|---------|-----------------|
| `cat()`/`print()` not `message()` | Lines ~539, 560, 767, 837 etc. use `cat()`/`print()` — CRAN requires `message()` |
| `dplyr` in data.table code | Line ~563 uses `%>% select() %>% distinct()` — inconsistent with data.table convention and adds a dependency |
| Stale commented-out code | Lines ~1702–1717, ~1786 — large blocks of dead commented-out code |
| `data.adult[, id := line]` | Line ~1662 overwrites user's `id` — may lose original id for adult rows |
| Hardcoded batch size | Batch size 2000 on line ~568 — should be a parameter |
| `print("Batching Begin AAA")` | Line ~560 — debug print that should be removed or gated behind `quietly` |

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
