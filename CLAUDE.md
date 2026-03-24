# CLAUDE.md — gc-github-latest (growthcleanr)

**Last updated:** 2026-03-24 (edits: pitfalls section, valid() warning, v3.0.0 callout, test as-of date, Known Issues rename)

## Overview

growthcleanr (gc) is an R package that identifies and flags
implausible anthropometric measurements (height, weight, head
circumference) in electronic health record data. It does not
remove any values — each row gets an `exclude` code:
`"Include"` or an exclusion code naming the step and reason.

Repository: https://github.com/carriedaymont/growthcleanr
Branch: `efficiency-updates` (v3.0.0, not yet merged to main)
CRAN: v2.2.0 (outdated — do not use for development)

This CLAUDE.md is for **Claude Code agents working on gc
code**. For pipeline/Qual-AD context, see
`__Pipeline/CLAUDE.md`.

---

## Algorithms

| Algorithm | Status | Default path | Notes |
|-----------|--------|--------------|-------|
| Child | Active, primary | `child_clean.R` | Default pediatric path in v3.0.0 |
| Adult | Active | `adult_clean.R` + `02a_support.R` | Separate logic from child |
| Legacy pediatric | Deprecated | `growth.R` | `use_legacy_algorithm = TRUE` |

The child algorithm replaces the legacy pediatric algorithm.
The legacy pediatric algorithm is maintained for backward
compatibility (`use_legacy_algorithm = TRUE`).
`prelim_infants` is deprecated — it triggers a warning and
maps to `use_legacy_algorithm` (verified; tested in
test-child-regression.R and test-cleangrowth.R). The
user-facing `prelim_infants` parameter has no other effect;
internal `read_anthro(prelim_infants = TRUE)` calls are
hardcoded and independent of the user-facing parameter.

**This CLAUDE.md currently covers the child algorithm only.**
Adult algorithm documentation will be added later.

---

## Code Structure

### Key files

| File | What it contains |
|------|------------------|
| `R/growth.R` | `cleangrowth()` entry point, exports; also contains legacy pediatric algorithm |
| `R/child_clean.R` | `cleanbatch_child()` main algorithm + support functions (`valid()`, `temporary_extraneous_infants()`, `calc_otl_evil_twins()`, `calc_and_recenter_z_scores()`, `ewma()`, `ewma_cache_init()`/`ewma_cache_update()`, `get_dop()`, `read_anthro()`) |
| `R/utils.R` | Shared utilities (does NOT contain `valid()` — see note below) |
| `R/adult_clean.R` | Adult algorithm main |
| `R/02a_support.R` | Adult support functions |
| `inst/extdata/` | Reference tables (growth charts, recentering file, velocity tables) |
| `tests/testthat/` | Test suite (see Testing section) |

### Preprocessing → main algorithm split

`cleangrowth()` (in `growth.R` / `child_clean.R`) handles:
- Input validation, data.table construction
- Imperial → metric conversion (HEIGHTIN, WEIGHTLBS)
- LENGTHCM → HEIGHTCM reclassification (param only,
  no measurement adjustment — known limitation)
- Age cutpoint split (pediatric vs. adult)
- Z-score calculation (CSD method, WHO/CDC blending)
- GA correction (Step 2b) for potcorr subjects
- Recentering → tbc.sd and ctbc.sd
- Missing value identification
- Batching and dispatch to `cleanbatch_child()`

`cleanbatch_child()` processes one batch through all
cleaning steps (see Complete List of Steps below).

### Batching

Two layered batching systems:

1. **Outer wrapper** (memory management): Divides subjects
   into batches of 2,000. Batch size is hard-coded — should
   be parameterized (future work).

2. **Inner batching** (parallelism): Subdivides for parallel
   processing. With `parallel = FALSE` (current default),
   this is a no-op.

---

## Input/Output Format

### Required input (individual vectors, not a dataframe)

| Parameter | Type | Description |
|-----------|------|-------------|
| `subjid` | any | Subject identifier |
| `param` | character | `"WEIGHTKG"`, `"WEIGHTLBS"`, `"HEIGHTCM"`, `"HEIGHTIN"`, `"LENGTHCM"`, or `"HEADCM"` |
| `agedays` | numeric | Age in days at measurement |
| `sex` | any | `0`/`"m"`/`"M"` = male; `1`/`"f"`/`"F"` = female |
| `measurement` | numeric | Recorded value in units specified by `param` |
| `id` | numeric/char (optional) | Unique row identifier; if NULL, generated as `1:N` |

### Input handling

- `measurement == 0` → replaced with `NaN` (treated as missing)
- Imperial converted to metric before processing
- LENGTHCM relabeled to HEIGHTCM (no measurement adjustment)
- Data split by `adult_cutpoint` (default 20 years)

### Sorting (critical for deterministic results)

```
setkey(data.df, subjid, param, agedays, id)
```

The `id` field breaks ties for same-day duplicates. Without
it in the sort key, SDE resolution order is undefined.

### Output

`cleangrowth()` returns a data.table with:

| Column | Description |
|--------|-------------|
| `id` | Row identifier |
| `exclude` | Exclusion code or `"Include"` |
| `param` | Growth parameter |
| `cf_rescued` | CF rescue reason (empty if not rescued) |
| `sd.orig_who` | WHO CSD z-score |
| `sd.orig_cdc` | CDC CSD z-score |
| `sd.orig` | Blended CSD z-score |
| `tbc.sd` | Recentered blended z-score |
| `ctbc.sd` | Recentered corrected z-score |
| `final_tbc` | ctbc.sd for potcorr, tbc.sd for others |

**Breaking change from v2.2.0:** v3.0.0 returns a data.table
with multiple columns (not a character vector of exclusion
codes like CRAN v2.2.0). Any code consuming `cleangrowth()`
output must expect a data.table, not a vector.

---

## Complete List of Exclusion Codes (child algorithm)

| Code | Step | Param | Description |
|------|------|-------|-------------|
| `Include` | — | All | Value passes all checks |
| `Missing` | Init | All | NA, NaN, agedays < 0; also HC ≥ 5y |
| `Not cleaned` | Init | HEADCM | HC with agedays > 3 × 365.25 |
| `Unit-Error-High` | 3 | HT, WT | Unit error detected (too high); corrected |
| `Unit-Error-Low` | 3 | HT, WT | Unit error detected (too low); corrected |
| `Swapped-Measurements` | 4 | HT, WT | HT/WT swapped; corrected |
| `Exclude-Carried-Forward` | 6 | All | Identical to prior-day value (not rescued) |
| `Exclude-Absolute-BIV` | 7 | All | Outside absolute biological limits |
| `Exclude-Standardized-BIV` | 7 | All | Z-score beyond age-dependent cutoffs |
| `Exclude-Evil-Twins` | 9 | All | Adjacent extreme value pair/group |
| `Exclude-EWMA1-Extreme` | 11 | All | Extreme EWMA outlier |
| `Exclude-SDE-Identical` | 13 | All | Same-day duplicate, identical value |
| `Exclude-SDE-All-Exclude` | 13 | All | All same-day values excluded |
| `Exclude-SDE-All-Extreme` | 13 | All | All same-day values extreme |
| `Exclude-SDE-EWMA` | 13 | All | SDE resolved by EWMA comparison |
| `Exclude-SDE-EWMA-All-Extreme` | 13 | All | All SDE values extreme by EWMA |
| `Exclude-SDE-One-Day` | 13 | All | All measurements on single day |
| `Exclude-Temporary-Extraneous-Same-Day` | 5 | All | Temp SDE (may persist) |
| `Exclude-EWMA2-middle` | 15 | All | Moderate EWMA: middle value |
| `Exclude-EWMA2-birth-WT` | 16 | WT | Birth weight EWMA outlier |
| `Exclude-EWMA2-birth-WT-ext` | 16 | WT | Birth weight EWMA (extended) |
| `Exclude-EWMA2-first` | 15 | All | Moderate EWMA: first value |
| `Exclude-EWMA2-first-ext` | 15 | All | Moderate EWMA: first (extended) |
| `Exclude-EWMA2-last` | 15 | All | Moderate EWMA: last value |
| `Exclude-EWMA2-last-high` | 15 | All | Moderate EWMA: last (high) |
| `Exclude-EWMA2-last-ext` | 15 | All | Moderate EWMA: last (extended) |
| `Exclude-EWMA2-last-ext-high` | 15 | All | Moderate EWMA: last (ext, high) |
| `Exclude-EWMA2-birth-HT-HC` | 16 | HT, HC | Birth HT/HC EWMA outlier |
| `Exclude-EWMA2-birth-HT-HC-ext` | 16 | HT, HC | Birth HT/HC EWMA (extended) |
| `Exclude-Min-diff` | 17 | HT, HC | Height/HC decreased > allowed |
| `Exclude-Max-diff` | 17 | HT, HC | Height/HC increased > allowed |
| `Exclude-2-meas->1-year` | 19 | All | 2 meas >1 yr apart, one excluded |
| `Exclude-2-meas-<1-year` | 19 | All | 2 meas <1 yr apart, one excluded |
| `Exclude-1-meas` | 19 | All | Single measurement excluded |
| `Exclude-Error-load` | 21 | All | Error ratio exceeds threshold |

The `exclude.levels` list includes many legacy codes from
the original pediatric algorithm for backward compatibility.
All exclusion codes will be renamed at the end of development.

### CF rescue reason codes (cf_rescued column)

| Code | Meaning |
|------|---------|
| `""` | Not CF, or CF not rescued |
| `1 CF deltaZ <0.05` | Single CF, z-score diff < 0.05 |
| `1 CF deltaZ <0.1 wholehalfimp` | Single CF, diff < 0.1, imperial |
| `Teen 2 plus CF deltaZ <0.05` | ≥2 CFs, teen, diff < 0.05 |
| `Teen 2 plus CF deltaZ <0.1 wholehalfimp` | ≥2 CFs, teen, diff < 0.1, imperial |

CF rescue thresholds may be tuned in future updates.

---

## Complete List of Steps (child algorithm)

Step numbers are not consecutive — they align with the Stata
implementation. Steps 3, 4, 8, 10, 12, 14, 18, 20 either
do not exist or are handled within other steps.

| Step | Name | Brief description |
|------|------|-------------------|
| (preprocessing) | Z-score calculation | WHO/CDC CSD z-scores, blending, rounding |
| (preprocessing) | Step 2b: GA correction | Fenton-based correction for potcorr subjects |
| (preprocessing) | Recentering | Subtract population medians → tbc.sd/ctbc.sd |
| Early 13 | SDE-Identicals | Remove same-day identical values before CF detection |
| 5 | Temporary SDE | Temporarily flag same-day duplicates |
| 6 | Carried Forwards | Identify and optionally rescue carried-forward values |
| 7 | BIV | Exclude biologically implausible values |
| 9 | Evil Twins | Exclude adjacent extreme values |
| 11 | EWMA1 (Extreme) | Exclude extreme EWMA outliers |
| 13 | Final SDE | Resolve remaining same-day duplicates |
| 15 | EWMA2 (Moderate) | Exclude moderate EWMA outliers |
| 16 | Birth HT/HC | EWMA2 variant for birth values |
| 17 | Height/HC Velocity | Exclude values exceeding velocity limits |
| 19 | Pairs and Singles | Evaluate subjects with 1–2 remaining measurements |
| 21 | Error Load | Exclude all if error ratio too high |
| 22 | Output | Assemble return columns |

---

## Configurable Parameters

| Parameter | Default | Used in | Description |
|-----------|---------|---------|-------------|
| `recover.unit.error` | FALSE | Step 3 | Attempt to identify and correct unit errors |
| `sd.extreme` | 25 | Step 7 | SD-score cutoff for extreme exclusion |
| `z.extreme` | 25 | Step 7 | Z-score cutoff for extreme exclusion |
| `height.tolerance.cm` | 2.5 | Step 17 | Maximum height decrease tolerated |
| `error.load.mincount` | 2 | Step 21 | Min exclusions before evaluating error load |
| `error.load.threshold` | 0.5 | Step 21 | Error ratio above this excludes all |
| `sd.recenter` | NA | Recentering | Child always uses fixed file |
| `include.carryforward` | FALSE | Step 6 | If TRUE, skip CF detection entirely |
| `ewma.exp` | -1.5 | EWMA steps | Exponent for EWMA weighting |
| `lt3.exclude.mode` | "default" | Step 19 | Mode for pairs/singles exclusion |
| `ewma_window` | 15 | EWMA steps | Max observations on each side for EWMA |
| `adult_cutpoint` | 20 | Preprocessing | Age (years) dividing pediatric/adult |
| `use_legacy_algorithm` | FALSE | Dispatch | If TRUE, use legacy pediatric algorithm |
| `quietly` | TRUE | All | Suppress progress messages |

---

## Key Concepts for Code Work

### The `valid()` function

**Two separate definitions exist:** `child_clean.R` (line
5044) and `pediatric_support_legacy.R` (line 10). They serve
different algorithm paths. **Do not create additional copies.**
A previous bug was caused by duplicate `valid()` definitions
where R's alphabetical file load order caused the wrong copy
to take precedence (see Fixed issues below).

Critical gatekeeper — returns a logical vector of which rows
are eligible for each step. The `include.*` flags control
whether temporarily excluded categories participate:

- `valid(df)` — only non-excluded, non-Missing, non-"Not cleaned"
- `valid(df, include.temporary.extraneous = TRUE)` — also temp SDEs
- `valid(df, include.extraneous = TRUE)` — also permanent SDEs
- `valid(df, include.carryforward = TRUE)` — also CFs

Works by checking the text of the `exclude` column, not
factor ordering.

### CSD z-scores (not LMS)

All z-scores use the Conditional Standard Deviation method:
```
If measurement < M:  sd = (measurement - M) / csd_neg
If measurement >= M: sd = (measurement - M) / csd_pos
```
This is intentionally more sensitive to extreme high values
than standard LMS z-scores. Intermediate z-score rounding
has been removed (previously used Stata-style rounding for
cross-platform validation).

### WHO/CDC age blending

| Age range | HT/WT formula | HC |
|-----------|---------------|-----|
| < 2 years | WHO only | WHO only |
| 2–5 years | `(cdc × (age-2) + who × (5-age)) / 3` | WHO only |
| > 5 years | CDC only | WHO only (through 5y) |

HC is WHO-only at all ages. HC > 3y is "Not cleaned"; HC ≥ 5y
is "Missing" (no reference data).

### Age-dependent id tiebreaking

At birth (agedays == 0): keep lowest id (earliest measurement,
before postnatal fluid shifts). All other ages: keep highest
id (later measurement, likely more careful). This differs from
the adult algorithm, which always keeps the highest id.

### Designated Other Parameter (DOP)

Weight's DOP is height; height's DOP is weight; HC's DOP is
height. Used in Steps 5, 6, 11, 13, 15, 19 for cross-parameter
plausibility checks.

### data.table reference semantics

The child algorithm keeps all rows in a single `data.df`
data.table — rows are never physically removed. Exclusions
work by setting the `exclude` column and using `valid()` to
filter. This differs from the adult algorithm, which copies
rows into separate shrinking dataframes.

`:=` modifies `data.df` in place. Be careful not to
accidentally modify copies or joined tables.

---

## Testing

### testthat suite

4 test files in `tests/testthat/`:

| File | Tests | Assertions | Coverage |
|------|-------|------------|----------|
| `test-cleangrowth.R` | 4 | 54 | Legacy pediatric, adult, edge cases |
| `test-child-regression.R` | 7 | 49 | Structural invariants, frozen counts, spot checks, cross-sample stability, Missing handling, HC, prelim_infants compat |
| `test-child-parameters.R` | 10 | 24 | All configurable parameters |
| `test-child-edge-cases.R` | 12 | 23 | Single subject, sparse data, all-NA, mixed NA, SDE-Identical, negative agedays, HEADCM >3yr, extreme values, density mix, CF, deterministic |

Total (as of 2026-03-24): 33 tests, 150 assertions. Existing
`test-cdc.R` and `test-utils.R` not modified.

Run with:
```bash
NOT_CRAN=true Rscript -e \
  'testthat::test_file("tests/testthat/test-child-regression.R")'
```

### Runtime benchmarks

| Dataset | Size | Time |
|---------|------|------|
| syn.csv (full) | 77,721 rows, 2,697 subjects | ~1.2–1.9 min |
| syn.csv (1,036 subj subset) | 28,434 rows | ~18 sec |
| 500-subject subsample | ~14K rows | ~9 sec |

Always use `parallel = FALSE`. Run in background from Claude
Code (`run_in_background: true`).

---

## Known Issues

### Open

- [ ] **Debug `stop()` in Step 5:** Lines 2783–2792 have a
  hard stop on duplicate Include values after temp SDE.
  Should be `warning()`. The check is useful but shouldn't
  halt production runs.
- [ ] **Parallel processing broken on macOS:** `parallel = TRUE`
  produces incorrect results. Likely `data.table` + forked
  process incompatibility. Potential fixes: PSOCK clusters,
  `setDTthreads()`, or `future`/`furrr` framework.
- [ ] **Batch size hard-coded at 2,000:** Should be a
  `batch_size` parameter on `cleangrowth()`.
- [ ] **Batch-invariant operations inside loop:**
  `exclude.levels` definition, Tanner/WHO velocity reads,
  and `read_anthro()` calls are repeated per batch. Move
  before the loop.
- [ ] **Exclusion code rename pending:** All codes will be
  renamed at end of development. Legacy codes in
  `exclude.levels` should be cleaned up at the same time.
- [ ] **CF rescue thresholds:** Current thresholds (0.05/0.1
  with wholehalfimp, teen age cutoffs) may be tuned in
  future updates.
- [ ] **LENGTHCM → HEIGHTCM:** No measurement adjustment for
  the ~0.5–0.7 cm supine/standing difference. Known
  limitation; future update planned.
- [ ] **Unused variables in Step 19:** `abs_tbd.sd`,
  `abs_ctbd.sd`, `med_dop`, `med_cdop` are computed but
  never used. Leftover from refactoring. Harmless but
  should be cleaned up.

### Fixed (recent)

- [x] **Outer batching wrapper defeated:** Lines 440–441
  overwrote batch filter with all subjects. Fixed.
- [x] **Missing data bug:** NA measurements coded as
  SDE-EWMA instead of Missing. Root cause: duplicate
  `valid()` definitions in two files; R's alphabetical
  load order caused the unfixed copy to overwrite the
  fixed one. Fixed in v3.0.0.
- [x] **BIV threshold for preterm weights:** Previous
  threshold excluded legitimate preterm weights (0.7–1.0
  kg). Changed to <0.2 kg for agedays ≤ 365 and <1 kg
  for agedays > 365.
- [x] **CF rescue re-inclusion:** Rescued CFs are now set
  back to `exclude = "Include"` so they participate in
  downstream steps.
- [x] **`calc_and_recenter_z_scores()` blending bug:**
  CDC-only cutoff was `>= 4 years` instead of `> 5 years`.
  Confirmed bug — confused WHO/CDC blending window with
  corrected/uncorrected smoothing window.

---

## Coding Conventions

- **Sex:** 0 = male, 1 = female
- **Param values:** `"WEIGHTKG"`, `"HEIGHTCM"`, `"HEADCM"`
  (after input conversion)
- **Use `data.table`** throughout — consistent with gc internals
- **Rounding:** Stata-style rounding (`round_stata()`) has
  been removed from the child algorithm z-score pipeline.
  Z-scores are no longer rounded to 0.001 at intermediate
  steps. The function was only needed for Stata comparison.
- **Factor levels for `exclude`:** Assigning an unlisted
  string to a factor silently produces NA — always verify
  new codes exist in `exclude.levels`
- **Sort order:** `setkey(data.df, subjid, param, agedays, id)`
  — assumed by many steps. Re-sort after any modification that
  could change order.
- **`by`-group correctness:** `by = subjid` vs.
  `by = .(subjid, param)` — wrong grouping silently produces
  wrong results. Check carefully.
- **C++/Rcpp ruled out** — CRAN target, Rtools not always
  available on enterprise systems
- **Report weight before height** in summaries and output

---

## Pitfalls (for Claude agents)

Things that have caused bugs before or fail silently:

- **Do not create a second `valid()` definition.** Two copies
  in different files caused the missing data bug (R's
  alphabetical load order made the wrong copy win). Child
  version is in `child_clean.R`; legacy version is in
  `pediatric_support_legacy.R`. Do not add a third.
- **Do not use `parallel = TRUE`.** Produces incorrect results
  on macOS due to `data.table` + forked process issues.
- **Do not modify sort order without re-sorting.** Many steps
  assume `setkey(data.df, subjid, param, agedays, id)`. If
  you add/modify rows, call `setkey()` again.
- **Do not assign unlisted strings to the `exclude` factor.**
  Unlisted values silently become NA. Always verify new
  exclusion codes exist in `exclude.levels` first.
- **Do not use `by = subjid` when you need
  `by = .(subjid, param)`.** Wrong grouping produces wrong
  results with no error or warning.
- **No intermediate z-score rounding.** `round_stata()` has
  been removed from the child algorithm. Do not reintroduce
  rounding at intermediate z-score steps.

---

## GitHub and Distribution

- **GitHub:** `carriedaymont/growthcleanr`, branch
  `efficiency-updates`
- **CRAN:** v2.2.0 (does not include child algorithm or
  v3.0.0 changes)
- **Install from local:**
  ```r
  devtools::install_local("gc-github-latest",
    force = TRUE, upgrade = "never")
  ```
- **Install from GitHub:**
  ```r
  devtools::install_github("carriedaymont/growthcleanr",
    ref = "efficiency-updates",
    force = TRUE, upgrade = "never")
  ```
- v3.0.0 has not been released on GitHub or CRAN yet

---

## Child Narrative Document

The technical narrative for the child algorithm is in
`child-gc-narrative-2026-03-18.md` (in gc-github-latest).
It documents what the current R code does, serving as both
a debugging aid and long-term reference.

### Status

**Front matter: COMPLETE.** Covers Key Concepts, Architecture,
Data Requirements, Z-Score Infrastructure, Working/Output
Dataframe, Variable Glossary, Output Format (with full
exclusion code tables), Configurable Parameters, and Complete
List of Steps.

**Step-by-step documentation: COMPLETE for all steps.**
Steps 5, 6, 7, 9, 11, 13, 15, 16, 17, 19, 21, 22 all
documented with summary tables, overview, logic details,
rationale, and code review checklist findings.

**Open questions and issues:** Tracked in the narrative's
"Open Questions and Issues to Investigate" section. Items
marked `[x]` are resolved; `[ ]` items are still open.

### Structure

Each step section includes:
- Summary table (scope, prior/next step, exclusion codes)
- Overview
- Key terms and variable names
- Logic and implementation
- Rationale for selected decisions
- Code review checklist findings (12-point checklist)

---

## Terms Reference

This is the **canonical** terms reference (primary source).
`__Pipeline/CLAUDE.md` has a convenience copy; update both
if adding terms.

| Term | Full name |
|------|-----------|
| ewma | Exponentially weighted moving average |
| dewma | Deviation from EWMA |
| absdewma | Absolute value of dewma |
| spa | Subject-parameter-ageday combination |
| sde | Same-day extraneous value |
| dop | Designated other parameter |
| tbc | To-be-cleaned (recentered z-score) |
| ctbc | Corrected to-be-cleaned (recentered corrected z-score) |
| potcorr | Potentially correctable (prematurity) |
| biv | Biologically implausible value |
| csd | Conditional standard deviation |
| nnte | No need to evaluate (legacy; always FALSE in current code) |
| otl | Out of line (used in Evil Twins step; formerly "oob"/"out of bounds") |

---

## Running growthcleanr

### Standalone (from gc-github-latest directory)

```r
library(growthcleanr)
library(data.table)

dt <- fread("path/to/data.csv")
result <- cleangrowth(
  subjid = dt$subjid,
  param = dt$param,
  agedays = dt$agedays,
  sex = dt$sex,
  measurement = dt$measurement,
  parallel = FALSE
)
```

### From Pipeline (error-impact)

See `__Pipeline/CLAUDE.md` for pipeline-specific wrapper
instructions. Working directory must be `error-impact/` so
`here::here()` resolves correctly.
